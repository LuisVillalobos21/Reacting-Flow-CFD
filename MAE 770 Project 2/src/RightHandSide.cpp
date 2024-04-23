#include "RightHandSide.hpp"

VariableVector::VariableVector(const SimParameters& params) {

	vec = Eigen::VectorXd::Zero(params.nvariables);
}

CellResiduals::CellResiduals(const SimParameters& params, const Mesh& mesh, const Species& species, const CellStateVars& state, const Chemistry& chem)
: params(params), mesh(mesh), state(state), species(species), chem(chem) {

	cell_res_vec.resize(mesh.jmax + 1, VariableVector(params));
	cell_flux_vec.resize(mesh.jmax + 1, VariableVector(params));
	cell_src_vec.resize(mesh.jmax + 1, VariableVector(params));
}

void CellResiduals::updateRHS() {

	for (int cell_idx = 0; cell_idx < mesh.jmax; ++cell_idx) {

		cell_flux_vec[cell_idx].vec = calcFaceFluxLDFSS(cell_idx);
	}

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		double area = mesh.getCellArea1D(cell_idx);

		cell_res_vec[cell_idx].vec = -calcFluxVec(cell_idx);

		cell_res_vec[cell_idx].vec(params.vel_idx) += calcQuasi_1DPressure(cell_idx);

		double temp = state.getTemp(cell_idx);

		if (temp > 1000) {
			cell_res_vec[cell_idx].vec(params.Tv_idx) += area * calcRelaxSrcTerm(cell_idx);
		}
	}
}

void CellResiduals::updateRHSChem() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		double temp = state.getTemp(cell_idx);
		double area = mesh.getCellArea1D(cell_idx);

		if (temp > 1000) {
			cell_res_vec[cell_idx].vec += area * calcChemSrcVec(cell_idx); // this needs to get multipled by area
		}
	}
}
 
Eigen::VectorXd CellResiduals::calcFluxVec(int cell_idx) {

	double dx = mesh.getCelldx1D(cell_idx);

	return (cell_flux_vec[cell_idx].vec - cell_flux_vec[cell_idx - 1].vec) / dx;
}

Eigen::VectorXd CellResiduals::calcFaceFluxLDFSS(int cell_idx) const {
	
	double ahalf = 0.5 * (state.getSoundSpeed(cell_idx) + state.getSoundSpeed(cell_idx + 1));

	double ul = state.getVel_components(cell_idx);
	double ur = state.getVel_components(cell_idx + 1);

	double xml = ul / ahalf;
	double xmr = ur / ahalf;

	double all = 0.5 * (1. + std::copysign(1.0, xml));
	double alr = 0.5 * (1. - std::copysign(1.0, xmr));

	double btl = -std::max(0.0, 1.0 - static_cast<double>(static_cast<int>(std::abs(xml))));
	double btr = -std::max(0.0, 1.0 - static_cast<double>(static_cast<int>(std::abs(xmr))));

	double xmml = 0.25 * (xml + 1.0) * (xml + 1.0); 
	double xmmr = -0.25 * (xmr - 1.0) * (xmr - 1.0); 

	double xmhalf = std::sqrt(0.5 * (xml * xml + xmr * xmr)); 
	double xmc = 0.25 * btl * btr * (xmhalf - 1.0) * (xmhalf - 1.0); 

	double Pl = state.getPressure(cell_idx);
	double Pr = state.getPressure(cell_idx + 1);
	double delp = Pl - Pr;
	double psum = Pl + Pr;

	double xmcp = xmc * std::max(0.0, (1.0 - (delp / psum + 2.0 * std::abs(delp) / Pl)));
	double xmcm = xmc * std::max(0.0, (1.0 + (delp / psum - 2.0 * std::abs(delp) / Pr)));
	double cvlp = all * (1.0 + btl) * xml - btl * xmml;
	double cvlm = alr * (1.0 + btr) * xmr - btr * xmmr;
	double cep = cvlp - xmcp;
	double cem = cvlm + xmcm;

	double fml = mesh.getFaceArea1D(cell_idx) * state.getRho(cell_idx) * ahalf * cep;
	double fmr = mesh.getFaceArea1D(cell_idx) * state.getRho(cell_idx + 1) * ahalf * cem;

	double ppl = 0.25 * (xml + 1.0) * (xml + 1.0) * (2.0 - xml);
	double ppr = 0.25 * (xmr - 1.0) * (xmr - 1.0) * (2.0 + xmr);

	double pnet = (all * (1.0 + btl) - btl * ppl) * Pl + (alr * (1.0 + btr) - btr * ppr) * Pr;

	Eigen::VectorXd flux = fmr * state.getFluxVars(cell_idx + 1) + fml * state.getFluxVars(cell_idx);
	flux(params.vel_idx) += mesh.getFaceArea1D(cell_idx) * pnet;

	return flux;
}

double CellResiduals::calcQuasi_1DPressure(int cell_idx) const {

	double pressure = state.getPressure(cell_idx);
	double face_area_l = mesh.getFaceArea1D(cell_idx - 1);
	double face_area_r = mesh.getFaceArea1D(cell_idx);
	double dx = mesh.getCelldx1D(cell_idx);

	return pressure * (face_area_r - face_area_l) / dx;
}

double CellResiduals::calcRelaxSrcTerm(int cell_idx) const {

	double value = 0.0;

	double temp_tr = state.getTemp(cell_idx);
	double temp_V = state.getTemp_V(cell_idx);

	for (int idx = 0; idx < params.nspecies_vib; ++idx) {

		int vib_idx = params.vib_idxs(idx);

		double rho_s = state.getRho_s(cell_idx, vib_idx);
		double ev_tr = species.getIntEnergyV(vib_idx, temp_tr);
		double ev_v = species.getIntEnergyV(vib_idx, temp_V);
		double tau = state.calcRelaxTime(cell_idx, vib_idx, idx);

		value += rho_s * (ev_tr - ev_v) / tau;
	}

	return value;
}

double CellResiduals::calcChemVibProd(int cell_idx) const {

	double value = 0.0;

	double temp_tr = state.getTemp(cell_idx);
	double temp_V = state.getTemp_V(cell_idx);
	Eigen::VectorXd species_con = state.cell_vec[cell_idx].concentrations;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double ev_v = species.getIntEnergyV(species_idx, temp_V);

		value += ev_v * chem.calcSpeciesProduction(species_con, temp_tr, temp_V, species_idx);
	}

	return value;
}

Eigen::VectorXd CellResiduals::calcChemSrcVec(int cell_idx) {

	Eigen::VectorXd src_vec = Eigen::VectorXd::Zero(params.nvariables);

	double temp_tr = state.getTemp(cell_idx);
	double temp_V = state.getTemp_V(cell_idx);
	Eigen::VectorXd species_con = state.cell_vec[cell_idx].concentrations;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		src_vec(species_idx) += chem.calcSpeciesProduction(species_con, temp_tr, temp_V, species_idx);
	}

	src_vec(params.Tv_idx) += calcChemVibProd(cell_idx);

	return src_vec;
}




