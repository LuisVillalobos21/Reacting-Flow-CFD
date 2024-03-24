#include "RightHandSide.hpp"


VariableVector::VariableVector(const SimParameters& params) {

	vec = Eigen::VectorXd::Zero(params.nvariables);
}

CellResiduals::CellResiduals(const SimParameters& params, const Mesh& mesh, const CellStateVars& state)
	: params(params), mesh(mesh), state(state) {

	cell_flux_vec.resize(mesh.jmax + 1, VariableVector(params));
	cell_res_vec.resize(mesh.jmax + 1, VariableVector(params));
}

void CellResiduals::updateRHS() {

	for (int cell_idx = 0; cell_idx < mesh.jmax; ++cell_idx) {

		cell_flux_vec[cell_idx].vec.setZero();
		cell_flux_vec[cell_idx].vec = calcFluxLDFSS(cell_idx);
	}

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {
		
		cell_res_vec[cell_idx].vec = calcResidual(cell_idx);
		cell_res_vec[cell_idx].vec(params.vel_idx) += calcQuasi_1DPressure(cell_idx);
	}
}

Eigen::VectorXd CellResiduals::calcResidual(int cell_idx) const {

	return mesh.getFaceArea1D(cell_idx) * cell_flux_vec[cell_idx].vec - mesh.getFaceArea1D(cell_idx - 1) * cell_flux_vec[cell_idx - 1].vec;
}

Eigen::VectorXd CellResiduals::calcFluxLDFSS(int cell_idx) const {
	
	double ahalf = 0.5 * (state.getSoundSpeed(cell_idx) + state.getSoundSpeed(cell_idx + 1));

	double ul = state.getVel_components(cell_idx)(0); // just the first velocity index for now
	double ur = state.getVel_components(cell_idx + 1)(0);

	double xml = ul / ahalf;
	double xmr = ur / ahalf;

	double all = 0.5 * (1. + std::copysign(1.0, xml));
	double alr = 0.5 * (1. - std::copysign(1.0, xmr));

	double btl = -std::max(0.0, 1.0 - static_cast<double>(static_cast<int>(std::abs(xml))));
	double btr = -std::max(0.0, 1.0 - static_cast<double>(static_cast<int>(std::abs(xmr))));

	double xmml = 0.25 * std::pow(xml + 1.0, 2);
	double xmmr = -0.25 * std::pow(xmr - 1.0, 2);

	double xmhalf = std::sqrt(0.5 * (xml * xml + xmr * xmr));
	double xmc = 0.25 * btl * btr * std::pow(xmhalf - 1.0, 2);

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

	double fml = state.getRho(cell_idx) * ahalf * cep;
	double fmr = state.getRho(cell_idx + 1) * ahalf * cem;

	double ppl = 0.25 * std::pow(xml + 1.0, 2) * (2.0 - xml);
	double ppr = 0.25 * std::pow(xmr - 1.0, 2) * (2.0 + xmr);

	double pnet = (all * (1.0 + btl) - btl * ppl) * Pl + (alr * (1.0 + btr) - btr * ppr) * Pr;

	Eigen::VectorXd left_vars = fml * state.getFluxVars(cell_idx); 
	Eigen::VectorXd right_vars = fmr * state.getFluxVars(cell_idx + 1);

	Eigen::VectorXd flux = right_vars - left_vars;
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