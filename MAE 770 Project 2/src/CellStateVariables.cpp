#include "CellStateVariables.hpp"

PrimVars::PrimVars(const SimParameters& params) 
	: var_vec(Eigen::VectorXd::Zero(params.nvariables)), mass_fracs(Eigen::VectorXd::Zero(params.nspecies)), Pressure(0.0), Rho(0.0) {}

CellStateVars::CellStateVars(const SimParameters& params, const Mesh& mesh, const Species& species)
	:params(params), mesh(mesh), species(species) {

	cell_vec.resize(mesh.jmax + 1, PrimVars(params));

	initializeFlowFieldShock();
	//initializeFlowFieldSuperSonic();
}

void CellStateVars::initializeFlowFieldShock() {

	//int middle = mesh.jmax + 1 - static_cast<int>((mesh.jmax + 1) / 2);
	int middle = static_cast<int>((mesh.jmax + 1) / 3);

	for (int cell_idx = 0; cell_idx < middle; ++cell_idx) {

		for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
			cell_vec[cell_idx].var_vec(species_idx) = params.ref_rho_s(species_idx);
			cell_vec[cell_idx].mass_fracs(species_idx) = params.ref_species_mass_frac(species_idx);
		}

		cell_vec[cell_idx].var_vec(params.vel_idx) = params.ref_velocity;
		cell_vec[cell_idx].var_vec(params.T_idx) = params.ref_temperature;
		cell_vec[cell_idx].var_vec(params.Tv_idx) = params.ref_temperature;

		cell_vec[cell_idx].Rho = calcRho(cell_idx);
		cell_vec[cell_idx].Pressure = calcPressure(cell_idx);
	}

	for (int cell_idx = middle; cell_idx < mesh.jmax + 1; ++cell_idx) {

		double temp_init = 4000;
		double mach = std::sqrt(1.4 * 287 * temp_init);

		cell_vec[cell_idx].var_vec(params.vel_idx) = mach * .75; // some subsonic velocity
		cell_vec[cell_idx].var_vec(params.T_idx) = temp_init;
		cell_vec[cell_idx].var_vec(params.Tv_idx) = temp_init;
		cell_vec[cell_idx].Pressure = params.Pback;
		cell_vec[cell_idx].Rho = params.Pback / (temp_init * 287);

		for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

			cell_vec[cell_idx].mass_fracs(species_idx) = params.ref_species_mass_frac(species_idx);
			cell_vec[cell_idx].var_vec(species_idx) = getRho(cell_idx) * getMassFrac(cell_idx, species_idx);
		}
	}
}

void CellStateVars::initializeFlowFieldSuperSonic() {

	for (int cell_idx = 0; cell_idx < mesh.jmax + 1; ++cell_idx) {

		for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
			cell_vec[cell_idx].var_vec(species_idx) = params.ref_rho_s(species_idx);
			cell_vec[cell_idx].mass_fracs(species_idx) = params.ref_species_mass_frac(species_idx);
		}

		cell_vec[cell_idx].var_vec(params.vel_idx) = params.ref_velocity;
		cell_vec[cell_idx].var_vec(params.T_idx) = params.ref_temperature;
		cell_vec[cell_idx].var_vec(params.Tv_idx) = params.ref_temperature;

		cell_vec[cell_idx].Rho = calcRho(cell_idx);
		cell_vec[cell_idx].Pressure = calcPressure(cell_idx);
	}

	cell_vec[mesh.jmax].Pressure = params.Pback;
}

double CellStateVars::getRho_s(int cell_idx, int species_idx) const {
	return cell_vec[cell_idx].var_vec[species_idx];
}

double CellStateVars::getRho(int cell_idx) const {

	return cell_vec[cell_idx].Rho;
}

double CellStateVars::calcRho(int cell_idx) const {

	double rho = 0.0;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		rho += getRho_s(cell_idx, species_idx);
	}
	return rho;
}

double CellStateVars::calcRhoOutFlow(int cell_idx) const {

	double rho = 0.0;
	double denom = 0.0;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		denom += (getMassFrac(cell_idx, species_idx) / (species.species_vector[species_idx].molecular_weight)) * R_U * getTemp(cell_idx);
	}

	return params.Pback / denom;
}

double CellStateVars::getMassFrac(int cell_idx, int species_idx) const {

	return cell_vec[cell_idx].mass_fracs[species_idx];
}

Eigen::VectorXd CellStateVars::calcMassfracs(int cell_idx) const {

	Eigen::VectorXd fracs = Eigen::VectorXd::Zero(params.nspecies);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		fracs(species_idx) = getRho_s(cell_idx, species_idx) / getRho(cell_idx);
	}

	return fracs;
}

double CellStateVars::calcPressure(int cell_idx) const {

	double P = 0.0;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		P += getRho_s(cell_idx, species_idx) * species.getR_s(species_idx) * getTemp(cell_idx);
	}

	return P; // no longer multipled by 1000
}

double CellStateVars::getPressure(int cell_idx) const {

	return cell_vec[cell_idx].Pressure;
}

double CellStateVars::calcTotalEnergy(int cell_idx) const {

	double e_t = 0.0;
	double temp = getTemp(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		e_t += getMassFrac(cell_idx, species_idx) * (species.getIntEnergyTR(species_idx, temp) + species.getR_s(species_idx) * temp);
		e_t += getMassFrac(cell_idx, species_idx) * species.getFormationEnergy(species_idx);
	}

	e_t += getKineticEnergy(cell_idx) / getRho(cell_idx);
	e_t += getIntEnergyVMix(cell_idx);

	return e_t;
}

double CellStateVars::getIntEnergyVMix(int cell_idx) const {

	double e_v = 0.0;
	double temp_V = getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		e_v += getMassFrac(cell_idx, species_idx) * species.getIntEnergyV(species_idx, temp_V);
	}

	return e_v;
}

Eigen::VectorXd CellStateVars::getFluxVars(int cell_idx) const {

	Eigen::VectorXd flux_vars = Eigen::VectorXd::Zero(params.nvariables);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		flux_vars(species_idx) = getMassFrac(cell_idx, species_idx);
	}

	flux_vars.segment(params.vel_idx, params.ndimension) = getVel_components(cell_idx);

	flux_vars(params.T_idx) = calcTotalEnergy(cell_idx);

	flux_vars(params.Tv_idx) = getIntEnergyVMix(cell_idx);

	return flux_vars;
}

Eigen::VectorXd CellStateVars::getVel_components(int cell_idx) const {
	return cell_vec[cell_idx].var_vec.segment(params.vel_idx, params.ndimension);
}

double CellStateVars::getKineticEnergy(int cell_idx) const {

	Eigen::VectorXd vel_vec = getVel_components(cell_idx);
	double rho = getRho(cell_idx);
	return (0.5 * rho * vel_vec.squaredNorm()); // no longer multiplied by 1000
}

double CellStateVars::getTemp(int cell_idx) const {
	return cell_vec[cell_idx].var_vec[params.T_idx];
}

double CellStateVars::getTemp_V(int cell_idx) const {
	return cell_vec[cell_idx].var_vec[params.Tv_idx];
}

double CellStateVars::getRhoCV_Mix(int cell_idx) const {

	double rhoCV = 0.0;
	double temp_tr = getTemp(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = getRho_s(cell_idx, species_idx);
		double CV = species.getCV_s(species_idx, temp_tr);

		rhoCV += rho_s * CV;
	}

	return rhoCV;
}

double CellStateVars::getRhoR_Mix(int cell_idx) const {

	double rhoRs = 0.0;
	double temp_tr = getTemp(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = getRho_s(cell_idx, species_idx);
		double R_s = species.getR_s(species_idx);

		rhoRs += rho_s * R_s;
	}

	return rhoRs;
}

double CellStateVars::getSoundSpeed(int cell_idx) const {

	double P = getPressure(cell_idx);
	double rho = getRho(cell_idx);
	double rhoR_mix = getRhoR_Mix(cell_idx);
	double rhoCV_mix = getRhoCV_Mix(cell_idx);

	double a2 = (P / rho) * (1 + rhoR_mix / rhoCV_mix); //rhoR and rhocv might be in J instead of KJ, units prolly cancel
	return sqrt(a2);
}

double CellStateVars::calcReducedMw(int cell_idx, int species_idx1, int species_idx2) const {

	double Wn = species.getMw(species_idx1);
	double Wm = species.getMw(species_idx2);

	return (Wm * Wn) / (Wm + Wn);
}

double CellStateVars::calcA_relax(int cell_idx, int species_idx1, int species_idx2, int idx) const {

	double mu = calcReducedMw(cell_idx, species_idx1, species_idx2);
	double mu_root = std::sqrt(mu);
	double theta_v_m = params.charact_temps_vib(idx);

	double theta_v_m_cbrt = std::cbrt(theta_v_m);
	double theta_v_m_4_3 = theta_v_m_cbrt * theta_v_m_cbrt * theta_v_m_cbrt * theta_v_m_cbrt;

	return 1.16e-3 * mu_root * theta_v_m_4_3;
}

double CellStateVars::calcB_relax(int cell_idx, int species_idx1, int species_idx2) const {

	double mu = calcReducedMw(cell_idx, species_idx1, species_idx2);
	double mu_root_fourth = std::sqrt(std::sqrt(mu));

	return 0.015 * mu_root_fourth;
}

double CellStateVars::calcRelaxTime(int cell_idx, int species_vib, int idx) const {

	double tau = 0.0;
	double YnWn = 0.0;
	double P_atm = getPressure(cell_idx) / 101325; 
	double temp_tr = getTemp(cell_idx);
	double term = 0.0;

	for (int species_n = 0; species_n < params.nspecies; ++species_n) {

		double Yn = getMassFrac(cell_idx, species_n);
		double Wn = species.getMw(species_n);
		YnWn += Yn / Wn;
	}

	for (int species_n = 0; species_n < params.nspecies; ++species_n) {

		double A = calcA_relax(cell_idx, species_n, species_vib, idx);
		double B = calcB_relax(cell_idx, species_n, species_vib);
		term += 1 / exp(A * (1 / std::cbrt(temp_tr) - B) - 18.42);
	}

	tau = YnWn / (P_atm * YnWn * term);

	return tau;
}





