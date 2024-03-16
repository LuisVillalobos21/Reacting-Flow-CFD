#include "CellStateVariables.hpp"

PrimVars::PrimVars(const SimParameters& params) : var_vec(Eigen::VectorXd::Zero(params.nvariables)) {}

CellStateVars::CellStateVars(const SimParameters& params, const Mesh& mesh, const Species& species)
	:params(params), species(species) {
		cell_vec.resize(mesh.ncells, PrimVars(params));


}

void CellStateVars::initializeFlowField() {

}

double CellStateVars::getRho_s(int cell_idx, int species_idx) const {
	return cell_vec[cell_idx].var_vec[species_idx];
}

double CellStateVars::getRho_tot(int cell_idx) const {

	double rho = 0.0;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		rho += getRho_s(cell_idx, species_idx);
	}

	return rho;
}

Eigen::VectorXd CellStateVars::getVel_components(int cell_idx) const {

	Eigen::VectorXd vel_components = Eigen::VectorXd::Zero(params.ndimension);

	for (int i = 0; i < params.ndimension; ++i) {
		vel_components[i] = cell_vec[cell_idx].var_vec[params.vel_idx + i];
	}

	return vel_components;
}

double CellStateVars::getKineticEnergy(int cell_idx) const {

	Eigen::VectorXd vel_vec = getVel_components(cell_idx);
	return 0.5 * vel_vec.squaredNorm();
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
