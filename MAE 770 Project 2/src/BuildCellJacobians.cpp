#include "BuildCellJacobians.hpp"


LHS::LHS(
	const SimParameters& params, 
	const Mesh& mesh, 
	const Species& species, 
	const CellStateVars& state)
	: params(params), mesh(mesh), species(species), state(state) {
	matrix = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
}

void LHS::updateLHS() {

	// jacobian of conservative variables wrt primitive variables
	matrix.setZero();
	matrix.diagonal().head(params.nspecies).setOnes();
	calcPartial_Et_RhoS();
	calcPartial_Et_u();
	calcPartial_Et_T_tr();
	calcPartial_Et_Tv();
	calcPartial_Ev_RhoS();
	calcPartial_Ev_Tv();
}

void LHS::calcPartial_Et_RhoS() {

	double temp_tr = state.getTemp(cell_idx);
	double temp_V = state.getTemp_V(cell_idx);
	double rhoCV = state.getRhoCV_Mix(cell_idx);
	double rhoRs = state.getRhoR_Mix(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int row_idx = params.T_idx;
		int col_idx = species_idx;

		double jacobian_value = 0.0;

		jacobian_value += species.getIntEnergyTR(species_idx, temp_tr);
		jacobian_value += species.getFormationEnergy(species_idx);
		jacobian_value += (rhoCV / rhoRs) * species.getR_s(species_idx) * temp_tr;
		jacobian_value += state.getKineticEnergy(cell_idx) / state.getRho_tot(cell_idx);
		jacobian_value += species.getIntEnergyV(species_idx, temp_V);

		matrix(row_idx, col_idx) += jacobian_value;
	}
}

void LHS::calcPartial_Et_u() {

	int row_idx = params.T_idx;

	Eigen::VectorXd vel_components = state.getVel_components(cell_idx);
	double rho = state.getRho_tot(cell_idx);

	for (int i = 0; i < params.ndimension; ++i) {

		int col_idx = params.vel_idx + i;

		matrix(row_idx, col_idx) += rho * vel_components(i);
	}
}

void LHS::calcPartial_Et_T_tr() {

	int row_idx = params.T_idx;
	int col_idx = params.T_idx;

	double jacobian_value = 0.0;
	double temp_tr = state.getTemp(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV = species.getCV_s(species_idx, temp_tr);

		jacobian_value += rho_s * CV;
	}

	matrix(row_idx, col_idx) += jacobian_value;
}

void LHS::calcPartial_Et_Tv() {

	int row_idx = params.T_idx;
	int col_idx = params.Tv_idx;

	double jacobian_value = 0.0;
	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV_V = species.getCV_V(species_idx, temp_V);

		jacobian_value += rho_s * CV_V;
	}

	matrix(row_idx, col_idx) += jacobian_value;
}

void LHS::calcPartial_Ev_RhoS() {

	int row_idx = params.Tv_idx;

	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int col_idx = species_idx;
		matrix(row_idx, col_idx) += species.getIntEnergyV(species_idx, temp_V);
	}
}

void LHS::calcPartial_Ev_Tv() {

	int row_idx = params.Tv_idx;
	int col_idx = params.Tv_idx;

	double jacobian_value = 0.0;
	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV_V = species.getCV_V(species_idx, temp_V);

		jacobian_value += rho_s * CV_V;
	}

	matrix(row_idx, col_idx) += jacobian_value;
}


CellJacobians::CellJacobians(
	const SimParameters& params, 
	const Mesh& mesh, 
	const Species& species, 
	const CellStateVars& state) {

	cell_vec.resize(mesh.ncells, LHS(params, mesh, species, state));

	for (int i = 0; i < mesh.ncells; ++i) { // assign a cell index to each LHS struct
		cell_vec[i].cell_idx = i;
	}
}
