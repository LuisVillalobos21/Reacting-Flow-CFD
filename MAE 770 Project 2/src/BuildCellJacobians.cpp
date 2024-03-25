#include "BuildCellJacobians.hpp"


LHS::LHS(const SimParameters& params) {

	matrix = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
}

CellJacobians::CellJacobians(

	const SimParameters& params,
	const Mesh& mesh,
	const Species& species,
	const CellStateVars& state)

	: params(params), mesh(mesh), species(species), state(state) {

	cell_vec.resize(mesh.jmax + 1, LHS(params));
}


void CellJacobians::updateLHS() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		cell_vec[cell_idx].matrix.setZero();
		cell_vec[cell_idx].matrix.diagonal().head(params.nspecies).setOnes();
		cell_vec[cell_idx].matrix(params.vel_idx, params.vel_idx) = state.getRho(cell_idx);

		for (int dim = 0; dim < params.ndimension; ++dim) {

			double vel = state.getVel_components(cell_idx)(dim);

			cell_vec[cell_idx].matrix.block(params.vel_idx, 0, 1, params.nspecies) =
				vel * Eigen::RowVectorXd::Ones(params.nspecies);
		}

		calcPrimVarJacobian(cell_idx);

		//std::cout << "matrix before multipled by cell area: " << '\n';
		//std::cout << cell_vec[cell_idx].matrix << '\n';

		cell_vec[cell_idx].matrix *= mesh.getCellArea1D(cell_idx);

		//std::cout << "cell area: " << mesh.getCellArea1D(cell_idx) << '\n';
		//std::cout << cell_vec[cell_idx].matrix << '\n';

	}
}

void CellJacobians::calcPrimVarJacobian(int cell_idx) {

	calcPartial_Et_RhoS(cell_idx);
	calcPartial_Et_u(cell_idx);
	calcPartial_Et_T_tr(cell_idx);
	calcPartial_Et_Tv(cell_idx);
	calcPartial_Ev_RhoS(cell_idx);
	calcPartial_Ev_Tv(cell_idx);
}

void CellJacobians::calcPartial_Et_RhoS(int cell_idx) {

	double temp_tr = state.getTemp(cell_idx);
	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int row_idx = params.T_idx;
		int col_idx = species_idx;

		double jacobian_value = 0.0;

		jacobian_value += species.getIntEnergyTR(species_idx, temp_tr);
		jacobian_value += species.getFormationEnergy(species_idx);
		jacobian_value += state.getKineticEnergy(cell_idx) / state.getRho(cell_idx);
		jacobian_value += species.getIntEnergyV(species_idx, temp_V);

		cell_vec[cell_idx].matrix(row_idx, col_idx) += jacobian_value;
	}
}

void CellJacobians::calcPartial_Et_u(int cell_idx) {

	int row_idx = params.T_idx;

	Eigen::VectorXd vel_components = state.getVel_components(cell_idx);
	double rho = state.getRho(cell_idx);

	for (int i = 0; i < params.ndimension; ++i) {

		int col_idx = params.vel_idx + i;

		cell_vec[cell_idx].matrix(row_idx, col_idx) += rho * vel_components(i); // no longer multipled by 1000
	}
}

void CellJacobians::calcPartial_Et_T_tr(int cell_idx) {

	int row_idx = params.T_idx;
	int col_idx = params.T_idx;

	double jacobian_value = 0.0;
	double temp_tr = state.getTemp(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV = species.getCV_s(species_idx, temp_tr);

		jacobian_value += rho_s * CV;
	}

	cell_vec[cell_idx].matrix(row_idx, col_idx) += jacobian_value;
}

void CellJacobians::calcPartial_Et_Tv(int cell_idx) {

	int row_idx = params.T_idx;
	int col_idx = params.Tv_idx;

	double jacobian_value = 0.0;
	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV_V = species.getCV_V(species_idx, temp_V);

		jacobian_value += rho_s * CV_V;
	}

	cell_vec[cell_idx].matrix(row_idx, col_idx) += jacobian_value;
}

void CellJacobians::calcPartial_Ev_RhoS(int cell_idx) {

	int row_idx = params.Tv_idx;

	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int col_idx = species_idx;
		cell_vec[cell_idx].matrix(row_idx, col_idx) += species.getIntEnergyV(species_idx, temp_V);
	}
}

void CellJacobians::calcPartial_Ev_Tv(int cell_idx) {

	int row_idx = params.Tv_idx;
	int col_idx = params.Tv_idx;

	double jacobian_value = 0.0;
	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		double rho_s = state.getRho_s(cell_idx, species_idx);
		double CV_V = species.getCV_V(species_idx, temp_V);

		jacobian_value += rho_s * CV_V;
	}

	cell_vec[cell_idx].matrix(row_idx, col_idx) += jacobian_value;
}


