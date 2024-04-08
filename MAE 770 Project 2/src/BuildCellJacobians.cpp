#include "BuildCellJacobians.hpp"


LHS::LHS(const SimParameters& params) {

	prim_var_jac = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
	src_term_jac = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
	quasi_1D_jac = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
	matrix = Eigen::MatrixXd::Zero(params.nvariables, params.nvariables);
}

CellJacobians::CellJacobians(

	const SimParameters& params,
	const Mesh& mesh,
	const Species& species,
	const CellStateVars& state,
	const Chemistry& chem)

	: params(params), mesh(mesh), species(species), state(state), chem(chem) {

	cell_vec.resize(mesh.jmax + 1, LHS(params));
}


void CellJacobians::updateLHS() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		calcPrimVarJacobian(cell_idx);

		cell_vec[cell_idx].prim_var_jac *= mesh.getCellArea1D(cell_idx);

		calcQuasi1DJacobian(cell_idx);

		calcNonEqSrcTermJacobian(cell_idx);

		cell_vec[cell_idx].src_term_jac *= mesh.getCellArea1D(cell_idx);
	}
}

void CellJacobians::updateLHSChem() {
	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		calcChemSrcTermJacobian(cell_idx);
	}
}


void CellJacobians::addJacobians() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		cell_vec[cell_idx].matrix = cell_vec[cell_idx].prim_var_jac - cell_vec[cell_idx].src_term_jac - 
			cell_vec[cell_idx].quasi_1D_jac;
	}
}

void CellJacobians::calcPrimVarJacobian(int cell_idx) {

	cell_vec[cell_idx].prim_var_jac.setZero();
	cell_vec[cell_idx].prim_var_jac.diagonal().head(params.nspecies).setOnes();
	cell_vec[cell_idx].prim_var_jac(params.vel_idx, params.vel_idx) = state.getRho(cell_idx);

	for (int dim = 0; dim < params.ndimension; ++dim) {

		double vel = state.getVel_components(cell_idx)(dim);

		cell_vec[cell_idx].prim_var_jac.block(params.vel_idx, 0, 1, params.nspecies) =
			vel * Eigen::RowVectorXd::Ones(params.nspecies);
	}

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

		cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += jacobian_value;
	}
}

void CellJacobians::calcPartial_Et_u(int cell_idx) {

	int row_idx = params.T_idx;

	Eigen::VectorXd vel_components = state.getVel_components(cell_idx);
	double rho = state.getRho(cell_idx);

	for (int i = 0; i < params.ndimension; ++i) {

		int col_idx = params.vel_idx + i;

		cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += rho * vel_components(i); // no longer multipled by 1000
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

	cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += jacobian_value;
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

	cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += jacobian_value;
}

void CellJacobians::calcPartial_Ev_RhoS(int cell_idx) {

	int row_idx = params.Tv_idx;

	double temp_V = state.getTemp_V(cell_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int col_idx = species_idx;
		cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += species.getIntEnergyV(species_idx, temp_V);
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

	cell_vec[cell_idx].prim_var_jac(row_idx, col_idx) += jacobian_value;
}

void CellJacobians::calcQuasi1DJacobian(int cell_idx) {

	cell_vec[cell_idx].quasi_1D_jac.setZero();
	calcPartial_Quasi_RhoS(cell_idx);
	calcPartial_Quasi_T_tr(cell_idx);
}

void CellJacobians::calcPartial_Quasi_RhoS(int cell_idx) {

	double temp = state.getTemp(cell_idx);

	double face_area_l = mesh.getFaceArea1D(cell_idx - 1);
	double face_area_r = mesh.getFaceArea1D(cell_idx);
	double dx = mesh.getCelldx1D(cell_idx);
	double dA_dx = (face_area_r - face_area_l) / dx;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		int row_idx = params.vel_idx;
		int col_idx = species_idx;

		double jacobian_value = 0.0;

		jacobian_value = species.getR_s(species_idx) * temp * dA_dx;

		cell_vec[cell_idx].quasi_1D_jac(row_idx, col_idx) = jacobian_value;
	}
}

void CellJacobians::calcPartial_Quasi_T_tr(int cell_idx) {

	int row_idx = params.vel_idx;
	int col_idx = params.T_idx;

	double face_area_l = mesh.getFaceArea1D(cell_idx - 1);
	double face_area_r = mesh.getFaceArea1D(cell_idx);
	double dx = mesh.getCelldx1D(cell_idx);
	double dA_dx = (face_area_r - face_area_l) / dx;

	double jacobian_value = 0.0;

	double rhos_Rs = state.getRhoR_Mix(cell_idx);
	jacobian_value = rhos_Rs * dA_dx;

	cell_vec[cell_idx].quasi_1D_jac(row_idx, col_idx) = jacobian_value;
}

void CellJacobians::calcNonEqSrcTermJacobian(int cell_idx) {

	calcPartial_Relax_T_tr(cell_idx);
	calcPartial_Relax_T_v(cell_idx);
}


void CellJacobians::calcPartial_Relax_T_tr(int cell_idx) {

	int row_idx = params.Tv_idx;
	int col_idx = params.T_idx;

	double jacobian_value = 0.0;
	double temp_tr = state.getTemp(cell_idx);

	for (int idx = 0; idx < params.nspecies_vib; ++idx) {

		int vib_idx = params.vib_idxs(idx);

		double rho_s = state.getRho_s(cell_idx, vib_idx);
		double tau = state.calcRelaxTime(cell_idx, vib_idx, idx);
		double CVe_T = species.getCV_V(vib_idx, temp_tr);

		jacobian_value += rho_s * CVe_T / tau;
	}

	cell_vec[cell_idx].src_term_jac(row_idx, col_idx) = jacobian_value;
}

void CellJacobians::calcPartial_Relax_T_v(int cell_idx) {

	int row_idx = params.Tv_idx;
	int col_idx = params.Tv_idx;

	double jacobian_value = 0.0;
	double temp_V = state.getTemp_V(cell_idx);

	for (int idx = 0; idx < params.nspecies_vib; ++idx) {

		int vib_idx = params.vib_idxs(idx);

		double rho_s = state.getRho_s(cell_idx, vib_idx);
		double tau = state.calcRelaxTime(cell_idx, vib_idx, idx);
		double CVe_T = species.getCV_V(vib_idx, temp_V);

		jacobian_value += rho_s * CVe_T / tau;
	}

	cell_vec[cell_idx].src_term_jac(row_idx, col_idx) = -jacobian_value;
}

void CellJacobians::calcChemSrcTermJacobian(int cell_idx) {

	double temp = state.getTemp(cell_idx);

	if (temp > 1000) {
		//calcPartialOmega_Rho_s_Diag(cell_idx);
		calcPartialOmega_Rho_s_OffDiag(cell_idx);
	}
}

void CellJacobians::calcPartialOmega_Rho_s_Diag(int cell_idx) {
	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		int row_idx = species_idx;
		int col_idx = species_idx;

		double temp_tr = state.getTemp(cell_idx);
		double temp_V = state.getTemp_V(cell_idx);

		Eigen::VectorXd rho_vec = state.cell_vec[cell_idx].var_vec.segment(0, params.nspecies);
		double d_rho = abs(rho_vec(species_idx)) * sqrt(2.22e-9);


		Eigen::VectorXd rho_vec_plus = rho_vec;
		Eigen::VectorXd rho_vec_minus = rho_vec;

		rho_vec_plus(species_idx) += d_rho;
		rho_vec_minus(species_idx) -= d_rho;

		double left_value = chem.calcSpeciesProduction(rho_vec_minus, temp_tr, temp_V, species_idx);
		double right_value = chem.calcSpeciesProduction(rho_vec_plus, temp_tr, temp_V, species_idx);

		double jacobian_value = (right_value - left_value) / (2 * d_rho);

		cell_vec[cell_idx].src_term_jac(row_idx, col_idx) = jacobian_value * mesh.getCellArea1D(cell_idx);
	}
}

void CellJacobians::calcPartialOmega_Rho_s_OffDiag(int cell_idx) {
	for (int row_species_idx = 0; row_species_idx < params.nspecies; ++row_species_idx) {
		for (int col_species_idx = 0; col_species_idx < params.nspecies; ++col_species_idx) {

			//if (row_species_idx == col_species_idx) {
			//	continue;
			//}

			double temp_tr = state.getTemp(cell_idx);
			double temp_V = state.getTemp_V(cell_idx);

			Eigen::VectorXd rho_vec = state.cell_vec[cell_idx].var_vec.segment(0, params.nspecies);

			double d_rho = abs(rho_vec(col_species_idx)) * sqrt(2.22e-9);

			Eigen::VectorXd rho_vec_plus = rho_vec;
			Eigen::VectorXd rho_vec_minus = rho_vec;

			rho_vec_plus(col_species_idx) += d_rho;
			rho_vec_minus(col_species_idx) -= d_rho;

			double left_value = chem.calcSpeciesProduction(rho_vec_minus, temp_tr, temp_V, row_species_idx);
			double right_value = chem.calcSpeciesProduction(rho_vec_plus, temp_tr, temp_V, row_species_idx);

			double jacobian_value = (right_value - left_value) / (2 * d_rho);

			cell_vec[cell_idx].src_term_jac(row_species_idx, col_species_idx) = jacobian_value * mesh.getCellArea1D(cell_idx);
		}
	}
}




