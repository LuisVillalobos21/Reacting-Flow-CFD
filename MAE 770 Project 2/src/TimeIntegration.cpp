#include "TimeIntegration.hpp"

TimeEvolveSolution::TimeEvolveSolution(const SimParameters& params,
	const Mesh& mesh,
	const Species& species,
	CellStateVars& state)

	: params(params), mesh(mesh), species(species), state(state), 
	SolutionJacobians(params, mesh, species, state),
	SolutionResiduals(params, mesh, state) {

	cell_vec.resize(mesh.jmax + 1, TimeEvolveCell());
}

double TimeEvolveSolution::computeTimeStep(int cell_idx) {

	double dx = mesh.getCelldx1D(cell_idx);
	double u = state.getVel_components(cell_idx).norm();
	double a = state.getSoundSpeed(cell_idx);

	return (params.CFL * dx) / (u + a);
}

void TimeEvolveSolution::updateTimeSteps() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		cell_vec[cell_idx].dt = computeTimeStep(cell_idx);
	}
}

void TimeEvolveSolution::LinSolveCells() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		Eigen::VectorXd vec = SolutionResiduals.cell_res_vec[cell_idx].vec;
		auto luDecomp = SolutionJacobians.cell_vec[cell_idx].matrix.partialPivLu();

		Eigen::VectorXd solution = luDecomp.solve(vec);

		state.cell_vec[cell_idx].var_vec += solution;
	}
}

void TimeEvolveSolution::evolveCells() {

	updateTimeSteps();

	SolutionJacobians.updateLHS();

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		SolutionJacobians.cell_vec[cell_idx].matrix /= cell_vec[cell_idx].dt;
	}

	SolutionResiduals.updateRHS();

	LinSolveCells();
}

void TimeEvolveSolution::updateDerivedVars() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		state.cell_vec[cell_idx].Pressure = state.calcPressure(cell_idx);
		state.cell_vec[cell_idx].Rho = state.calcRho(cell_idx);
		state.cell_vec[cell_idx].mass_fracs = state.calcMassfracs(cell_idx);
	}
}

void TimeEvolveSolution::updateOutFlowBC() {

	state.cell_vec[mesh.jmax].var_vec = state.cell_vec[mesh.jmax - 1].var_vec;
	state.cell_vec[mesh.jmax].Pressure = state.cell_vec[mesh.jmax - 1].Pressure;
	state.cell_vec[mesh.jmax].Rho = state.cell_vec[mesh.jmax - 1].Rho;
	state.cell_vec[mesh.jmax].mass_fracs = state.cell_vec[mesh.jmax - 1].mass_fracs;
}


void TimeEvolveSolution::solve() {

	for (int step = 0; step < params.num_time_steps; ++step) {

		evolveCells();
		updateDerivedVars();
		updateOutFlowBC();

		//std::cout << "Current time step: " << step + 1 << '\n';
	}
}