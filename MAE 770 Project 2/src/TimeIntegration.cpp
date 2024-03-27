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

		if (cell_idx == mesh.jmax - 2) {
			Eigen::VectorXd vec = -SolutionResiduals.cell_res_vec[cell_idx].vec;
			auto luDecomp = SolutionJacobians.cell_vec[cell_idx].matrix.partialPivLu();

			Eigen::VectorXd solution = luDecomp.solve(vec);

			//std::cout << solution << '\n' << '\n';

			state.cell_vec[cell_idx].var_vec += solution;
		}

		Eigen::VectorXd vec = -SolutionResiduals.cell_res_vec[cell_idx].vec;
		auto luDecomp = SolutionJacobians.cell_vec[cell_idx].matrix.partialPivLu();

		Eigen::VectorXd solution = luDecomp.solve(vec);

		//std::cout << solution << '\n' << '\n';

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

		state.cell_vec[cell_idx].Rho = state.calcRho(cell_idx);
		state.cell_vec[cell_idx].mass_fracs = state.calcMassfracs(cell_idx);
		state.cell_vec[cell_idx].Pressure = state.calcPressure(cell_idx);
	}
}

void TimeEvolveSolution::updatePressureBoundary(int cell_idx) {

	int ghost_idx = cell_idx;
	int adjacent_idx = cell_idx - 1;

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

		state.cell_vec[ghost_idx].mass_fracs(species_idx) = state.getMassFrac(adjacent_idx, species_idx);
	}

	state.cell_vec[ghost_idx].var_vec(params.vel_idx) = state.cell_vec[adjacent_idx].var_vec(params.vel_idx);
	state.cell_vec[ghost_idx].var_vec(params.T_idx) = state.cell_vec[adjacent_idx].var_vec(params.T_idx);
	state.cell_vec[ghost_idx].var_vec(params.Tv_idx) = state.cell_vec[adjacent_idx].var_vec(params.Tv_idx);

	state.cell_vec[mesh.jmax].Rho = state.calcRhoOutFlow(ghost_idx);

	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		state.cell_vec[ghost_idx].var_vec(species_idx) = state.getRho(ghost_idx) * state.getMassFrac(ghost_idx, species_idx);
	}
}


void TimeEvolveSolution::solve() {

	for (int step = 0; step < params.num_time_steps; ++step) {

		evolveCells();
		updateDerivedVars();
		updatePressureBoundary(mesh.jmax);

		//std::cout << "Current time step: " << step + 1 << '\n';
	}
}