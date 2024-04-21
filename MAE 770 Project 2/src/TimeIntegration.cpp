#include "TimeIntegration.hpp"

TimeEvolveSolution::TimeEvolveSolution(const SimParameters& params,
	const Mesh& mesh,
	const Species& species,
	CellStateVars& state,
	const Chemistry& chem)

	: params(params), mesh(mesh), species(species), state(state), chem(chem),
	SolutionJacobians(params, mesh, species, state, chem),
	SolutionResiduals(params, mesh, species, state, chem) {

	cell_vec.resize(mesh.jmax + 1, TimeEvolveCell());

	res_vec_momentum.resize(mesh.jmax + 1);
	res_vec_momentum.setZero();

	res_vec_energy.resize(mesh.jmax + 1);
	res_vec_energy.setZero();

	res_vec_vibe.resize(mesh.jmax + 1);
	res_vec_vibe.setZero();

	inverted_M_w.resize(params.nspecies);
	for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
		double M_w = species.getMw(species_idx);
		inverted_M_w(species_idx) = 1.0 / M_w;
	}
}

double TimeEvolveSolution::computeTimeStep(int cell_idx) {

	double dx = mesh.getCelldx1D(cell_idx);
	double u = state.getVel_components(cell_idx);
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

	//if (chem_switch == 1) {
	//	SolutionJacobians.updateLHSChem();
	//}

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		SolutionJacobians.cell_vec[cell_idx].prim_var_jac /= cell_vec[cell_idx].dt;
	}

	SolutionJacobians.addJacobians();

	SolutionResiduals.updateRHS();

	//if (chem_switch == 1) {
	//	SolutionResiduals.updateRHSChem();
	//}

	LinSolveCells();
}

void TimeEvolveSolution::updateDerivedVars() {

	for (int cell_idx = 1; cell_idx < mesh.jmax; ++cell_idx) {

		state.cell_vec[cell_idx].Rho = state.calcRho(cell_idx);
		state.cell_vec[cell_idx].mass_fracs = state.calcMassfracs(cell_idx);
		state.cell_vec[cell_idx].concentrations = state.cell_vec[cell_idx].var_vec.segment(0, params.nspecies).cwiseProduct(inverted_M_w);
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

void TimeEvolveSolution::calcConvergence(int step) {

	for (int i = 1; i < mesh.jmax; ++i) {

		res_vec_momentum(i) = SolutionResiduals.cell_res_vec[i].vec(params.vel_idx);
		res_vec_energy(i) = SolutionResiduals.cell_res_vec[i].vec(params.T_idx);
		res_vec_vibe(i) = SolutionResiduals.cell_res_vec[i].vec(params.Tv_idx);
	}

	if (step == 0) {
		resnorm_momentum_0 = res_vec_momentum.norm();
		resnorm_energy_0 = res_vec_energy.norm();
		resnorm_vibe_0 = res_vec_vibe.norm();

		resnorm_momentum = res_vec_momentum.norm() / resnorm_momentum_0;
		resnorm_energy = res_vec_energy.norm() / resnorm_energy_0;
		resnorm_vibe = res_vec_vibe.norm() / resnorm_vibe_0;

		residual_history_momentum.emplace_back(resnorm_momentum);
		residual_history_energy.emplace_back(resnorm_energy);
		residual_history_vibe.emplace_back(resnorm_vibe);

		return;
	}

	resnorm_momentum = res_vec_momentum.norm() / resnorm_momentum_0;
	resnorm_energy = res_vec_energy.norm() / resnorm_energy_0;
	resnorm_vibe = res_vec_vibe.norm() / resnorm_vibe_0;

	residual_history_momentum.emplace_back(resnorm_momentum);
	residual_history_energy.emplace_back(resnorm_energy);
	residual_history_vibe.emplace_back(resnorm_vibe);
}



void TimeEvolveSolution::solve() {

	bool chem_operations_done = false;

	auto start = std::chrono::high_resolution_clock::now();

	for (int step = 0; step < params.num_time_steps; ++step) {

		if (resnorm_momentum < params.rel_tol && resnorm_energy < params.rel_tol && resnorm_vibe < params.rel_tol) {
			break;
		}

		if (chem_switch == 0 && resnorm_momentum < chem_tol && resnorm_energy < chem_tol && resnorm_vibe < chem_tol) {
			chem_switch = 1;
		}

		evolveCells();
		calcConvergence(step);
		updateDerivedVars();
		updatePressureBoundary(mesh.jmax);

		if (step % 100 == 0 || step == params.num_time_steps - 1) {
			std::cout << std::scientific; 
			std::cout << "Time step: " << step + 1;
			std::cout << ", Mom. eqn. resnorm: " << resnorm_momentum << ", Ener. eqn resnorm: " <<
				resnorm_energy << ", Vibe eqn resnorm: " << resnorm_vibe << '\n';
			std::cout << std::defaultfloat; 
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = end - start;

	std::cout << '\n';

	if (resnorm_momentum < params.rel_tol && resnorm_energy < params.rel_tol && resnorm_vibe < params.rel_tol) {

		std::cout << std::scientific;
		std::cout << "Residuals converged below specified tolerance: " << params.rel_tol << "\n\n";
		std::cout << std::defaultfloat;
	}
	else {

		std::cout << std::scientific;
		std::cout << "Residuals did not fall below specified tolerance: " << params.rel_tol << "\n\n";
		std::cout << std::defaultfloat;
	}

	std::cout << "Solving process complete!" << '\n' << '\n';
	std::cout << "Solve duration: " << std::fixed << std::setprecision(4) << elapsed.count() << " seconds" << '\n' << '\n';
}