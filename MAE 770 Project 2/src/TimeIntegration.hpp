#pragma once

#include "ProjectIncludes.hpp"
#include "SimParameters.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"
#include "BuildCellJacobians.hpp"
#include "RightHandSide.hpp"
#include "Chemistry.hpp"

struct TimeEvolveCell {

	double dt;
};

struct TimeEvolveSolution {
	
	const SimParameters& params;
	const Mesh& mesh;
	const Species& species;
	CellStateVars& state;
	const Chemistry& chem;

	CellJacobians SolutionJacobians;
	CellResiduals SolutionResiduals;

	std::vector<TimeEvolveCell> cell_vec;

	int chem_switch = 0;
	double chem_tol = 1e-2;

	double resnorm_momentum_0;
	double resnorm_energy_0;
	double resnorm_vibe_0;

	double resnorm_momentum = 1.0;
	double resnorm_energy = 1.0;
	double resnorm_vibe = 1.0;

	Eigen::VectorXd res_vec_momentum;
	Eigen::VectorXd res_vec_energy;
	Eigen::VectorXd res_vec_vibe;

	Eigen::VectorXd residual_history_momentum;
	Eigen::VectorXd residual_history_energy;
	Eigen::VectorXd residual_history_vibe;

	TimeEvolveSolution(const SimParameters& params,
		const Mesh& mesh,
		const Species& species,
		CellStateVars& state,
		const Chemistry& chem);

	double computeTimeStep(int cell_idx);

	void updateTimeSteps();

	void LinSolveCells();

	void evolveCells();

	void updateDerivedVars();

	void updatePressureBoundary(int cell_idx);

	void calcConvergence(int step);

	void solve();
};