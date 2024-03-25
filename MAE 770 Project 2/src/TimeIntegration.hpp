#pragma once

#include "SimParameters.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"
#include "BuildCellJacobians.hpp"
#include "RightHandSide.hpp"

struct TimeEvolveCell {

	double dt;
};

struct TimeEvolveSolution {
	
	const SimParameters& params;
	const Mesh& mesh;
	const Species& species;
	CellStateVars& state;

	CellJacobians SolutionJacobians;
	CellResiduals SolutionResiduals;

	std::vector<TimeEvolveCell> cell_vec;

	TimeEvolveSolution(const SimParameters& params,
		const Mesh& mesh,
		const Species& species,
		CellStateVars& state);

	double computeTimeStep(int cell_idx);

	void updateTimeSteps();

	void LinSolveCells();

	void evolveCells();

	void updateDerivedVars();

	void updateGhostCellsPressureBoundary();

	void solve();
};