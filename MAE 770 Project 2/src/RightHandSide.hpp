#pragma once
#include "CellStateVariables.hpp"

struct VariableVector {

	Eigen::VectorXd vec;

	VariableVector(const SimParameters& params);
};

struct CellResiduals {

	std::vector<VariableVector> cell_flux_vec;
	std::vector<VariableVector> cell_res_vec;

	const SimParameters& params;
	const CellStateVars& state;
	const Mesh& mesh;

	CellResiduals(const SimParameters& params, const Mesh& mesh, const CellStateVars& state);

	void updateRHS();

	Eigen::VectorXd calcResidual(int cell_idx) const;

	Eigen::VectorXd calcFluxLDFSS(int cell_idx) const;

	double calcQuasi_1DPressure(int cell_idx) const;
};
 