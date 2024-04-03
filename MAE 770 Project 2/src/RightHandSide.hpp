#pragma once
#include "CellStateVariables.hpp"

struct VariableVector {

	Eigen::VectorXd vec;

	VariableVector(const SimParameters& params);
};

struct CellResiduals {

	std::vector<VariableVector> cell_flux_vec;
	std::vector<VariableVector> cell_src_vec;
	std::vector<VariableVector> cell_q1D_vec;
	std::vector<VariableVector> cell_res_vec;

	const SimParameters& params;
	const CellStateVars& state;
	const Mesh& mesh;
	const Species& species;

	CellResiduals(const SimParameters& params, const Mesh& mesh, const Species& species, const CellStateVars& state);

	void updateRHS();

	Eigen::VectorXd calcFluxVec(int cell_idx);

	Eigen::VectorXd calcQ1DVec(int cell_idx);

	Eigen::VectorXd calcSrcVec(int cell_idx);

	Eigen::VectorXd calcFaceFluxLDFSS(int cell_idx) const;

	double calcRelaxSrcTerm(int cell_idx) const;

	double calcQuasi_1DPressure(int cell_idx) const;


};
 