#pragma once
#include "SimParameters.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"

struct LHS {

    Eigen::MatrixXd prim_var_jac;
    Eigen::MatrixXd src_term_jac;
    Eigen::MatrixXd quasi_1D_jac;

    Eigen::MatrixXd matrix;

    LHS(const SimParameters& params);
};

struct CellJacobians {

    std::vector<LHS> cell_vec;

    const SimParameters& params;
    const Mesh& mesh;
    const Species& species;
    const CellStateVars& state;

    CellJacobians(
        const SimParameters& params, 
        const Mesh& mesh,
        const Species& species, 
        const CellStateVars& state);

    void updateLHS();

    void addJacobians();

    void calcPrimVarJacobian(int cell_idx);

    void calcPartial_Et_RhoS(int cell_idx);

    void calcPartial_Et_u(int cell_idx);

    void calcPartial_Et_T_tr(int cell_idx);

    void calcPartial_Et_Tv(int cell_idx);

    void calcPartial_Ev_RhoS(int cell_idx);

    void calcPartial_Ev_Tv(int cell_idx);

    void calcQuasi1DJacobian(int cell_idx);

    void calcPartial_Quasi_RhoS(int cell_idx);

    void calcPartial_Quasi_T_tr(int cell_idx);

    void calcPartial_Relax_T_tr(int cell_idx);

    void calcPartial_Relax_T_v(int cell_idx);

    void calcSrcTermJacobian(int cell_idx);
};