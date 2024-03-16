#pragma once
#include "SimParameters.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"

struct LHS {
    int cell_idx;
    Eigen::MatrixXd matrix;
    const SimParameters& params;
    const Mesh& mesh;
    const Species& species;
    const CellStateVars& state;

    LHS(
        const SimParameters& params, 
        const Mesh& mesh, 
        const Species& species, 
        const CellStateVars& state);

    void updateLHS();

    void calcPartial_Et_RhoS();

    void calcPartial_Et_u();

    void calcPartial_Et_T_tr();

    void calcPartial_Et_Tv();

    void calcPartial_Ev_RhoS();

    void calcPartial_Ev_Tv();
};

struct CellJacobians {
    std::vector<LHS> cell_vec;

    CellJacobians(
        const SimParameters& params, 
        const Mesh& mesh,
        const Species& species, 
        const CellStateVars& state);
};