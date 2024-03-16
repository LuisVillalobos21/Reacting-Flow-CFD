#pragma once
#include "SimParameters.hpp"
#include "SpeciesThermo.hpp"

struct PrimVars {
    Eigen::VectorXd var_vec;

    PrimVars(const SimParameters& params);
};

struct CellStateVars {

    const SimParameters& params;
    const Species& species;

    std::vector<PrimVars> cell_vec;

    CellStateVars(const SimParameters& params, const Mesh& mesh, const Species& species);

    void initializeFlowField();

    double getRho_s(int cell_idx, int species_idx) const;

    double getRho_tot(int cell_idx) const;

    Eigen::VectorXd getVel_components(int cell_idx) const;

    double getKineticEnergy(int cell_idx) const;

    double getTemp(int cell_idx) const;

    double getTemp_V(int cell_idx) const;

    double getRhoCV_Mix(int cell_idx) const;

    double getRhoR_Mix(int cell_idx) const;
};