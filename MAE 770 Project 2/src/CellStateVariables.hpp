#pragma once

#include "ProjectIncludes.hpp"
#include "SimParameters.hpp"
#include "SpeciesThermo.hpp"
#include "MeshProccesing.hpp"

struct PrimVars {
    Eigen::VectorXd var_vec;

    Eigen::VectorXd mass_fracs;

    double Rho;

    double Pressure;

    PrimVars(const SimParameters& params);
};

struct CellStateVars {

    const SimParameters& params;
    const Species& species;
    const Mesh& mesh;

    std::vector<PrimVars> cell_vec;

    CellStateVars(const SimParameters& params, const Mesh& mesh, const Species& species);

    void initializeFlowFieldShock();

    void initializeFlowFieldSuperSonic();

    double getRho_s(int cell_idx, int species_idx) const;

    double getMassFrac(int cell_idx, int species_idx) const;

    Eigen::VectorXd calcMassfracs(int cell_idx) const;

    double calcRho(int cell_idx) const;

    double getRho(int cell_idx) const;

    double calcRhoOutFlow(int cell_idx) const;

    double calcPressure(int cell_idx) const;

    double getPressure(int cell_idx) const;

    Eigen::VectorXd getVel_components(int cell_idx) const;

    double getKineticEnergy(int cell_idx) const;

    double getTemp(int cell_idx) const;

    double getTemp_V(int cell_idx) const;

    double getIntEnergyVMix(int cell_idx) const;

    double getRhoCV_Mix(int cell_idx) const;

    double getRhoR_Mix(int cell_idx) const;

    double getSoundSpeed(int cell_idx) const;

    double calcTotalEnergy(int cell_idx) const;

    Eigen::VectorXd getFluxVars(int cell_idx) const;

    double calcReducedMw(int cell_idx, int species_idx1, int species_idx2) const;

    double calcA_relax(int cell_idx, int species_idx1, int species_idx2, int idx) const;

    double calcB_relax(int cell_idx, int species_idx1, int species_idx2) const;

    double calcRelaxTime(int cell_idx, int species_vib, int idx) const;
};