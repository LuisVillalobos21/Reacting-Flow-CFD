#pragma once
#include "ProjectIncludes.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"

struct SimParameters {

    int ndimension;
    int nspecies;
    int nvariables;
    int vel_idx;
    int T_idx;
    int Tv_idx;

    std::string gridPath;
    std::string speciesThermoDataPath;

    double ref_velocity;
    double ref_temperature;
    double ref_mixture_rho;
    Eigen::VectorXd ref_species_mass_frac;
    Eigen::VectorXd ref_rho_s;

    SimParameters(const std::string& filename, Mesh& mesh, Species& species);

    void readInputFile(const std::string& filename);
};