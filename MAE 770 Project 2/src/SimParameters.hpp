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

    int nspecies_vib;
    Eigen::VectorXi vib_idxs;
    Eigen::VectorXd charact_temps_vib;

    std::string gridPath;
    std::string speciesThermoDataPath;
    std::string reactionPath;

    int num_time_steps;
    double CFL;
    double rel_tol;
    double ref_velocity;
    double ref_temperature;
    double ref_mixture_rho;
    double Pback;
    Eigen::VectorXd ref_species_mass_frac;
    Eigen::VectorXd ref_rho_s;

    SimParameters(const std::string& filename);

    void readInputFile(const std::string& filename);
};