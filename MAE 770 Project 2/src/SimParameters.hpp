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

    SimParameters(const std::string& filename, Mesh& mesh, Species& species);

    void readInputFile(const std::string& filename);
};