#include "SimParameters.hpp"

SimParameters::SimParameters(const std::string& filename, Mesh& mesh, Species& species)
	: ndimension(0), nspecies(0), nvariables(0) {

	readInputFile(filename);
	mesh.readMesh(gridPath);
	species.readSpeciesInput(speciesThermoDataPath);

	nvariables = nspecies + ndimension + 2;
	vel_idx = nspecies;
	T_idx = nspecies + ndimension;
	Tv_idx = T_idx + 1;
}

void SimParameters::readInputFile(const std::string& filename) {
	std::ifstream inputFile(filename);
	if (!inputFile.is_open()) {
		std::cerr << "ERROR: Could not open " << filename << std::endl;
		return;
	}

	std::string line;
	while (getline(inputFile, line)) {
		std::istringstream iss(line);
		std::string keyword;
		iss >> keyword;

		if (keyword == "DIMENSION") {
			if (!(iss >> ndimension)) {
				std::cerr << "ERROR: Invalid DIMENSION value" << std::endl;
			}
		}
		else if (keyword == "GRID_PATH") {
			getline(inputFile, gridPath); // Assuming the path is on the next line
		}
		else if (keyword == "NUMBER_SPECIES") {
			if (!(iss >> nspecies)) {
				std::cerr << "ERROR: Invalid NUMBER_SPECIES value" << std::endl;
			}
		}
		else if (keyword == "SPECIES_THERMO_DATA") {
			getline(inputFile, speciesThermoDataPath); // Assuming the path is on the next line
		}
	}

	inputFile.close();
}