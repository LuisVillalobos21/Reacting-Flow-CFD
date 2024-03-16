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
			ref_species_mass_frac.resize(nspecies);

			if (getline(inputFile, line)) {
				std::istringstream massFractionStream(line);
				for (int i = 0; i < nspecies; ++i) {
					if (!(massFractionStream >> ref_species_mass_frac(i))) {
						std::cerr << "ERROR: Invalid mass fraction value" << std::endl;
						break;
					}
				}
			}
		}
		else if (keyword == "SPECIES_THERMO_DATA") {
			getline(inputFile, speciesThermoDataPath); // Assuming the path is on the next line
		}
		else if (keyword == "REFERENCE_VELOCITY") {
			if (!(iss >> ref_velocity)) {
				std::cerr << "ERROR: Invalid REFERENCE_VELOCITY value" << std::endl;
			}
		}
		else if (keyword == "REFERENCE_TEMPERATURE") {
			if (!(iss >> ref_temperature)) {
				std::cerr << "ERROR: Invalid REFERENCE_TEMPERATURE value" << std::endl;
			}
		}
		else if (keyword == "REFERENCE_MIXTURE_DENSITY") {
			if (!(iss >> ref_mixture_rho)) {
				std::cerr << "ERROR: Invalid REFERENCE_DENSITY value" << std::endl;
			}
		}
	}
	inputFile.close();

	ref_rho_s.resize(nspecies);
	for (int i = 0; i < nspecies; ++i) {
		ref_rho_s(i) = ref_mixture_rho * ref_species_mass_frac(i);
	}

}