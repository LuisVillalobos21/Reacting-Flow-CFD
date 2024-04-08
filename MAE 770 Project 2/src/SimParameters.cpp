#include "SimParameters.hpp"

SimParameters::SimParameters(const std::string& filename)
	: ndimension(0), nspecies(0), nvariables(0) {

	readInputFile(filename);

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
		else if (keyword == "NUMBER_VIB_SPECIES") {
			if (!(iss >> nspecies_vib)) {
				std::cerr << "ERROR: Invalid NUMBER_VIB_SPECIES value" << std::endl;
				continue;
			}
		}
		else if (keyword == "VIB_SPECIES_IDX") {
			vib_idxs.resize(nspecies_vib); 
			for (int i = 0; i < nspecies_vib; ++i) {
				if (!(iss >> vib_idxs(i))) {
					std::cerr << "ERROR: Invalid VIB_SPECIES_IDX value" << std::endl;
					break;
				}
				vib_idxs(i) -= 1;
			}
		}
		else if (keyword == "CHARACT_TEMPS_VIB") {
			charact_temps_vib.resize(nspecies_vib); 
			for (int i = 0; i < nspecies_vib; ++i) {
				if (!(iss >> charact_temps_vib(i))) {
					std::cerr << "ERROR: Invalid CHARACT_TEMPS_VIB value" << std::endl;
					break;
				}
			}
		}

		else if (keyword == "SPECIES_THERMO_DATA") {
			getline(inputFile, speciesThermoDataPath); // Assuming the path is on the next line
		}

		else if (keyword == "REACTIONS_PATH") {
			getline(inputFile, reactionPath); // Assuming the path is on the next line
		}

		else if (keyword == "NUM_REACTIONS") {
			if (!(iss >> num_reactions)) {
				std::cerr << "ERROR: Invalid NUM_REACTIONS value" << std::endl;
			}
		}

		else if (keyword == "SPECIES_PROD_COMBO") {
			spec_react_comp.resize(num_reactions); 
			for (int r = 0; r < num_reactions; ++r) {
				if (!getline(inputFile, line)) break; 

				std::istringstream lineStream(line);
				std::vector<int> tempCombo;
				int val;
				while (lineStream >> val) {
					tempCombo.push_back(val);
				}

				Eigen::VectorXi combo(tempCombo.size());
				for (size_t i = 0; i < tempCombo.size(); ++i) {
					combo[i] = tempCombo[i] - 1;
				}
				spec_react_comp[r].type = combo; 
			}
		}
		else if (keyword == "SPECIES_PROD_COEFF") {
			spec_react_coeff.resize(num_reactions); 
			for (int r = 0; r < num_reactions; ++r) {
				if (!getline(inputFile, line)) break; 

				std::istringstream lineStream(line);
				std::vector<double> tempCoeff;
				double val;
				while (lineStream >> val) {
					tempCoeff.push_back(val);
				}

				Eigen::VectorXd coeff(tempCoeff.size());
				for (size_t i = 0; i < tempCoeff.size(); ++i) {
					coeff[i] = tempCoeff[i];
				}
				spec_react_coeff[r].coeff = coeff; 
			}
		}

		else if (keyword == "CFL") {
			if (!(iss >> CFL)) {
				std::cerr << "ERROR: Invalid CFL value" << std::endl;
			}
		}
		else if (keyword == "RELATIVE_TOLERANCE") {
			if (!(iss >> rel_tol)) {
				std::cerr << "ERROR: Invalid RELATIVE_TOLERANCE value" << std::endl;
			}
		}
		else if (keyword == "NUM_TIME_STEPS") {
			if (!(iss >> num_time_steps)) {
				std::cerr << "ERROR: Invalid NUM_TIME_STEPS value" << std::endl;
			}
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
		else if (keyword == "PRESSURE_OUTLET") {
			if (!(iss >> Pback)) {
				std::cerr << "ERROR: Invalid PRESSURE_OUTLET value" << std::endl;
			}
		}
	}
	inputFile.close();

	ref_rho_s.resize(nspecies);
	for (int i = 0; i < nspecies; ++i) {
		ref_rho_s(i) = ref_mixture_rho * ref_species_mass_frac(i);
	}

}