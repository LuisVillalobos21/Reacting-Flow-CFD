#include "Chemistry.hpp"

Reaction::Reaction(const std::string& filepath, const SimParameters& params)
    :params(params) {

    readReactionFile(filepath);
}

Eigen::VectorXi Reaction::readVector(std::ifstream& inputFile, int count) {
    Eigen::VectorXi values(count);
    for (int i = 0; i < count; ++i) {
        int tempValue; // Temporary variable to store the input
        if (!(inputFile >> tempValue)) {
            std::cerr << "ERROR: Invalid coefficient/index value" << std::endl;
            break;
        }
        values(i) = tempValue - 1; // Assign the read value to the Eigen vector
    }
    return values;
}

void Reaction::readReactionFile(const std::string& filepath) {
    std::ifstream inputFile(filepath);
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: Could not open " << filepath << std::endl;
        return;
    }

    std::string line;
    bool isTB_Dissoc = false;
    while (getline(inputFile, line)) {
        std::istringstream iss(line);
        std::string keyword;
        iss >> keyword;

        if (keyword == "TYPE") {
            std::string value;
            iss >> value;
            if (value == "TB_DISSOCIATION") {

                type = ReactionType::Dissociation;
                isTB_Dissoc = true;
                TB_coeff.resize(params.nspecies);
                alpha = 1000.0; // hardcocde alpha 
            }
            else {
                type = ReactionType::Exchange;
                alpha = 1.0;
            }
        }

        else if (keyword == "NUM_REACTANTS") {
            iss >> num_reactants;
            product_coeffs.resize(num_reactants);
            product_idxs.resize(num_reactants);
        }
        else if (keyword == "NUM_PRODUCTS") {
            iss >> num_products;
            reactant_coeffs.resize(num_products);
            reactant_idxs.resize(num_products);
        }
        else if (keyword == "REACTANT_IDXS") {
            reactant_idxs = readVector(inputFile, num_reactants);
        }
        else if (keyword == "PRODUCT_IDXS") {
            product_idxs = readVector(inputFile, num_products);
        }
        else if (keyword == "REACTANT_COEFFS") {
            reactant_coeffs = readVector(inputFile, num_reactants);
        }
        else if (keyword == "PRODUCT_COEFFS") {
            product_coeffs = readVector(inputFile, num_products);
        }

        else if (keyword == "C_f") {
            iss >> C_f;
        }
        else if (keyword == "ETA") {
            iss >> eta;
        }
        else if (keyword == "THETA_d") {
            iss >> theta_d;
        }
        else if (keyword == "K_EQ_COEFF") {
            if (getline(inputFile, line)) {
                std::istringstream k_eq_coeffStream(line);
                K_eq_coeff.resize(5);
                for (int i = 0; i < 5; ++i) { // hard coding 5 k_eq coefficients
                    if (!(k_eq_coeffStream >> K_eq_coeff(i))) {
                        std::cerr << "ERROR: Invalid K_eq coefficient" << std::endl;
                        break;
                    }
                }
            }
        }
        else if (keyword == "THIRD_BODY_COEFF" && isTB_Dissoc) {
            if (getline(inputFile, line)) {
                std::istringstream thirdbodyStream(line);
                for (int i = 0; i < params.nspecies; ++i) {
                    if (!(thirdbodyStream >> TB_coeff(i))) {
                        std::cerr << "ERROR: Invalid third body coefficient" << std::endl;
                        break;
                    }
                }
            }
        }
    }

    A1 = K_eq_coeff(0);
    A2 = K_eq_coeff(1);
    A3 = K_eq_coeff(2);
    A4 = K_eq_coeff(3);
    A5 = K_eq_coeff(4);

    inputFile.close();
}

double Reaction::calcForwardRate(double temp, double temp_V) const {

    if (type != ReactionType::Dissociation) {
        return C_f * pow(temp, eta) * exp(-theta_d / temp);
    }

    double temp_combo = sqrt(temp * temp_V);

    return C_f * pow(temp_combo, eta) * exp(-theta_d / temp_combo);
}

double Reaction::calcEqRate(double temp) const {

    double z = 10000 / temp;
    double z2 = z * z;
    double z3 = z2 * z;
    double z4 = z3 * z;

    return alpha * exp(A1 +
        A2 * z +
        A3 * z2 +
        A4 * z3 +
        A5 * z4);
}

double Reaction::calcBackwardRate(double temp, double temp_V) const {

    double k_f = calcForwardRate(temp, temp_V);
    double k_eq = calcEqRate(temp);

    return k_f / k_eq;
}





Chemistry::Chemistry(
    const std::string& inputListFile, 
    const SimParameters& params, 
    const Species& species, 
    const CellStateVars& state)
    :params(params), species(species), state(state) {

    inverted_M_w.resize(params.nspecies);
    for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {
        double M_w = species.getMw(species_idx);
        inverted_M_w(species_idx) = 1.0 / M_w;
    }

    readReactionInput(inputListFile);
}

void Chemistry::readReactionInput(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: could not open " << filename << '\n';
        return;
    }

    std::string line;
    while (getline(inputFile, line)) {
        Reaction temp_reaction(line, params);
        reaction_vec.push_back(temp_reaction);
    }

    inputFile.close();
}

Eigen::VectorXd Chemistry::calcConcentrations(const Eigen::VectorXd& species_Rho) const {
    return species_Rho.cwiseProduct(inverted_M_w);
}

double Chemistry::calcTBFactor(const Reaction& reaction, const Eigen::VectorXd& concentrations) const {
    if (reaction.type != ReactionType::Dissociation) {
        return 1.0;
    }

    return reaction.TB_coeff.dot(concentrations);
}


double Chemistry::calcLawMassAction(
    const Reaction& reaction,
    const Eigen::VectorXd& species_Rho,
    double temp,
    double temp_V) const {

    Eigen::VectorXd concentrations = calcConcentrations(species_Rho);
    double TB_factor = calcTBFactor(reaction, concentrations);
    double k_f = reaction.calcForwardRate(temp, temp_V);
    double k_b = reaction.calcBackwardRate(temp, temp_V);

    double forward_concentration = (reaction.num_reactants > 0) ? 1.0 : 0.0;
    double backward_concentration = (reaction.num_products > 0) ? 1.0 : 0.0;

    for (int idx = 0; idx < reaction.num_reactants; ++idx) {
        int reactant_idx = reaction.reactant_idxs(idx);
        forward_concentration *= concentrations(reactant_idx);
    }

    for (int idx = 0; idx < reaction.num_products; ++idx) {
        int product_idx = reaction.product_idxs(idx);
        backward_concentration *= concentrations(product_idx);
    }

    return (k_f * forward_concentration - k_b * backward_concentration) * TB_factor;
}

double Chemistry::calcSpeciesProduction(
    Eigen::VectorXd& rho_vec,
    double temp_tr,
    double temp_V,
    int species_idx) const {

    Eigen::VectorXi idxs = params.spec_react_comp[species_idx].type;

    double omega_dot = 0;

    for (int i = 0; i < idxs.size(); ++i) {

        double coeff = params.spec_react_coeff[species_idx].coeff(i);

        omega_dot += coeff * calcLawMassAction(reaction_vec[idxs(i)], rho_vec, temp_tr, temp_V);
    }

    return species.getMw(species_idx) * omega_dot;
}