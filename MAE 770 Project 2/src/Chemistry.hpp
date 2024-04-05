#pragma once
#include "SpeciesThermo.hpp"
 

enum class ReactionType {
    Collision,
    Exchange,
};

struct Reaction {

    const SimParameters& params;

    ReactionType type;

    double alpha;
    double C_f;
    double eta;
    double theta_d;

    Eigen::VectorXd K_eq_coeff;
    double A1, A2, A3, A4, A5;

    Eigen::VectorXd TB_coeff;

    int num_reactants;
    int num_products;
    Eigen::VectorXi reactant_idxs;
    Eigen::VectorXi product_idxs;
    Eigen::VectorXi reactant_coeffs;
    Eigen::VectorXi product_coeffs;

    Reaction(const std::string& filepath, const SimParameters& params) 
    :params(params){

        readReactionFile(filepath);
    }

    Eigen::VectorXi readVector(std::ifstream& inputFile, int count) {
        Eigen::VectorXi values(count);
        for (int i = 0; i < count; ++i) {
            double tempValue; // Temporary variable to store the input
            if (!(inputFile >> tempValue)) {
                std::cerr << "ERROR: Invalid coefficient/index value" << std::endl;
                break;
            }
            values(i) = tempValue; // Assign the read value to the Eigen vector
        }
        return values;
    }


    void readReactionFile(const std::string& filepath) {
        std::ifstream inputFile(filepath);
        if (!inputFile.is_open()) {
            std::cerr << "ERROR: Could not open " << filepath << std::endl;
            return;
        }

        std::string line;
        bool isTB_Collision = false; 
        while (getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string keyword;
            iss >> keyword;

            if (keyword == "TYPE") {
                std::string value;
                iss >> value;
                if (value == "TB_COLLISION") {

                    type = ReactionType::Collision;
                    isTB_Collision = true; 
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
            else if (keyword == "THIRD_BODY_COEFF" && isTB_Collision) {
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

    double calcK_f(double temp) {

        return C_f * pow(temp, eta) * exp(-theta_d / temp);
    }

    double calcK_eq(double temp) {

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

    double calcK_b(double temp) {

        double k_f = calcK_f(temp);
        double k_eq = calcK_eq(temp);

        return k_f / k_eq;
    }
};


struct Chemistry {

    std::vector<Reaction> reaction_vec;

    const SimParameters& params;

    Chemistry(const std::string& inputListFile, const SimParameters& params) 
    :params(params){

        readReactionInput(inputListFile);
    }

    void readReactionInput(const std::string& filename) {
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
};