#pragma once
#include "SpeciesThermo.hpp"
#include "SimParameters.hpp"

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

    Reaction(const std::string& filepath, const SimParameters& params);

    Eigen::VectorXi readVector(std::ifstream& inputFile, int count);

    void readReactionFile(const std::string& filepath);

    double calcForwardRate(double temp) const;

    double calcEqRate(double temp) const;

    double calcBackwardRate(double temp) const;
};

struct Chemistry {

    std::vector<Reaction> reaction_vec;

    const SimParameters& params;
    const Species& species;

    Chemistry(const std::string& inputListFile, const SimParameters& params, const Species& species);

    void readReactionInput(const std::string& filename);

    Eigen::VectorXd calcConcentrations(const Eigen::VectorXd& species_Rho);

    double Chemistry::calcTBFactor(const Reaction& reaction, const Eigen::VectorXd& concentrations);

    double calcLawMassAction(
        const Reaction& reaction,
        const Eigen::VectorXd& species_Rho,
        double temp);
};

