#pragma once
#include "ProjectIncludes.hpp"
#include "GlobalConstants.hpp"

struct SpeciesThermo {

    Eigen::VectorXd coeff_high;
    Eigen::VectorXd coeff_low;
    double A_low, B_low, C_low, D_low, E_low, F_low, G_low;
    double A_high, B_high, C_high, D_high, E_high, F_high, G_high;
    double molecular_weight;
    double h_ref;
    double lowT;
    double highT;
    double temp_ref = 298.15;
    double h_formation;
    double e_formation;
    double R_s;

    SpeciesThermo(const std::string& filename);

    void readSpeciesData(const std::string& filename);

    // Enthalpy Functions
    double calcPolyEnthalpy(const Eigen::MatrixXd& coeffs, double temp) const;

    double calcPolyEnthalpyTR(double temp) const;

    double calcPolyEntropy(const Eigen::MatrixXd& coeffs, double temp) const;

    double calcEnthalpyEQ(double temp) const;

    double calcEnthalpyTransRot(double temp) const;

    double calcEnthalpyVib(double temp) const;

    double calcEnthalpyNonEQ(double temp_tr, double temp_V) const;

    // Internal Energy Functions
    double calcIntEnergyEQ(double temp) const;

    double calcIntEnergyTransRot(double temp) const;

    double calcIntEnergyVib(double temp) const;

    double calcIntEnergyNonEQ(double temp_tr, double temp_V) const;

    double calcCV_tr(double temp) const;

    double calcPolyCV_V(const Eigen::MatrixXd& coeffs, double temp) const;

    double calcCV_V(double temp) const;

    // Entropy Functions 
    double calcEntropy(double temp) const;

    // Rest
    double calcGibbs(double temp) const;

    void printSpeciesThermoData();
};

struct Species {

    std::vector<SpeciesThermo> species_vector;

    void readSpeciesInput(const std::string& filename);

    double getCV_s(int species_idx, double temp) const;

    double getCV_V(int species_idx, double temp) const;

    double getR_s(int species_idx) const;

    double getMw(int species_idx) const;

    double getIntEnergyTR(int species_idx, double temp_tr) const;

    double getIntEnergyV(int species_idx, double temp_V) const;

    double getFormationEnergy(int species_idx) const;
};