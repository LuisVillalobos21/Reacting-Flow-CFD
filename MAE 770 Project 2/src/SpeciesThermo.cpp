#include "SpeciesThermo.hpp"

// SPECIESTHERMO STRUCT
SpeciesThermo::SpeciesThermo(const std::string& filename) {

    readSpeciesData(filename);

    h_formation = calcEnthalpyEQ(temp_ref);
    e_formation = calcIntEnergyEQ(temp_ref);
    R_s = R_U / molecular_weight;
}

void SpeciesThermo::readSpeciesData(const std::string& filename) {

    coeff_high.resize(7);
    coeff_low.resize(7);

    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: could not open " << filename << '\n';
        return;
    }

    std::string line;
    if (!std::getline(inputFile, line)) {
        std::cerr << "ERROR: could not read from " << filename << '\n';
        return;
    }

    std::istringstream iss(line); // read in the molecular weight, temp ranges, starting from end of line
    std::vector<std::string> tokens;
    std::string token;

    while (iss >> token) {
        tokens.push_back(token);
    }

    molecular_weight = std::stod(tokens[tokens.size() - 2]);
    highT = std::stod(tokens[tokens.size() - 3]);
    lowT = std::stod(tokens[tokens.size() - 4]);

    double temp; // reading in all thermo coefficients
    for (int i = 0; i < 5; ++i) {
        inputFile >> coeff_high(i);
    }

    inputFile >> temp;

    for (int i = 5; i < 7; ++i) {
        inputFile >> coeff_high(i);
    }

    for (int i = 0; i < 3; ++i) {
        inputFile >> coeff_low(i);
    }

    inputFile >> temp;

    for (int i = 3; i < 7; ++i) {
        inputFile >> coeff_low(i);
    }

    inputFile >> h_ref;

    inputFile.close();
}

double SpeciesThermo::calcPolyEnthalpy(const Eigen::MatrixXd& coeffs, double temp) const {
    double temp2 = temp * temp;   
    double temp3 = temp2 * temp; 
    double temp4 = temp3 * temp;  
    double temp5 = temp4 * temp; 

    return (R_U / molecular_weight) * (coeffs(0) * temp +
        0.5 * coeffs(1) * temp2 +
        (1.0 / 3.0) * coeffs(2) * temp3 +
        0.25 * coeffs(3) * temp4 +
        0.2 * coeffs(4) * temp5 +
        coeffs(5));
}

double SpeciesThermo::calcPolyEnthalpyTR(double temp) const {

    if (temp < lowT) {
        throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
        return 0.0;
    }

    //if (temp > highT) {
    //    throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
    //    return 0.0;
    //}

    if (temp < 1000)
        return (R_U / molecular_weight) * (coeff_low(0) * temp);
    else
        return (R_U / molecular_weight) * (coeff_high(0) * temp);
}

double SpeciesThermo::calcEnthalpyEQ(double temp)  const {

    if (temp < lowT) {
        throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
        return 0.0;
    }


    //if (temp > highT) {
    //    throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
    //    return 0.0;
    //}

    if (temp < 1000)
        return calcPolyEnthalpy(coeff_low, temp);
    else
        return calcPolyEntropy(coeff_high, temp);
}

double SpeciesThermo::calcEnthalpyTransRot(double temp)  const {

    return calcPolyEnthalpyTR(temp) - calcPolyEnthalpyTR(temp_ref) + h_formation;
}

double SpeciesThermo::calcEnthalpyVib(double temp)  const {

    double h_tr = calcPolyEnthalpyTR(temp) - calcPolyEnthalpyTR(temp_ref);
    double h = calcEnthalpyEQ(temp);
    return h - h_tr - h_formation; // refer to vulcan manual for formula
}

double SpeciesThermo::calcEnthalpyNonEQ(double temp_tr, double temp_V)  const {

    // see vulcan manual for formula

    return calcEnthalpyVib(temp_V) + calcEnthalpyTransRot(temp_tr);
}

double SpeciesThermo::calcIntEnergyEQ(double temp)  const {

    return calcEnthalpyEQ(temp) - (R_U / molecular_weight) * temp;
}

double SpeciesThermo::calcIntEnergyTransRot(double temp)  const {

    return calcEnthalpyTransRot(temp) - (R_U / molecular_weight) * temp;
}

double SpeciesThermo::calcIntEnergyNonEQ(double temp_tr, double temp_V)  const {

    return calcIntEnergyTransRot(temp_tr) + calcEnthalpyVib(temp_V);
}

double SpeciesThermo::calcCV_tr(double temp) const{

	if (temp < 1000)
		return (R_U / molecular_weight) * coeff_low(0) - (R_U / molecular_weight);
	else
		return (R_U / molecular_weight) * coeff_high(0) - (R_U / molecular_weight);
}

double SpeciesThermo::calcPolyCV_V(const Eigen::MatrixXd& coeffs, double temp) const {
    double temp2 = temp * temp;   // temp^2
    double temp3 = temp2 * temp;  // temp^3
    double temp4 = temp3 * temp;  // temp^4

    return (R_U / molecular_weight) * (coeffs(1) * temp +
        coeffs(2) * temp2 +
        coeffs(3) * temp3 +
        coeffs(4) * temp4);
}

double SpeciesThermo::calcCV_V(double temp) const {

    if (temp < 1000) 
        return calcPolyCV_V(coeff_low, temp);
    else 
        return calcPolyCV_V(coeff_high, temp);
}

double SpeciesThermo::calcPolyEntropy(const Eigen::MatrixXd& coeffs, double temp)  const {

    return R_U * (coeffs(0) * log(temp) + coeffs(1) * temp + 0.5 * coeffs(2) * pow(temp, 2.0) +
        (1.0 / 3.0) * coeffs(3) * pow(temp, 3.0) + 0.25 * coeffs(4) * pow(temp, 4.0) + coeffs(6)); // probably need to divide my M_w
}

double SpeciesThermo::calcEntropy(double temp)  const {

    if (temp < lowT) {
        throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
        return 0.0;
    }


    if (temp > highT) {
        throw std::runtime_error("ERROR: temperature of species exceeds curve fit bounds.");
        return 0.0;
    }

    if (temp < 1000)
        return calcPolyEntropy(coeff_low, temp);
    else
        return calcPolyEntropy(coeff_high, temp);
}

double SpeciesThermo::calcGibbs(double temp)  const {

    return calcEnthalpyEQ(temp) - temp * calcEntropy(temp);
}

void SpeciesThermo::printSpeciesThermoData() {

    std::cout << "Lower Temperature Bound: " << lowT << '\n';
    std::cout << "Upper Temperature Bound: " << highT << '\n';
    std::cout << "Molecular Weight: " << molecular_weight << '\n';
    std::cout << "Reference Enthalpy: " << h_ref << '\n';
    std::cout << "Enthalpy of Formation : " << h_formation << '\n';
    std::cout << "High Fit Coefficients:\n" << coeff_high << '\n';
    std::cout << "Low Fit Coefficients:\n" << coeff_low << '\n';
}

// SPECIES STRUCT

void Species::readSpeciesInput(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile.is_open()) {
        std::cerr << "ERROR: could not open " << filename << '\n';
        return;
    }

    std::string line;
    while (getline(inputFile, line)) {
        SpeciesThermo tempThermo(line);
        species_vector.push_back(tempThermo);
    }

    inputFile.close();
}

double Species::getCV_s(int species_idx, double temp) const {
    return species_vector[species_idx].calcCV_tr(temp);
}

double Species::getCV_V(int species_idx, double temp) const {
    return species_vector[species_idx].calcCV_V(temp);
}

double Species::getR_s(int species_idx) const {
    return species_vector[species_idx].R_s;
}

double Species::getIntEnergyTR(int species_idx, double temp_tr) const {

    return species_vector[species_idx].calcIntEnergyTransRot(temp_tr);
}

double Species::getIntEnergyV(int species_idx, double temp_V) const {
    return species_vector[species_idx].calcEnthalpyVib(temp_V);
}

double Species::getFormationEnergy(int species_idx) const {
    return species_vector[species_idx].e_formation;
}

