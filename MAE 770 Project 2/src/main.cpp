#include "ProjectIncludes.hpp"
#include "GlobalConstants.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"
#include "BuildCellJacobians.hpp"

int main() {

    Mesh mesh;
    Species species;

    std::string inputfilepath = "C:\\Users\\luis2\\Documents\\MAE 770\\Project 2\\inputfile.dat";
    SimParameters params(inputfilepath, mesh, species);

    CellStateVars state(params, mesh, species);

    CellJacobians LHS(params, mesh, species, state);

    for (int i = 0; i < mesh.ncells; ++i) {
        LHS.cell_vec[i].updateLHS();
    }

    std::cout << LHS.cell_vec[3].matrix << '\n';

    mesh.printMeshData(3);

	double temp_tr = 298.15;
	double temp_V = temp_tr;
	std::cout << "Noneq enthalpy is: " << species.species_vector[2].calcEnthalpyNonEQ(temp_tr, temp_V) << '\n';
	std::cout << "Eq enthalpy is: " << species.species_vector[2].calcEnthalpyEQ(temp_tr) << '\n';

	std::cout << "Eq internal energy is: " << species.species_vector[2].calcIntEnergyEQ(temp_tr) << '\n';
	std::cout << "Noneq internal energy is: " << species.species_vector[2].calcIntEnergyNonEQ(temp_tr, temp_V) << '\n';

    return 0;
}
