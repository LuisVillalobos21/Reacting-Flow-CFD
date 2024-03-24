#include "ProjectIncludes.hpp"
#include "GlobalConstants.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"
#include "BuildCellJacobians.hpp"
#include "RightHandSide.hpp"
#include "TimeIntegration.hpp"


int main() {

    Mesh mesh;

    Species species;

    std::string inputfilepath = "C:\\Users\\luis2\\Documents\\MAE 770\\Project 2\\inputfile.dat";

    SimParameters params(inputfilepath, mesh, species);

    CellStateVars state(params, mesh, species);

    TimeEvolveSolution evolve(params, mesh, species, state);

    evolve.solve();

    std::cout << "Final temps" << '\n';
    for (int i = 0; i < mesh.jmax + 1; ++i) {

        std::cout << std::fixed << std::setprecision(9) << state.cell_vec[i].var_vec(params.T_idx) << '\n';
    }

    //for (int i = 0; i < mesh.jmax + 1; ++i) {

    //    std::cout << std::fixed << std::setprecision(9) << state.cell_vec[i].var_vec(params.T_idx) << '\n';
    //}

    //for (int i = 0; i < mesh.jmax + 1; ++i) {

    //    std::cout << std::fixed << std::setprecision(9) << state.cell_vec[i].var_vec(params.T_idx) << '\n';
    //}

    return 0;
}
