#include "ProjectIncludes.hpp"
#include "GlobalConstants.hpp"
#include "MeshProccesing.hpp"
#include "SpeciesThermo.hpp"
#include "CellStateVariables.hpp"
#include "BuildCellJacobians.hpp"
#include "RightHandSide.hpp"
#include "TimeIntegration.hpp"


struct PostProcess {

    const SimParameters& params;
    const Mesh& mesh;
    const CellStateVars& state;

    PostProcess(const SimParameters& params, const Mesh& mesh, CellStateVars& state)
    : params(params), mesh(mesh), state(state) {
    }

    void writeDataToFile(int flag, const std::string& filename) {
        std::ofstream outFile(filename);
        if (!outFile.is_open()) {
            std::cerr << "Failed to open file for writing: " << filename << '\n';
            return;
        }

        if (flag == -1) {
            for (int cell_idx = 0; cell_idx < mesh.jmax + 1; ++cell_idx) {
                outFile << state.getPressure(cell_idx) << '\n';
            }
            std::cout << "Data written to " << filename << '\n';
            return;
        }

        if (flag == -2) {
            for (int cell_idx = 0; cell_idx < mesh.jmax + 1; ++cell_idx) {
                outFile << state.getRho(cell_idx) << '\n';
            }
            std::cout << "Data written to " << filename << '\n';
            return;
        }

        if (flag == -3) {
            for (int cell_idx = 0; cell_idx < mesh.jmax + 1; ++cell_idx) {
                for (int species_idx = 0; species_idx < params.nspecies; ++species_idx) {

                    outFile << state.cell_vec[cell_idx].var_vec(species_idx);

                    if (species_idx < params.nspecies - 1) {
                        outFile << ", ";
                    }
                    else {
                        outFile << '\n'; 
                    }
                }
            }
            std::cout << "Data written to " << filename << '\n';
            return;
        }

        outFile << std::fixed << std::setprecision(9);
        for (int cell_idx = 0; cell_idx < mesh.jmax + 1; ++cell_idx) {
            outFile << state.cell_vec[cell_idx].var_vec(flag) << '\n';
        }

        outFile.close();
        std::cout << "Data written to " << filename << '\n';
    }
};

int main() {

    Mesh mesh;

    Species species;

    std::string inputfilepath = "C:\\Users\\luis2\\Documents\\MAE 770\\Project 2\\simple_input.dat";

    SimParameters params(inputfilepath, mesh, species);

    CellStateVars state(params, mesh, species);

    TimeEvolveSolution evolve(params, mesh, species, state);

    evolve.solve();

    PostProcess PP(params, mesh, state);

    PP.writeDataToFile(params.vel_idx, "soln_vel.dat");
    PP.writeDataToFile(params.T_idx, "soln_temp.dat");
    PP.writeDataToFile(params.Tv_idx, "soln_tempV.dat");
    PP.writeDataToFile(-1, "soln_pressure.dat");
    PP.writeDataToFile(-2, "soln_density.dat");
    PP.writeDataToFile(-3, "soln_species_densities.dat");

    return 0;
}
