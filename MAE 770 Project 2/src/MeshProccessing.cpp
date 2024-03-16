#include "ProjectIncludes.hpp"
#include "MeshProccesing.hpp"

void Mesh::calcCellVolumes() {
    double xa, ya, xb, yb, xc, yc, xd, yd;
    double da, db, dc, dap, dbp, dcp;
    double s1, AT1, s2, AT2;

    for (int j = 1; j < jmax; ++j) {
        for (int i = 1; i < imax; ++i) {

            xa = x_coords(i - 1, j);
            ya = y_coords(i - 1, j);
            xb = x_coords(i, j);
            yb = y_coords(i, j);
            xc = x_coords(i, j - 1);
            yc = y_coords(i, j - 1);
            xd = x_coords(i - 1, j - 1);
            yd = y_coords(i - 1, j - 1);

            da = sqrt(pow(xb - xa, 2) + pow(yb - ya, 2));
            db = sqrt(pow(xd - xb, 2) + pow(yd - yb, 2));
            dc = sqrt(pow(xa - xd, 2) + pow(ya - yd, 2));

            s1 = (da + db + dc) / 2;
            AT1 = sqrt(s1 * (s1 - da) * (s1 - db) * (s1 - dc));

            dap = sqrt(pow(xb - xc, 2) + pow(yb - yc, 2));
            dbp = sqrt(pow(xd - xb, 2) + pow(yd - yb, 2));
            dcp = sqrt(pow(xc - xd, 2) + pow(yc - yd, 2));

            s2 = (dap + dbp + dcp) / 2;
            AT2 = sqrt(s2 * (s2 - dap) * (s2 - dbp) * (s2 - dcp));

            Data[i][j].volume = AT1 + AT2;
        }
    }
}

void Mesh::calcCellCentroids() {
    for (int i = 1; i < imax; ++i) {
        for (int j = 1; j < jmax; ++j) {
            double x_centroid = 0.25 * (x_coords(i, j) + x_coords(i - 1, j) + x_coords(i, j - 1) + x_coords(i - 1, j - 1));
            double y_centroid = 0.25 * (y_coords(i, j) + y_coords(i - 1, j) + y_coords(i, j - 1) + y_coords(i - 1, j - 1));
            Data[i][j].centroid(0) = x_centroid; 
            Data[i][j].centroid(1) = y_centroid; 
        }
    }

    for (int i = 1; i < imax; ++i) {
        Data[i][jmax].centroid = 2 * Data[i][jmax - 1].centroid - Data[i][jmax - 2].centroid;
        Data[i][0].centroid = 2 * Data[i][1].centroid - Data[i][2].centroid;
    }

    for (int j = 1; j < jmax; ++j) {
        Data[imax][j].centroid = 2 * Data[imax - 1][j].centroid - Data[imax - 2][j].centroid;
        Data[0][j].centroid = 2 * Data[1][j].centroid - Data[2][j].centroid;
    }
}


void Mesh::calcAreaVectors() {
    for (int i = 1; i < imax; ++i) {
        for (int j = 0; j < jmax; ++j) {
            Data[i][j].i_surface_norms(0) = y_coords(i, j) - y_coords(i - 1, j);
            Data[i][j].i_surface_norms(1) = -(x_coords(i, j) - x_coords(i - 1, j));
        }
    }

    for (int j = 1; j < jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            Data[i][j].j_surface_norms(0) = y_coords(i, j - 1) - y_coords(i, j);
            Data[i][j].j_surface_norms(1) = -(x_coords(i, j - 1) - x_coords(i, j));
        }
    }
}

void Mesh::printMeshData(int flag) {
    if (flag == 1) {
        std::cout << "Cell Centroid Vectors:" << '\n';
        for (int i = 0; i < Data.size(); ++i) {
            for (int j = 0; j < Data[i].size(); ++j) {
                std::cout << "(" << Data[i][j].centroid(0) << ", " << Data[i][j].centroid(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 2) {
        std::cout << "Cell Volumes:" << '\n';
        for (int i = 0; i < Data.size(); ++i) {
            for (int j = 0; j < Data[i].size(); ++j) {
                std::cout << Data[i][j].volume << " ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 3) { 
        std::cout << "i Surface Normal Vectors:" << '\n';
        for (int i = 0; i < Data.size(); ++i) {
            for (int j = 0; j < Data[i].size(); ++j) {
                std::cout << "(" << Data[i][j].i_surface_norms(0) << ", " << Data[i][j].i_surface_norms(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 4) { 
        std::cout << "j Suface Normal Vectors:" << '\n';
        for (int i = 0; i < Data.size(); ++i) {
            for (int j = 0; j < Data[i].size(); ++j) {
                std::cout << "(" << Data[i][j].j_surface_norms(0) << ", " << Data[i][j].j_surface_norms(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
}

void Mesh::readMesh(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    file >> jmax >> imax;

    ncells = (imax - 1) * (jmax - 1);
    x_coords = Eigen::MatrixXd(imax, jmax);
    y_coords = Eigen::MatrixXd(imax, jmax);
    Data.resize(imax + 1, std::vector<Cell>(jmax + 1));

    double value;

    for (int i = 0; i < imax; ++i) {
        for (int j = 0; j < jmax; ++j) {
            file >> value;
            x_coords(i, j) = value;
        }
    }

    for (int i = 0; i < imax; ++i) {
        for (int j = 0; j < jmax; ++j) {
            file >> value;
            y_coords(i, j) = value;
        }
    }

    file.close();

    calcCellVolumes();
    calcCellCentroids();
    calcAreaVectors();
}



