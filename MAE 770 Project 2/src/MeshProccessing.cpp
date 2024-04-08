#include "MeshProccesing.hpp"

Mesh::Mesh(const std::string& filename) {
    processMesh(filename);
}

void Mesh::processMesh(const std::string& filename) {

    readMesh(filename);
    calcCellVolumes1D();
    calcCellCentroids();
    calcAreaVectors();
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
    cell_data.resize(imax + 1, std::vector<Cell>(jmax + 1));
    face_data.resize(imax + 1, std::vector<Face>(jmax + 1));

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
}

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

            cell_data[i][j].volume = AT1 + AT2;
        }
    }
}

void Mesh::calcDx() {

    for (int j = 1; j < jmax; ++j) {
        for (int i = 1; i < imax; ++i) {

            cell_data[i][j].dx = x_coords(i, j) - x_coords(i, j - 1);
        }
    }
}

void Mesh::calcCellVolumes1D() {

    calcDx();
    calcCellVolumes();

	for (int j = 1; j < jmax; ++j) {
		for (int i = 1; i < imax; ++i) {

            cell_data[i][j].volume /= cell_data[i][j].dx;
		}
	}
}

void Mesh::calcCellCentroids() {
    for (int i = 1; i < imax; ++i) {
        for (int j = 1; j < jmax; ++j) {
            double x_centroid = 0.25 * (x_coords(i, j) + x_coords(i - 1, j) + x_coords(i, j - 1) + x_coords(i - 1, j - 1));
            double y_centroid = 0.25 * (y_coords(i, j) + y_coords(i - 1, j) + y_coords(i, j - 1) + y_coords(i - 1, j - 1));
            cell_data[i][j].centroid(0) = x_centroid; 
            cell_data[i][j].centroid(1) = y_centroid; 
        }
    }

    for (int i = 1; i < imax; ++i) {
        cell_data[i][jmax].centroid = 2 * cell_data[i][jmax - 1].centroid - cell_data[i][jmax - 2].centroid;
        cell_data[i][0].centroid = 2 * cell_data[i][1].centroid - cell_data[i][2].centroid;
    }

    for (int j = 1; j < jmax; ++j) {
        cell_data[imax][j].centroid = 2 * cell_data[imax - 1][j].centroid - cell_data[imax - 2][j].centroid;
        cell_data[0][j].centroid = 2 * cell_data[1][j].centroid - cell_data[2][j].centroid;
    }
}

void Mesh::calcAreaVectors() {
    for (int i = 1; i < imax; ++i) {
        for (int j = 0; j < jmax; ++j) {
            face_data[i][j].i_surface_norms(0) = y_coords(i, j) - y_coords(i - 1, j);
            face_data[i][j].i_surface_norms(1) = -(x_coords(i, j) - x_coords(i - 1, j));
        }
    }

    for (int j = 1; j < jmax; ++j) {
        for (int i = 0; i < imax; ++i) {
            face_data[i][j].j_surface_norms(0) = y_coords(i, j - 1) - y_coords(i, j);
            face_data[i][j].j_surface_norms(1) = -(x_coords(i, j - 1) - x_coords(i, j));
        }
    }
}


double Mesh::getCellArea1D(int cell_idx) const {

    return cell_data[1][cell_idx].volume;
}

double Mesh::getCelldx1D(int cell_idx) const {

    return cell_data[1][cell_idx].dx;
}

double Mesh::getFaceArea1D(int cell_idx) const {

    return face_data[1][cell_idx].i_surface_norms(0);
}

Eigen::Vector2d Mesh::getCellCentroid(int cell_idx) const {

    return cell_data[1][cell_idx].centroid;
}

void Mesh::printMeshData(int flag) {
    if (flag == 1) {
        std::cout << "Cell Centroid Vectors:" << '\n';
        for (int i = 0; i < cell_data.size(); ++i) {
            for (int j = 0; j < cell_data[i].size(); ++j) {
                std::cout << "(" << cell_data[i][j].centroid(0) << ", " << cell_data[i][j].centroid(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 2) {
        std::cout << "Cell Volumes:" << '\n';
        for (int i = 0; i < cell_data.size(); ++i) {
            for (int j = 0; j < cell_data[i].size(); ++j) {
                std::cout << cell_data[i][j].volume << " ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 3) { 
        std::cout << "i Surface Normal Vectors:" << '\n';
        for (int i = 0; i < cell_data.size(); ++i) {
            for (int j = 0; j < cell_data[i].size(); ++j) {
                std::cout << "(" << face_data[i][j].i_surface_norms(0) << ", " << face_data[i][j].i_surface_norms(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
    else if (flag == 4) { 
        std::cout << "j Suface Normal Vectors:" << '\n';
        for (int i = 0; i < cell_data.size(); ++i) {
            for (int j = 0; j < cell_data[i].size(); ++j) {
                std::cout << "(" << face_data[i][j].j_surface_norms(0) << ", " << face_data[i][j].j_surface_norms(1) << ") ";
            }
            std::cout << std::endl;
        }
    }
}




