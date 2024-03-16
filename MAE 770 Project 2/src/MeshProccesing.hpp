#pragma once
#include "ProjectIncludes.hpp"

struct Cell {
    double volume;
    Eigen::Vector2d centroid;
    Eigen::Vector2d i_surface_norms;
    Eigen::Vector2d j_surface_norms;
};

struct Mesh {
    int imax;
    int jmax;
    double ncells;
    Eigen::MatrixXd x_coords;
    Eigen::MatrixXd y_coords;
    Eigen::MatrixXd cell_vols;
    Eigen::MatrixXd x_centroids;
    Eigen::MatrixXd y_centroids;
    std::vector<std::vector<Cell>> Data;

    void readMesh(const std::string& filename);

    void calcCellVolumes();

    void calcCellCentroids();

    void calcAreaVectors();

    void printMeshData(int flag);
};


