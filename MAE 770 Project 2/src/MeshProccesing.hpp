#pragma once
#include "ProjectIncludes.hpp"

struct Cell {
    double volume;
    double dx;
    Eigen::Vector2d centroid;
};

struct Face {
    Eigen::Vector2d i_surface_norms;
    Eigen::Vector2d j_surface_norms;
};

struct Mesh {
    int imax;
    int jmax;
    int ncells;

    Eigen::MatrixXd x_coords;
    Eigen::MatrixXd y_coords;
    Eigen::MatrixXd cell_vols;
    Eigen::MatrixXd x_centroids;
    Eigen::MatrixXd y_centroids;

    std::vector<std::vector<Cell>> cell_data;
    std::vector<std::vector<Face>> face_data;

    Mesh(const std::string& filename);

    void processMesh(const std::string& filename);

    void readMesh(const std::string& filename);

    void calcCellVolumes();

    void calcDx();

    void calcCellVolumes1D();

    void calcCellCentroids();

    void calcAreaVectors();

    void printMeshData(int flag);

    double getCellArea1D(int cell_idx) const;

    double getCelldx1D(int cell_idx) const;

    double getFaceArea1D(int cell_idx) const;

    Eigen::Vector2d getCellCentroid(int cell_idx) const;
};


