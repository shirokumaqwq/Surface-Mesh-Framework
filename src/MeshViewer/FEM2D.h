#pragma once

#include <vector>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>

#include "MeshDefinition.h"

#define	Pi 3.1415926535898

typedef Eigen::Triplet<double> Tri;

class FEM2D {
public:
	FEM2D(const Mesh& _mesh);
	~FEM2D() { area.clear(); };

	Eigen::SparseMatrix<double> ConstructMatrix() const;
	Eigen::VectorXd Solve(Eigen::SparseMatrix<double> M) const;

private:
	Mesh mesh;
	std::vector<double> area;
	std::vector<double> varea;
};