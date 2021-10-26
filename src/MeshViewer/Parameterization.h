#pragma once

#include <vector>
#include <set>
#include <Eigen/Core>
#include "MeshDefinition.h"

class Parameterization {
public:
	Parameterization(Mesh* mesh);
	~Parameterization() {
		mesh = nullptr;
		Boundary_idx.clear();
	};
	void FindBoundary();
	std::vector<OpenMesh::Vec2d> CircleBoundaryCoordinate();

protected:
	Mesh* mesh;

	std::vector<int> Boundary_idx;
};