#pragma once

#include <vector>

#include "MeshDefinition.h"

#define	Pi 3.1415926535898

class LoopSubdivision {
public:
	LoopSubdivision(Mesh* p_mesh) :mesh(p_mesh){}
	~LoopSubdivision() {
		mesh = nullptr;
	}
	
	std::vector<OpenMesh::Vec3d> AddEdgePoints() const;
	void UpdateVertices();
	
private:
	inline double Beta(int n) const;
protected:
	Mesh* mesh;
};