#pragma once

#include <vector>
#include "MeshDefinition.h"

class Morphing2D {
public:
	Morphing2D(Mesh& init,Mesh& final);
	~Morphing2D();

	void Local();

private:

	Mesh mesh_init;
	Mesh mesh_final;
	double delta_t;
};
