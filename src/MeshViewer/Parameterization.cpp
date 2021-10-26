#include <iostream>
#include "Parameterization.h"

Parameterization::Parameterization(Mesh* mesh)
{
	this->mesh = mesh;
}

void Parameterization::FindBoundary()
{
	if (mesh == nullptr)
		throw"mesh is null!";
	Boundary_idx.clear();

	OpenMesh::SmartHalfedgeHandle bhe;
	int start_idx = -1;
	for (auto& heh : mesh->halfedges())
	{
		if (heh.is_boundary())
		{
			bhe = heh;
			Boundary_idx.push_back(bhe.from().idx());
			start_idx = bhe.from().idx();

			break;
		}
	}
	if (start_idx == -1)
		throw"Boundary is empty!";
	while (bhe.to().idx() != start_idx)
	{
		bhe = bhe.next();
		Boundary_idx.push_back(bhe.from().idx());
	}
}


std::vector<OpenMesh::Vec2d> Parameterization::CircleBoundaryCoordinate()
{
	if (Boundary_idx.empty())
		FindBoundary();

	int v_num = Boundary_idx.size();
	std::vector<OpenMesh::Vec2d> boundary_coordinate(v_num);
	double pi = 3.1415926535;
	for (int n = 0; n < v_num; n++)
	{
		boundary_coordinate[n] = { cos(2 * pi / v_num * n),sin(2 * pi / v_num * n) };
	}
	return boundary_coordinate;
}
