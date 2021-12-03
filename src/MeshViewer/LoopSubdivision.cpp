#include "LoopSubdivision.h"

std::vector<OpenMesh::Vec3d> LoopSubdivision::AddEdgePoints() const {
	if (mesh == nullptr)
		throw"MESH IS NULL!";

	std::vector<OpenMesh::Vec3d> new_edge_points;
	new_edge_points.reserve(mesh->n_edges());
	for (auto& eh: mesh->edges())
	{
		if (eh.is_boundary())
		{
			auto p = 0.5 * (mesh->point(eh.v0()) + mesh->point(eh.v1()));
			new_edge_points.push_back(p);
		}
		else
		{
			auto p = 0.375 * (mesh->point(eh.v0()) + mesh->point(eh.v1()));
			p += 0.125 * mesh->point(eh.h0().next().to());
			p += 0.125 * mesh->point(eh.h1().next().to());
			new_edge_points.push_back(p);
		}
	}

	return new_edge_points;
}

void LoopSubdivision::UpdateVertices()
{
	for (auto& vh : mesh->vertices())
	{
		OpenMesh::Vec3d new_point;
		if (vh.is_boundary())
		{
			new_point = 0.75 * mesh->point(vh);
			for (auto& vvh : mesh->vv_range(vh))
			{
				if (vvh.is_boundary())
					new_point += 0.125 * mesh->point(vvh);
				else
					continue;
			}
		}
		else
		{
			new_point = OpenMesh::Vec3d(0);
			int n = 0;
			for (auto& hh : vh.outgoing_halfedges())
			{
				new_point += mesh->point(hh.to());
				n++;
			}
			double beta = Beta(n);
			new_point = beta * new_point + (1 - n * beta) * mesh->point(vh);
		}
		mesh->set_point(vh, new_point);
	}
}

inline double LoopSubdivision::Beta(int n) const
{
	return  (0.625 - pow((0.375 + 0.25 * cos(2 * Pi / n)), 2)) / n;
}
