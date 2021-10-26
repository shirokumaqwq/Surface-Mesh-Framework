#include "MVCParameterization.h"

#include <Eigen/Sparse>
#include <iostream>
void MVCParameterization::Parameterize()
{
	if(Boundary_idx.empty())
		FindBoundary();
	//std::cout << "bnum:" << Boundary_idx.size() << std::endl;
	TanWeight();
	std::vector<OpenMesh::Vec2d> bound_coor = CircleBoundaryCoordinate();
	for (auto iter : bound_coor)
		std::cout << iter << std::endl;
	std::vector<T> triplelist;
	triplelist.reserve(mesh->n_halfedges());
	int vnum = mesh->n_vertices();

	for (auto& vh:mesh->vertices())
	{
		int vidx = vh.idx();
		triplelist.push_back(T(vidx, vidx, 1));

		for (auto& hvh : vh.outgoing_halfedges())
		{
			auto v2 = hvh.to();
			if (!v2.is_boundary())
			{
				triplelist.push_back(T(vidx, v2.idx(), -tan_weight[hvh.idx()]));
			}
			
		}
	}

	Eigen::SparseMatrix<double> MVC(vnum, vnum);
	MVC.setFromTriplets(triplelist.begin(), triplelist.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
	solver.compute(MVC);

	Eigen::MatrixX2d RH = Eigen::MatrixX2d::Zero(vnum, 2);

	for (int i = 0; i < Boundary_idx.size(); i++)
	{
		int bv_idx = Boundary_idx[i];
		RH(bv_idx, 0) = bound_coor[i][0];
		RH(bv_idx, 1) = bound_coor[i][1];

		for (auto& vhh : mesh->vih_range(mesh->vertex_handle(bv_idx)))
		{
			int v2_idx = vhh.from().idx();
			RH(v2_idx, 0) += tan_weight[vhh.idx()] * bound_coor[i][0];
			RH(v2_idx, 1) += tan_weight[vhh.idx()] * bound_coor[i][1];
		}
	}
	Eigen::VectorXd solu = solver.solve(RH.col(0));
	Eigen::VectorXd solv = solver.solve(RH.col(1));
	
	for (auto& vh : mesh->vertices())
	{
		mesh->set_point(vh, OpenMesh::Vec3d(solu[vh.idx()], solv[vh.idx()], 0));
		if (vh.is_boundary())
			std::cout << vh.idx() << ":" << solu[vh.idx()] << "," << solv[vh.idx()] << std::endl;
	}
}

void MVCParameterization::TanWeight()
{
	tan_weight.clear();
	int hnum = mesh->n_halfedges();
	tan_weight.reserve(hnum);

	for (auto& hh:mesh->halfedges())
	{
		if (hh.from().is_boundary())
		{
			tan_weight.push_back(0);
			continue;
		}
		auto v1 = mesh->calc_edge_vector(hh);
		auto v2 = mesh->calc_edge_vector(hh.prev().opp());
		//°ë½Ç¹«Ê½£¬tan(a/2)=sin(a)/(1+cos(a))=|v1 x v2|/(|v1||v2| + v1 * v2)
		double tan1 = v1.cross(v2).norm() / (v1.norm() * v2.norm() + v1.dot(v2));

		v2 = mesh->calc_edge_vector(hh.opp().next());
		double tan2 = v1.cross(v2).norm() / (v1.norm() * v2.norm() + v1.dot(v2));

		tan_weight.push_back((tan1 + tan2) / v1.norm());
	}

	for (auto& vh : mesh->vertices())
	{
		if (!vh.is_boundary())
		{
			double qwq = 0;
			for (auto& vhh : vh.outgoing_halfedges())
			{
				qwq += tan_weight[vhh.idx()];
			}
			for (auto& vhh : vh.outgoing_halfedges())
			{
				tan_weight[vhh.idx()] /= qwq;
			}
		}
	}
}

