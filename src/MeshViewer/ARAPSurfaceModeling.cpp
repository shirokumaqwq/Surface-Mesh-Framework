#include "ARAPSurfaceModeling.h"

#include <Eigen/SVD>
#include <Eigen/dense>

#include <omp.h>
#include <iostream>
using namespace std;

ARAPSurfaceModeling::ARAPSurfaceModeling(Mesh* input_mesh)
{
	mesh = input_mesh;
}

ARAPSurfaceModeling::~ARAPSurfaceModeling()
{
	mesh = nullptr;
	delta_h.clear();
	fixed_vertex.clear();
}

void ARAPSurfaceModeling::SetMovingVertex(int point_idx, OpenMesh::Vec3d point_target_pos)
{
	this->moving_point = point_idx;
	this->target_pos = point_target_pos;
	//cout << moving_point << ":" << target_pos << endl;
}


void ARAPSurfaceModeling::SetFixedVertex(std::vector<int> selectedVertex)
{
	for (auto iter : selectedVertex)
		fixed_vertex.insert(iter);
	
	std::cout << "Fixed point id:";
	for (auto iter : fixed_vertex)
		std:cout << iter << ",";
	cout << endl;
}

void ARAPSurfaceModeling::InitDeltah()
{
	delta_h.reserve(mesh->n_halfedges());
	for (auto& heh : mesh->halfedges())
	{
		/*if (heh.to().idx() == moving_point)
		{
			auto p = mesh->point(heh.from());
			delta_h.push_back({ target_pos[0] - p[0],target_pos[1] - p[1],target_pos[2] - p[2] });
		}
		else if(heh.from().idx() == moving_point)
		{
			auto p = mesh->point(heh.from());
			delta_h.push_back({ -target_pos[0] + p[0],-target_pos[1] + p[1],-target_pos[2] + p[2] });
		}
		else
		{*/
			auto v = mesh->calc_edge_vector(heh);
			delta_h.push_back({ v[0],v[1],v[2] });
		//}
	}
}

void ARAPSurfaceModeling::InitUniformWeight()
{
	weight.clear();
	weight.resize(mesh->n_halfedges());
	for (int i = 0; i < weight.size(); i++)
		weight[i] = 1;
}

void ARAPSurfaceModeling::InitCotWeight()
{
	weight.clear();
	weight.reserve(mesh->n_halfedges());;
	for (auto& heh : mesh->halfedges())
	{
		if (heh.is_boundary())
		{
			weight.push_back(0);
			continue;
		}
		auto he1 = mesh->calc_edge_vector(heh.prev());
		auto he2 = mesh->calc_edge_vector(heh.next().opp());
		weight.push_back(he1.dot(he2) / he1.cross(he2).norm());
	}
}

void ARAPSurfaceModeling::PreComputeMatrix()
{
	std::vector<T> triplist;
	triplist.reserve(mesh->n_halfedges());

	for (auto& vh : mesh->vertices())
	{
		if (fixed_vertex.count(vh.idx()) > 0)
		{
			triplist.push_back(T(vh.idx(), vh.idx(), 1));
			continue;
		}
		double qwq_weight = 0;
		for (auto& hvh : vh.outgoing_halfedges())
		{
			double temp_weight = weight[hvh.idx()] + weight[hvh.opp().idx()];
			qwq_weight += temp_weight;
			triplist.push_back(T(vh.idx(), hvh.to().idx(), -temp_weight));
		}
		//cout << vh.idx() << ":" << qwq_weight << endl;
		triplist.push_back(T(vh.idx(), vh.idx(), qwq_weight));
	}
	Eigen::SparseMatrix<double> S;
	S.resize(mesh->n_vertices(), mesh->n_vertices());
	S.setFromTriplets(triplist.begin(),triplist.end());

	solver.compute(S);
}

void ARAPSurfaceModeling::Init(std::vector<int> selectedVertex)
{
	SetFixedVertex(selectedVertex);
	InitDeltah();
	//InitCotWeight();
	InitUniformWeight();
	PreComputeMatrix();
}


std::vector<Eigen::Matrix3d>  ARAPSurfaceModeling::Local()
{

	std::vector<Eigen::Matrix3d> rotation_m;
	rotation_m.resize(mesh->n_vertices());

#pragma omp parallel for
	for(int n = 0;n < mesh->n_vertices(); n++)
	{
		auto vh = mesh->vertex_handle(n);
		Eigen::Matrix3d S = Eigen::Matrix3d::Zero();
		for (auto& heh : mesh->voh_range(vh))
		{
			int hidx = heh.idx();
			auto deltax = mesh->calc_edge_vector(mesh->halfedge_handle(hidx));
			S += 0.5*(weight[hidx]+weight[heh.opp().idx()]) * Eigen::Vector3d({ deltax[0],deltax[1],deltax[2] }) * delta_h[hidx].transpose();
		}
		Eigen::JacobiSVD<Eigen::Matrix3d> svd(S, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix3d U = svd.matrixU();
		Eigen::Matrix3d V = svd.matrixV();
		S = V * U.transpose();
		if (S.determinant() < 0)
		{
			U(0, 2) *= -1;
			U(1, 2) *= -1;
			U(2, 2) *= -1;
			S = V * U.transpose();
		}
		//std::cout << S << std::endl;
		rotation_m[n] = S;
	}

	return rotation_m;
}

Eigen::MatrixX3d ARAPSurfaceModeling::Global(std::vector<Eigen::Matrix3d> rotation)
{
	//Eigen::VectorXd X = Eigen::VectorXd::Zero(mesh->n_vertices());
	int vnum = mesh->n_vertices();
	Eigen::MatrixX3d RH = Eigen::MatrixX3d::Zero(vnum, 3);
	Eigen::MatrixX3d sol = Eigen::MatrixX3d::Zero(vnum, 3);

#pragma omp parallel for
	for (int n = 0; n < vnum; n++)
	{
		auto vh = mesh->vertex_handle(n);
		Eigen::Vector3d qwq = { 0,0,0 };
		for (auto& heh : mesh->voh_range(vh))
		{
			auto vh2 = heh.to();
			auto hv = mesh->calc_edge_vector(heh);
			qwq += 0.5 * (weight[heh.idx()] + weight[heh.opp().idx()]) * (rotation[vh.idx()] + rotation[vh2.idx()]) * (Eigen::Vector3d(-hv[0], -hv[1], -hv[2]));
		}
		for (int i = 0; i < 3; i++)
			RH(n, i) = qwq[i];
	}

	//set fixed vertices' position
	for (auto iter : fixed_vertex)
	{
		OpenMesh::Vec3d pv = mesh->point(mesh->vertex_handle(iter));
		for (int i = 0; i < 3; i++)
			RH(iter, i) = pv[i];
	}

	for (int i = 0; i < 3; i++)
		RH(moving_point, i) = target_pos[i];

	for (int i = 0; i < 3; i++)
		sol.col(i) = solver.solve(RH.col(i));

#pragma omp parallel for
	for (int n = 0; n < mesh->n_halfedges(); n++)
	{
		auto hh = mesh->halfedge_handle(n);
		int v1 = mesh->from_vertex_handle(hh).idx();
		int v2 = mesh->to_vertex_handle(hh).idx();
		delta_h[n] = (sol.row(v2) - sol.row(v1));
	}
	return sol;
}

void ARAPSurfaceModeling::DoARAP(int n)
{
	Eigen::MatrixX3d sol;
	for (int i = 0; i < n-1; i++)
	{
			Global(Local());
	}
	sol = Global(Local());
#pragma omp parallel for
	for (int n = 0; n < mesh->n_vertices(); n++)
		mesh->set_point(mesh->vertex_handle(n), MeshTraits::Point(sol(n, 0), sol(n, 1), sol(n, 2)));
}