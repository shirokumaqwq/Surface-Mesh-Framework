#include "Morphing2D.h"
#include <omp.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/dense>


#include <iostream>

Morphing2D::Morphing2D(Mesh& init, Mesh & final)
{
	mesh_init = init;
	mesh_final = final;
}

Morphing2D::~Morphing2D()
{
	mesh_init.clear();
	mesh_final.clear();
}

void Morphing2D::Local()
{
	if (mesh_init.n_faces() != mesh_final.n_faces())
	{
		std::cout << "初始曲面与终曲面无对应" << std::endl;
	}
	int nface = mesh_init.n_faces();
	for (int n = 0; n < nface; n++)
	{
		auto fh1 = mesh_init.face_handle(n);
		auto fh2 = mesh_final.face_handle(n);

		Eigen::Matrix2d L;
		std::vector<OpenMesh::Vec3d> tempe1,tempe2;
		for (auto& hh1 : mesh_init.fh_range(fh1))
		{
			tempe1.push_back(mesh_init.calc_edge_vector(hh1));
			//std::cout << tempe1.back() << " ";
		}
		//std::cout << std::endl;
		for (auto& hh2 : mesh_init.fh_range(fh2))
		{
			tempe2.push_back(mesh_final.calc_edge_vector(hh2));
			//std::cout << tempe2.back() << " ";

		}		//std::cout << std::endl;

		Eigen::Matrix2d uv, xy;
		uv << tempe2[0][0], -tempe2[2][0], tempe2[0][1], -tempe2[2][1];
		xy << tempe1[0][0], -tempe1[2][0], tempe1[0][1], -tempe1[2][1];
		std::cout << "cross:" << tempe2[0].cross(-tempe2[2])[2] << "  " << tempe1[0].cross(-tempe1[2])[2] << std::endl;
		//std::cout << uv << '\n' << xy << std::endl;
		L = uv * xy.inverse();
		std::cout << "J1=" << L << std::endl;

		Eigen::JacobiSVD<Eigen::Matrix2d> svd(L, Eigen::ComputeFullU | Eigen::ComputeFullV);
		Eigen::Matrix2d U = svd.matrixU();
		Eigen::Matrix2d V = svd.matrixV();
		Eigen::Matrix2d sigma;
		sigma(0, 1) = sigma(1, 0) = 0;
		sigma(0, 0) = svd.singularValues()[0];
		sigma(1, 1) = svd.singularValues()[1];
		std::cout <<"rotation:\n"<< U << std::endl;
		Eigen::Matrix2d R = U * V.transpose();
		Eigen::Matrix2d S = V * sigma* V.transpose();

	}
}