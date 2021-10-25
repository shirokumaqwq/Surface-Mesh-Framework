#pragma once

#include <vector>
#include <unordered_set>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "MeshDefinition.h"

typedef Eigen::Triplet<double> T;

class ARAPSurfaceModeling {
public:
	ARAPSurfaceModeling(Mesh *input_mesh);
	~ARAPSurfaceModeling();


	void SetMovingVertex(int point_idx, OpenMesh::Vec3d point_target_pos);
	void SetFixedVertex(std::vector<int> selectedVertex);

	void DoARAP(int n);
	void Init(std::vector<int> selectedVertex);

protected:
	std::vector<Eigen::Matrix3d> Local();
	Eigen::MatrixX3d Global(std::vector<Eigen::Matrix3d> rotation);
	
	Mesh* mesh;
	std::vector<Eigen::Vector3d> delta_h;  //�洢�ƶ�֮���λ��
	std::vector<double> weight;

	OpenMesh::Vec3d target_pos;
	int moving_point = 0;
	std::unordered_set<int> fixed_vertex;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

private:
	void InitDeltah();
	void InitUniformWeight();
	void InitCotWeight();
	void PreComputeMatrix();
};