#pragma once
#include "Parameterization.h"
#include <Eigen/Sparse>

typedef Eigen::Triplet<double> T;
class MVCParameterization :public Parameterization {
public:
	MVCParameterization(Mesh* mesh) : Parameterization(mesh) {
	}
	~MVCParameterization() {
		tan_weight.clear();
	}

	void Parameterize();

protected:
	void TanWeight();

	std::vector<double> tan_weight;
};