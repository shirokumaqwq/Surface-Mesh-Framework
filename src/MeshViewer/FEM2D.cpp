#include "FEM2D.h"

inline double f(OpenMesh::Vec2d p)
{
	double x = p[0]; double y = p[1];
	//return 2*Pi * Pi * sin(Pi * p[0]) * sin(Pi * p[1]);
	return -10 * (-2 * cos(1 - y) * cos(y) * sin(1 - x) * sin(x) - 2 * cos(1 - x) * cos(x) * sin(1 - y) * sin(y) - 4 * sin(1 - x) * sin(x) * sin(1 - y) * sin(y));
	//return 2 * sin(Pi * p[0]) + Pi * Pi * (1 - p[1]) * p[1] * sin(Pi * p[0]);
	//return p[0] * p[1];
}

inline std::vector<double> phi1(const OpenMesh::Vec2d& p1, const OpenMesh::Vec2d& p2, const OpenMesh::Vec2d& p3)
{
	std::vector<double> coff(3);

	double temp = -p2[0] * p1[1] + p3[0] * p1[1] + p1[0] * p2[1] - p3[0] * p2[1] - p1[0] * p3[1] + p2[0] * p3[1];
	coff[0] = (p2[1] - p3[1]) / temp;
	coff[1] = -(p2[0] - p3[0]) / temp;
	coff[2] = -(p3[0] * p2[1] - p2[0] * p3[1]) / temp;
	return coff;
}

inline std::vector<double> phi2(const OpenMesh::Vec2d& p1, const OpenMesh::Vec2d& p2, const OpenMesh::Vec2d& p3)
{
	return phi1(p2, p3, p1);
}

inline std::vector<double> phi3(const OpenMesh::Vec2d& p1, const OpenMesh::Vec2d& p2, const OpenMesh::Vec2d& p3)
{
	return phi1(p3, p1, p2);
}

//不用这个，瞎写的
double Gauss2(double f(OpenMesh::Vec2d), OpenMesh::Vec2d p1, OpenMesh::Vec2d p2, OpenMesh::Vec2d p3)
{
	//std::cout << "f=" << f(OpenMesh::Vec2d(0, 0)) << std::endl;
	//三角形中没有钝角三角形，否则p3应当为钝角顶点
	OpenMesh::Vec2d v1 = (p2 - p1).normalize();  //以一条边为轴
	double left_d = (p3 - p1).dot(v1);  //左侧直角三角形底边长
	double right_d = (p2 - p3).dot(v1);  //右侧直角三角形底边长
	OpenMesh::Vec2d v2 = OpenMesh::Vec2d(-v1[1], v1[0]);  //高所在轴
	double height = (p3 - p1).dot(v2);

	double sum = 0;
	std::vector<OpenMesh::Vec2d> int_p;
	std::vector<double> coff;
	//Gauss2积分在[0,1]积分节点为(3-sqrt(3))/3, sqrt(3)/3
	const double c1 = (3.0 - sqrt(3)) / 3; const double c2 = sqrt(3) / 3;

	//left triangle:
	OpenMesh::Vec2d Lor = p1;
	int_p.push_back(Lor + c1 * left_d * v1 + c1 * c2 * height * v2);
	coff.push_back(0.5 * left_d * 0.5 * c1 * height);
	int_p.push_back(Lor + c1 * left_d * v1 + c1 * c1 * height * v2);
	coff.push_back(0.5 * left_d * 0.5 * c1 * height);
	int_p.push_back(Lor + c2 * left_d * v1 + c2 * c2 * height * v2);
	coff.push_back(0.5 * left_d * 0.5 * c2 * height);
	int_p.push_back(Lor + c2 * left_d * v1 + c2 * c1 * height * v2);
	coff.push_back(0.5 * left_d * 0.5 * c2 * height);

	//right triangle:
	OpenMesh::Vec2d Ror = p1 + left_d * v1;
	int_p.push_back(Ror + c1 * right_d * v1 + c1 * c2 * height * v2);
	coff.push_back(0.5 * right_d * 0.5 * height * c2);
	int_p.push_back(Ror + c1 * right_d * v1 + c2 * c2 * height * v2);
	coff.push_back(0.5 * right_d * 0.5 * height * c2);
	int_p.push_back(Ror + c2 * right_d * v1 + c1 * c1 * height * v2);
	coff.push_back(0.5 * right_d * 0.5 * height * c1);
	int_p.push_back(Ror + c2 * right_d * v1 + c2 * c1 * height * v2);
	coff.push_back(0.5 * right_d * 0.5 * height * c1);
	for (size_t n = 0; n < coff.size(); n++)
		sum += coff[n] * f(int_p[n]);

	return sum;
}

inline double Simpson(double f(OpenMesh::Vec2d), double area, const OpenMesh::Vec2d& p1, const OpenMesh::Vec2d& p2, const OpenMesh::Vec2d& p3)
{
	return area / 3.0f * (f((p1 + p2) / 2) + f((p2 + p3) / 2) + f((p3 + p1) / 2));
}

inline double phi(const OpenMesh::Vec2d& x, const std::vector<double>& coff)
{
	return coff[0] * x[0] + coff[1] * x[1] + coff[2];
}

inline double Simpson_phi(double f(OpenMesh::Vec2d), double area, const OpenMesh::Vec2d& p1, const OpenMesh::Vec2d& p2, const OpenMesh::Vec2d& p3)
{
	auto coff = phi1(p1, p2, p3);
	auto e1 = (p1 + p2) / 2;
	auto e2 = (p2 + p3) / 2;
	auto e3 = (p3 + p1) / 2;

	return area / 3.0f * (f(e1) * phi(e1, coff) + f(e2) * phi(e2, coff) + f(e3) * phi(e3, coff));
}

FEM2D::FEM2D(const Mesh& _mesh) :mesh(_mesh) {
	varea = std::vector<double>(mesh.n_vertices());
	//初始化三角形面积
	double m_diam = 0;
	double avg_diam = 0;
	for (auto& hh : mesh.halfedges())
	{
		m_diam = m_diam > mesh.calc_edge_length(hh) ? m_diam : mesh.calc_edge_length(hh);
		avg_diam += mesh.calc_edge_length(hh);
	}
	for (auto& fh : mesh.faces())
	{
		auto hh1 = fh.halfedge();
		auto hh2 = hh1.next();

		double a = 0.5 * mesh.calc_edge_vector(hh1).cross(mesh.calc_edge_vector(hh2)).norm();
		area.push_back(a);
		varea[hh1.from().idx()] += a;
		varea[hh1.to().idx()] += a;
		varea[hh2.to().idx()] += a;
	}
	std::cout << "Init FEM2d solver successfully!\n Diam(mesh) = " << m_diam << "avg diam:" << avg_diam / mesh.n_halfedges() << std::endl;


};

Eigen::SparseMatrix<double> FEM2D::ConstructMatrix() const{
	if (mesh.vertices_empty())
		throw"mesh is empty!";

	std::vector<Tri> triplet;
	int IntV_num = 0;  //内部点个数
	for (auto& vh : mesh.vertices())
	{
		if (!vh.is_boundary())
		{
			double selfgrad = 0;
			for (auto& hh : vh.outgoing_halfedges())
			{
				OpenMesh::Vec2d vv1(mesh.point(vh)[0], mesh.point(vh)[1]);
				auto v2 = mesh.point(hh.to()); auto v3 = mesh.point(hh.next().to()); auto v4 = mesh.point(hh.opp().next().to());
				OpenMesh::Vec2d vv2(v2[0], v2[1]);
				OpenMesh::Vec2d vv3(v3[0], v3[1]);
				OpenMesh::Vec2d vv4(v4[0], v4[1]);

				//上方三角形的梯度计算积分
				auto coff1 = phi1(vv1, vv2, vv3);
				auto coff2 = phi2(vv1, vv2, vv3);
				double grad1 = coff1[0] * coff2[0] + coff1[1] * coff2[1];
				//std::cout << "faceself:" << hh.face().idx() << (coff1[0] * coff1[0] + coff1[1] * coff1[1]) * area[hh.face().idx()] << std::endl;
				selfgrad += (coff1[0] * coff1[0] + coff1[1] * coff1[1]) * area[hh.face().idx()];
				if (!hh.to().is_boundary())
				{
					//下方三角形的梯度计算积分
					coff1 = phi1(vv1, vv2, vv4);
					coff2 = phi2(vv1, vv2, vv4);
					double grad2 = coff1[0] * coff2[0] + coff1[1] * coff2[1];

					double qwq = area[hh.face().idx()] * grad1 + area[hh.opp().face().idx()] * grad2;
					triplet.push_back(Tri(vh.idx(), hh.to().idx(), qwq));

					std::cout << vh.idx() << ", " << hh.to().idx() << ": " << qwq << std::endl;
				}
			}
			triplet.push_back(Tri(vh.idx(), vh.idx(), selfgrad));
			std::cout << vh.idx() << ", " << vh.idx() << ": " << selfgrad << std::endl;
		}
		else
		{
			triplet.push_back(Tri(vh.idx(), vh.idx(), 1));
			//std::cout << vh.idx() << ", " << vh.idx() << ": " << 1 << std::endl;
		}
	}
	int vnum = mesh.n_vertices();
	Eigen::SparseMatrix<double> M(vnum, vnum);
	M.setFromTriplets(triplet.begin(), triplet.end());

	return M;
}

Eigen::VectorXd FEM2D::Solve(Eigen::SparseMatrix<double> M) const
{
	Eigen::VectorXd f_phi(mesh.n_vertices());
	for (auto& vh : mesh.vertices())
	{
		if (!vh.is_boundary())
		{
			double qwq = 0;
			for (auto& hh : vh.outgoing_halfedges())
			{
				OpenMesh::Vec2d vv1(mesh.point(vh)[0], mesh.point(vh)[1]);
				auto v2 = mesh.point(hh.to()); auto v3 = mesh.point(hh.next().to());
				OpenMesh::Vec2d vv2(v2[0], v2[1]);
				OpenMesh::Vec2d vv3(v3[0], v3[1]);

				qwq += Simpson_phi(f, area[hh.face().idx()], vv1, vv2, vv3);
			}
			f_phi(vh.idx()) = qwq;
		}
		else
		{
			f_phi(vh.idx()) = 0;
		}
	}
	std::cout << "f:\n" << f_phi << std::endl;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
	solver.compute(M);
	Eigen::VectorXd u = solver.solve(f_phi);
	std::cout << "solve:\n" << u << std::endl;
	double delta_0 = 0;
	double delta_1 = 0;
	/*for (int n = 0; n < mesh.n_vertices(); n++)
	{
		double x = mesh.point(mesh.vertex_handle(n))[0];
		double y = mesh.point(mesh.vertex_handle(n))[1];

		//double delta = abs(sin(Pi * x) * sin(Pi * y) - u[n]);
		double delta = abs(y * (1 - y) * sin(Pi * x) - u[n]);
		//std::cout<< sin(Pi*x)*sin(Pi*y) <<" - "<<u[n] << " d:" << delta << std::endl;
		delta_0 = delta > delta_0 ? delta : delta_0;
		delta_1 += delta * varea[n];
	}*/

	for (auto& fh : mesh.faces())
	{
		auto p = mesh.calc_centroid(fh);
		//double x = p[0]; double y = p[1];

		auto h = fh.halfedge();
		auto& v1 = mesh.point(h.from()); auto& v2 = mesh.point(h.to()); auto& v3 = mesh.point(h.next().to());
		OpenMesh::Vec2d p1(v1[0], v1[1]); OpenMesh::Vec2d p2(v2[0], v2[1]); OpenMesh::Vec2d p3(v3[0], v3[1]);

		double qwq = 0; std::vector<double> c;
		c = phi1(p1, p2, p3); qwq += (c[0] * p[0] + c[1] * p[1] + c[2]) * u[h.from().idx()];
		c = phi2(p1, p2, p3); qwq += (c[0] * p[0] + c[1] * p[1] + c[2]) * u[h.to().idx()];
		c = phi3(p1, p2, p3); qwq += (c[0] * p[0] + c[1] * p[1] + c[2]) * u[h.next().to().idx()];

		///double delta = abs(qwq - sin(Pi * p[0]) * sin(Pi * p[1]));
		double delta = abs(qwq - p[1] * (1 - p[1]) * sin(Pi * p[0]));
		delta_0 = delta > delta_0 ? delta : delta_0;
		delta_1 += delta * area[fh.idx()];
	}
	std::cout << "& " << delta_0 << " &  & " << delta_1 << "&" << std::endl;
	return u;
}
