#include <iostream>
#include <algorithm>
#include <cmath>

#include <qapplication.h>

#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <Eigen/Core>

#include "MeshViewerWidget.h"
#include "../Common/CommonDefinitions.h"

#include "MVCParameterization.h"
#include "Morphing2D.h"
#include "FEM2D.h"
#include "LoopSubdivision.h"

using namespace Qt;

MeshViewerWidget::MeshViewerWidget(QWidget* parent)
 : QGLViewerWidget(parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::MeshViewerWidget(QGLFormat& _fmt, QWidget* _parent)
: QGLViewerWidget(_fmt, _parent)
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();

	mesh_vector.clear(); mesh_vector_index = -1;

	draw_BBox_OK = false;
	draw_mesh_boundary_ok = false;
	first_init = true;
}

MeshViewerWidget::~MeshViewerWidget()
{
}

void MeshViewerWidget::updateMeshCenter()
{
	typedef Mesh::Point Point;
	Mesh::VertexIter vIt = mesh.vertices_begin();
	Mesh::VertexIter vEnd = mesh.vertices_end();
	bbMin = bbMax = OpenMesh::vector_cast<OpenMesh::Vec3d>(mesh.point(vIt));

	size_t count = 0;
	for(; vIt != vEnd; ++vIt, ++count)
	{
		bbMin.minimize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
		bbMax.maximize( OpenMesh::vector_cast<OpenMesh::Vec3d>( mesh.point(vIt) ) );
	}

	Mesh::EdgeIter e_it = mesh.edges_begin();
	Mesh::EdgeIter e_end = mesh.edges_end();
	double aveLen = 0.0; double maxLen = 0.0; double minLen = mesh.calc_edge_length(e_it);
	double e_len = 0.0;
	for(; e_it != e_end; ++e_it)
	{
		double e_len = mesh.calc_edge_length(e_it);
		if( e_len > maxLen )
		{
			maxLen = e_len;
		}
		else if(e_len < minLen )
		{
			minLen = e_len;
		}
		aveLen += e_len;
	}

	// set center and radius and box's radius.
	/*OpenMesh::Vec3d c = (bbMin + bbMax)*0.5;
	for (vIt = mesh.vertices_begin(); vIt != vEnd; ++vIt, ++count)
	{
		OpenMesh::Vec3d p = mesh.point(vIt) - c;
		mesh.set_point(vIt, p);
	}*/

	if( first_init )
	{
		set_scene_pos( (bbMin+bbMax)*0.5, (bbMin-bbMax).norm()*0.5 );
		first_init = false;
	}
	else
	{
		set_scene_pos((bbMin + bbMax)*0.5, (bbMin - bbMax).norm()*0.5);
	}
	
	
	printf("BoundingBox:\nX : [ %f , %f ]\n",bbMin[0],bbMax[0]);
	printf("Y : [ %f , %f ]\n",bbMin[1],bbMax[1]);
	printf("Z : [ %f , %f ]\n",bbMin[2],bbMax[2]);
	printf("Diag length of BBox : %f\n", (bbMax-bbMin).norm() );
	printf("Edge Length : Max : %f; Min : %f; AVG : %f\n",maxLen,minLen, aveLen/mesh.n_edges() );
}

void MeshViewerWidget::updateMeshNormals()
{
	mesh.update_face_normals();
	mesh.update_vertex_normals();
}

bool MeshViewerWidget::openMesh(const char* filename)
{
	clearAllMesh();
	

	bool read_OK = OpenMesh::IO::read_mesh( mesh, filename );

	printf("%s\n", filename);
	if ( read_OK )
	{
		initMesh();
		/*Mesh::EdgeHandle eh0 = mesh.edge_handle(960);
		flip_openmesh(eh0, mesh);
		Mesh::EdgeHandle eh1 = mesh.edge_handle(85);
		flip_openmesh(eh1, mesh);*/
		/*FILE* f_de = fopen("b.de", "w");
		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			if (mesh.is_boundary(v_it))
			{
				OpenMesh::Vec3d& p = mesh.point(v_it);
				fprintf(f_de, "%d %15.14f %15.14f\n", v_it.handle().idx(), p[0], p[1]);
			}
		}
		fclose(f_de);*/
		// loading done
		mesh_vector.push_back(mesh); mesh_vector_index++;
		return true;
	}
	return false;
}

void MeshViewerWidget::initMesh()
{
	mesh.request_vertex_status();
	mesh.request_edge_status();
	mesh.request_face_status();

	mesh.request_face_normals();
	mesh.request_vertex_normals();
	printBasicMeshInfo();
	updateMesh();
}

void MeshViewerWidget::printBasicMeshInfo()
{
	if (mesh.n_vertices() == 0)
		printf("No Mesh\n");

	checkMeshMode();

	QString meshType;
	switch(meshMode())
	{
	case TRIANGLE:
		printf("Triangle Mesh.\n");
		break;
	case QUAD:
		printf("Quadrilateral Mesh.\n");
		break;
	default:
		printf("General Mesh.\n");
		break;
	}

	printf("Information of the input mesh:\nVertex : %d;\nFace : %d;\nEdge : %d, HalfEdge : %d\n",
		mesh.n_vertices(),mesh.n_faces(),mesh.n_edges(),mesh.n_halfedges());;

	//
	/*for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d& p = mesh.point(v_it);
		p[2] = 0.0;
	}*/

	//save tet mesh
	/*int nf = mesh.n_faces(); int nv = mesh.n_vertices(); int ne = mesh.n_edges();
	std::vector<OpenMesh::Vec3d> face_n(nf);
	std::vector<double> edge_len(ne);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<OpenMesh::Vec3d> one_face;
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it); fv_it; ++fv_it)
		{
			one_face.push_back(mesh.point(fv_it));
		}
		OpenMesh::Vec3d n = OpenMesh::cross(one_face[1] - one_face[0], one_face[2] - one_face[0]);
		face_n[f_it->idx()] = n.normalize();
	}
	for (Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		edge_len[e_it->idx()] = mesh.calc_edge_length(e_it);
	}
	std::vector<OpenMesh::Vec3d> v_v(nv);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		double ave_edge_len = 0.0; OpenMesh::Vec3d vn(0, 0, 0); int edge_count = 0;
		for (Mesh::VertexOHalfedgeIter voh_it = mesh.voh_iter(v_it); voh_it; ++voh_it)
		{
			ave_edge_len += edge_len[mesh.edge_handle(voh_it).idx()]; edge_count += 1;
			Mesh::FaceHandle fh = mesh.face_handle(voh_it);
			if (fh != Mesh::InvalidFaceHandle)
			{
				vn += face_n[fh.idx()];
			}
		}
		vn.normalize();
		v_v[v_it->idx()] = mesh.point(v_it) - vn*ave_edge_len / edge_count;
	}
	std::vector<std::vector<int>> all_tets(nf*6);
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<int> one_face; 
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(f_it); fv_it; ++fv_it)
		{
			one_face.push_back( fv_it->idx() );
		}
		int f_id = f_it->idx();
		std::vector<int> one_tet(4);
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[0] + nv;
		all_tets[6 * f_id + 0] = one_tet;
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[1] + nv;
		all_tets[6 * f_id + 1] = one_tet;
		one_tet[0] = one_face[0]; one_tet[1] = one_face[2]; one_tet[2] = one_face[1]; one_tet[3] = one_face[2] + nv;
		all_tets[6 * f_id + 2] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[0];
		all_tets[6 * f_id + 3] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[1];
		all_tets[6 * f_id + 4] = one_tet;
		one_tet[0] = one_face[0] + nv; one_tet[1] = one_face[1] + nv; one_tet[2] = one_face[2] + nv; one_tet[3] = one_face[2];
		all_tets[6 * f_id + 5] = one_tet;
	}

	FILE* f_tet = fopen("tet.tet", "w");
	fprintf(f_tet, "Vertices %d", nv+nv);
	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		OpenMesh::Vec3d p = mesh.point(v_it);
		fprintf(f_tet, "\n%.19f %.19f %.19f", p[0], p[1], p[2]);
	}
	for (int i = 0; i < nv;++i)
	{
		fprintf(f_tet, "\n%.19f %.19f %.19f", v_v[i][0], v_v[i][1], v_v[i][2]);
	}
	fprintf(f_tet, "\nTetrahedra %d", nf*6);
	for (int i = 0; i < all_tets.size(); ++i)
	{
		fprintf(f_tet, "\n%d %d %d %d", all_tets[i][0] + 1, all_tets[i][1] + 1, all_tets[i][2] + 1, all_tets[i][3] + 1);
	}
	fclose(f_tet);*/
}


void MeshViewerWidget::solveLaplaceEquation()
{
	//����FEM�����
	FEM2D femtool(mesh);
	Eigen::VectorXd u = femtool.Solve(femtool.ConstructMatrix());

	for (auto& vh : mesh.vertices())
		mesh.set_point(vh, OpenMesh::Vec3d(mesh.point(vh)[0], mesh.point(vh)[1], u[vh.idx()]));

	updateMesh();
	update();
}

void MeshViewerWidget::domvcparameterization()
{
	std::cout << "MVC Parameterization" << std::endl;
	MVCParameterization para(&mesh);
	para.Parameterize();
	updateMesh();
	update();
}

void MeshViewerWidget::setInitialMesh()
{
	mesh_init.clear();
	mesh_init.assign(mesh);
	std::cout <<"n:"<< mesh_init.n_vertices() << std::endl;
	std::cout << "Set initial mesh successfully" << std::endl;
}

void MeshViewerWidget::setFinalMesh()
{
	mesh_final.clear();
	mesh_final.assign(mesh);

	std::cout << "Set final mesh successfully" << std::endl;
}

void MeshViewerWidget::doMorphing2D()
{
	std::cout << "Morphing" << std::endl;

	mesh.clear();
	mesh = mesh_init;
	initMesh();
	Morphing2D morphing(mesh, mesh_final);
	morphing.Local();
	update();
}

void MeshViewerWidget::doLoopSubdivision()
{
	std::cout << "Loop Subdivision..." << std::endl;
	LoopSubdivision subdivision(&mesh);
	std::vector<OpenMesh::Vec3d> EP = subdivision.AddEdgePoints();  //compute new egde points.
	subdivision.UpdateVertices();  //update vertices

	mesh.request_face_status();
	mesh.request_edge_status();
	mesh.request_vertex_status();
	std::vector<OpenMesh::SmartVertexHandle> EPH; 
	EPH.reserve(mesh.n_edges());
	for (auto& piter : EP)  //add new edge points to mesh.
	{
		EPH.push_back(mesh.add_vertex(piter));
	}

	for (auto& fh : mesh.faces())  //add new faces.
	{
		auto h1 = fh.halfedge();
		auto p1 = h1.from(); auto p2 = h1.to(); auto p3 = h1.next().to();
		OpenMesh::SmartVertexHandle ep1 = EPH[h1.edge().idx()];
		OpenMesh::SmartVertexHandle ep2 = EPH[h1.next().edge().idx()];
		OpenMesh::SmartVertexHandle ep3 = EPH[h1.prev().edge().idx()];
		//std::cout << mesh.point(p1) << std::endl;
		//std::cout << mesh.point(p2) << std::endl;
		//std::cout << mesh.point(p3) << std::endl;
		//std::cout << mesh.point(ep1) << std::endl;
		//std::cout << mesh.point(ep2) << std::endl;
		//std::cout << mesh.point(ep3) << std::endl;
		std::vector< OpenMesh::SmartVertexHandle> f1 = { p1,ep1,ep3 };
		std::vector< OpenMesh::SmartVertexHandle> f2 = { ep1,ep2,ep3 };
		std::vector< OpenMesh::SmartVertexHandle> f3 = { ep1,p2,ep2 };
		std::vector< OpenMesh::SmartVertexHandle> f4 = { ep2,p3,ep3 };
		mesh.delete_face(mesh.face_handle(fh.idx()),false);  //delete old faces.
		mesh.add_face(f1);
		mesh.add_face(f2);
		mesh.add_face(f3);		
		mesh.add_face(f4);
	}
	mesh.garbage_collection();
	updateMesh();
	update();

	std::cout << "v num:" << mesh.n_vertices() << "\n f num:" << mesh.n_faces() << std::endl;

}

bool MeshViewerWidget::saveMesh(const char* filename)
{
	return OpenMesh::IO::write_mesh(mesh, filename);
}
bool MeshViewerWidget::saveScreen(const char* filePath)
{
	QImage image = grabFrameBuffer();
	image.save(filePath);
	return true;
}

void MeshViewerWidget::checkMeshMode()
{
	Mesh::FaceIter fIt = mesh.faces_begin();
	Mesh::FaceIter fEnd = mesh.faces_end();
	Mesh::FaceEdgeIter fe_it;
	int count = 1;
	int meshType[3] = {0};
	for(fIt; fIt != fEnd; ++fIt)
	{
		fe_it = mesh.fe_iter(fIt);
		while(--fe_it)
		{
			++count;
		}
		if(count == 4)
		{
			meshType[1]++;
		}
		else if(count == 3)
		{
			meshType[0]++;
		}
		else
		{
			meshType[2]++;
		}
		count = 1;
	}
	int faceNum = mesh.n_faces();
	if(meshType[0] == faceNum)//triangle
	{
		setMeshMode(TRIANGLE);
	}
	else if(meshType[1] == faceNum)//no 
	{
		setMeshMode(QUAD);
	}
	else
	{
		setMeshMode(N_MESH_MODES);
	}
}

void MeshViewerWidget::updateIndices()
{
	//Mesh::ConstFaceIter f_it(mesh.faces_sbegin()), 
	//	f_end(mesh.faces_end());
	//Mesh::ConstFaceVertexIter  fv_it;

	////Indices.clear(); VIndices.clear();
	//int PolygonClass = 0 ;
	//switch(meshMode())
	//{
	//case TRIANGLE:
	//	PolygonClass = 3;
	//	break;
	//case QUAD:
	//	PolygonClass = 4;
	//	break;
	//default:
	//	return;
	//}
	//
	//Indices.reserve(mesh.n_faces()*PolygonClass);
	//for (; f_it!=f_end; ++f_it)
	//{
	//	fv_it = mesh.cfv_iter(f_it); 
	//	for( fv_it; fv_it; ++fv_it)
	//	{
	//		Indices.push_back((fv_it).handle().idx());
	//	}
	//}

	//VIndices.resize(mesh.n_vertices());
	//for(int i=0;i<mesh.n_vertices();++i)
	//{
	//	VIndices[i] = i;
	//}
}

void MeshViewerWidget::draw_scene(int drawmode)
{
	QFont Text_Font("Courier", 12);
	glViewport ( 0,0, width(),height());
	glMatrixMode( GL_PROJECTION );
	glLoadMatrixd( &ProjectionMatrix[0] );
	glMatrixMode( GL_MODELVIEW );
	glLoadMatrixd( &ModelViewMatrix[0] );

	/*double r = (bbMin-bbMax).norm()*0.5;
	OpenMesh::Vec3d c = (bbMin+ bbMax)*0.5;
	OpenMesh::Vec3d x(1.0*r,0,0); OpenMesh::Vec3d y(0,1.0*r,0); OpenMesh::Vec3d z(0,0,1.0*r);
	OpenMesh::Vec3d temp;
	glDisable(GL_LIGHTING);
	glLineWidth(2.0);
	glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + x);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0] + 0.01*r, temp[1], temp[2] , "X", Text_Font );

	glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + y);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0], temp[1] + 0.01*r, temp[2] , "Y", Text_Font );

	glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3dv( c.data() );
	temp = (c + z);
	glVertex3dv( temp.data() );
	glEnd();
	renderText(temp[0], temp[1], temp[2] + 0.01*r , "Z", Text_Font );*/

	if(draw_BBox_OK)
	{
		OpenMesh::Vec3d temp0 = bbMin;
		OpenMesh::Vec3d temp1;
		glLineWidth(2.0);
		glColor3f(1.0, 1.0, 0.0);
		glBegin(GL_LINES);
		temp1 = bbMin; temp1[0] = bbMax[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[1] = bbMax[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMin; temp1[2] = bbMax[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[0] = bbMax[0];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[1] = bbMax[1];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMin; temp0[2] = bbMax[2];
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );

		temp0 = bbMax;
		temp1 = bbMax; temp1[0] = bbMin[0];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[1] = bbMin[1];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		temp1 = bbMax; temp1[2] = bbMin[2];
		glVertex3dv( temp0.data() );
		glVertex3dv( temp1.data() );
		glEnd();
	}

	if(draw_mesh_boundary_ok)
	{
		glLineWidth(2.0);
		glColor3f(1.0, 0.5, 0.0);
		glBegin(GL_LINES);
		for(Mesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
		{
			if( mesh.is_boundary(e_it ) )
			{
				Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle(e_it, 0);
				glVertex3dv( mesh.point(mesh.to_vertex_handle(heh0)).data() );
				glVertex3dv( mesh.point(mesh.from_vertex_handle(heh0)).data() );
			}
		}
		glEnd();

		/*FILE* f_bde = fopen("bde.de", "w");
		for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
		{
			if (mesh.is_boundary(v_it))
			{
				OpenMesh::Vec3d& np = mesh.point(v_it);
				fprintf(f_bde, "%d %20.19f %20.19f 0.0\n", v_it.handle().idx(), np[0], np[1]);
			}
		}
		fclose(f_bde);*/
	}

	draw_scene_mesh(drawmode);
}

void MeshViewerWidget::draw_scene_mesh(int drawmode)
{
	if(mesh.n_vertices() == 0) { return; }

	switch (drawmode)
	{
	case WIRE_FRAME:
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case HIDDEN_LINES:
		glDisable(GL_LIGHTING);
		draw_mesh_hidden_lines();
		break;
	case SOLID_FLAT:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case FLAT_POINTS:
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.5f,2.0f);
		glEnable(GL_LIGHTING);
		glShadeModel(GL_FLAT);
		draw_mesh_solidflat();
		glDisable(GL_POLYGON_OFFSET_FILL);
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		draw_mesh_wireframe();
		//draw_meshpointset();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		break;
	case SOLID_SMOOTH:
		glEnable(GL_LIGHTING);
		glShadeModel(GL_SMOOTH);
		draw_mesh_solidsmooth();
		//draw_meshpointset();
		glDisable(GL_LIGHTING);
		break;
	case POINT_SET:
		glDisable(GL_LIGHTING);
		draw_mesh_pointset();
		break;
	case CURVATURE:
		
		break;
	default:
		break;
	}
}


void MeshViewerWidget::draw_mesh_wireframe()
{
	glLineWidth(1);
	//glColor3f(0.753, 0.753, 0.753);
	//glColor3f(0.0, 0.0, 0.0);
	glColor3f(0.0, 0.0, 0.25);
	//if(meshMode() != TRIANGLE && meshMode() != QUAD)
	{
		Mesh::ConstFaceIter fIt(mesh.faces_begin()),
			fEnd(mesh.faces_end());
		Mesh::ConstFaceVertexIter fvIt;
		for (; fIt != fEnd; ++fIt)
		{
			//GL::glNormal(dualMesh.normal(f_it));
			fvIt = mesh.cfv_iter(fIt); 
			glBegin(GL_POLYGON);
			for( fvIt; fvIt; ++fvIt )
			{
				glVertex3dv( mesh.point(fvIt).data() );
			}
			glEnd();
		}
	}

	/*OpenMesh::Vec3d pos = mesh.point(mesh.vertex_handle(0));
	GLdouble  winX, winY, winZ;
	GLint     viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluProject(pos[0],pos[1],pos[2],&ModelViewMatrix[0][0],&ProjectionMatrix[0][0],viewport,&winX,&winY,&winZ);
	int x = (long)winX;
	int y = viewport[3]-(long)winY;
	QString str = "fxum";
	render_text(x,y,str);*/
	/*else
	{
		glEnableClientState(GL_VERTEX_ARRAY);
		glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
		switch(meshMode())
		{
		case TRIANGLE:
			glDrawElements(GL_TRIANGLES,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		case QUAD:
			glDrawElements(GL_QUADS,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
			break;
		default:
			break;
		}

		glDisableClientState(GL_VERTEX_ARRAY);
	}*/

}


void MeshViewerWidget::draw_mesh_hidden_lines() const
{
	Mesh::ConstFaceIter f_it(mesh.faces_begin());
	Mesh::ConstFaceIter f_end(mesh.faces_end());
	Mesh::ConstFaceVertexIter fv_it;
	glLineWidth(2.0);
	glColor3f(0.0, 1.0, 1.0);

	glDrawBuffer(GL_NONE);
	glDepthRange(0.01, 1.0);
	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}

		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDrawBuffer(GL_BACK);
	glDepthRange(0.0, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glDepthFunc(GL_LEQUAL);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (f_it = mesh.faces_begin(); f_it!=f_end; ++f_it)
		{
			fv_it = mesh.cfv_iter(f_it); 
			glBegin(GL_POLYGON);
			for(fv_it;fv_it;++fv_it)
			{
				glVertex3dv(&mesh.point(fv_it)[0]);
			}
			glEnd();
		}
		break;
	}

	glDepthFunc(GL_LESS);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void MeshViewerWidget::draw_mesh_solidflat() const
{
	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	GLfloat mat_a[] = { 0.7f, 0.7f, 0.7f, 1.0f };
	GLfloat mat_d[] = { 0.88f, 0.84f, 0.76f, 1.0f };
	GLfloat mat_s[] = { 0.4f, 0.4f, 0.4f, 1.0f };
	GLfloat shine[] = { 120.0f };

	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_a);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_d);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_s);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);

	switch(meshMode())
	{
	case TRIANGLE:
		glBegin(GL_TRIANGLES);
		for (fIt; fIt != fEnd; ++fIt) 
		{
			glNormal3dv(mesh.normal(fIt).data());
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(mesh.point(fvIt).data());
			}
		}
		glEnd();
		break;
	case QUAD:
		glBegin(GL_QUADS);
		for (;fIt != fEnd; ++fIt) 
		{
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
		}
		glEnd();
		break;
	default:
		for (; fIt != fEnd; ++fIt) 
		{
			glBegin(GL_POLYGON);
			glNormal3dv(&mesh.normal(fIt)[0]);
			fvIt = mesh.cfv_iter(fIt.handle());
			for(fvIt; fvIt; ++fvIt)
			{
				glVertex3dv(&mesh.point(fvIt)[0]);
			}
			glEnd();
		}
		break;
	}
}

void MeshViewerWidget::draw_mesh_solidsmooth() const 
{
	bool drawOK = false;
	glLoadName(mesh.n_vertices());

	//glEnableClientState(GL_VERTEX_ARRAY);
	//glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	//glEnableClientState(GL_NORMAL_ARRAY);
	//glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());
	//switch(meshMode())
	//{
	//case TRIANGLE:
	//	glDrawElements(GL_TRIANGLES,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
	//	drawOK = true;
	//	break;
	//case QUAD:
	//	glDrawElements(GL_QUADS,Indices.size(),GL_UNSIGNED_INT,&Indices[0]);
	//	drawOK = true;
	//	break;
	//default:
	//	break;
	//}
	//glDisableClientState(GL_VERTEX_ARRAY);
	//glDisableClientState(GL_NORMAL_ARRAY);

	//if(drawOK) return;

	Mesh::ConstFaceIter fIt(mesh.faces_begin()),
		fEnd(mesh.faces_end());
	Mesh::ConstFaceVertexIter fvIt;

	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_DOUBLE, 0, mesh.vertex_normals());

	for (; fIt != fEnd; ++fIt)
	{
		glBegin(GL_POLYGON);
		fvIt = mesh.cfv_iter(fIt.handle());
		for(fvIt; fvIt ; --fvIt)
		{
			glArrayElement(fvIt.handle().idx());
		}
		glEnd();
	}

	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
}

void MeshViewerWidget::draw_mesh_pointset() const 
{
	glPolygonMode(GL_FRONT_AND_BACK, GL_POINTS);

	glColor3f(0.0, 1.0, 1.0);
	glPointSize(10);
	Mesh::VertexIter v_it = mesh.vertices_begin();
	glBegin(GL_POINTS);
	for(v_it;v_it != mesh.vertices_end();++v_it)
	{
		glVertex3dv(mesh.point(v_it).data());
	}
	glEnd();
	/*glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_DOUBLE, 0, mesh.points());
	glDrawElements(GL_POINTS,mesh.n_vertices(),GL_UNSIGNED_INT,&VIndices[0]);
	glDisableClientState(GL_VERTEX_ARRAY);*/

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}
