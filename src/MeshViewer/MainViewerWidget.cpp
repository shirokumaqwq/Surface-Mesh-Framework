#include "MainViewerWidget.h"

MainViewerWidget::MainViewerWidget(QWidget* _parent/* =0 */)
{
	initViewerWindow();
	LoadMeshSuccess = false;
}
MainViewerWidget::~MainViewerWidget()
{
};

void MainViewerWidget::initViewerWindow()
{
	createParamIDialog();
	createViewerDialog();

	//this->setOrientation( Qt::Horizontal );

	//this->addWidget(debugDialog);
	//OpenGL mesh viewer
	/*this->addWidget(MeshParam);
	this->addWidget(MeshViewer);

	//set the splitter line color
	this->setStyleSheet("QSplitter::handle { background-color: green }");
	QSplitterHandle* splitterHandle = this->handle(1);
	splitterHandle->setDisabled(true);*/

	QHBoxLayout* main_layout = new QHBoxLayout();
	main_layout->addWidget(MeshParam, 1);
	main_layout->addWidget(MeshViewer, 5);
	this->setLayout(main_layout);

	connect(MeshViewer,SIGNAL(setMouseMode_signal(int)),SIGNAL(setMouseMode_signal_main(int)));
	connect(MeshViewer,SIGNAL(setDrawMode_signal(int)),SIGNAL(setDrawMode_signal_main(int)));
	connect(MeshViewer,SIGNAL(set_edit_undo_enable_viewer_signal(bool)),SIGNAL(set_edit_undo_enable_signal(bool)));
	connect(MeshViewer,SIGNAL(set_edit_redo_enable_viewer_signal(bool)),SIGNAL(set_edit_redo_enable_signal(bool)));

	connect(MeshParam, SIGNAL(print_info_signal()), SLOT(print_info()));
	connect(MeshParam, SIGNAL(do_mvc_parameterization()), SLOT(mvc_parameterization()));
	connect(MeshParam, SIGNAL(set_initial_mesh()), SLOT(initial_mesh()));
	connect(MeshParam, SIGNAL(set_final_mesh()), SLOT(final_mesh()));
	connect(MeshParam, SIGNAL(do_morphing2d()), SLOT(morphing_2d()));
	connect(MeshParam, SIGNAL(solve_laplace1()), SLOT(solve_Laplace_equation()));
	connect(MeshParam, SIGNAL(loop_subdivision()), SLOT(do_loop_subdivision()));

}

void MainViewerWidget::createParamIDialog()
{
	MeshParam = new MeshParamDialog();
}

void MainViewerWidget::createViewerDialog()
{
	QGLFormat glFormat;
	glFormat.setSampleBuffers(true);
	glFormat.setSamples(16);

	MeshViewer = new InteractiveViewerWidget(glFormat, NULL);
	MeshViewer->setAcceptDrops(true);
	connect(MeshViewer,SIGNAL(loadMeshOK(bool,QString)), this, SLOT(LoadMeshFromInner(bool,QString)) );
}

void MainViewerWidget::open_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->openMesh(fname.toLocal8Bit())) 
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
	else
	{
		LoadMeshSuccess = true;
		MeshViewer->setDrawMode(InteractiveViewerWidget::FLAT_POINTS);
		MeshViewer->setMouseMode(InteractiveViewerWidget::TRANS);
		if(LoadMeshSuccess)
		{
			SetMeshForALL();
		}
		emit(haveLoadMesh(fname));
	}
}

void MainViewerWidget::save_mesh_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveMesh(fname.toLocal8Bit()))
	{
		QString msg = "Cannot read mesh from file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::save_screen_gui(QString fname)
{
	if (fname.isEmpty() || !MeshViewer->saveScreen(fname.toLocal8Bit()))
	{
		QString msg = "Cannot save image to file:\n '";
		msg += fname;
		msg += "'";
		QMessageBox::critical(NULL, windowTitle(), msg);
	}
}

void MainViewerWidget::print_info()
{
	MeshViewer->printBasicMeshInfo();
}

void MainViewerWidget::mvc_parameterization()
{
	MeshViewer->domvcparameterization();
}

void MainViewerWidget::initial_mesh()
{
	MeshViewer->setInitialMesh();
}

void MainViewerWidget::final_mesh()
{
	MeshViewer->setFinalMesh();
}

void MainViewerWidget::morphing_2d()
{
	MeshViewer->doMorphing2D();
}

void MainViewerWidget::solve_Laplace_equation()
{
	MeshViewer->solveLaplaceEquation();
}

void MainViewerWidget::do_loop_subdivision()
{
	MeshViewer->doLoopSubdivision();
}