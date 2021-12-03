#include "MeshParamDialog.h"
#include <QApplication>
#include <QDesktopWidget>

#include "../Common/CommonDefinitions.h"

MeshParamDialog::MeshParamDialog(QWidget* parent /* = 0 */)
	:QDialog(parent)
{
	initDialog();
}

MeshParamDialog::~MeshParamDialog()
{
}

void MeshParamDialog::initDialog()
{
	createWidget();
	createLayout();
}

void MeshParamDialog::createWidget()
{
	create_Basic_Operation_Information_Widget();
}

void MeshParamDialog::createLayout()
{
	tabWidget = new QTabWidget();
	tabWidget->addTab(view_BOI, "QP");

	QGridLayout *layout = new QGridLayout();
	layout->addWidget(tabWidget, 0, 0, 1, 1);
	setLayout(layout);
}

void MeshParamDialog::create_Basic_Operation_Information_Widget()
{
	print_info = new QPushButton("Print Mesh Information");
	MVC_parameterization = new QPushButton("MVC Parameterization");
	initial_mesh = new QPushButton("Set Initial Mesh");
	final_mesh = new QPushButton("Set Final Mesh");
	morphing2d = new QPushButton("Morphing 2D");
	solve_laplace = new QPushButton("Solve Laplace equation");
	loop_subdivision_button = new QPushButton("Loop subdivision");
	leftLabel_BOI = new QLabel("");

	QGridLayout* mainLayout = new QGridLayout(); int main_index = 0;
	mainLayout->addWidget(print_info, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(MVC_parameterization, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(initial_mesh, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(final_mesh, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(morphing2d, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(solve_laplace, main_index, 0, 1, 2); main_index += 1;
	mainLayout->addWidget(loop_subdivision_button, main_index, 0, 1, 2); main_index += 1;

	mainLayout->addWidget(leftLabel_BOI, main_index, 0, 1, 40);

	Basic_Operation_And_Information = new QWidget();
	Basic_Operation_And_Information->setLayout(mainLayout);

	view_BOI = new QScrollArea;
	view_BOI->setFocusPolicy(Qt::NoFocus);
	view_BOI->setFrameStyle(QFrame::NoFrame);
	view_BOI->setWidget(Basic_Operation_And_Information);
	view_BOI->setWidgetResizable(true);

	connect(print_info, SIGNAL(clicked()), SIGNAL(print_info_signal()));
	connect(MVC_parameterization, SIGNAL(clicked()), SIGNAL(do_mvc_parameterization()));
	connect(initial_mesh, SIGNAL(clicked()), SIGNAL(set_initial_mesh()));
	connect(final_mesh, SIGNAL(clicked()), SIGNAL(set_final_mesh()));
	connect(morphing2d, SIGNAL(clicked()), SIGNAL(do_morphing2d()));
	connect(solve_laplace, SIGNAL(clicked()), SIGNAL(solve_laplace1()));
	connect(loop_subdivision_button, SIGNAL(clicked()), SIGNAL(loop_subdivision()));

	
}
