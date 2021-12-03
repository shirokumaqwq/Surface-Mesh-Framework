#ifndef MESHPROCESSING_MESHPARAMDIALOG_H
#define MESHPROCESSING_MESHPARAMDIALOG_H

#include <QDialog>
#include <QtGui>
#include <QtWidgets>

class MeshParamDialog : public QDialog
{
	Q_OBJECT
public:
	MeshParamDialog(QWidget* parent=0);
	~MeshParamDialog();

	QSize sizeHint()
	{
		QRect rect = QApplication::desktop()->screenGeometry();
		return QSize( int( rect.width()*0.15), rect.height() );
	}

private:
	QTabWidget* tabWidget;

signals:
	void print_info_signal();
	void do_mvc_parameterization();
	void set_initial_mesh();
	void set_final_mesh();
	void do_morphing2d();
	void solve_laplace1();
	void loop_subdivision();
private:
	QWidget* Basic_Operation_And_Information;
	QScrollArea *view_BOI;

	QLabel* leftLabel_BOI;
	QPushButton* print_info;
	QPushButton* MVC_parameterization;
	QPushButton* initial_mesh;
	QPushButton* final_mesh;
	QPushButton* morphing2d;
	QPushButton* solve_laplace;
	QPushButton* loop_subdivision_button;

private:
	void create_Basic_Operation_Information_Widget();

private:
	void initDialog();
	void createWidget();
	void createLayout();

};

#endif