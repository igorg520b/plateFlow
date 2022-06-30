#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

#include <QSizePolicy>
#include <QPushButton>
#include <QSplitter>
#include <QLabel>
#include <QVBoxLayout>
#include <QTreeWidget>
#include <QProgressBar>
#include <QMenu>
#include <QList>
#include <QDebug>
#include <QComboBox>
#include <QMetaEnum>
#include <QDir>
#include <QString>
#include <QCheckBox>
#include <QFile>
#include <QTextStream>
#include <QIODevice>

#include <QtCharts>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>

#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>

#include "objectpropertybrowser.h"
#include "backgroundworker.h"
#include "meshrepresentation.h"
#include "model.h"


QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT
private:
    Ui::MainWindow *ui;

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

    void showEvent( QShowEvent* event ) override;
    void closeEvent( QCloseEvent* event ) override;
    void paintEvent(QPaintEvent *event) override;

private Q_SLOTS:
    void background_worker_paused();    // updates the GUI when background worker pauses
    void updateGUI();   // when simulation is started/stopped or when a step is advanced

    void quit_triggered();

    void sliderValueChanged(int val);
    void comboboxIndexChanged_visualizations(int index);

    void single_step_triggered();
    void start_triggered(bool checked);
    void gotoStep0_triggered();
    void trim_triggered();

    void cameraReset_triggered();

    void showPlots_toggled(bool arg1);
    void importMesh_triggered();
    void save_triggered();
    void open_triggered();

    void screenshot_triggered();

private:
    icy::Model model;
    icy::MeshRepresentation meshRepresentation;
    BackgroundWorker *worker;

    QString settingsFileName;       // includes current dir
    QLabel *statusLabel;                    // statusbar
    QLabel *labelElapsedTime;
    QLabel *statusLabelStepFactor;             // attempt # at advancing a time step
    QLabel *labelStepCount;
    QLabel *labelNodeCount;
    QLabel *labelElemCount;
    QComboBox *comboBox_visualizations;
    QSlider *slider1;

    // splitter and the right side window
    QSplitter *splitter_main;

    // properties
    ObjectPropertyBrowser *pbrowser;

    // VTK
    vtkNew<vtkGenericOpenGLRenderWindow> renderWindow;
    QVTKOpenGLNativeWidget *qt_vtk_widget;
    vtkNew<vtkRenderer> renderer;
    vtkNew<vtkScalarBarActor> scalarBar;

    // layout for the main part of the window
    QHBoxLayout *right_side_layout;
    QWidget *right_side_container;

    // layout for charts
    QWidget *container_charts;
    QVBoxLayout *layout_charts;

    // chart for indenter force
    QChartView *cvForce;
    QChart *chart_line_indenter_force;
    QLineSeries *series_indenter_force;

    // screenshot
    vtkNew<vtkWindowToImageFilter> windowToImageFilter;
    vtkNew<vtkPNGWriter> writerPNG;

    // other
    void OpenH5(QString fileName);
    void updateGUIFileInfo();       // update file name and nubmer of nodes,elems,czs
    void GeneratePlots();
};
#endif // MAINWINDOW_H
