#include <QFileDialog>
#include <QList>
#include <QPointF>
#include <algorithm>
#include <cmath>
#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <gmsh.h>

MainWindow::~MainWindow() {delete ui;}

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    // chart for indenter force
    chart_line_indenter_force = new QChart;
    series_indenter_force = new QLineSeries;
    series_indenter_force->setName("Indenter Force");
    chart_line_indenter_force->addSeries(series_indenter_force);
    chart_line_indenter_force->createDefaultAxes();

    cvForce = new QChartView;
    cvForce->setChart(chart_line_indenter_force);

    // layout for charts
    container_charts = new QWidget;
    layout_charts = new QVBoxLayout;
    container_charts->setLayout(layout_charts);
    layout_charts->addWidget(cvForce);
    container_charts->hide();

    // property browser
    pbrowser = new ObjectPropertyBrowser(this);

    // VTK
    qt_vtk_widget = new QVTKOpenGLNativeWidget();
    qt_vtk_widget->setRenderWindow(renderWindow);

    renderer->SetBackground(1.0,1.0,1.0);
    renderWindow->AddRenderer(renderer);

    // VTK - scalar bar
    renderer->AddActor(scalarBar);
    scalarBar->SetMaximumWidthInPixels(130);
    scalarBar->SetBarRatio(0.07);
    scalarBar->SetMaximumHeightInPixels(400);
    scalarBar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedDisplay();
    scalarBar->GetPositionCoordinate()->SetValue(0.01,0.015, 0.0);
    scalarBar->SetLabelFormat("%.1e");
    scalarBar->GetLabelTextProperty()->BoldOff();
    scalarBar->GetLabelTextProperty()->ItalicOff();
    scalarBar->GetLabelTextProperty()->ShadowOff();
    scalarBar->GetLabelTextProperty()->SetColor(0.1,0.1,0.1);


    // right frame
    right_side_container = new QWidget;
    right_side_layout = new QHBoxLayout;
    right_side_container->setLayout(right_side_layout);
    right_side_layout->setSpacing(0);
    right_side_layout->setMargin(0);
    right_side_layout->addWidget(qt_vtk_widget,1);
    right_side_layout->addWidget(container_charts, 1);

    // splitter
    splitter_main = new QSplitter(Qt::Orientation::Horizontal);
    splitter_main->addWidget(pbrowser);
    splitter_main->addWidget(right_side_container);
    splitter_main->setSizes(QList<int>({100, 500}));
    setCentralWidget(splitter_main);

    // toolbar - combobox
    comboBox_visualizations = new QComboBox();
    ui->toolBar->addWidget(comboBox_visualizations);

    // slider
    ui->toolBar->addSeparator();
    slider1 = new QSlider(Qt::Horizontal);
    ui->toolBar->addWidget(slider1);
    slider1->setTracking(true);
    slider1->setMinimum(0);
    slider1->setMaximum(0);
    connect(slider1, SIGNAL(valueChanged(int)), this, SLOT(sliderValueChanged(int)));

    // statusbar
    statusLabel = new QLabel();
    labelElapsedTime = new QLabel();
    labelStepCount = new QLabel();
    statusLabelStepFactor = new QLabel();
    labelNodeCount = new QLabel();
    labelElemCount = new QLabel();

    QSizePolicy sp;
    const int status_width = 60;
    sp.setHorizontalPolicy(QSizePolicy::Fixed);
    statusLabelStepFactor->setSizePolicy(sp);
    statusLabelStepFactor->setFixedWidth(status_width);
    labelStepCount->setSizePolicy(sp);
    labelStepCount->setFixedWidth(status_width);
    labelElapsedTime->setSizePolicy(sp);
    labelElapsedTime->setFixedWidth(status_width);
    labelNodeCount->setSizePolicy(sp);
    labelNodeCount->setFixedWidth(status_width);
    labelElemCount->setSizePolicy(sp);
    labelElemCount->setFixedWidth(status_width);

    ui->statusbar->addWidget(statusLabel);
    ui->statusbar->addPermanentWidget(labelElapsedTime);
    ui->statusbar->addPermanentWidget(labelStepCount);
    ui->statusbar->addPermanentWidget(statusLabelStepFactor);
    ui->statusbar->addPermanentWidget(labelNodeCount);
    ui->statusbar->addPermanentWidget(labelElemCount);

// anything that includes the Model
    meshRepresentation.mesh = &model.mesh;
    meshRepresentation.model = &model;

    pwrapper.p = &model.p;
    pbrowser->setActiveObject(&pwrapper);
    scalarBar->SetLookupTable(meshRepresentation.hueLut);
    renderer->AddActor(meshRepresentation.actor_mesh_deformable);
    renderer->AddActor(meshRepresentation.actor_czs);
    renderer->AddActor(meshRepresentation.cylinderActor);
    meshRepresentation.actor_mesh_deformable->SetVisibility(false);

    worker = new BackgroundWorker(&model);
    connect(worker, SIGNAL(workerPaused()), SLOT(background_worker_paused()));
    connect(worker, SIGNAL(stepCompleted()), SLOT(updateGUI()));

    // populate combobox
    QMetaEnum qme = QMetaEnum::fromType<icy::MeshRepresentation::VisOpt>();
    for(int i=0;i<qme.keyCount();i++) comboBox_visualizations->addItem(qme.key(i));

    // VTK - screenshot
    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn(); // read from the back buffer
    writerPNG->SetInputConnection(windowToImageFilter->GetOutputPort());

    connect(comboBox_visualizations, QOverload<int>::of(&QComboBox::currentIndexChanged),
            [=](int index){ comboboxIndexChanged_visualizations(index); });

    gmsh::initialize();

    // read/restore saved settings
    settingsFileName = QDir::currentPath() + "/cm.ini";
    QFileInfo fi(settingsFileName);

    if(fi.exists())
    {
        QSettings settings(settingsFileName,QSettings::IniFormat);
        QVariant var;

        splitter_main->restoreState(settings.value("splitter_main").toByteArray());
        comboBox_visualizations->setCurrentIndex(settings.value("vis_option").toInt());

        vtkCamera* camera = renderer->GetActiveCamera();
        renderer->ResetCamera();
        camera->ParallelProjectionOff();

        var = settings.value("camData");
        if(!var.isNull())
        {
            double *vec = (double*)var.toByteArray().constData();
            camera->SetPosition(vec[0],vec[1],vec[2]);
            camera->SetFocalPoint(vec[3],vec[4],vec[5]);
            camera->SetViewUp(vec[6],vec[7],vec[8]);
            camera->SetViewAngle(vec[9]);
        }
        // open last file
        var = settings.value("lastFile");
        if(!var.isNull()) OpenH5(var.toString());
    }


    updateGUI();

    connect(ui->action_quit, &QAction::triggered, this, &MainWindow::quit_triggered);

    connect(ui->action_simulation_single_step, &QAction::triggered, this, &MainWindow::single_step_triggered);
    connect(ui->action_simulation_start, &QAction::triggered, this, &MainWindow::start_triggered);
    connect(ui->action_GotoStep0, &QAction::triggered, this, &MainWindow::gotoStep0_triggered);
    connect(ui->actionTrim, &QAction::triggered, this, &MainWindow::trim_triggered);
    connect(ui->action_camera_reset, &QAction::triggered, this, &MainWindow::cameraReset_triggered);

    connect(ui->actionShow_Plots, &QAction::triggered, this, &MainWindow::showPlots_toggled);
    connect(ui->action_Import_Mesh, &QAction::triggered, this, &MainWindow::importMesh_triggered);

    connect(ui->actionSave, &QAction::triggered, this, &MainWindow::save_triggered);
    connect(ui->actionOpen, &QAction::triggered, this, &MainWindow::open_triggered);

    connect(ui->actionScreenshot, &QAction::triggered, this, &MainWindow::screenshot_triggered);
}

void MainWindow::OpenH5(QString fileName)
{
    QFileInfo fi(fileName);
    if (fileName.isEmpty() || fi.suffix()!="h5" || !fi.exists()) return;
    model.HDF5Open(fileName.toStdString());
    meshRepresentation.UnsafeSynchronizeTopology();

    updateGUI();
    updateGUIFileInfo();

    QSettings settings(settingsFileName,QSettings::IniFormat);
    settings.setValue("lastFile",fileName);
    qDebug() << "opened file " << fileName;
}

void MainWindow::updateGUI()
{
    slider1->blockSignals(true);
    slider1->setMaximum(model.stepHistory.size()-1);
    slider1->setValue(model.currentStep.stepNumber);
    slider1->setEnabled(!worker->running);
    slider1->blockSignals(false);

    ui->action_simulation_single_step->setEnabled(!worker->running);
    labelStepCount->setText(QString{"step: %1"}.arg(model.currentStep.stepNumber));
    statusLabelStepFactor->setText(QString{"ltf: %1"}.arg(-std::log10(model.currentStep.currentStepFactor), 5, 'f', 3, '0'));
    labelElapsedTime->setText(QString{"t: %1"}.arg(model.currentStep.time, 5, 'f', 3, '0'));

    meshRepresentation.UnsafeSynchronizeValues();

    renderWindow->Render();
}






void MainWindow::showEvent(QShowEvent*) {}

void MainWindow::quit_triggered() { this->close(); }

void MainWindow::closeEvent( QCloseEvent* event )
{
    // save settings and stop simulation
    QSettings settings(settingsFileName,QSettings::IniFormat);
    qDebug() << "MainWindow: closing main window; " << settings.fileName();

    settings.setValue("splitter_main", splitter_main->saveState());

    double data[10];
    renderer->GetActiveCamera()->GetPosition(&data[0]);
    renderer->GetActiveCamera()->GetFocalPoint(&data[3]);
    renderer->GetActiveCamera()->GetViewUp(&data[6]);
    data[9]=renderer->GetActiveCamera()->GetViewAngle();

    QByteArray arr((char*)&data[0], sizeof(double)*10);
    settings.setValue("camData", arr);
    settings.setValue("vis_option", comboBox_visualizations->currentIndex());

    // kill backgroundworker
    worker->Finalize();
    event->accept();
}

void MainWindow::comboboxIndexChanged_visualizations(int index)
{
    meshRepresentation.ChangeVisualizationOption(index);
    scalarBar->SetVisibility(index != 0 && index != 1);
    renderWindow->Render();
}

void MainWindow::cameraReset_triggered()
{
    renderer->ResetCamera();
    renderWindow->Render();
}

void MainWindow::single_step_triggered()
{
    model.Prepare();
    model.Step();
    updateGUI();
}


void MainWindow::start_triggered(bool checked)
{
    if(!worker->running && checked)
    {
        statusLabel->setText("starting simulation");
        worker->Resume();
    }
    else if(worker->running && !checked)
    {
        statusLabel->setText("pausing simulation");
        worker->Pause();
        ui->action_simulation_start->setEnabled(false);
    }
}

void MainWindow::background_worker_paused()
{
    // enable the "Start" button
    ui->action_simulation_start->blockSignals(true);
    ui->action_simulation_start->setEnabled(true);
    ui->action_simulation_start->setChecked(false);
    ui->action_simulation_start->blockSignals(false);
    statusLabel->setText("stopped");
    updateGUI();
}

void MainWindow::sliderValueChanged(int val)
{
    model.GoToStep(val);
    updateGUI();
}

void MainWindow::gotoStep0_triggered()
{
    model.GoToStep(0);
    updateGUI();
}

void MainWindow::trim_triggered()
{
    model.Trim();
    updateGUI();
}

void MainWindow::showPlots_toggled(bool arg1)
{
    container_charts->setVisible(arg1);
    qt_vtk_widget->setHidden(arg1);
    if(arg1)
        GeneratePlots();
    else
        renderWindow->Render();
}


void MainWindow::GeneratePlots()
{

}


void MainWindow::importMesh_triggered()
{
    if(worker->running) return;
    QString fileName = QFileDialog::getOpenFileName(this, "Open STL", QDir::current().path(), "STL Files (*.stl)");
    if (fileName.isEmpty() || !QFile::exists(fileName)) return;
    model.Import(fileName.toStdString());
    meshRepresentation.UnsafeSynchronizeTopology();
    updateGUI();
    updateGUIFileInfo();
}


void MainWindow::save_triggered()
{
    if(worker->running) return;
    QString fileName = QFileDialog::getSaveFileName(this, "Save H5", QDir::current().path(), "H5 Files (*.h5)");
    if (fileName.isEmpty()) return;
    model.HDF5SaveNew(fileName.toStdString());
    QSettings settings(settingsFileName,QSettings::IniFormat);
    settings.setValue("lastFile",fileName);
    updateGUIFileInfo();
}


void MainWindow::open_triggered()
{
    if(worker->running) return;
    QString fileName = QFileDialog::getOpenFileName(this, "Open H5", QDir::current().path(), "Mesh Files (*.h5)");
    OpenH5(fileName);
}



void MainWindow::screenshot_triggered()
{
    QFileInfo fi(QString::fromStdString(model.fileName));
    QString bn = fi.baseName();
    QDir dir(QDir::currentPath()+ "/" + bn);
    if(!dir.exists()) dir.mkdir(QDir::currentPath()+ "/" + bn);
    QString outputPath = dir.path() + "/" + QString::number(model.currentStep.stepNumber) + ".png";

    windowToImageFilter->SetInput(renderWindow);
    windowToImageFilter->SetScale(1); // image quality
    windowToImageFilter->SetInputBufferTypeToRGBA(); //also record the alpha (transparency) channel
    windowToImageFilter->ReadFrontBufferOn();
    windowToImageFilter->Update();
    windowToImageFilter->Modified();

    writerPNG->Modified();
    writerPNG->SetFileName(outputPath.toUtf8().constData());
    writerPNG->Write();
    qDebug() << "screenshot: " << outputPath;
}

void MainWindow::updateGUIFileInfo()
{
    this->setWindowFilePath(QString::fromStdString(model.fileName));
    qDebug() << "setWindowFilePath " << QString::fromStdString(model.fileName);
    labelNodeCount->setText(QString{"nds: %1"}.arg(model.mesh.nodes.size()));
    labelElemCount->setText(QString{"els: %1"}.arg(model.mesh.elems.size()));
}



void MainWindow::paintEvent(QPaintEvent *event)
{
//    if(model.replayMode)
//    {
//        renderWindow->WaitForCompletion();
//        model.cv_.notify_one();
//        qDebug() << "paintEvent taking screenshot " << model.currentStep;
        //on_actionScreenshot_triggered();
//    }
}

