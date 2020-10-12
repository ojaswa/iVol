#include "mainwindow.h"
#include "ui_mainwindow.h"

//Qt includes
#include <QFileDialog>

// VTK includes
#include <vtkSmartPointer.h>
// -- I/O
#include <vtkNrrdReader.h>
#include <vtkDICOMImageReader.h>
// --3D display
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkVolumeProperty.h>
#include <vtkPiecewiseFunction.h>
#include <vtkColorTransferFunction.h>
// --2D chart
#include <vtkContextView.h>
#include <vtkChartXY.h>
#include <vtkTable.h>
#include <vtkContextScene.h>
#include <vtkPlot.h>
#include <vtkFloatArray.h>
#include <vtkPlotPoints.h>
#include <vtkChartHistogram2D.h>
#include <vtkColorTransferFunction.h>

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingVolumeOpenGL2);
VTK_MODULE_INIT(vtkRenderingContextOpenGL2);
VTK_MODULE_INIT(vtkRenderingFreeType);

#ifdef Q_OS_OSX
#include "osxHelper.h"
#endif

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    this->setUnifiedTitleAndToolBarOnMac(true);
#ifdef Q_OS_OSX
    disableGLHiDPI(ui->qvtkWidget->winId());
#endif

    ui->qvtkWidget->show();
    volumeManager = new VolumeManager;
    transferFunctionWindow = new TransferFunctionWindow;
    transferFunctionWindow->setWindowFlags(Qt::Tool|Qt::WindowSystemMenuHint|Qt::WindowMinMaxButtonsHint);
}

MainWindow::~MainWindow()
{
    delete ui;
    if (volumeManager) delete volumeManager;
    if(transferFunctionWindow) delete transferFunctionWindow;

}

void MainWindow::on_action_Open_triggered()
{
    QFileDialog dialog(this); //declare a dialog box to help assist selection
    dialog.setNameFilter(
                tr("Supported formats (*.dcm *.nhdr *.nrrd *.vtk);;DICOM (*.dcm);;NRRD (*.nhdr *.nrrd);;VTK (*.vtk)"));

    dialog.setViewMode(QFileDialog::Detail); // set the view mode
    QStringList list; // List to store the selections
    if(dialog.exec()) //execute the dialog box
    {
        list = dialog.selectedFiles(); // capture list
    }
    readVolume(list.first());
    displayVolume();

}

void MainWindow::readVolume(QString filename)
{
    QFileInfo fi(filename);
    QString filetype = fi.suffix();

    fprintf(stderr, "Reading %s from disk...\n", filename.toStdString().c_str());
    if (filetype.compare("nhdr")==0 | filetype.compare("nrrd")==0) //Load NRRD file
    {
        vtkSmartPointer<vtkNrrdReader> reader = vtkSmartPointer<vtkNrrdReader>::New();
        reader->SetFileName(filename.toStdString().c_str());
        reader->Update();
        volumeManager->setData(reader->GetOutput());
    } else if (filetype.compare("dcm") ==0) //Load DICOM file
    {
        vtkSmartPointer<vtkDICOMImageReader> reader = vtkSmartPointer<vtkDICOMImageReader>::New();
        //Get parent folder from the selected dcm file
        QString folder  = fi.path();
        reader->SetDirectoryName(folder.toStdString().c_str());
        reader->Update();
        volumeManager->setData(reader->GetOutput());
    } else if (filetype.compare("vtk") == 0) // Load VTK file
    {
        //TODO:
    }

    int *dims  = volumeManager->getData()->GetDimensions();
    int size  = volumeManager->getData()->GetScalarSize();
    const char* type  = volumeManager->getData()->GetScalarTypeAsString();
    fprintf(stderr, "Read %dx%dx%d voxels, of %s of size %d bytes\n", dims[0], dims[1], dims[2], type, size);
}

void MainWindow::displayVolume()
{
    vtkRenderWindow *renWin = ui->qvtkWidget->GetRenderWindow();
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renWin->AddRenderer(renderer);
    renWin->Render();

    vtkSmartPointer<vtkSmartVolumeMapper> volumeMapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();
    volumeMapper->SetBlendModeToComposite(); // composite first
    volumeMapper->SetInputData(volumeManager->getData());
    volumeMapper->SetRequestedRenderModeToGPU();
    vtkSmartPointer<vtkVolumeProperty> volumeProperty = vtkSmartPointer<vtkVolumeProperty>::New();
    volumeProperty->ShadeOff();
    volumeProperty->SetInterpolationType(VTK_LINEAR_INTERPOLATION);

    vtkSmartPointer<vtkPiecewiseFunction> compositeOpacity = vtkSmartPointer<vtkPiecewiseFunction>::New();
    compositeOpacity->AddPoint(0.0,0.0);
    compositeOpacity->AddPoint(80.0,1.0);
    compositeOpacity->AddPoint(80.1,0.0);
    compositeOpacity->AddPoint(255.0,0.0);
    volumeProperty->SetScalarOpacity(compositeOpacity); // composite first.

    vtkSmartPointer<vtkColorTransferFunction> color = vtkSmartPointer<vtkColorTransferFunction>::New();
    color->AddRGBPoint(0.0, 0.0, 1.0, 1.0);
    color->AddRGBPoint(40.0, 1.0, 0.0, 0.0);
    color->AddRGBPoint(255.0, 1.0, 1.0, 1.0);
    volumeProperty->SetColor(color);

    vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
    volume->SetMapper(volumeMapper);
    volume->SetProperty(volumeProperty);
    renderer->AddViewProp(volume);
    renderer->ResetCamera();
    renWin->Render();
}

template <class T>
void MainWindow::displayLHScatter()
{
    QVTKWidget *vtkWidget = transferFunctionWindow->getQVTKWidget();
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->SetRenderWindow(vtkWidget->GetRenderWindow());

    //Create a table
    vtkSmartPointer<vtkTable> table = vtkSmartPointer<vtkTable>::New();
    vtkSmartPointer<vtkFloatArray> arrX = vtkSmartPointer<vtkFloatArray>::New();
    arrX->SetName("L");
    table->AddColumn( arrX.GetPointer() );
    vtkSmartPointer<vtkFloatArray> arrY = vtkSmartPointer<vtkFloatArray>::New();
    arrY->SetName("H");
    table->AddColumn( arrY.GetPointer() );

    int dim[3];
    vtkImageData *lhData = volumeManager->getLH();
    T *imageLHPtr = static_cast<T*> (lhData->GetScalarPointer());
    lhData->GetDimensions(dim);
    int size = dim[0]*dim[1]*dim[2];
    table->SetNumberOfRows(size);
    for(long i=0; i<size; i++) {
        table->SetValue(i, 0, *(imageLHPtr + i*2));
        table->SetValue(i, 1, *(imageLHPtr + i*2 + 1));
    }

    //Create a chart
    vtkSmartPointer<vtkChartXY> chart = vtkSmartPointer<vtkChartXY>::New();
    chart->SetTitle("LH plot");
    //chart->SetAutoAxes(true);
    //chart->SetAutoSize(true);
    //chart->SetLayoutStrategy(vtkChartXY::FILL_SCENE);

    view->GetScene()->AddItem(chart);
    vtkPlot *scatterPlot = chart->AddPlot(vtkChart::POINTS);
    scatterPlot->SetInputData(table.GetPointer(), 0, 1);
    scatterPlot->SetColor(0, 0, 0, 255);
    scatterPlot->SetWidth(1.0);
    scatterPlot->SetTooltipLabelFormat("(%x, %y)");

    vtkPlotPoints::SafeDownCast(scatterPlot)->SetMarkerStyle(vtkPlotPoints::CIRCLE);
    vtkPlotPoints::SafeDownCast(scatterPlot)->SetMarkerSize(0.5);

    vtkWidget->GetRenderWindow()->SetMultiSamples(0);
    vtkWidget->GetRenderWindow()->Render();
    transferFunctionWindow->show();

}

void MainWindow::displayLHDensity()
{
    //Display density image as a plot
    QVTKWidget *vtkWidget = transferFunctionWindow->getQVTKWidget();
    vtkSmartPointer<vtkContextView> view = vtkSmartPointer<vtkContextView>::New();
    view->SetRenderWindow(vtkWidget->GetRenderWindow());

    //Create a chart
    vtkSmartPointer<vtkChartHistogram2D> chart = vtkSmartPointer<vtkChartHistogram2D>::New();
    chart->SetTitle("LH log-density plot");
    view->GetScene()->AddItem(chart);
    chart->SetInputData(volumeManager->getLHDensity());
    vtkSmartPointer<vtkColorTransferFunction> transferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
    transferFunction->SetColorSpaceToDiverging();
    transferFunction->AddRGBPoint(volumeManager->getLHDensity()->GetScalarRange()[0], 0.436, 0.308, 0.631);
    transferFunction->AddRGBPoint(volumeManager->getLHDensity()->GetScalarRange()[1], 0.759, 0.334, 0.046);
    transferFunction->SetScaleToLog10();
    transferFunction->Build();
    chart->SetTransferFunction(transferFunction);

    vtkWidget->GetRenderWindow()->SetMultiSamples(0);
    vtkWidget->GetRenderWindow()->Render();
    transferFunctionWindow->show();
}

void MainWindow::on_action_Quit_triggered()
{
    QApplication::quit();
}

void MainWindow::on_actionCompute_LH_triggered()
{
    //Compute LH histogram
    switch (volumeManager->getData()->GetScalarType())
    {
    case VTK_UNSIGNED_CHAR:
        volumeManager->computeLH<unsigned char>();
        //displayLHScatter<unsigned char>();
        break;
    case VTK_UNSIGNED_SHORT:
        volumeManager->computeLH<unsigned short>();
        //displayLHScatter<unsigned short>();
        break;
    }
    displayLHDensity();
}

void MainWindow::on_actionSave_LH_triggered()
{
    volumeManager->writeLH();
}


void MainWindow::on_actionK_means_triggered()
{
       volumeManager->clusteringKMeans();
}
