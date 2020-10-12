#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QToolButton>

#include "algorithm/volumemanager.h"
#include "ui/transferfunctionwindow.h"

#include <vtkImageData.h>

namespace Ui {
    class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();
    void readVolume(QString filename);
    void displayVolume();
    template<class T> void displayLHScatter();
    void displayLHDensity();

private slots:
    void on_action_Open_triggered();
    void on_action_Quit_triggered();
    void on_actionCompute_LH_triggered();
    void on_actionSave_LH_triggered();
    void on_actionK_means_triggered();

private:
    Ui::MainWindow *ui;
    VolumeManager *volumeManager;
    TransferFunctionWindow *transferFunctionWindow;
};

#endif // MAINWINDOW_H
