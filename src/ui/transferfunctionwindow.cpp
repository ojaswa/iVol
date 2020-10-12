#include "transferfunctionwindow.h"
#include "ui_transferfunctionwindow.h"

#ifdef Q_OS_OSX
#include "osxHelper.h"
#endif

TransferFunctionWindow::TransferFunctionWindow(QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TransferFunctionWindow)
{
    ui->setupUi(this);
#ifdef Q_OS_OSX
    disableGLHiDPI(ui->qvtkWidget->winId());
#endif
}

TransferFunctionWindow::~TransferFunctionWindow()
{
    delete ui;
}

 QVTKWidget* TransferFunctionWindow::getQVTKWidget()
 {
     return ui->qvtkWidget;
 }
