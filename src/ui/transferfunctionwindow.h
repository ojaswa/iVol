#ifndef TRANSFERFUNCTIONWINDOW_H
#define TRANSFERFUNCTIONWINDOW_H

#include <QWidget>
class QVTKWidget;

namespace Ui {
class TransferFunctionWindow;
}

class TransferFunctionWindow : public QWidget
{
    Q_OBJECT

public:
    explicit TransferFunctionWindow(QWidget *parent = 0);
    ~TransferFunctionWindow();
    QVTKWidget* getQVTKWidget();

private:
    Ui::TransferFunctionWindow *ui;
};

#endif // TRANSFERFUNCTIONWINDOW_H
