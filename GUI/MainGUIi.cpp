#include "MainGUI.hpp"

MainGUI::MainGUI(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainGUI)
{

    ui->setupUi(this);

}

MainGUI::~MainGUI()
{
    delete ui;
}
