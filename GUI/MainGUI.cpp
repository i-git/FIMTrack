/*****************************************************************************
 * Copyright (c) 2011-2014 The FIMTrack Team as listed in CREDITS.txt        *
 * http://fim.uni-muenster.de                                             	 *
 *                                                                           *
 * This file is part of FIMTrack.                                            *
 * FIMTrack is available under multiple licenses.                            *
 * The different licenses are subject to terms and condition as provided     *
 * in the files specifying the license. See "LICENSE.txt" for details        *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * FIMTrack is free software: you can redistribute it and/or modify          *
 * it under the terms of the GNU General Public License as published by      *
 * the Free Software Foundation, either version 3 of the License, or         *
 * (at your option) any later version. See "LICENSE-gpl.txt" for details.    *
 *                                                                           *
 * FIMTrack is distributed in the hope that it will be useful,               *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *
 * GNU General Public License for more details.                              *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * For non-commercial academic use see the license specified in the file     *
 * "LICENSE-academic.txt".                                                   *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 * If you are interested in other licensing models, including a commercial-  *
 * license, please contact the author at fim@uni-muenster.de      			 *
 *                                                                           *
 *****************************************************************************/

#include "MainGUI.hpp"
#include "ui_MainGUI.h"

using namespace cv;
using std::vector;
using std::string;

MainGUI::MainGUI(QWidget *parent) : 
    QMainWindow(parent), 
    ui(new Ui::MainGUI)
{
    
    try
    {
        ui->setupUi(this);
//        this->readParameters();
        this->mScene = this->ui->graphicsView->getScene();
        this->ui->labelCameraParameter->setStyleSheet("QLabel {color : red; }");
        this->ui->labelCameraParameter->setText("No Cameraparameter");
        
        this->ui->treeView->setColumnWidth(0, 120);
        
        /* initialice pointer */
        this->mPreferencesDialogWindow = NULL;
        this->mTrackingThread = NULL;
        this->mLogWindow = NULL;
        
        this->mPreferencesDialogWindow = new PreferencesDialogWindow(this);
        this->mLogWindow = new LOGWindow(this);
        
        this->mTrackingThread = new QThread();
        
        dynamic_cast<TrackerGraphicsView*>(this->ui->graphicsView)->enableRIOContextMenu();
        dynamic_cast<TrackerGraphicsView*>(this->ui->graphicsView)->enableContextMenu(true);
        
        qRegisterMetaType<cv::Mat>("cv::Mat");
        qRegisterMetaType<std::vector<std::vector<std::string> > >("std::vector<std::vector<std::string> >");
        qRegisterMetaType<Undistorter>("Undistorter");
//        qRegisterMetaType<RegionOfInterestContainer>("RegionOfInterestContainer");
        
        connect(this,SIGNAL(startTrackingSignal(std::vector<std::vector<std::string> >,bool,Undistorter, const RegionOfInterestContainer*)),
                &mTracker, SLOT(startTrackingSlot(std::vector<std::vector<std::string> >,bool,Undistorter, const RegionOfInterestContainer*)));
        
        connect(&mTracker, SIGNAL(previewTrackingImageSignal(cv::Mat)),
                this, SLOT(previewTrackingImageSlot(cv::Mat)));
        
        connect(&mTracker, SIGNAL(progressBarChangeSignal(int)),
                this,SLOT(on_progressBar_valueChanged(int)));
        
        connect(&mTracker, SIGNAL(trackingDoneSignal()), this, SLOT(trackingDoneSlot()));
        
        connect(ui->actionQuit, SIGNAL(triggered()), qApp, SLOT(quit()));
        
        connect(this, SIGNAL(logMessage(QString,LOGLEVEL)), Logger::getInstance(), SLOT(handleLogMessage(QString,LOGLEVEL)));
        connect(Logger::getInstance(), SIGNAL(newMessageForMainGUI(QString)), 
                this, SLOT(showMessage(QString)));
        
        /* Preferencesdialog Connections */
        connect(this->mPreferencesDialogWindow,      SIGNAL(loadCameraParameter(QString)),                   this,                               SLOT(initUndistorer(QString)));
        connect(this,                               SIGNAL(loadCameraParameter(QString)),                   this->mPreferencesDialogWindow,      SLOT(updateCameraParameterFile(QString)));
        connect(this,                               SIGNAL(updatePreferenceParameter()),                    this->mPreferencesDialogWindow,      SLOT(setParameters ()));
        connect(this,                               SIGNAL(setBackgroundSpinboxesValues(unsigned int)),     this->mPreferencesDialogWindow,      SLOT(updateBackgroungSpinboxesValues(unsigned int)));
        
//        connect(this->ui->treeView, SIGNAL(itemSelectionChanged()), this, SLOT(updateCurrentJobID()));
        
        
        this->setupBaseGUIElements(false);
        this->showMaximized();
        emit logMessage("FIMTrack started!", INFO);
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::~MainGUI"), FATAL);
    }
}

MainGUI::~MainGUI()
{
    try
    {
        emit logMessage("Program terminated!", DEBUG);
        
        delete this->ui;
        
        if(this->mTrackingThread)
        {
            this->mTrackingThread->quit();
            this->mTrackingThread->wait();
            delete this->mTrackingThread;
        }
        
        if(this->mPreferencesDialogWindow)
        {
            delete this->mPreferencesDialogWindow;
        }
        
        if(this->mLogWindow)
        {
            delete this->mLogWindow;
        }
        
        if(GeneralParameters::bSaveLog)
        {
            Logger::getInstance()->saveLog();
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::~MainGUI"), FATAL);
    }
}

void MainGUI::on_btnLoad_clicked()
{
    try
    {
        QStringList fileNames = QFileDialog::getOpenFileNames(this, tr("Open Image Files"), "", tr("Image Files (*.tif *.tiff)"));
        if(!fileNames.isEmpty())
        {
            this->mFileNames = fileNames;
            this->updateTreeViewer();
            this->setupBaseGUIElements(true);
            this->showImage(this->mFileNames.first());
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_btnLoad_clicked"), FATAL);
    }
    
}

void MainGUI::showImage(QString path/*const QModelIndex &index*/)
{
    try
    {
        //        if(index.isValid())
        //        {
        this->readParameters();
        
        // get the selected image fileNames list
        //            cv::Mat img = cv::imread(QtOpencvCore::qstr2str(this->fileNames.at(index.row())), 0);
        cv::Mat img = cv::imread(QtOpencvCore::qstr2str(path), 0);
        
        // undistor the inputimage
        if(this->mUndistorter.isReady())
        {
            cv::Mat dst;
            this->mUndistorter.getUndistortImage(img, dst);
            dst.copyTo(img);
        }
        
        // generate a contours container
        contoursType contours;
        
        contoursType collidedContours;
        
        // preprocess the image for preview
        Preprocessor::preprocessPreview2(img, 
                                         contours, 
                                         collidedContours, 
                                         /*this->ui->spinBox_graythresh->value()*/GeneralParameters::iGrayThreshold, 
                                         /*this->ui->spinBox_minSizeThresh->value()*/GeneralParameters::iMinLarvaeArea, 
                                         /*this->ui->spinBox_maxSizeThresh->value()*/GeneralParameters::iMaxLarvaeArea);
        
        cv::Mat grayImg = img.clone();
        
        // set the color to BGR to draw contours and contour sizes into the image
        cv::cvtColor(img, img, cv::COLOR_GRAY2BGR);
        
        // draw contours into the image
        cv::drawContours(img,contours,-1, Scalar(255,0,0), 6);
        cv::drawContours(img,collidedContours,-1,Scalar(0,0,255),8);
        
        // draw contour sizes into the image
        for(contoursType::iterator it = contours.begin(); it != contours.end(); ++it)
        {
            // get the current contour size
            int currentContourSize = cv::contourArea((*it));
            
            // convert int to string
            std::stringstream ss;
            ss << currentContourSize;
            
            // draw the contour size into the image
            cv::putText(img, ss.str(), (*it).at(0), cv::FONT_HERSHEY_PLAIN, 3, Scalar(255,255,0), 4);
            
            // draw raw larval informations
            // generate a RawLarva Object
            RawLarva rawLarva(*it,grayImg);
            
            // get the discrete spine
            vector<Point> discreteSpine = rawLarva.getDiscreteSpine();
            
            vector<double> larvalRadii = rawLarva.getLarvalRadii();
            
            bool isCoiled = rawLarva.getIsCoiledIndicator();
            Scalar color;
            if(isCoiled)
                color = Scalar(100,0,180);
            else
                color = Scalar(0,255,255);
            
            // draw the sharpest contour point
            cv::circle(img,discreteSpine.at(0),2,Scalar(0,0,255),2);
            // draw the spine points
            for (unsigned int i=1; i<discreteSpine.size()-1; ++i)
            {
                cv::circle(img,discreteSpine.at(i),larvalRadii.at(i),color,2);
            }
            // draw the second sharpest contour point
            cv::circle(img,discreteSpine.at(discreteSpine.size()-1), 2, Scalar(0,255,0),2);
        }
        
        QImage qimg = QtOpencvCore::img2qimg(img);
        
        // convert the opencv image to a QPixmap (to show in a QLabel)
        QPixmap pixMap = QPixmap::fromImage(qimg);
        // scale pixMap image to fit into the QLabel
        //        pixMap = pixMap.scaled(this->ui->graphicsView->size(), Qt::KeepAspectRatio);
        
        this->mScene->setPixmap(pixMap);
        
        // plot the image in the QLabel of the main gui (imgLabel)
        //        this->ui->imgLabel->setPixmap(pixMap);
        //        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::showImage"), FATAL);
    }
}

// button "Delete Job"
void MainGUI::on_btnDelete_clicked()
{   
    int index = -1;
    
    if(this->ui->treeView->model()->hasChildren(this->ui->treeView->currentIndex()))
    {
        // index stores the top level index of the selected job or is -1 if the selected item is not a top-level index
        index = this->ui->treeView->indexOfTopLevelItem(this->ui->treeView->currentItem());
    }
    else
    {
        if(this->ui->treeView->model()->parent(this->ui->treeView->currentIndex()).isValid())
        {
            index = this->ui->treeView->indexOfTopLevelItem(this->ui->treeView->currentItem()->parent());
        }
    }
    
    // Index of the selected item in the job list is -1 if it is not a top-level item, thus the
    // job is only removed if a top level item (i.e. a job) is selected
    if(index > -1)
    {
        QTreeWidgetItem *item = this->ui->treeView->currentItem();
        this->ui->treeView->takeTopLevelItem(index);
        delete item;
        this->mJobs.removeAt(index);
        QTreeWidgetItemIterator it = QTreeWidgetItemIterator(this->ui->treeView, QTreeWidgetItemIterator::HasChildren);
        
        if(*it)
        {
            while(*it)
            {
                int idx = this->ui->treeView->indexOfTopLevelItem((*it));
                if(idx > -1 && idx >= index)
                {
                    QString text = (*it)->text(0);
                    int jobNumber = text.split(" ").at(1).toInt();
                    (*it)->setText(0, QString("Job ").append(QString::number(--jobNumber)));
                }
                ++it;
            }
        }
        else
        {
            this->mScene->clearScene();
        }
    }
    
    if(this->ui->treeView->topLevelItemCount() == 0)
    {
        this->setupBaseGUIElements(false);
        this->mScene->clearScene();
    }
}

// button "Reset All Jobs"
void MainGUI::on_btnReset_clicked()
{
    this->resetListViewe();
}

void MainGUI::on_btnTrack_clicked()
{
    try
    {
        if(this->ui->btnTrack->text() == "Track")
        {
            this->readParameters();

            QString msg;
            msg.append("Start Tracking with parameters: GrayThresh:");
            msg.append(QString::number(/*this->ui->spinBox_graythresh->value()*/GeneralParameters::iGrayThreshold));
            msg.append(" MinSizeThresh:");
            msg.append(QString::number(/*this->ui->spinBox_minSizeThresh->value()*/GeneralParameters::iMinLarvaeArea));
            msg.append(" MaxSizeThresh:");
            msg.append(QString::number(/*this->ui->spinBox_maxSizeThresh->value()*/GeneralParameters::iMaxLarvaeArea));
            emit logMessage(msg,DEBUG);

            this->mTracker.moveToThread(mTrackingThread);
            this->mTrackingThread->start();

            this->ui->btnTrack->setText(QString("Stop"));
            this->ui->btnPreview->setEnabled(false);
            this->ui->btnReset->setEnabled(false);
            this->ui->btnDelete->setEnabled(false);
            this->ui->btnLoad->setEnabled(false);
            this->ui->spinBox_graythresh->setEnabled(false);
            this->ui->spinBox_maxSizeThresh->setEnabled(false);
            this->ui->spinBox_minSizeThresh->setEnabled(false);
            this->ui->progressBar->setEnabled(true);
            this->ui->treeView->setEnabled(false);

            std::vector<std::vector<std::string> > multiImgPaths;
            std::vector<std::string> strList;
            for(int i = 0; i < this->ui->treeView->topLevelItemCount(); i++)
            {
                for(int j = 0; j < this->ui->treeView->topLevelItem(i)->childCount(); j++)
                {
                    QString path = this->ui->treeView->topLevelItem(i)->text(1);
                    path.append("/");
                    path.append(this->ui->treeView->topLevelItem(i)->child(j)->text(1));
                    strList.push_back(QtOpencvCore::qstr2str(path));
                }
                multiImgPaths.push_back(strList);
                strList.clear();
            }
            
            if(this->mScene->hasROIContainer())
            {
                emit startTrackingSignal(multiImgPaths, this->ui->checkBoxShowTrackingProgress->isChecked(), this->mUndistorter, this->mScene->getROIContainer());
            }
            else
            {
                emit startTrackingSignal(multiImgPaths, this->ui->checkBoxShowTrackingProgress->isChecked(),this->mUndistorter, NULL);
            }
        }
        
        else if(this->ui->btnTrack->text() == "Stop")
        {
            this->mTrackingThread->terminate();
            this->mTrackingThread->wait();

            this->mFileNames.clear();
            this->ui->btnTrack->setText(QString("Track"));
            this->ui->btnPreview->setEnabled(true);
            this->ui->btnReset->setEnabled(true);
            this->ui->btnDelete->setEnabled(true);
            this->ui->btnLoad->setEnabled(true);
            this-> ui->spinBox_graythresh->setEnabled(true);
            this-> ui->spinBox_maxSizeThresh->setEnabled(true);
            this->ui->spinBox_minSizeThresh->setEnabled(true);

            this->ui->treeView->setEnabled(true);
            this->ui->progressBar->setValue(0);
            this->ui->progressBar->setEnabled(false);
            this->mScene->setPixmap(QPixmap());

        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_btnTrack_clicked"), FATAL);
    }
}

void MainGUI::readParameters()
{
    try
    {
        GeneralParameters::iGrayThreshold = this->ui->spinBox_graythresh->value();
        GeneralParameters::iMinLarvaeArea = this->ui->spinBox_minSizeThresh->value();
        GeneralParameters::iMaxLarvaeArea = this->ui->spinBox_maxSizeThresh->value();
        GeneralParameters::bShowTrackingProgress = this->ui->checkBoxShowTrackingProgress->isChecked();
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::readParameters"), FATAL);
    }
}

void MainGUI::setupBaseGUIElements(bool enable)
{
    this->ui->btnDelete->setEnabled(enable);
    this->ui->btnPreview->setEnabled(enable);
    this->ui->btnReset->setEnabled(enable);
    this->ui->btnTrack->setEnabled(enable);
    
    this->ui->spinBox_graythresh->setEnabled(enable);
    this->ui->spinBox_graythresh->setValue(GeneralParameters::iGrayThreshold);
    this->ui->spinBox_maxSizeThresh->setEnabled(enable);
    this->ui->spinBox_maxSizeThresh->setValue(GeneralParameters::iMaxLarvaeArea);
    this->ui->spinBox_minSizeThresh->setEnabled(enable);
    this->ui->spinBox_minSizeThresh->setValue(GeneralParameters::iMinLarvaeArea);
    
    this->ui->checkBoxShowTrackingProgress->setEnabled(enable);
    this->ui->checkBoxShowTrackingProgress->setChecked(GeneralParameters::bShowTrackingProgress);
    
    this->ui->treeView->setEnabled(enable);
    this->ui->graphicsView->setEnabled(enable);
    
    this->ui->progressBar->setEnabled(enable);
}

void MainGUI::on_btnPreview_clicked()
{
    if(this->ui->treeView->currentItem()->parent())
    {
        QString path = this->ui->treeView->currentItem()->parent()->text(this->ui->treeView->currentColumn());
        path.append("/");
        path.append(this->ui->treeView->currentItem()->text(this->ui->treeView->currentColumn()));
        this->showImage(path);
    }
}

void MainGUI::previewTrackingImageSlot(cv::Mat img)
{
    try
    {
        QImage qimg = QtOpencvCore::img2qimg(img);
        
        // convert the opencv image to a QPixmap (to show in a QLabel)
        QPixmap pixMap = QPixmap::fromImage(qimg);     
        
        this->mScene->setPixmap(pixMap);
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::previewTrackingImageSlot"), FATAL);
    }
}

void MainGUI::trackingDoneSlot()
{
    try
    {
        //        this->ui->imgLabel->clear();
        this->mFileNames.clear();
        this->ui->btnTrack->setText(QString("Track"));
        this->ui->btnPreview->setEnabled(true);
        this->ui->btnReset->setEnabled(true);
        this->ui->btnDelete->setEnabled(true);
        this->ui->btnLoad->setEnabled(true);
        this-> ui->spinBox_graythresh->setEnabled(true);
        this-> ui->spinBox_maxSizeThresh->setEnabled(true);
        this->ui->spinBox_minSizeThresh->setEnabled(true);
        this->ui->treeView->setEnabled(true);
        this->ui->progressBar->setEnabled(false);
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::trackingDoneSlot"), FATAL);
    }
}

void MainGUI::initUndistorer(const QString &file)
{
    try
    {
        this->mUndistorter.setPath(QtOpencvCore::qstr2str(file));
        if(this->mUndistorter.isReady())
        {
            this->ui->labelCameraParameter->setStyleSheet("QLabel {color : green; }");
            this->ui->labelCameraParameter->setText("Cameraparameter Loaded");
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::initUndistorer"), FATAL);
    }
}

void MainGUI::updateTreeViewer()
{
    try
    {
        if(!this->mFileNames.isEmpty())
        {
            this->mJobs.push_back(mFileNames);
            QFileInfo fInfo(mFileNames.at(0));
            
            QTreeWidgetItem *item = new QTreeWidgetItem(this->ui->treeView);
            item->setText(0, QString("Job ").append(QString::number(this->mJobs.size())));
            item->setText(1, fInfo.absolutePath());
            this->ui->treeView->addTopLevelItem(item);
            
            for (int i = 0; i < this->mFileNames.size(); i++) 
            {
                QFileInfo fi = QFileInfo(this->mFileNames.at(i));
                QTreeWidgetItem *childItem = new QTreeWidgetItem();
                childItem->setText(1, QString(fi.fileName()));
                item->addChild(childItem);
            }
            
            emit logMessage("Load files...", DEBUG);
            emit logMessage(fInfo.absolutePath(), DEBUG);
            QString msg;
            msg.append(QString::number(mFileNames.size()));
            msg.append(" images selected!");
            emit logMessage(msg, DEBUG);
        }
        else
        {
            emit logMessage("Process \"Load files...\" was aborted", INFO);
        }
        
        emit setBackgroundSpinboxesValues(mFileNames.size());
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_btnLoad_clicked"), FATAL);
    }
}

void MainGUI::resetListViewe()
{
    try
    {
        this->mScene->clearScene();
        this->mFileNames.clear();
        this->ui->treeView->clear();
        this->mJobs.clear();
        
        this->setupBaseGUIElements(false);
        
        emit setBackgroundSpinboxesValues(0);
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_btnReset_clicked\n"), FATAL);
    }
}

void MainGUI::showMessage(QString msg)
{
    this->ui->statusBar->showMessage(msg);
}

void MainGUI::updateCurrentJobID()
{
    QTreeWidgetItem* item = this->ui->treeView->currentItem();
    if(item)
    {
        if(item->isSelected())
        {
            // we has a toplevelItem
            if(item->childCount() > 0)
            {
                // index stores the top level index of the selected job or is -1 if the selected item is not a top-level index
                int index = this->ui->treeView->indexOfTopLevelItem(item);
                
                // Index of the selected item in the job list is -1 if it is not a top-level item, thus the
                // job is only removed if a top level item (i.e. a job) is selected
                if(index > -1)
                {
//                    qDebug() << item->text(0);
                }
            }
            // we have a cildItem (Filenameitem)
            else
            {
                QTreeWidgetItem* parentItem = item->parent();
                if(parentItem)
                {
                    // index stores the top level index of the selected job or is -1 if the selected item is not a top-level index
                    int index = this->ui->treeView->indexOfTopLevelItem(parentItem);
                    
                    // Index of the selected item in the job list is -1 if it is not a top-level item, thus the
                    // job is only removed if a top level item (i.e. a job) is selected
                    if(index > -1)
                    {
//                        qDebug() << parentItem->text(0);
                    }
                }
            }
        }
    }
}

void MainGUI::on_actionPreferences_triggered()
{
    try
    {
        emit updatePreferenceParameter();
        this->mPreferencesDialogWindow->show();
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionPreferences_triggered"), FATAL);
    }
}

void MainGUI::on_actionOpen_LOG_Window_triggered()
{
    try
    {
        this->mLogWindow->show();
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionOpen_LOG_Window_triggered"), FATAL);
    }
}

void MainGUI::on_progressBar_valueChanged(int value)
{
    try
    {
        ui->progressBar->setValue(value);
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_progressBar_valueChanged"), FATAL);
    }
}

void MainGUI::on_actionLoad_Camera_Setting_triggered()
{
    try
    {
        QString file = QFileDialog::getOpenFileName(this, tr("Open Camera Parameter"), "", tr("Camera Parameter File (*.yaml)"));
        if(!file.isNull() && !file.isEmpty())
        {
            this->mUndistorter.setPath(QtOpencvCore::qstr2str(file));
            if(this->mUndistorter.isReady())
            {
                this->ui->labelCameraParameter->setStyleSheet("QLabel {color : green; }");
                this->ui->labelCameraParameter->setText("Cameraparameter Loaded");
            }
            emit loadCameraParameter(file);
        }
        
        CameraParameter::File = file;
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionLoad_Camera_Setting_triggered"), FATAL);
    }
}

void MainGUI::on_actionSave_triggered()
{
    try
    {
        if(this->mConfigurationFile.isNull() || this->mConfigurationFile.isEmpty())
        {
            this->mConfigurationFile = QFileDialog::getSaveFileName(this, "Save Configuration", "", "FIMTrack Configuration File (*.conf)");
            if(!this->mConfigurationFile.isNull() && !this->mConfigurationFile.isEmpty())
            {
                OutputGenerator::saveConfiguration(QtOpencvCore::qstr2str(this->mConfigurationFile));
            }
        }
        else
        {
            OutputGenerator::saveConfiguration(QtOpencvCore::qstr2str(this->mConfigurationFile));
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionSave_triggered"), FATAL);
    }
}

void MainGUI::on_actionSave_As_triggered()
{
    try
    {
        this->mConfigurationFile = QFileDialog::getSaveFileName(this, "Save Configuration As", "", "FIMTrack Configuration File (*.conf)");
        if(!this->mConfigurationFile.isNull() && !this->mConfigurationFile.isEmpty())
        {
            OutputGenerator::saveConfiguration(QtOpencvCore::qstr2str(this->mConfigurationFile));
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionSave_As_triggered"), FATAL);
    }
}

void MainGUI::on_actionOpen_triggered()
{
    try
    {
        this->mConfigurationFile = QFileDialog::getOpenFileName(this, tr("Open Configuration File"), "", "FIMTrack Configuration File (*.conf)");
        if(!this->mConfigurationFile.isNull() && !this->mConfigurationFile.isEmpty())
        {
            InputGenerator::loadConfiguration(QtOpencvCore::qstr2str(this->mConfigurationFile));

            /* Restore Cameraparameters if there are some */
            if(!CameraParameter::File.isNull() && !CameraParameter::File.isEmpty())
            {
                this->initUndistorer(CameraParameter::File);
            }
            else
            {
                this->ui->labelCameraParameter->setStyleSheet("QLabel {color : red; }");
                this->ui->labelCameraParameter->setText("No Cameraparameter");
            }
            this->setupBaseGUIElements(true);
        }
    }
    catch(...)
    {
        emit logMessage(QString("Fatal Error in MainGUI::on_actionLoad_Camera_Setting_triggered"), FATAL);
    }
}

void MainGUI::on_actionNew_triggered()
{
    TrackerConfig::reset();
    this->resetListViewe();
    this->mConfigurationFile = "";
    this->mUndistorter.setReady(false);
    this->ui->labelCameraParameter->setStyleSheet("QLabel {color : red; }");
    this->ui->labelCameraParameter->setText("No Cameraparameter");
    this->setupBaseGUIElements(false);
}

void MainGUI::on_actionResults_Viewer_triggered()
{
    this->mResultsViewer = new ResultsViewer(this);
    mResultsViewer->setAttribute(Qt::WA_DeleteOnClose);
    mResultsViewer->show();
}

void MainGUI::on_treeView_itemClicked(QTreeWidgetItem *item, int column)
{
    if(item->parent())
    {
        QString path = item->parent()->text(column);
        path.append("/");
        path.append(item->text(column));
        this->showImage(path);
    }
}
