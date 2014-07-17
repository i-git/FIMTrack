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

#include "ResultsViewer.hpp"
#include "ui_ResultsViewer.h"

ResultsViewer::ResultsViewer(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ResultsViewer)
{
    ui->setupUi(this);
    
    this->mDrawOpenCV       = false;
    
    this->mTimer            = new QTimer(this);
    this->mZoomFactor       = 0;
    this->mPlayingModeOn    = false;
    
    this->mScene            = this->ui->graphicsView->getScene();    
    this->mBtnAddTab        = new QToolButton(this->ui->larvaTabWidget);
    this->mBtnAddTab->setIcon(QIcon::fromTheme("list-add"));
    this->mBtnAddTab->setText("+");
    this->mBtnAddTab->setCursor(Qt::ArrowCursor);
    this->mBtnAddTab->setAutoRaise(true);
    this->ui->larvaTabWidget->setCornerWidget(this->mBtnAddTab, Qt::TopRightCorner);
    
    dynamic_cast<TrackerGraphicsView*>(this->ui->graphicsView)->enableLandmarkContextMenu();
    dynamic_cast<TrackerGraphicsView*>(this->ui->graphicsView)->enableContextMenu(true);
    this->ui->graphicsView->setEnabled(false);
    this->ui->pbtOneStepPrevTime->setEnabled(false);
    this->ui->pbtPlayPause->setEnabled(false);
    this->ui->pbtOneStepNextTime->setEnabled(false);
    this->ui->larvaTabWidget->setEnabled(false);
    this->ui->horizontalSlider_images->setEnabled(false);
    this->ui->tableView->setEnabled(false);
    this->ui->tabWidget_images->setEnabled(false);
    
    this->mPlottingTabVisible = false;
    this->mLarvaIDForCropping = 0;
    
    this->ui->tab_3->setLarvaeContainerPointer(&this->mLarvaeContainer);
    
    this->setupConnections();
    
    this->setWindowFlags(Qt::Window);
    this->showMaximized();
}

ResultsViewer::~ResultsViewer()
{
    delete ui;
}

void ResultsViewer::showImage(int const index)
{
    if (!this->mFileNames.empty())
    {
        /* get the selected image fileNames list */
        cv::Mat img = cv::imread(QtOpencvCore::qstr2str(this->mFileNames.at(index)), 0);
        
        this->mImageSize.setWidth(img.size().width);
        this->mImageSize.setHeight(img.size().height);
        
        emit sendNewImageSize(this->mImageSize);
        
        if(mUndistorer.isReady())
        {
            cv::Mat tmpImg;
            mUndistorer.getUndistortImage(img,tmpImg);
            tmpImg.copyTo(img);
        }
        
        cv::Mat cImg = cv::Mat::zeros(img.size(), CV_8UC3);
        cv::cvtColor(img,cImg,CV_GRAY2BGR);
        
        //        if(!this->mLoadedLarvae.empty())
        //        {
        //            //            this->mScene->updateLarvae(index);
        
        //            if(this->mDrawOpenCV)
        //            {
        //                for (auto it = this->mLoadedLarvae.begin(); it != this->mLoadedLarvae.end(); ++it)
        //                {
        //                    spineType spine;
        //                    if(it->getSpineAt(index,spine))
        //                    {
        //                        cv::Point head = spine.at(0);
        //                        cv::Point tail = spine.at(it->getNSpinePoints()-1);
        //                        cv::Point mid = spine.at(it->getSpineMidPointIndex());
        
        //                        bool isCoiled;
        //                        cv::Scalar color;
        //                        it->getIsCoiledIndicatorAt(index,isCoiled);
        //                        if(!isCoiled)
        //                            color = cv::Scalar(0,255,255);
        //                        else
        //                            color = cv::Scalar(255,0,255);
        
        //                        cv::circle(cImg,head,2,cv::Scalar(0,0,255),2);
        //                        cv::circle(cImg,tail,2,cv::Scalar(255,0,0),2);
        //                        cv::circle(cImg,mid,2,color,2);
        
        //                        cv::line(cImg,head,mid,color);
        //                        cv::line(cImg,mid,tail,color);
        
        //                        bool isLeftBended;
        //                        bool isRightBended;
        
        //                        cv::Scalar bendingTextColor(255,255,255);
        
        //                        if (it->getLeftBendingIndicatorAt(index,isLeftBended) &&
        //                                it->getRightBendingIndicatorAt(index,isRightBended))
        //                        {
        //                            if(isLeftBended)
        //                                bendingTextColor = cv::Scalar(0,0,255);
        //                            else if(isRightBended)
        //                                bendingTextColor = cv::Scalar(0,255,0);
        //                        }
        
        //                        std::stringstream ss;
        //                        ss << it->getID();
        //                        std::string text = "";
        
        //                        int goPhase;
        //                        if(it->getGoPhaseIndicatorAt(index,goPhase))
        //                        {
        //                            text = ":stop";
        
        //                            if(goPhase == 1)
        //                                text = ":go";
        
        //                            ss << text;
        
        //                        }
        //                        cv::putText(cImg,ss.str(),head,cv::FONT_HERSHEY_PLAIN, 1, bendingTextColor, 1);
        
        //                    } 
        //                }
        //            }
        //        }
        
        QImage qimg = QtOpencvCore::img2qimg(cImg);
        
        /* convert the opencv image to a QPixmap (to show in a QLabel) */
        QPixmap pixMap = QPixmap::fromImage(qimg);
        
        this->mScene->setPixmap(pixMap);
    }
}

void ResultsViewer::showTable(const int index)
{
    if (!this->mFileNames.empty())
    {
        if(this->mLarvaeContainer.hasLoadedLarvae())
        {
            unsigned int nLarvae = this->mLarvaeContainer.getNumberOfLarvae();
            Larva l;
            int nLandmarks = 0;
            QStringList landmarkNames;
            if( mScene->hasLandmarkContainer())
            {
                landmarkNames = mScene->getLandmarkContainer()->getAllLandmarkNames();
                nLandmarks = landmarkNames.size();
            }

            unsigned int nParameters = 16 +2*nLandmarks;

            
            QStandardItemModel *tableModel = new QStandardItemModel(nParameters,nLarvae,this);
            for (unsigned int cols = 0; cols < nLarvae; ++cols)
            {
                if(this->mLarvaeContainer.getLarva(cols, l))
                {
                    QString hHeader("larva(");
                    hHeader.append(QString::number(l.getID()));
                    hHeader.append(")");
                    tableModel->setHorizontalHeaderItem(cols, new QStandardItem(hHeader));
                }
            }

            tableModel->setVerticalHeaderItem(0, new QStandardItem(QString("area:")));
            tableModel->setVerticalHeaderItem(1, new QStandardItem(QString("mom_x:")));
            tableModel->setVerticalHeaderItem(2, new QStandardItem(QString("mom_y:")));
            tableModel->setVerticalHeaderItem(3, new QStandardItem(QString("bodyBending:")));
            tableModel->setVerticalHeaderItem(4, new QStandardItem(QString("isCoiled:")));
            tableModel->setVerticalHeaderItem(5, new QStandardItem(QString("isLeftBended:")));
            tableModel->setVerticalHeaderItem(6, new QStandardItem(QString("isRightBended:")));
            tableModel->setVerticalHeaderItem(7, new QStandardItem(QString("isGoPhase:")));
            tableModel->setVerticalHeaderItem(8, new QStandardItem(QString("movDirection:")));
            tableModel->setVerticalHeaderItem(9, new QStandardItem(QString("spineLength:")));
            tableModel->setVerticalHeaderItem(10, new QStandardItem(QString("perimeter:")));
            tableModel->setVerticalHeaderItem(11, new QStandardItem(QString("momDist:")));
            tableModel->setVerticalHeaderItem(12, new QStandardItem(QString("accDist:")));
            tableModel->setVerticalHeaderItem(13, new QStandardItem(QString("distToOrigin:")));
            tableModel->setVerticalHeaderItem(14, new QStandardItem(QString("velocity:")));
            tableModel->setVerticalHeaderItem(15, new QStandardItem(QString("acceleration:")));

            int landmarkIndex = 0;
            for (int i=16; i<nParameters; i+=2)
            {
                QString curLandmarkName = landmarkNames[landmarkIndex];
                tableModel->setVerticalHeaderItem(i, new QStandardItem(curLandmarkName+QString(" dist")));
                tableModel->setVerticalHeaderItem(i+1, new QStandardItem(curLandmarkName+QString(" angle")));
                ++landmarkIndex;
            }
            
            cv::Point mom;
            double area;
            double bodyBending;
            double spineLength;
            double perimeter;
            double momDist;
            double accDist;
            double distToOrigin;
            bool isCoiled;
            bool isLeftBended;
            bool isRightBended;
            int goPhaseIndicator;
            double movementDirection;
            int isCoiledInt;
            int isLeftBendedInt;
            int isRightBendedInt;
            double velocity;
            double acceleration;
            double distToLandmark;
            double bearingAngle;

            
            for (unsigned int cols = 0; cols < nLarvae; ++cols)
            {
                if(this->mLarvaeContainer.getLarva(cols, l))
                {
                    if (l.getMomentumAt(index,mom))
                    {
                        l.getAreaAt(index,area);
                        l.getMainBodyBendingAngleAt(index,bodyBending);
                        l.getSpineLengthAt(index,spineLength);
                        l.getPerimeterAt(index,perimeter);
                        l.getMomentumDistAt(index,momDist);
                        l.getAccDistAt(index,accDist);
                        l.getDistToOriginAt(index,distToOrigin);                        
                        
                        tableModel->setItem(0,cols, new QStandardItem(QString::number(area)));
                        tableModel->setItem(1,cols, new QStandardItem(QString::number(mom.x)));
                        tableModel->setItem(2,cols, new QStandardItem(QString::number(mom.y)));
                        tableModel->setItem(3,cols, new QStandardItem(QString::number(bodyBending)));
                        if(l.getIsCoiledIndicatorAt(index,isCoiled))
                        {
                            isCoiledInt = isCoiled ? 1 : 0;
                            tableModel->setItem(4,cols, new QStandardItem(QString::number(isCoiledInt)));
                        }
                        
                        if(l.getLeftBendingIndicatorAt(index,isLeftBended))
                        {
                            isLeftBendedInt = isLeftBended ? 1 : 0;
                            tableModel->setItem(5,cols, new QStandardItem(QString::number(isLeftBendedInt)));
                        }
                        
                        if(l.getRightBendingIndicatorAt(index,isRightBended))
                        {
                            isRightBendedInt = isRightBended ? 1 : 0;
                            tableModel->setItem(6,cols, new QStandardItem(QString::number(isRightBendedInt)));
                        }
                        
                        if(l.getGoPhaseIndicatorAt(index,goPhaseIndicator))
                        {
                            tableModel->setItem(7,cols, new QStandardItem(QString::number(goPhaseIndicator)));
                        }
                        
                        if(l.getMovementDirectionAt(index,movementDirection))
                        {
                            tableModel->setItem(8,cols, new QStandardItem(QString::number(movementDirection)));
                        }

                        l.getVelosityAt(index, velocity);
                        l.getAccelerationAt(index, acceleration);
                        
                        tableModel->setItem(9,cols, new QStandardItem(QString::number(spineLength)));
                        tableModel->setItem(10,cols, new QStandardItem(QString::number(perimeter)));
                        tableModel->setItem(11,cols, new QStandardItem(QString::number(momDist)));
                        tableModel->setItem(12,cols, new QStandardItem(QString::number(accDist)));
                        tableModel->setItem(13,cols, new QStandardItem(QString::number(distToOrigin)));
                        tableModel->setItem(14, cols, new QStandardItem(QString::number(velocity)));
                        tableModel->setItem(15, cols, new QStandardItem(QString::number(acceleration)));

                        int landmarkIndex = 0;
                        for (int i=16; i<nParameters; i+=2)
                        {
                            QString curLandmarkName = landmarkNames[landmarkIndex];
                            if(l.getDistanceToLandmark(index, curLandmarkName.toStdString(), distToLandmark))
                            {
                                l.getBearingAngleToLandmark(index, curLandmarkName.toStdString(), bearingAngle);
                                tableModel->setItem(i, cols, new QStandardItem(QString::number(distToLandmark)));
                                tableModel->setItem(i+1, cols, new QStandardItem(QString::number(bearingAngle)));
                            }
                            ++landmarkIndex;
                        }
                    }
                    
                }
            }
            
            this->ui->tableView->setModel(tableModel);
        }
    }
}

void ResultsViewer::setupConnections()
{
    connect(this->mTimer,                               SIGNAL(timeout()),                              this,                                               SLOT(moveSliderSlot()));
    connect(this->ui->pbtPlayPause,                     SIGNAL(clicked()),                              this,                                               SLOT(playPause()));
    connect(this->ui->pushButton_loadAllResults,        SIGNAL(clicked()),                              &this->mLarvaeContainer,                            SLOT(removeAllLarvae()));
    connect(this->ui->pushButton_SaveResults,           SIGNAL(clicked()),                              this,                                               SLOT(saveResultLarvae()));
    
    connect(this->ui->pbtOneStepPrevTime,               SIGNAL(clicked()),                              this,                                               SLOT(goOneTimeStepPrev()));
    connect(this->ui->pbtOneStepNextTime,               SIGNAL(clicked()),                              this,                                               SLOT(goOneTimeStepNext()));
    
    //    connect(this->mBtnAddTab,                       SIGNAL(clicked()),                        this,                   SLOT(addTab()));
    
    connect(&this->mLarvaeContainer,                    SIGNAL(sendLarvaModelDeleted()),                this,                                               SLOT(updateLarvaeTabs()));
    connect(&this->mLarvaeContainer,                    SIGNAL(sendRemovedResultLarvaID(uint)),         this,                                               SLOT(removeTabByLarvaID(uint)));
    connect(&this->mLarvaeContainer,                    SIGNAL(reset()),                                this,                                               SLOT(resetAndClear()));
    connect(&this->mLarvaeContainer,                    SIGNAL(sendUpdatedResultLarvaID(uint)),         this,                                               SLOT(updateLarvaePath(uint)));
    
    connect(this->mScene->getLandmarkContainer(),       SIGNAL(landmarkAdded(const Landmark*)),         &this->mLarvaeContainer,                            SLOT(updateLandmark(const Landmark*)));
    connect(this->mScene->getLandmarkContainer(),       SIGNAL(landmarkRemoved(QString)),               &this->mLarvaeContainer,                            SLOT(removeLandmark(QString)));
    
    connect(this->ui->pbtRemoveAllWhichAreShorterThen,  SIGNAL(clicked()),                              this,                                               SLOT(removeAllLarvaeModelsWhichAreShorterThen()));
    connect(this,                                       SIGNAL(sendShortestLarvaeTrackLength(uint)),    &this->mLarvaeContainer,                            SLOT(removeShortTracks(uint)));
    
    connect(this->ui->tabWidget_images,                 SIGNAL(currentChanged(int)),                    this,                                               SLOT(adjustTimeSlider(int)));
    
    connect(this->ui->tab_3,                            SIGNAL(sendCurrentPlottingParameter(QString,int)), this,                                            SLOT(adjustPlottingValues(QString,int)));
    connect(this,                                       SIGNAL(sendAvailableLarvaIDs(QStringList)),     this->ui->tab_3,                                    SLOT(setAvailableLarvaIDs(QStringList)));
    connect(this,                                       SIGNAL(sendPlottingTimeStemp(int)),             this->ui->tab_3,                                    SLOT(bookmarkCroppedTimeStemp(int)));
    connect(this,                                       SIGNAL(sendCroppedImage(QImage)),               this->ui->tab_3,                                    SLOT(showCroppedImage(QImage)));
    connect(this,                                       SIGNAL(sendNewImageSize(QSize)),                this->ui->tab_3,                                    SLOT(setImageSize(QSize)));
    connect(this->mScene->getLandmarkContainer(),       SIGNAL(sendAllLandmarkNames(QStringList)),      this->ui->tab_3,                                    SLOT(setAvailableLandmarkNames(QStringList)));
}

void ResultsViewer::connectTab(int index)
{
    // setup Larvatab-connections
    this->ui->larvaTabWidget->widget(index)->blockSignals(false);
    
    connect(this->ui->larvaTabWidget->widget(index),    SIGNAL(requestCenterOn(qreal, qreal)),      this,                                           SLOT(bringLarvaIntoFocus(qreal, qreal)));
    connect(this->ui->larvaTabWidget->widget(index),    SIGNAL(sendAmbiguityTimepoint(int)),        this,                                           SLOT(on_horizontalSlider_images_valueChanged(int)));       
    connect(this,                                       SIGNAL(newTimeStep(uint)),                  this->ui->larvaTabWidget->widget(index),        SLOT(update(uint)));
}

void ResultsViewer::disconnectTab(int index)
{
    this->ui->larvaTabWidget->widget(index)->blockSignals(true);
    disconnect(this->ui->larvaTabWidget->widget(index),    SIGNAL(requestCenterOn(qreal, qreal)),      this,                                        SLOT(bringLarvaIntoFocus(qreal, qreal)));
    disconnect(this->ui->larvaTabWidget->widget(index),    SIGNAL(sendAmbiguityTimepoint(int)),        this,                                        SLOT(on_horizontalSlider_images_valueChanged(int)));       
    disconnect(this,                                       SIGNAL(newTimeStep(uint)),                  this->ui->larvaTabWidget->widget(index),     SLOT(update(uint)));    
}

void ResultsViewer::resetAndClear()
{
    while(this->ui->larvaTabWidget->count() > 0)
    {
        this->removeTabByIndex(0);
    }
    
    this->mCurrentTimestep = 0;
    
    this->mFileNames.clear();
    this->mYmlFileName.clear();
    this->mUndistFileName.clear();
    
    this->mScene->clearScene();
    
    this->mImgPaths.clear();
    
    this->mNumberOfImages = 0;
    
    this->mPlayingModeOn = false;
    this->mTimer->stop();
    
    this->mUndistorer.reset();
    
    this->loadAllResults();
}

void ResultsViewer::removeAllLarvaeModelsWhichAreShorterThen()
{
    QString result = QInputDialog::getText(0, "Remove Lavaemodels", "Minimum Tracklength (in Timesteps):");
    if(!result.isNull() && !result.isEmpty())
    {
        bool ok;
        uint minTimeSteps = result.toUInt(&ok);
        if(ok)
        {
            emit sendShortestLarvaeTrackLength(minTimeSteps);
        }
    }   
}

void ResultsViewer::adjustTimeSlider(int tabID)
{
    if(tabID == 2)
    {
        this->mPlottingTabVisible = true;
        
        this->ui->horizontalSlider_images->setValue(this->mCurrentTimePointForPlotting + 1);
        this->ui->horizontalSlider_images->setMinimum(this->mTimeIntervalForPlotting.first);
        this->ui->horizontalSlider_images->setMaximum(this->mTimeIntervalForPlotting.second);
    }
    else
    {
        this->mPlottingTabVisible = false;
        this->ui->horizontalSlider_images->setMinimum(1);
        this->ui->horizontalSlider_images->setValue(this->mCurrentTimestep + 1);
        this->ui->horizontalSlider_images->setMaximum(this->mNumberOfImages);
    }
}

void ResultsViewer::adjustPlottingValues(QString larvaID, int currentTimeStep)
{    
    bool ok;
    int lID = larvaID.toUInt(&ok);
    if(ok)
    {
        this->mTimeIntervalForPlotting          = this->mLarvaeContainer.getStartEndTimesteps(lID);
        this->mCurrentTimePointForPlotting      = currentTimeStep;//this->mTimeIntervalForPlotting.first;
        this->mTimeIntervalForPlotting.first    += 1;
        this->mTimeIntervalForPlotting.second   += 1;
        
        QString qstr;
        qstr.append(QString::number(this->mTimeIntervalForPlotting.first));
        qstr.append("/");
        qstr.append(QString::number(this->mTimeIntervalForPlotting.second));
        this->ui->label_imageNumber->setText(qstr);
        
        this->ui->horizontalSlider_images->setMinimum(this->mTimeIntervalForPlotting.first);
        this->ui->horizontalSlider_images->setMaximum(this->mTimeIntervalForPlotting.second);
        this->ui->horizontalSlider_images->setValue(this->mCurrentTimePointForPlotting+1);
        
        this->mLarvaIDForCropping = lID;
        this->cropImage(this->mLarvaIDForCropping);
    }
}

void ResultsViewer::moveSliderSlot(bool forward)
{
    int position;
    if(forward)
    {
        position = (this->ui->horizontalSlider_images->value() < this->mNumberOfImages)
                ? this->ui->horizontalSlider_images->value() : 0;
        
        this->ui->horizontalSlider_images->setValue(++position);
    }
    else
    {
        position = (this->ui->horizontalSlider_images->value() > 1)
                ? this->ui->horizontalSlider_images->value() : this->mNumberOfImages + 1;
        
        this->ui->horizontalSlider_images->setValue(--position);
    }
}

void ResultsViewer::bringLarvaIntoFocus(qreal x, qreal y)
{
    this->ui->graphicsView->centerOn(x,y);
}

void ResultsViewer::playPause()
{
    if(this->mPlayingModeOn)
    {
        this->mTimer->stop();
        this->ui->pbtPlayPause->setIcon(QIcon(":/ResultViewer/Icons/Icons/media-playback-start.png"));
        this->ui->pbtOneStepPrevTime->setEnabled(true);
        this->ui->pbtOneStepNextTime->setEnabled(true);
    }
    else
    {
        this->mTimer->start(1000/10);
        this->ui->pbtPlayPause->setIcon(QIcon(":/ResultViewer/Icons/Icons/media-playback-pause.png"));
        this->ui->pbtOneStepPrevTime->setEnabled(false);
        this->ui->pbtOneStepNextTime->setEnabled(false);
    }
    
    this->mPlayingModeOn = !this->mPlayingModeOn;
}

void ResultsViewer::resetView()
{
    // Reset all lables
    this->ui->label_imagePath->setText(QString("No image path selected"));
    this->ui->label_trackingResults->setText(QString("No tracking result selected"));
    this->ui->label_cameraMatrix->setText(QString("No camera matrix selected"));
    this->ui->label_imageNumber->setText(QString("0/0"));
    this->ui->horizontalSlider_images->setValue(1);
    this->ui->horizontalSlider_images->setMaximum(1);
}

void ResultsViewer::loadAllResults()
{
    this->resetView();
    
    if(this->loadImageFiles()) 
    {
        this->mNumberOfImages = mFileNames.size();
        
        if(this->loadYmlFile())
        {
            this->mScene->loadROIContainer(this->mYmlFileName);
            this->mScene->loadLandmarkContainer(this->mYmlFileName);
            this->initUndistorer();
            this->setupBaseGUIElements();
            this->setupLarvaeTabs();
            
            this->showImage(0);
            this->showTable(0);
        }
    }
}

bool ResultsViewer::loadImageFiles()
{
    // 1. Load image files...
    this->mFileNames.clear();
    this->mYmlFileName.clear();
    this->mUseUndist = false;
    this->mImgPaths.clear();
    this->mScene->clearScene();
    
    this->mUndistFileName.clear();
    
    this->mFileNames = QFileDialog::getOpenFileNames(this, tr("Open Image Files"), "", tr("Image Files (*.tif *.tiff)"));
    
    return !this->mFileNames.isEmpty();
}

bool ResultsViewer::loadYmlFile()
{
    this->mYmlFileName = QFileDialog::getOpenFileName(this, tr("Open Output File"), "", tr("YML File (*.yml)"));
    
    if (!this->mYmlFileName.isEmpty())
    {
        this->mLarvaeContainer.readLarvae(this->mYmlFileName,
                                          this->mImgPaths,
                                          this->mUseUndist);
        
        if( this->mFileNames.size() != this->mImgPaths.size())
        {
            QMessageBox messageBox;
            messageBox.critical(0,"Error", "Number of loaded images do not fit with output file!");
            messageBox.setFixedSize(500,200);
            messageBox.show();
            this->mLarvaeContainer.removeAllLarvae();
            return false;
        }
        else
        {
            return true;
        }
    }
    else 
    {
        return false;
    }
}

void ResultsViewer::initUndistorer()
{
    if(this->mUseUndist)
    {
        this->mUndistFileName = QFileDialog::getOpenFileName(this, tr("Open Distortion File"), "", tr("YAML File (*.yaml)"));
        
        if (!this->mUndistFileName.isEmpty())
        {
            QFileInfo fileInfo = QFileInfo(this->mUndistFileName);
            this->ui->label_cameraMatrix->setText(fileInfo.path());
            this->mUndistorer.setPath(QtOpencvCore::qstr2str(this->mUndistFileName));
            
        }
    }
}

void ResultsViewer::setupBaseGUIElements()
{
    QFileInfo fileInfo = QFileInfo(this->mFileNames.at(0));
    this->ui->label_imagePath->setText(fileInfo.path());
    fileInfo = QFileInfo(this->mYmlFileName);
    this->ui->label_trackingResults->setText(fileInfo.path());
    
    QString qstr;
    qstr.append(QString::number(1));
    qstr.append("/");
    qstr.append(QString::number(this->mNumberOfImages));
    this->ui->label_imageNumber->setText(qstr);
    
    this->ui->horizontalSlider_images->setValue(1);
    this->ui->horizontalSlider_images->setMaximum(this->mNumberOfImages);
    
    this->ui->graphicsView->setEnabled(true);
    this->ui->pbtOneStepPrevTime->setEnabled(true);
    this->ui->pbtPlayPause->setEnabled(true);
    this->ui->pbtOneStepNextTime->setEnabled(true);
    this->ui->larvaTabWidget->setEnabled(true);
    this->ui->horizontalSlider_images->setEnabled(true);
    this->ui->tableView->setEnabled(true);
    this->ui->tabWidget_images->setEnabled(true);
    
    emit sendAvailableLarvaIDs(this->mLarvaeContainer.getAllLarvaeIDs());
}

void ResultsViewer::setupLarvaeTabs()
{
    // create one tab for each Larva
    for (int i = 0; i < this->mLarvaeContainer.getNumberOfLarvae(); ++i)
    {
        qsrand(QTime::currentTime().msec());
        TrackerSceneLarva* tLarva = this->mScene->addLarva(
                    this->mLarvaeContainer.getLarvaPointer(i), 
                    QColor(qrand() % 256, qrand() % 256, qrand() % 256));
        
        this->ui->larvaTabWidget->addTab(
                    new LarvaTab(tLarva, 
                                 this->mFileNames.size(),
                                 &this->mLarvaeContainer,
                                 this->ui->larvaTabWidget),
                    tLarva->getID());
        
        if(this->mLarvaeContainer.getAllTimesteps(this->mLarvaeContainer.getLarvaPointer(i)->getID()).front() > 0)
        {
            this->ui->larvaTabWidget->setTabEnabled(i, false);
        }
        
        this->connectTab(i);
    }
}

void ResultsViewer::addTab()
{
    //    qsrand(QTime::currentTime().msec());
    
    //    TrackerSceneLarva* tLarva = this->mScene->addLarva(
    //                this->mLarvaeContainer.createDefaultLarva(this->mCurrentTimestep), 
    //                QColor(qrand() % 256, qrand() % 256, qrand() % 256));
    
    //    this->ui->larvaTabWidget->addTab(
    //                new LarvaTab(tLarva, 
    //                             this->mFileNames.size(),
    //                             &this->mLarvaeContainer,
    //                             this->ui->larvaTabWidget),
    //                tLarva->getID());
    
    //    if(this->mLarvaeContainer.getAllTimeSteps(tLarva->getID()).front() > 0)
    //    {
    //        this->ui->larvaTabWidget->setTabEnabled(this->ui->larvaTabWidget->count() - 1, false);
    //    }
    
    //    this->connectTab(this->ui->larvaTabWidget->count() - 1);
}

void ResultsViewer::updateLarvaeTabs()
{
    QPair<QVector<uint>, QVector<uint> > visibilities = this->mLarvaeContainer.getVisibleLarvaID(this->mCurrentTimestep);
    int tabIndex;
    
    foreach(uint i, visibilities.first)
    {
        if(this->findLarvaTab(i, tabIndex))
        {
            this->ui->larvaTabWidget->setTabEnabled(tabIndex, true);
        }
    }
    
    foreach(uint i, visibilities.second)
    {
        if(this->findLarvaTab(i, tabIndex))
        {
            this->ui->larvaTabWidget->setTabEnabled(tabIndex, false);
        }
    }
}

void ResultsViewer::updateLarvaePath(uint larvaID)
{
    int i;
    if(this->findLarvaTab(larvaID, i))
    {
        dynamic_cast<LarvaTab*>(this->ui->larvaTabWidget->widget(i))->redrawPaths();
    }
}

void ResultsViewer::removeTabByLarvaID(uint larvaID)
{
    int i;
    if(this->findLarvaTab(larvaID, i))
    {
        dynamic_cast<LarvaTab*>(this->ui->larvaTabWidget->widget(i))->getSceneLarva()->removeFromScene();
        this->ui->larvaTabWidget->widget(i)->deleteLater();
        this->disconnectTab(i);
        this->ui->larvaTabWidget->widget(i)->close();
        this->ui->larvaTabWidget->removeTab(i);
    }
    this->updateLarvaPointer();
}

void ResultsViewer::removeTabByIndex(int tabIndex)
{
    if(tabIndex < this->ui->larvaTabWidget->count())
    {
        dynamic_cast<LarvaTab*>(this->ui->larvaTabWidget->widget(tabIndex))->getSceneLarva()->removeFromScene();
        this->ui->larvaTabWidget->widget(tabIndex)->deleteLater();
        this->disconnectTab(tabIndex);
        this->ui->larvaTabWidget->widget(tabIndex)->close();
        this->ui->larvaTabWidget->removeTab(tabIndex);
    }
}

bool ResultsViewer::findLarvaTab(uint larvaID, int& index)
{
    QStringList l;
    for(int i = 0; i < this->ui->larvaTabWidget->count(); ++i)
    {
        l = this->ui->larvaTabWidget->tabText(i).split(" ");
        if(larvaID == l.value(l.size()-1).toUInt())
        {
            index = i;
            return true;
        }
    }
    
    return false;
}

void ResultsViewer::updateLarvaPointer()
{
    int tabIndex;
    Larva* l;
    for(int i = 0; i < this->mLarvaeContainer.getNumberOfLarvae(); ++i)
    {
        l = this->mLarvaeContainer.getLarvaPointer(i);
        if(this->findLarvaTab(l->getID(), tabIndex))
        {
            dynamic_cast<LarvaTab *>(this->ui->larvaTabWidget->widget(tabIndex))->setLarvaPointer(l);
        }
    }
}

void ResultsViewer::saveResultLarvae()
{    
    
    QVector<TrackerSceneLarva*> larvaVector;
    for(int i = 0; i < this->ui->larvaTabWidget->count(); ++i)
    {
        larvaVector.push_back(dynamic_cast<LarvaTab*>(this->ui->larvaTabWidget->widget(i))->getSceneLarva());
    }
    
    this->mLarvaeContainer.saveResultLarvae(this->mImgPaths, this->mScene->getSaveImage(larvaVector), this->mUndistorer.isReady(), this->mScene->getROIContainer(), this->mScene->getLandmarkContainer());

}

void ResultsViewer::on_horizontalSlider_images_valueChanged(int value)
{
    if(!this->mPlottingTabVisible)
    {
        QString qstr;
        qstr.append(QString::number(value));
        qstr.append("/");
        qstr.append(QString::number(this->mNumberOfImages));
        this->ui->label_imageNumber->setText(qstr);
        //        this->ui->horizontalSlider_images->setValue(value);
        
        this->showImage(value - 1);
        this->showTable(value - 1);
        
        this->mCurrentTimestep = value - 1;
        this->updateLarvaeTabs();
        emit newTimeStep(this->mCurrentTimestep);
    }
    else
    {
        QString qstr;
        qstr.append(QString::number(value));
        qstr.append("/");
        qstr.append(QString::number(this->mTimeIntervalForPlotting.second));
        this->ui->label_imageNumber->setText(qstr);
        
        this->mCurrentTimePointForPlotting = value - 1;
               
        this->cropImage(this->mLarvaIDForCropping);
        
        emit sendPlottingTimeStemp(this->mCurrentTimePointForPlotting);
    }
}

void ResultsViewer::goOneTimeStepPrev()
{
    this->moveSliderSlot(false);
}

void ResultsViewer::goOneTimeStepNext()
{
    this->moveSliderSlot(true);
}

void ResultsViewer::cropImage(uint larvaID)
{
    Larva l;
    if (!this->mFileNames.empty() && this->mCurrentTimePointForPlotting < this->mFileNames.size() && this->mLarvaeContainer.getLarva(larvaID, l))
    {
        /* get the selected image fileNames list */
        cv::Mat img = cv::imread(QtOpencvCore::qstr2str(this->mFileNames.at(this->mCurrentTimePointForPlotting)), 0);
        
        this->mImageSize.setWidth(img.size().width);
        this->mImageSize.setHeight(img.size().height);
        emit sendNewImageSize(this->mImageSize);
        
        if(mUndistorer.isReady())
        {
            cv::Mat tmpImg;
            mUndistorer.getUndistortImage(img,tmpImg);
            tmpImg.copyTo(img);
        }
        
        cv::Point mom;
        double spineLength = this->mLarvaeContainer.getMaxSpineLength();
        if(spineLength > 0.0 && l.getMomentumAt(this->mCurrentTimePointForPlotting, mom))
        {
            cv::Rect box;
            
            box.x = mom.x - spineLength;
            if(box.x < 0)
                box.x = 0;
            
            box.y = mom.y - spineLength;
            if(box.y < 0)
                box.y = 0;
            
            box.width = 2*spineLength;
            if(box.width > img.cols)
                box.width = img.cols-1;
            
            box.height = 2*spineLength;
            if(box.height > img.rows)
                box.height = img.rows-1;
            
            cv::Mat cImg = cv::Mat::zeros(img.size(), CV_8UC3);
            cv::cvtColor(img,cImg,CV_GRAY2BGR);
            
            cImg = cImg(box);
            
            QImage qimg = QtOpencvCore::img2qimg(cImg);
            emit sendCroppedImage(qimg);
        }
        
    }
}
