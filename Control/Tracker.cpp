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

#include "Tracker.hpp"

using namespace cv;
using std::vector;
using std::string;
using std::map;

Tracker::Tracker(QObject *parent) : QObject(parent)
{
    qRegisterMetaType<LOGLEVEL>("LOGLEVEL");
    connect(this,SIGNAL(logMessageSignal(QString,LOGLEVEL)), Logger::getInstance(), SLOT(handleLogMessage(QString,LOGLEVEL)));
}

void Tracker::startTrackingSlot(const std::vector<std::vector<std::string> > &multiImgPaths,
                                bool showProgress,
                                const Undistorter &undist,
                                RegionOfInterestContainer const* ROIContainer)
{
    this->multiImgPaths.clear();
    
    
    this->multiImgPaths = multiImgPaths;
    this->showTrackingProgress = showProgress;
    
    emit logMessageSignal("Start Tracking!", INFO);
    
    for (vector<vector<string> >::const_iterator pathsIt = this->multiImgPaths.begin(); pathsIt != this->multiImgPaths.end(); ++pathsIt)
    {
        this->larvaID = 0;
        
        curRawLarvae.clear();
        this->mLarvaeContainer.removeAllLarvae();
        vector<string> imgPaths = *pathsIt;
        
        Backgroundsubtractor bs(imgPaths,undist);
        
        track(imgPaths, bs, undist, ROIContainer);
        
        this->mLarvaeContainer.interplolateLarvae();
        
        std::string imgPath = imgPaths.at(0);
        
        QFileInfo fi(QtOpencvCore::str2qstr(imgPath));
        QString absPath = fi.absolutePath();
        absPath.append("/output_");
        
        QTime time = QTime::currentTime();
        QDate date = QDate::currentDate();
        QString strDate = date.toString("yyyy-MM-dd");
        QString strTime = time.toString("hh-mm-ss");
        absPath.append(strDate);
        absPath.append("_");
        absPath.append(strTime);
        
        QDir dir = QDir::root();
        dir.mkpath(absPath);
        
        QString tablePath = absPath;
        tablePath.append("/table.csv");
        OutputGenerator::writeLarvaeInverted(QtOpencvCore::qstr2str(tablePath),this->mLarvaeContainer.getAllLarvae(),imgPaths.size());
        
        QString ymlPath = absPath;
        ymlPath.append("/output.yml");
        OutputGenerator::writeOutputLarva(QtOpencvCore::qstr2str(ymlPath),this->mLarvaeContainer.getAllLarvae(), imgPaths, undist.isReady(), ROIContainer, NULL);
        
        QString trackImgPath = absPath;
        trackImgPath.append("/tracks.tif");
        OutputGenerator::drawTrackingResults(QtOpencvCore::qstr2str(trackImgPath),imgPaths,this->mLarvaeContainer.getAllLarvae());
        
        QString trackImgNoNumbersPath = absPath;
        trackImgNoNumbersPath.append("/tracksNoNumbers.tif");
        OutputGenerator::drawTrackingResultsNoNumbers(QtOpencvCore::qstr2str(trackImgNoNumbersPath),imgPaths,this->mLarvaeContainer.getAllLarvae());
        
    }
    
    emit trackingDoneSignal();
    
    emit logMessageSignal("Tracking done!", INFO);
}

void Tracker::track(const std::vector<std::string> &imgPaths, const Backgroundsubtractor &bs, const Undistorter &undist, RegionOfInterestContainer const* ROIContainer)
{   
    unsigned int timePoint = 0;
    
    for (vector<string>::const_iterator cIt = imgPaths.begin(); cIt != imgPaths.end(); ++cIt)
    {
        emit logMessageSignal(QString("Process Image: ").append(QtOpencvCore::str2qstr(*cIt)), INFO);
        cv::Mat img;
        cv::Mat imgTmp;
        if(undist.isReady())
        {
            imgTmp = imread(*cIt,0);
            emit logMessageSignal(QString("Undistortion used!"), DEBUG);
            undist.getUndistortImage(imgTmp, img);
        }
        else
        {
            img = imread(*cIt,0);
        }
        
        if (ROIContainer != NULL)
        {
            cv::Mat mask = cv::Mat::zeros(img.size(), img.type());

            for(int i = 0; i < ROIContainer->getRegionOfInterests().size(); ++i)
            {
                switch(ROIContainer->getRegionOfInterests().at(i).getType())
                {
                case RegionOfInterest::RECTANGLE:
                    mask(QtOpencvCore::qRect2Rect(ROIContainer->getRegionOfInterests().at(i).getBoundingBox())) = 255;
                    break;
                case RegionOfInterest::ELLIPSE:
                    cv::Mat ellipseMask = cv::Mat::zeros(mask.size(), mask.type());
                    cv::ellipse(ellipseMask, QtOpencvCore::qRect2RotatedRect(ROIContainer->getRegionOfInterests().at(i).getBoundingBox()), cv::Scalar(255,255,255), CV_FILLED);
                    mask |= ellipseMask;
                    break;
                }
            }
            img &= mask;
        }
        
        cv::Mat previewImg;
        
        if(showTrackingProgress)
        {
            cvtColor(img, previewImg, CV_GRAY2BGR);
        }
        
        extractRawLarvae(img, bs, &previewImg);
        //        this->assign(timePoint, minDist);
        this->assignByHungarian(timePoint);
        
        if(showTrackingProgress)
        {
            std::vector<Larva> larvae = this->mLarvaeContainer.getAllLarvae();
            for (vector<Larva>::iterator larvaIt = larvae.begin(); larvaIt != larvae.end(); ++larvaIt)
            {
                spineType spine;
                string goText ="";
                std::stringstream ss;
                if(larvaIt->getSpineAt(timePoint,spine))
                {
                    bool isCoiled;
                    larvaIt->getIsCoiledIndicatorAt(timePoint,isCoiled);
                    Scalar color;
                    if(isCoiled)
                    {
                        color = Scalar(100,0,180);
                    }
                    else
                    {
                        color = Scalar(0,255,255);
                    }
                    
                    for(vector<Point>::iterator sIt = spine.begin(); sIt != spine.end(); ++sIt)
                    {
                        circle(previewImg,(*sIt),2,color,2);
                    }
                    circle(previewImg,spine.at(0),3,Scalar(0,0,255),3);
                    circle(previewImg,spine.at(spine.size()-1),3,Scalar(255,0,0),3);
                    
                    int goPhaseIndicator;
                    if(larvaIt->getGoPhaseIndicatorAt(timePoint,goPhaseIndicator))
                    {
                        if(goPhaseIndicator == 1)
                            goText = "go";
                        else if(goPhaseIndicator == 0)
                        {
                            bool leftBended=false;
                            bool rightBended=false;
                            
                            larvaIt->getLeftBendingIndicatorAt(timePoint, leftBended);
                            larvaIt->getRightBendingIndicatorAt(timePoint, rightBended);
                            if(leftBended)
                            {
                                goText = "stop:left";
                            }
                            else if (rightBended)
                            {
                                goText = "stop:right";
                            }
                            else
                            {
                                goText = "stop";
                            }
                        }
                    }
                    ss << larvaIt->getID();
                    ss << ":" << goText;
                    putText(previewImg,ss.str(),spine.at(0),cv::FONT_HERSHEY_PLAIN, 2, Scalar(255,255,255), 2);
                }
            }
            emit previewTrackingImageSignal(previewImg);
        }
        
        ++timePoint;
        
        int progressbarValue = (int) ((timePoint*100) / imgPaths.size());
        emit progressBarChangeSignal(progressbarValue);
    }
}

void Tracker::extractRawLarvae(Mat const & img, Backgroundsubtractor const & bs, Mat * previewImg)
{
    contoursType contours;
    contoursType collidedContours;
    Preprocessor::preprocessTracking2(img,
                                      contours,
                                      collidedContours,
                                      GeneralParameters::iGrayThreshold,
                                      GeneralParameters::iMinLarvaeArea,
                                      GeneralParameters::iMaxLarvaeArea,bs);
    
    if(showTrackingProgress)
    {
        drawContours(*previewImg,contours,-1,Scalar(130,200,80),3);
        drawContours(*previewImg,collidedContours,-1,Scalar(0,0,255),8);
    }
    
    curRawLarvae.clear();
    curRawLarvae.reserve(contours.size());
    
    for (contoursType::const_iterator cIt = contours.begin(); cIt != contours.end(); ++cIt)
    {
        RawLarva rl((*cIt), img);
        curRawLarvae.push_back(rl);
    }
}

void Tracker::assignByHungarian(unsigned int timePoint)
{
    if(this->mLarvaeContainer.isEmpty())
    {
        for (vector<RawLarva>::const_iterator cIt = curRawLarvae.begin(); cIt != curRawLarvae.end(); ++cIt)
        {
            this->mLarvaeContainer.createNewLarva(timePoint, *cIt, larvaID);
            ++larvaID;
        }
    }
    else
    {
        cv::Point rawLarvaMomentum;
        cv::Point curLarvaMomentum;
        
        //        std::vector<Larva> larva    = this->mLarvaeContainer.getAllLarvae(timePoint-1);
        std::vector<int> validLarvaeIDs = this->mLarvaeContainer.getAllValidLarvaeIDS(timePoint-1);
        cv::Mat costMatrix              = cv::Mat::zeros(curRawLarvae.size(), validLarvaeIDs.size()/*this->mLarvaeContainer.getNumberOfLarvaAtTimeStep(timePoint-1)*/, CV_64F);
        
        for (size_t i = 0; i < curRawLarvae.size(); ++i)
        {
            rawLarvaMomentum = this->curRawLarvae.at(i).getMomentum();
            
            for (size_t j = 0; j < validLarvaeIDs.size(); ++j)
            {
                if(this->mLarvaeContainer.getMomentumAt(validLarvaeIDs.at(j), timePoint-1, curLarvaMomentum))
                {
                    double curDist = Calc::eucledianDist(curLarvaMomentum, rawLarvaMomentum);
                    
                    costMatrix.at<double>(i,j) = curDist;
                }
                else
                {
                    costMatrix.at<double>(i,j) = std::numeric_limits<double>::max();
                }
            }
        }
        Algorithms::Hungarian hs = Algorithms::Hungarian(costMatrix, Algorithms::HUNGARIAN_MODE_MINIMIZE_COST);
        cv::Mat assigments = hs.getAssignmentAsMatrix();
        bool foundLarvaForRawLarva = false;
        
        for(int i = 0; i < assigments.rows; ++i)
        {
            foundLarvaForRawLarva = false;
            for(int j = 0; j < assigments.cols; ++j)
            {
                Larva l_tmp;
                mLarvaeContainer.getLarva(validLarvaeIDs.at(j), l_tmp);
                if(assigments.at<uchar>(i,j) == 1 && this->larvaHasPointInRawLarva(timePoint, l_tmp, curRawLarvae.at(i)))
                {
                    foundLarvaForRawLarva = true;
                    this->mLarvaeContainer.insertRawLarva(validLarvaeIDs.at(j), timePoint, curRawLarvae.at(i));
                    break;
                }
            }
            
            if(!foundLarvaForRawLarva)
            {
                this->mLarvaeContainer.createNewLarva(timePoint, curRawLarvae.at(i), this->larvaID);
                ++this->larvaID;
            }   
        }
    }
}



bool Tracker::larvaHasPointInRawLarva(unsigned int const timePoint, Larva const & larva, RawLarva const& rawLarva)
{
    bool retBool = false;
    
    std::vector<cv::Point> curContour = rawLarva.getContour();
    std::vector<cv::Point> spine;
    cv::Point momentum;
    
    if (larva.getSpineAt(timePoint-1, spine))
    {
        larva.getMomentumAt(timePoint-1, momentum);
        double momInContour = cv::pointPolygonTest(curContour, momentum, false);
        // momentum in contour
        if (momInContour > 0)
        {
            retBool = true;
        }
        else
        {
            //            for(auto spinePoint : spine)
            foreach(auto spinePoint, spine)
            {
                double spinePointInContour = cv::pointPolygonTest(curContour, spinePoint, false);
                if (spinePointInContour > 0)
                {
                    retBool = true;
                    break;
                }
            }
        }
    }
    
    return retBool;
}
