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

#include "LarvaeContainer.hpp"

struct SplineSet
{
    double a;
    double b;
    double c;
    double d;
    double x;
};

LarvaeContainer::LarvaeContainer(QObject *parent) :
    QObject(parent)
{
    this->mMaxSpineLength = 0.0;
}

Larva* LarvaeContainer::createDefaultLarva(uint timeStep)
{
    Q_UNUSED(timeStep);
    return NULL;
}

void LarvaeContainer::readLarvae(const QString &ymlFileName, 
                                 std::vector<std::string> &imgPaths, 
                                 bool useUndist)
{
    InputGenerator::readOutputLarvae(QtOpencvCore::qstr2str(ymlFileName),
                                     this->mLarvae,
                                     imgPaths,
                                     useUndist);
    double spineLength;
    for(Larva l : this->mLarvae)
    {
        for(uint t : l.getAllTimeSteps())
        {
            if(l.getSpineLengthAt(t, spineLength))
            {
                this->mMaxSpineLength = std::max(spineLength, this->mMaxSpineLength);
            }
        }
    }
    
    for(int i = 0; i < this->mLarvae.size(); ++i)
        this->recalculateLarvaVelocityAndAcceleration(i);
}

void LarvaeContainer::removeAllLarvae()
{
    this->mLarvae.clear();
    emit reset();
}

void LarvaeContainer::updateLandmark(Landmark const* l)
{
    QPointF p;
    QLineF lin;
    QRectF r;
    QString name = l->getName();
    
    switch(l->getType())
    {
    case Landmark::POINT:
        p = l->getPoint();
        this->calcLandmarkParameter(name, p);
        break;
    case Landmark::LINE:
        lin = l->getLine();
        this->calcLandmarkParameter(name, lin);
        break;
    case Landmark::RECTANGLE:
        r = l->getRectagle();
        this->calcLandmarkParameter(name, r, false);
        break;
    case Landmark::ELLIPSE:
        r = l->getEllipse();
        this->calcLandmarkParameter(name, r, true);
        break;
    }
}

void LarvaeContainer::removeLandmark(QString name)
{
    std::vector<uint> timeSteps;
    std::string stdName = QtOpencvCore::qstr2str(name);
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        timeSteps = this->mLarvae.at(i).getAllTimeSteps();
        foreach(uint t, timeSteps)
        {
            this->mLarvae[i].parameters[t].distanceToLandmark.erase(stdName);
            this->mLarvae[i].parameters[t].isInLandmark.erase(stdName);
        }
    }
}

void LarvaeContainer::removeShortTracks(uint minTrackLenght)
{
    QVector<uint> removedLarvae;   
    this->mLarvae.erase(
                std::remove_if(this->mLarvae.begin(), this->mLarvae.end(),
                               
                               [&minTrackLenght, &removedLarvae](const Larva & l) -> bool
    {
        if(l.getAllTimeSteps().size() <= minTrackLenght)
        {
            removedLarvae.push_back(l.getID());
            return true;
        }
        else
        {
            return false;
        }
    }),this->mLarvae.end());
    
    foreach (uint i, removedLarvae) 
    {
        emit sendRemovedResultLarvaID(i);   
    }
}

QPair<QVector<uint>, QVector<uint> > LarvaeContainer::getVisibleLarvaID(uint time)
{
    QVector<uint> invisibleLarvae;
    QVector<uint> visibleLarvae;
    std::vector<uint> allTimeSteps;
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        allTimeSteps = this->mLarvae.at(i).getAllTimeSteps();
        if(std::find(allTimeSteps.begin(), allTimeSteps.end(), time) != allTimeSteps.end())
        {
            visibleLarvae << this->mLarvae.at(i).getID();
        }
        else
        {
            invisibleLarvae << this->mLarvae.at(i).getID();
        }
        
    }
    
    return QPair<QVector<uint>, QVector<uint> >(visibleLarvae, invisibleLarvae);
}

QStringList LarvaeContainer::getAllTimestepsBefore(uint id, uint time)
{
    QStringList result;
    size_t larvaIndex;
    
    if(this->getIndexOfLarva(id, larvaIndex))
    {
        std::vector<uint> timeSteps = this->mLarvae.at(larvaIndex).getAllTimeSteps();
        for(size_t i = 0; i < timeSteps.size(); ++i)
        {
            if(timeSteps.at(i) >= time)
            {
                std::reverse(result.begin(), result.end());
                return result;
            }
            else
            {
                result << QString::number(timeSteps.at(i)+1);
            }
        }
    }
    
    std::reverse(result.begin(), result.end());
    return result;
}

QStringList LarvaeContainer::getAllTimestepsAfter(uint id, uint time)
{
    QStringList result;
    size_t larvaIndex;
    
    if(this->getIndexOfLarva(id, larvaIndex))
    {
        std::vector<uint> timeSteps = this->mLarvae.at(larvaIndex).getAllTimeSteps();
        for(int i = static_cast<int>(timeSteps.size()) - 1; i >= 0; --i)
        {
            if(timeSteps.at(i) <= time)
            {
                std::reverse(result.begin(), result.end());
                return result;
            }
            else
            {
                result << QString::number(timeSteps.at(i)+1);
            }
        }
    }
    
    std::reverse(result.begin(), result.end());
    return result;
}

QStringList LarvaeContainer::getAllContemplableLarvaeIDsForAttach(uint id)
{
    QStringList result;
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        if(id != this->mLarvae.at(i).getID())
        {
            result << QString("Larva %1").arg(this->mLarvae.at(i).getID());
        }
    }
    
    return result;
}

void LarvaeContainer::invertLarva(uint larvaID, uint currentTime, uint toTime)
{    
    /* past */
    if(toTime < currentTime)
    {
        size_t larvaIndex;
        if(this->getIndexOfLarva(larvaID,larvaIndex))
        {
            for(int t = static_cast<int>(currentTime); t >= static_cast<int>(toTime); --t)
            {
                this->mLarvae.at(larvaIndex).invert(t);
            }
        }
    }
    
    /* now */
    if(currentTime == toTime)
    {
        size_t larvaIndex;
        if(this->getIndexOfLarva(larvaID,larvaIndex))
        {
            this->mLarvae.at(larvaIndex).invert(currentTime);
        }
    }
    
    /* future */
    else
    {
        size_t larvaIndex;
        if(this->getIndexOfLarva(larvaID,larvaIndex))
        {
            for(uint t = currentTime; t <= toTime; ++t)
            {
                this->mLarvae.at(larvaIndex).invert(t);
            }
        }
    }
}

void LarvaeContainer::attachToLarva(uint toLarvaID, uint fromLarvaID)
{
    size_t indexToLarva;
    size_t indexFromLarva;
    
    if(this->getIndexOfLarva(toLarvaID, indexToLarva) && this->getIndexOfLarva(fromLarvaID, indexFromLarva))
    {
        if((indexFromLarva != indexToLarva))
        {
            int minTime                             = qMin(this->mLarvae[indexToLarva].getAllTimeSteps().front(), 
                                                           this->mLarvae[indexFromLarva].getAllTimeSteps().front());
            
            int maxTime                             = qMax(this->mLarvae[indexToLarva].getAllTimeSteps().back(), 
                                                           this->mLarvae[indexFromLarva].getAllTimeSteps().back());
            
            for(int i = minTime; i <= maxTime; ++i)
            {
                if(this->mLarvae[indexFromLarva].parameters.find(i) == this->mLarvae[indexFromLarva].parameters.end())
                {
                    continue;
                }
                else
                {
                    this->mLarvae[indexToLarva].parameters[i] = this->mLarvae[indexFromLarva].parameters[i];
                }
            }
            
            std::vector<uint> timeSteps = this->mLarvae.at(indexToLarva).getAllTimeSteps();
            for(size_t i = 0; i < timeSteps.size(); ++i)
            {
                this->updateGoPhaseIndicator(indexToLarva, timeSteps.at(i));
                this->updateMovementDirection(indexToLarva, timeSteps.at(i));
            }
            
            this->recalculateLarvaDistanceParameter(toLarvaID); 
            this->recalculateLarvaVelocityAndAcceleration(toLarvaID);
            this->mLarvae.erase(this->mLarvae.begin() + indexFromLarva);
            emit sendRemovedResultLarvaID(fromLarvaID);
        }
    }
}

bool LarvaeContainer::copyModel2PrevTimeStep(uint larvaID, uint currentTime)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        this->mLarvae[larvaIndex].parameters[currentTime - 1] = this->mLarvae[larvaIndex].parameters[currentTime];
        
        std::vector<uint> timeSteps = this->mLarvae.at(larvaIndex).getAllTimeSteps();
        for(size_t i = 0; i < timeSteps.size(); ++i)
        {
            if(timeSteps.at(i) < currentTime)
            {
                continue;
            }
            else
            {
                this->updateGoPhaseIndicator(larvaID, timeSteps.at(i));
                this->updateMovementDirection(larvaID, timeSteps.at(i));
            }
        }
        this->recalculateLarvaDistanceParameter(larvaID);
        this->recalculateLarvaVelocityAndAcceleration(larvaIndex);
        return true;
    }
    return false;
}

bool LarvaeContainer::copyModel2NextTimeStep(uint larvaID, uint currentTime)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        this->mLarvae[larvaIndex].parameters[currentTime + 1] = this->mLarvae[larvaIndex].parameters[currentTime];
        this->updateGoPhaseIndicator(larvaID, currentTime + 1);
        this->updateMovementDirection(larvaID, currentTime + 1);
        this->recalculateLarvaDistanceParameter(larvaID);
        this->recalculateLarvaVelocityAndAcceleration(larvaIndex);
        return true;
    }
    return false;
}

void LarvaeContainer::recalculateLarvaDistanceParameter(uint larvaID)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        this->recalculateLarvaDistanceToOrigin(larvaIndex);
        this->recalculateLarvaMomentumDistance(larvaIndex);
    }
}

void LarvaeContainer::recalculateLarvaDistanceToOrigin(size_t larvaIndex)
{
    std::vector<uint> timeSteps = this->mLarvae.at(larvaIndex).getAllTimeSteps();
    cv::Point p0 = this->mLarvae.at(larvaIndex).getOrigin();
    cv::Point pt;
    foreach(uint t, timeSteps)
    {
        if(this->mLarvae.at(larvaIndex).getMomentumAt(t, pt))
        {
            this->mLarvae[larvaIndex].parameters[t].distToOrigin = Calc::eucledianDist(p0, pt);
        }
    }
}

void LarvaeContainer::recalculateLarvaMomentumDistance(size_t larvaIndex)
{
    std::vector<uint> timeSteps = this->mLarvae.at(larvaIndex).getAllTimeSteps();
    cv::Point p1;
    cv::Point p2;
    double accDist = 0.0;
    double momentumDist = 0.0;
    
    for(size_t i = 1; i < timeSteps.size(); ++i)
    {
        if(this->mLarvae.at(larvaIndex).getMomentumAt(timeSteps.at(i-1), p1) && this->mLarvae.at(larvaIndex).getMomentumAt(timeSteps.at(i), p2))
        {
            momentumDist = Calc::eucledianDist(p1,p2);
            accDist += momentumDist;
            
            this->mLarvae[larvaIndex].parameters[timeSteps.at(i)].momentumDist = momentumDist;
            this->mLarvae[larvaIndex].parameters[timeSteps.at(i)].accDist = accDist;
        }
    }
}

void LarvaeContainer::recalculateLarvaVelocityAndAcceleration(size_t larvaIndex)
{
    std::vector<uint> time = this->mLarvae.at(larvaIndex).getAllTimeSteps();
    cv::Point p0;
    cv::Point p1;
    double velosity = 0.0;
    int frameWindow = static_cast<int>(std::round(CameraParameter::dFPS));
    
    for(int i = time.size()-1; i >= frameWindow; --i)
    {
        int t0 = time.at(i - frameWindow);
        int t1 = time.at(i);
        
        if(std::abs(t0-t1) != frameWindow)
        {
            this->mLarvae[larvaIndex].parameters[t0].velosity = std::numeric_limits<double>::min();
        }
        else
        {            
            if (this->mLarvae.at(larvaIndex).getSpinePointAt(t0,mLarvae.at(larvaIndex).getNSpinePoints()-1, p0) &&
                    this->mLarvae.at(larvaIndex).getSpinePointAt(t1,mLarvae.at(larvaIndex).getNSpinePoints()-1, p1))
            {
                velosity     = Calc::eucledianDist(p0, p1) / static_cast<double>(std::abs(t1-t0));//(static_cast<double>(t0) - static_cast<double>(t1));
                velosity     *= CameraParameter::dScaleFactor;  // cm per second
                velosity     *= 10.0;                           // mm per second
                this->mLarvae[larvaIndex].parameters[t1].velosity = velosity;
            }
        }
    }
    
    for(int i = 0; i < time.size() && i < frameWindow; ++i)
    {
        int t = time.at(i);
        this->mLarvae[larvaIndex].parameters[t].velosity = std::numeric_limits<double>::min();
    }
    
    this->recalculateLarvaAcceleration(larvaIndex);
}

void LarvaeContainer::recalculateLarvaAcceleration(size_t larvaIndex)
{
    double v0, v1;
    double acceleration;
    std::vector<uint> time = this->mLarvae.at(larvaIndex).getAllTimeSteps();
    
    if (!time.empty())
    {
        // there is no acceleration value for the first timestep
        this->mLarvae[larvaIndex].parameters[time.front()].acceleration = std::numeric_limits<double>::min();

        for(int i = time.size()-1; i >= 1; --i)
        {
            int t0 = time.at(i-1);
            int t1 = time.at(i);

            if(this->mLarvae.at(larvaIndex).getVelosityAt(t0, v0) && this->mLarvae.at(larvaIndex).getVelosityAt(t1, v1))
            {
                acceleration = (v1 - v0) / static_cast<double>(std::abs(t1-t0));
                this->mLarvae[larvaIndex].parameters[t0].acceleration = acceleration;
            }
            else
            {
                this->mLarvae[larvaIndex].parameters[t0].acceleration = std::numeric_limits<double>::min();
            }
            
            this->mLarvae[larvaIndex].parameters[time.front()].acceleration = std::numeric_limits<double>::min();
        }
    }
}

void LarvaeContainer::calcLandmarkParameter(QString const& name, QPointF const& p)
{
    std::vector<uint> timeSteps;
    std::string stdName = QtOpencvCore::qstr2str(name);
    cv::Point momentum;
    double dist = 0.0;
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        timeSteps = this->mLarvae.at(i).getAllTimeSteps();
        foreach(uint t, timeSteps)
        {
            if(this->mLarvae.at(i).getMomentumAt(t, momentum))
            {
                dist = Calc::eucledianDist(momentum, p); 
                this->mLarvae[i].parameters[t].distanceToLandmark[stdName] = dist;
                this->mLarvae[i].parameters[t].isInLandmark[stdName] = false;
            }
        }
    }
    
    this->calcBearinAngle(name, p);
}

void LarvaeContainer::calcLandmarkParameter(QString const& name, QLineF const& l)
{
    std::vector<uint> timeSteps;
    std::string stdName = QtOpencvCore::qstr2str(name);
    cv::Point momentum;
    double dist = 0.0;
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        timeSteps = this->mLarvae.at(i).getAllTimeSteps();
        foreach(uint t, timeSteps)
        {
            if(this->mLarvae.at(i).getMomentumAt(t, momentum))
            {
                dist = Calc::eucledianDist(QtOpencvCore::point2qpoint(momentum), l); 
                this->mLarvae[i].parameters[t].distanceToLandmark[stdName] = dist;
                this->mLarvae[i].parameters[t].isInLandmark[stdName] = false;
            }
        }
    }
    
    this->calcBearinAngle(name, QPointF(l.p1().x() + l.dx()/2, l.p1().y() + l.dy()/2));
}

void LarvaeContainer::calcLandmarkParameter(const QString &name, const QRectF &r, const bool ellipse)
{
    std::vector<uint> timeSteps;
    std::string stdName = QtOpencvCore::qstr2str(name);
    cv::Point momentum;
    double dist = 0.0;
    QPointF p;
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        timeSteps = this->mLarvae.at(i).getAllTimeSteps();
        foreach(uint t, timeSteps)
        {
            if(this->mLarvae.at(i).getMomentumAt(t, momentum))
            {
                p = QtOpencvCore::point2qpoint(momentum);
                if(!r.contains(p))
                {
                    dist = Calc::eucledianDist(p, r, ellipse); 
                    this->mLarvae[i].parameters[t].distanceToLandmark[stdName] = dist;
                    this->mLarvae[i].parameters[t].isInLandmark[stdName] = false;
                }
                else
                {
                    this->mLarvae[i].parameters[t].distanceToLandmark[stdName] = 0.0;
                    this->mLarvae[i].parameters[t].isInLandmark[stdName] = true;
                }
                
            }
        }
    }
    
    this->calcBearinAngle(name, r.center());
}

void LarvaeContainer::calcBearinAngle(const QString &name, const QPointF &landmarkPoint)
{
    std::vector<uint> timeSteps;
    std::string stdName = QtOpencvCore::qstr2str(name);
    cv::Point momentum;
    cv::Point cvLandmarkPoint;
    double movementDirection;
    double offset = 10.0;
    double bearinAngle;
    
    cvLandmarkPoint.x = landmarkPoint.x();
    cvLandmarkPoint.y = landmarkPoint.y();
    
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        //        qDebug() << "LarvaID: " << this->mLarvae.at(i).getID() << "\n";
        timeSteps = this->mLarvae.at(i).getAllTimeSteps();
        foreach(uint t, timeSteps)
        {
            if(this->mLarvae.at(i).getMomentumAt(t, momentum) && this->mLarvae.at(i).getMovementDirectionAt(t, movementDirection))
            {
                cv::Point pRotated;
                movementDirection    = Calc::angleToRadian(movementDirection);
                
                double sinVal   = std::sin(movementDirection);
                double cosVal   = std::cos(movementDirection);
                
                /**
                 * because the angles are calculeted relativly to the y-axis 
                 * the meaning of the sin- and cos-values is switched. So if you look at the 
                 * circle-table you have to be carefully and keep this fact in mind.
                 */
                if(cosVal > 0.0 && sinVal > 0.0)
                {
                    pRotated.x = momentum.x + (offset * sinVal);
                    pRotated.y = momentum.y - (offset * cosVal);
                }
                else if(cosVal > 0.0 && sinVal < 0.0)
                {
                    pRotated.x = momentum.x + (offset * sinVal);
                    pRotated.y = momentum.y - (offset * cosVal);
                }
                else if(cosVal < 0.0 && sinVal > 0.0)
                {
                    pRotated.x = momentum.x + (offset * sinVal);
                    pRotated.y = momentum.y - (offset * cosVal);
                }
                else
                {
                    pRotated.x = momentum.x + (offset * sinVal);
                    pRotated.y = momentum.y - (offset * cosVal);
                }
                
                bearinAngle = Calc::calcSmallestAngle(momentum, pRotated, cvLandmarkPoint);
                
                //                QString debugString = QString("Time: %1,\t Momentum: [%2, %3],\t LandmarkPoint: [%4, %5],\t pRotated: [%6 %7],\t movementDirection(Radian): %8,\t bearinAngle: %9.")
                //                        .arg(t).arg(momentum.x).arg(momentum.y).arg(cvLandmarkPoint.x).arg(cvLandmarkPoint.y).arg(pRotated.x).arg(pRotated.y).arg(movementDirection).arg(bearinAngle);
                //                qDebug() << debugString;
                
                this->mLarvae[i].parameters[t].bearinAngle[stdName] = bearinAngle;
            }
        }
        //        qDebug() << "\n\n";
    }
}

bool LarvaeContainer::eraseLarvaAt(uint larvaID, uint time)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        this->eraseAt(larvaIndex, time);
        this->recalculateLarvaVelocityAndAcceleration(larvaIndex);
        emit sendLarvaModelDeleted();
        return true;
    }
    return false;
}

bool LarvaeContainer::eraseLarva(uint larvaID)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        this->mLarvae.erase(this->mLarvae.begin() + larvaIndex);
        emit sendRemovedResultLarvaID(larvaID);
        return true;
    }
    return false;
}

void LarvaeContainer::saveResultLarvae(const std::vector<std::string> &imgPaths, const QImage &img, 
                                       const bool useUndist, 
                                       const RegionOfInterestContainer *ROIContainer, 
                                       const LandmarkContainer *landmarkContainer)
{
    QString saveAs = QFileDialog::getSaveFileName(NULL, QString("Save (modified) Larvae As..."), NULL, tr("YAML-File (*.yml)"));
    if(!saveAs.isNull() && !saveAs.isEmpty())
    {
        OutputGenerator::writeOutputLarva(saveAs.toStdString(), this->mLarvae, imgPaths, useUndist, ROIContainer, landmarkContainer);
    }
    
    saveAs = QFileDialog::getSaveFileName(NULL, QString("Save (modified) Larvae As..."), NULL, tr("CSV-File (*.csv)"));
    if(!saveAs.isNull() && !saveAs.isEmpty())
    {
        OutputGenerator::writeLarvaeInverted(saveAs.toStdString(), this->mLarvae,imgPaths.size(), landmarkContainer);
    }
    
    saveAs = QFileDialog::getSaveFileName(NULL, QString("Save (modified) Larvae As..."), NULL, tr("TIF-File (*.tif)"));
    if(!saveAs.isNull() && !saveAs.isEmpty())
    {
        OutputGenerator::saveResultImage(saveAs, img);
    }
}

QStringList LarvaeContainer::getAllLarvaeIDs() const
{
    QStringList ids;
    for(Larva const& l : this->mLarvae)
    {
        ids << QString::number(l.getID());
    }
    
    return ids;
}

bool LarvaeContainer::getLarva(const uint index, Larva &l)
{
    if(index < this->mLarvae.size())
    {
        l = this->mLarvae.at(index);
        return true;
    }
    return false;
}

std::vector<uint> LarvaeContainer::getAllTimesteps(uint larvaID)
{
    size_t index;
    if(this->getIndexOfLarva(larvaID, index))
    {
        return this->mLarvae.at(index).getAllTimeSteps();
    }
    
    return std::vector<uint>();
}

QVector<double> LarvaeContainer::getAllTimestepsForPlotting(uint larvaID)
{
    QVector<double> res;
    size_t index;
    if(this->getIndexOfLarva(larvaID, index))
    {
        for(uint t : this->mLarvae.at(index).getAllTimeSteps())
        {
            res << static_cast<double>(t+1);
        }
    }
    
    return res;
}

QVector<QVector<double> > LarvaeContainer::getAllTimestepGaps(uint larvaID)
{
    QVector<double> allTimeSteps = this->getAllTimestepsForPlotting(larvaID);
    QVector<QVector<double> > res;
    
    int i = 0;
    
    while(i < allTimeSteps.size())
    {
        QVector<double> currentTimegap;
        while((i < allTimeSteps.size()-1) && (allTimeSteps.at(i) == allTimeSteps.at(i+1)-1))
        {
            currentTimegap << allTimeSteps.at(i);
            ++i;
        }
        
        currentTimegap << allTimeSteps.at(i);
        res << currentTimegap;
        ++i;
    }
    
    return res;
}

QPair<int, int> LarvaeContainer::getStartEndTimesteps(uint larvaID)
{
    QPair<int, int> res;
    size_t index;
    if(this->getIndexOfLarva(larvaID, index))
    {
        res.first   = this->mLarvae.at(index).getAllTimeSteps().front();
        res.second  = this->mLarvae.at(index).getAllTimeSteps().back();
    }
    
    return res;
}

std::vector<int> LarvaeContainer::getAllValidLarvaeIDS(uint timePoint)
{
    std::vector<int> res;
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        if(this->isAssignedAt(this->mLarvae.at(i).getID(), timePoint))
        {
            res.push_back(this->mLarvae.at(i).getID());
        }
    }
    return res;
}

bool LarvaeContainer::getSpineMidPointIndex(uint larvaID, uint& index) const
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))
    {
        index = this->mLarvae.at(i).getSpineMidPointIndex();
        return true;
    }
    
    return false;
}

bool LarvaeContainer::getSpinePointAt(const uint larvaID, const uint timePoint, const uint index, cv::Point &spinePoint) const
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))
    {
        return this->mLarvae.at(i).getSpinePointAt(timePoint, index, spinePoint);
    }
    
    return false;
}

bool LarvaeContainer::getAccDistAt(const uint larvaID, const uint timePoint, double &retAccDist)
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))
    {
        return this->mLarvae.at(i).getAccDistAt(timePoint, retAccDist);
    }
    
    return false;
}

bool LarvaeContainer::getDistToOriginAt(const uint larvaID, const uint timePoint, double &retDistToOrigin)
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))
    {
        return this->mLarvae.at(i).getDistToOriginAt(timePoint, retDistToOrigin);
    }
    
    return false;
}

bool LarvaeContainer::getIsCoiledIndicatorAt(const uint larvaID, const uint timePoint, bool &retIsCoiledIndicator)
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))
    {
        return this->mLarvae.at(i).getIsCoiledIndicatorAt(timePoint, retIsCoiledIndicator);
    }
    
    return false;
}

bool LarvaeContainer::getMomentumAt(const uint larvaID, const uint timePoint, cv::Point &retMomentum) const
{
    size_t i;
    if(this->getIndexOfLarva(larvaID, i))  
    {    
        this->mLarvae.at(i).getMomentumAt(timePoint, retMomentum);
        return true;    
    }

    return false;
}

QVector<cv::Point> LarvaeContainer::getAllMomentumValues(uint larvaID)
{
    QVector<cv::Point> res;
    size_t i;
    cv::Point val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getMomentumAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllAreaValues(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getAreaAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllMainBodybendingAngle(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getMainBodyBendingAngleAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllCoiledIndicator(uint larvaID)
{
    QVector<double> res;
    size_t i;
    bool val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getIsCoiledIndicatorAt(t, val))
            {
                res << static_cast<double>(val);
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllPerimeter(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getPerimeterAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllDistanceToOrigin(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getDistToOriginAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllMomentumDistance(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getMomentumDistAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllAccumulatedDistance(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getAccDistAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllGoPhaseIndicator(uint larvaID)
{
    QVector<double> res;
    size_t i;
    int val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getGoPhaseIndicatorAt(t, val))
            {
                res << static_cast<double>(val);
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllLeftBendingIndicator(uint larvaID)
{
    QVector<double> res;
    size_t i;
    bool val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getLeftBendingIndicatorAt(t, val))
            {
                res << static_cast<double>(val);
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllRightBendingIndicator(uint larvaID)
{
    QVector<double> res;
    size_t i;
    bool val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getRightBendingIndicatorAt(t, val))
            {
                res << static_cast<double>(val);
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllMovementDirection(uint larvaID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getMovementDirectionAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllDistancesToLandmark(uint larvaID, const QString &landmarkID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getDistanceToLandmark(t, landmarkID.toStdString(), val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAllBearingAnglesToLandmark(uint larvaID, const QString &landmarkID)
{
    QVector<double> res;
    size_t i;
    double val;
    if(this->getIndexOfLarva(larvaID, i))
    {
        for(uint t : this->mLarvae.at(i).getAllTimeSteps())
        {
            if(this->mLarvae.at(i).getBearingAngleToLandmark(t, landmarkID.toStdString(), val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

QVector<double> LarvaeContainer::getVelocity(uint larvaID)
{
    QVector<double> x;
    QVector<double> y;
    
    QVector<double> xReal;
    QVector<double> yReal;
    
    QVector<double> res;
    size_t larvaIndex;
    double val;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        
        /*
        std::vector<uint> time = this->mLarvae.at(larvaIndex).getAllTimeSteps();
        
        for(uint t : time)
        {
            if(this->mLarvae.at(larvaIndex).getMomentumDistAt(t, val))
            {
                xReal << static_cast<double>(t);
                yReal << val;
            }
        }
        
        for(int i = 0; i < time.size(); i += CameraParameter::dFPS)
        {
            uint t = time.at(i);
            if(this->mLarvae.at(larvaIndex).getMomentumDistAt(t, val))
            {
                x << static_cast<double>(t);
                y << val;
            }
        }
        int n = x.size()-1;
        
        QVector<double> a = y;
        QVector<double> b = QVector<double>(n, 0.0);
        QVector<double> c;
        QVector<double> d = QVector<double>(n, 0.0);
        QVector<double> h;
        QVector<double> e = QVector<double>(n, 0.0);;
        QVector<double> matrixA, matrixB, matrixC, matrixD;
        
        for(int k = 0; k <= n-1; ++k)
            h << (x[k+1] - x[k]);
            
        for(int k = 1; k <= n-1; ++k)
            e[k] = 3*((a[k+1]-a[k])/h[k] - (a[k]-a[k-1])/h[k-1]); 
            
        matrixA << 0.0;
        for(int k = 1; k < n-1; ++k)
            matrixA << h[k];
            
        for(int i = 0; i < n-1; ++i)
            matrixB << (2*(h[i]+h[i+1]));
            
        for(int i = 1; i < n-1; ++i)
            matrixC << h[i];
        matrixC << 0.0;
        
        for(int i = 0; i < n-1; ++i)
            matrixD << ( 3*(a[i+2]-a[i+1])/h[i+1] - 3*(a[i+1]-a[i])/h[i] );
            
        matrixC[0] /= matrixB[0];
        matrixD[0] /= matrixB[0];
        
        for(int i = 1; i < n-1; ++i)
        {
            matrixC[i] = (matrixC[i] / (matrixB[i]-matrixA[i]*matrixC[i-1]));
            matrixD[i] = ((matrixD[i] - matrixA[i]*matrixD[i-1])/(matrixB[i]-matrixA[i]*matrixC[i-1]));
        }
        
        matrixD[n-1] = ((matrixD[n-1] - matrixA[n-1]*matrixD[n-1-1])/(matrixB[n-1]-matrixA[n-1]*matrixC[n-1-1]));
        

        
        QVector<double> c = QVector<double>(n+1, 0.0);
        QVector<double> l = QVector<double>(n+1, 0.0);
        QVector<double> mu = QVector<double>(n+1, 0.0);
        QVector<double> z = QVector<double>(n+1, 0.0);
        
        l[0] = 1.0;
        mu[0] = z[0] = 0.0;
        
        for(int i = 1; i < n; ++i)
        {
            l[i]    = (2*(x[i+1]-x[i-1])) - (h[i-1]*mu[i-1]);
            mu[i]   = h[i] / l[i];
            z[i]    = (r[i] - (h[i-1]*z[i-1]))/l[i];
        }
        
        l[n] = 1.0;
        z[n] = c[n] = 0.0;
        
        for(int j = n-1; j >= 0; --j)
        {
            c[j] = z[j] - (mu[j]*c[j+1]);
            b[j] = ((a[j+1]-a[j])/h[j]) - ((h[j]*(c[j+1]+(2*c[j])))/3);
            d[j] = (c[j+1] - c[j]) / (3*h[j]);
        }

        QVector<SplineSet> outputSet;
        for(int i = 0; i < n; ++i)
        {
            SplineSet spSet;
            spSet.a = a[i];
            spSet.b = b[i];
            spSet.c = c[i];
            spSet.d = d[i];
            spSet.x = x[i];
            
            for(int j = i*CameraParameter::dFPS; j < i*CameraParameter::dFPS + CameraParameter::dFPS; ++j)
            {
                double velosity         = spSet.b + 2*spSet.c*(xReal[j] - spSet.x) + 3*spSet.d*(xReal[j] - spSet.x);
                double acceleration     = 2*spSet.c + 6*spSet.d*(xReal[j] - spSet.x);
                
                QString debugString = QString("Time: %1,\t Accumulated Distance: %2,\t Velosity: %3,\t Acceleration: %4")
                        .arg(xReal[j]).arg(yReal[j]).arg(velosity).arg(acceleration);
                        
                qDebug() << debugString;
                
                res << velosity;
            }
            
            
            //            qDebug() << debugString;
            
            //            outputSet << spSet;
            
            //            res << spSet.b;
        }*/
        
        for(uint t : this->mLarvae.at(larvaIndex).getAllTimeSteps())
        {
            if(this->mLarvae.at(larvaIndex).getVelosityAt(t, val))
            {
                res << val;
            }
        }
        
        
        
        /*
        QVector<QVector<double> > timeGaps = this->getAllTimestepGaps(larvaID);
        
        for(QVector<double> const& time : timeGaps)
        {
            if(time.size() <= frameWindowSize)
            {
                for(int t = 0; t < time.size(); ++i)
                    res << 0.0;
            }
            else
            {

                // Forwardpass
                for(int t = 0; t < time.size() - frameWindowSize; ++t)
                {
                    int t0 = time.at(t);
                    int t1 = time.at(t+frameWindowSize);
                    
                    double s0;
                    double s1;
                    
                    if(this->mLarvae.at(i).getAccDistAt(t0, s0) && this->mLarvae.at(i).getAccDistAt(t1, s1))
                    {
                        val = (s1 - s0) / (static_cast<double>(t1) - static_cast<double>(t0));
                        val *= CameraParameter::dScaleFactor; // cm per second
                        val *= 10.0; // mm per second
                        res << val;
                    }
                }
                
                // Backwardpass
                for(int t = time.size()-1; t >= frameWindowSize; --t)
                {
                    int t0 = time.at(t);
                    int t1 = time.at(t-frameWindowSize);
                    
                    double s0;
                    double s1;
                    
                    if(this->mLarvae.at(i).getAccDistAt(t0, s0) && this->mLarvae.at(i).getAccDistAt(t1, s1))
                    {
                        val = (s1 - s0) / (static_cast<double>(t1) - static_cast<double>(t0));
                        val *= CameraParameter::dScaleFactor; // cm per second
                        val *= 10.0; // mm per second
                        res << val;
                    }
                }
            }
        }*/
    }
    
    return res;
}

QVector<double> LarvaeContainer::getAcceleration(uint larvaID)
{
    QVector<double> res;
    size_t larvaIndex;
    double val;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        for(uint t : this->mLarvae.at(larvaIndex).getAllTimeSteps())
        {
            if(this->mLarvae.at(larvaIndex).getAccelerationAt(t, val))
            {
                res << val;
            }
        }
    }
    
    return res;
}

void LarvaeContainer::updateLarvaValues(TrackerSceneLarva const* tLarva)
{
    int time = tLarva->getCurrentTimePoint();    
    size_t larvaIndex;
    
    if(this->getIndexOfLarva(tLarva->getBaseID(), larvaIndex))
    {      
        this->updateLarvaSpine(                 larvaIndex,         time,   tLarva->getSpine(), tLarva->getSpineRadii());
        this->updateLarvaMomentum(              larvaIndex,         time,   tLarva->getPoligon());
        this->updateLarvaArea(                  larvaIndex,         time,   tLarva->getPoligon());
        this->updateLarvaPerimeter(             larvaIndex,         time,   tLarva->getPoligon());
        
        this->updateIsCoiledIndicator(          larvaIndex,         time);
        this->updateGoPhaseIndicator(           larvaIndex,         time);
        this->updateTurnIndicator(              larvaIndex,         time);
        this->updateMovementDirection(          larvaIndex,         time);
        this->updateLarvaAccumulatedDistance(   larvaIndex);
        this->updateLarvaDistance2Origin(       larvaIndex);
        
        if(dynamic_cast<TrackerScene*>(tLarva->getScene())->hasLandmarkContainer())
        {
            LandmarkContainer const* landmarkContainer = dynamic_cast<TrackerScene*>(tLarva->getScene())->getLandmarkContainer();
            for(int i = 0; i < landmarkContainer->getSize(); ++i)
            {
                this->updateLandmark(landmarkContainer->getLandmarkAt(i));
            }
            
        }
        this->recalculateLarvaVelocityAndAcceleration(larvaIndex);
        emit sendUpdatedResultLarvaID(larvaIndex);
    }
}

spineType LarvaeContainer::calcSpine(QPainterPath const& spinePath)
{
    spineType spine;
    QPainterPath::Element e;
    
    for(int i = 0; i < spinePath.elementCount(); ++i) 
    {
        e = spinePath.elementAt(i);
        spine.push_back(cv::Point(e.x, e.y));
    }
    
    return spine;
}

void LarvaeContainer::updateLarvaSpine(int index, 
                                       uint time, 
                                       QPainterPath const& paintSpine, 
                                       std::vector<double> const& radii)
{
    spineType spine                                                         = this->calcSpine(paintSpine);
    double spineLength                                                      = Calc::calcSpineLength(spine);
    double mainBodyBendingAngle                                             = Calc::calcAngle(spine.at((spine.size() - 1) / 2), 
                                                                                              spine.at(0), 
                                                                                              spine.at((spine.size() - 1)));    
    this->mLarvae[index].parameters[time].spine                       = spine;
    this->mLarvae[index].parameters[time].spineLength                 = spineLength;
    this->mLarvae[index].parameters[time].mainBodyBendingAngle        = mainBodyBendingAngle;
    this->mLarvae[index].parameters[time].spineRadii                  = radii;
}

void LarvaeContainer::updateLarvaMomentum(int index, 
                                          uint time, 
                                          QPolygonF const& paintPolygon)
{
    cv::Point consecutiveMomentum;
    cv::Point momentum      = Calc::calcPolygonCenterOfMass(paintPolygon);
    double momentumDist     = 0.0;
    
    if(this->mLarvae[index].getMomentumAt(time-1, consecutiveMomentum))
    {
        momentumDist = Calc::eucledianDist(momentum, consecutiveMomentum);
    }
    this->mLarvae[index].parameters[time].momentum        = momentum;
    this->mLarvae[index].parameters[time].momentumDist    = momentumDist;
}

void LarvaeContainer::updateLarvaArea(int index, 
                                      uint time, 
                                      QPolygonF const& paintPolygon)
{
    double area                                         = Calc::calcPolygonArea(paintPolygon);
    this->mLarvae[index].parameters[time].area    = area;
}

void LarvaeContainer::updateLarvaPerimeter(int index, 
                                           uint time, 
                                           QPolygonF const& paintPolygon)
{
    double perimeter                                        = Calc::calcPerimeter(paintPolygon);
    this->mLarvae[index].parameters[time].perimeter   = perimeter;
}

void LarvaeContainer::updateLarvaDistance2Origin(int index)
{
    cv::Point momentum;
    double distToOrigin = 0.0;
    std::vector<unsigned int> timeSteps = this->mLarvae.at(index).getAllTimeSteps();
    
    for(unsigned int i = 0; i < timeSteps.size(); ++i) 
    {
        if(this->mLarvae.at(index).getMomentumAt(i, momentum))
        {
            distToOrigin = Calc::eucledianDist(this->mLarvae.at(index).getOrigin(), momentum);
            this->mLarvae[index].parameters[timeSteps.at(i)].distToOrigin = distToOrigin;
        }
    }
}

void LarvaeContainer::updateLarvaAccumulatedDistance(int index)
{
    cv::Point p1, p2;
    double accDist = 0.0;
    std::vector<unsigned int> timeSteps = this->mLarvae.at(index).getAllTimeSteps();
    for(unsigned int i = 1; i < timeSteps.size(); ++i) 
    {
        if(this->mLarvae.at(index).getMomentumAt(i-1, p1) && this->mLarvae.at(index).getMomentumAt(i, p2))
        {
            accDist += Calc::eucledianDist(p1, p2);
            this->mLarvae[index].parameters[timeSteps.at(i)].accDist = accDist;
        }
    }
    
}

void LarvaeContainer::updateIsCoiledIndicator(int index, uint time)
{
    double perimeter;
    double spineLength;
    double peri2spineLengthRatio;
    
    uint midPointIndex;
    double midPointRadius;
    double midCirclePeri2PeriRatio;
    
    int nPoints;
    
    std::vector<double> radii;
    
    if(     this->mLarvae.at(index).getPerimeterAt(time, perimeter) && 
            this->mLarvae.at(index).getSpineLengthAt(time, spineLength) &&
            this->mLarvae.at(index).getSpineRadiiAt(time, radii))
    {
        double peri2spineLengthThresh = 2.6;
        double midCirclePeri2PeriThresh = 0.5;
        
        if(!LarvaeExtractionParameters::bUseDefault)
        {
            nPoints = LarvaeExtractionParameters::iNumerOfSpinePoints;
            
            peri2spineLengthThresh = LarvaeExtractionParameters::CoiledRecognitionParameters::dPerimeterToSpinelenghtThreshold;
            midCirclePeri2PeriThresh = LarvaeExtractionParameters::CoiledRecognitionParameters::dMidcirclePerimeterToPerimeterThreshold;
            
        }
        
        midPointIndex = this->mLarvae.at(index).getSpineMidPointIndex();
        peri2spineLengthRatio = perimeter / spineLength;
        
        if(midPointIndex < radii.size())
        {
            midPointRadius = radii.at(midPointIndex);
            midCirclePeri2PeriRatio = (2*midPointRadius*CV_PI) / perimeter;
            
            if( peri2spineLengthRatio >= peri2spineLengthThresh &&
                    midCirclePeri2PeriRatio >= midCirclePeri2PeriThresh)
            {
                this->mLarvae[index].parameters[time].isCoiled = true;
            }
            else
            {
                this->mLarvae[index].parameters[time].isCoiled = false;
            }
        }
    }
}

void LarvaeContainer::updateGoPhaseIndicator(int index, uint time)
{
    // body bending must be +/- 30 degree deviation from 180 degree
    int angleThresh = 30;
    // traveled distance is calculated between one second (i.e. framesForSpeedCalc = FPS)
    int framesForSpeedCalc = static_cast<int>(CameraParameter::dFPS);
    
    double speedThreshTmp = static_cast<double>(GeneralParameters::iMaxLarvaeArea);
    // larval length is approximatly 2*sqrt(max_area)
    speedThreshTmp = 2 * std::sqrt(speedThreshTmp);
    // assuming a framesForSpeedCalc value equal to the frame rate: speed is calculated for one second
    // a further assumption, that a crawling larva is travaleing 5% of its body length per second (during a go) leads to:
    int speedThresh = static_cast<int>((speedThreshTmp * 0.05));
    
    if(!LarvaeExtractionParameters::StopAndGoCalculation::bUseDynamicStopAndGoParameterCalculation)
    {
        angleThresh = LarvaeExtractionParameters::StopAndGoCalculation::iAngleThreshold;
        framesForSpeedCalc = LarvaeExtractionParameters::StopAndGoCalculation::iFramesForSpeedCalculation;
        speedThresh = LarvaeExtractionParameters::StopAndGoCalculation::iSpeedThreshold;
    }
    
    cv::Point pastMomentum;
    cv::Point curMomentum;
    int phaseIndicator = -1;
    double curBending;
    
    if(
            this->mLarvae.at(index).getMomentumAt(time - framesForSpeedCalc, pastMomentum) && 
            this->mLarvae.at(index).getMomentumAt(time, curMomentum) && 
            this->mLarvae.at(index).getMainBodyBendingAngleAt(time, curBending))
    {
        double angleLowThresh = 180 - angleThresh;
        double angleHighThresh = 180 + angleThresh;
        double distance = Calc::eucledianDist(curMomentum, pastMomentum);
        
        if(distance >= speedThresh && curBending >= angleLowThresh && curBending <= angleHighThresh)
        {
            phaseIndicator = 1;
        }
        else
        {
            phaseIndicator = 0;
        }
    }
    this->mLarvae[index].parameters[time].goPhase =  phaseIndicator;
}

void LarvaeContainer::updateTurnIndicator(int index, uint time)
{
    int bendingAngleThresh = 30;
    double mainBodyBendingAngle;
    if(!LarvaeExtractionParameters::BodyBendingParameters::bUseDynamicBodyBendingCalculation)
    {
        bendingAngleThresh = static_cast<int>(LarvaeExtractionParameters::BodyBendingParameters::dAngleThreshold);
    }
    
    if(this->mLarvae.at(index).getMainBodyBendingAngleAt(time, mainBodyBendingAngle))
    {
        this->mLarvae[index].parameters[time].isLeftBended    = this->calcLeftTurnIndicator(mainBodyBendingAngle, bendingAngleThresh);
        this->mLarvae[index].parameters[time].isRightBended   = this->calcRightTurnIndicator(mainBodyBendingAngle, bendingAngleThresh);
    }
}

bool LarvaeContainer::calcLeftTurnIndicator(const double curBending, const uint bendingThresh)
{
    return (curBending > (180+bendingThresh));
}

bool LarvaeContainer::calcRightTurnIndicator(const double curBending, const uint bendingThresh)
{
    return (curBending < (180-bendingThresh));
}

void LarvaeContainer::updateMovementDirection(int index, uint time)
{
    int framesForMovementDirectionCalc = static_cast<int>(CameraParameter::dFPS);
    if(!LarvaeExtractionParameters::MovementDirectionParameters::bUseDynamicMovementDirectionParameterCalculation)
    {
        framesForMovementDirectionCalc = LarvaeExtractionParameters::MovementDirectionParameters::iFramesForMovementDirectionCalculation;
    }
    double movementDirection = -1;
    cv::Point prevMomentum;
    cv::Point curMomentum;
    if(
            this->mLarvae.at(index).getMomentumAt(time - framesForMovementDirectionCalc, prevMomentum) &&
            this->mLarvae.at(index).getMomentumAt(time, curMomentum))
    {
        // if there is no movement at all angle calculation will fail.
        // Thus, we asure that the distance between the two points is not zero!
        if((Calc::normL2(curMomentum - prevMomentum)) != 0)
        {
            movementDirection = Calc::calcAngleToYAxes(prevMomentum,curMomentum);
        }
    }
    
    this->mLarvae[index].parameters[time].movementDirection = movementDirection;
}

bool LarvaeContainer::getIndexOfLarva(uint id, size_t& index) const
{
    for(unsigned int i = 0; i < this->mLarvae.size(); ++i)
    {
        if(this->mLarvae.at(i).getID() == id)
        {
            index = i;
            return true;
        }
    }
    return false;
}

void LarvaeContainer::insertRawLarva(uint larvaID, 
                                     uint timePoint, 
                                     const RawLarva &rawLarva)
{
    size_t larvaIndex;
    if(this->getIndexOfLarva(larvaID, larvaIndex))
    {
        spineType spine;
        spine.reserve(rawLarva.getDiscreteSpine().size());
        
        std::vector<double> radii;
        radii.reserve(rawLarva.getLarvalRadii().size());
        
        bool invertParameters = this->changeDirectionality(larvaIndex, timePoint,rawLarva);
        
        spine = (invertParameters) ? this->reverseVec(rawLarva.getDiscreteSpine())  : rawLarva.getDiscreteSpine();
        radii = (invertParameters) ? this->reverseVec(rawLarva.getLarvalRadii())    : rawLarva.getLarvalRadii();
        
        this->mLarvae[larvaIndex].values.spine          = spine;
        this->mLarvae[larvaIndex].values.momentum       = rawLarva.getMomentum();
        this->mLarvae[larvaIndex].values.area           = rawLarva.getArea();
        this->mLarvae[larvaIndex].values.spineRadii     = radii;
        
        int midIndex = static_cast<int>((this->mLarvae[larvaIndex].getNSpinePoints() - 1) / 2);
        this->mLarvae[larvaIndex].values.mainBodyBendingAngle = Calc::calcAngle(spine.at(midIndex), spine.at(0), spine.at((this->mLarvae[larvaIndex].getNSpinePoints()-1)));
        
        this->mLarvae[larvaIndex].values.spineLength    = rawLarva.getSpineLength();
        this->mMaxSpineLength                           = std::max(this->mLarvae[larvaIndex].values.spineLength, this->mMaxSpineLength);
        this->mLarvae[larvaIndex].values.perimeter      = rawLarva.getContourPerimeter();
        
        this->mLarvae[larvaIndex]. values.isCoiled      = rawLarva.getIsCoiledIndicator();
        
        this->mLarvae[larvaIndex].values.distToOrigin   = calcDistToOrigin(larvaIndex, rawLarva.getMomentum());
        
        this->mLarvae[larvaIndex].values.momentumDist   = calcMomentumDist(larvaIndex, timePoint, rawLarva.getMomentum());
        this->mLarvae[larvaIndex].values.accDist        += this->mLarvae[larvaIndex].values.momentumDist;
        
        // body bending must be +/- 30 degree deviation from 180 degree
        int angleThresh = 30;
        // traveled distance is calculated between one second (i.e. framesForSpeedCalc = FPS)
        int framesForSpeedCalc = static_cast<int>(CameraParameter::dFPS);
        // movement (in pixel) must be more than 15% of the spine length
        //int speedThresh = (int) (rawLarva.getSpineLength() * 15 / 100);
        
        double speedThreshTmp = static_cast<double>(GeneralParameters::iMaxLarvaeArea);
        // larval length is approximatly 2*sqrt(max_area)
        speedThreshTmp = 2 * std::sqrt(speedThreshTmp);
        // assuming a framesForSpeedCalc value equal to the frame rate: speed is calculated for one second
        // a further assumption, that a crawling larva is travaleing 5% of its body length per second (during a go) leads to:
        int speedThresh = static_cast<int>(speedThreshTmp * 0.05);
        
        if(!LarvaeExtractionParameters::StopAndGoCalculation::bUseDynamicStopAndGoParameterCalculation)
        {
            angleThresh = LarvaeExtractionParameters::StopAndGoCalculation::iAngleThreshold;
            framesForSpeedCalc = LarvaeExtractionParameters::StopAndGoCalculation::iFramesForSpeedCalculation;
            speedThresh = LarvaeExtractionParameters::StopAndGoCalculation::iSpeedThreshold;
        }
        
        this->mLarvae[larvaIndex].values.goPhase = calcGoPhaseIndicator(larvaIndex,
                                                                        timePoint,
                                                                        rawLarva.getMomentum(),
                                                                        this->mLarvae[larvaIndex].values.mainBodyBendingAngle,
                                                                        speedThresh,
                                                                        angleThresh,
                                                                        framesForSpeedCalc);
        
        int bendingAngleThresh = 30;
        if(!LarvaeExtractionParameters::BodyBendingParameters::bUseDynamicBodyBendingCalculation)
        {
            bendingAngleThresh = static_cast<int>(LarvaeExtractionParameters::BodyBendingParameters::dAngleThreshold);
        }
        
        this->mLarvae[larvaIndex].values.isLeftBended  = this->calcLeftTurnIndicator(this->mLarvae[larvaIndex].values.mainBodyBendingAngle, bendingAngleThresh);
        this->mLarvae[larvaIndex].values.isRightBended = this->calcRightTurnIndicator(this->mLarvae[larvaIndex].values.mainBodyBendingAngle, bendingAngleThresh);
        
        int framesForMovementDirectionCalc = static_cast<int>(CameraParameter::dFPS);
        if(!LarvaeExtractionParameters::MovementDirectionParameters::bUseDynamicMovementDirectionParameterCalculation)
        {
            framesForMovementDirectionCalc = LarvaeExtractionParameters::MovementDirectionParameters::iFramesForMovementDirectionCalculation;
        }
        
        this->mLarvae[larvaIndex].values.movementDirection = calcMovementDirection(larvaIndex, timePoint,framesForMovementDirectionCalc,rawLarva.getMomentum());
        
        this->recalculateLarvaVelocityAndAcceleration(larvaIndex);
        
        this->mLarvae[larvaIndex].parameters.insert(std::pair<unsigned int, Larva::ValuesType>(timePoint, this->mLarvae[larvaIndex].values));
    }
}

bool LarvaeContainer::isAssignedAt(uint larvaID, uint timePoint) const
{
    size_t index;
    if(this->getIndexOfLarva(larvaID, index))
    {
        return (this->mLarvae.at(index).parameters.count(timePoint) > 0);
    }
    else
    {
        return false;
    }
}

void LarvaeContainer::interplolateLarvae()
{
    // change Head-Tail position if it is not consistent over time
    this->interpolateHeadTailOverTime();
    // fill sampling gaps caused by time windows (e.g. 10 fps are oversampled for movement direction etc.)
    this->fillTimeSamplingGaps();
}

void LarvaeContainer::interpolateHeadTailOverTime(uint larvaIndex)
{
    int changeDirectionalityCounter = 0;
    for (std::map<unsigned int, Larva::ValuesType>::iterator it = this->mLarvae[larvaIndex].parameters.begin(); it != this->mLarvae[larvaIndex].parameters.end(); ++it)
    {
        double movementDirection;
        if (this->mLarvae.at(larvaIndex).getMovementDirectionAt(it->first,movementDirection))
        {
            cv::Point head = it->second.spine.at(0);
            cv::Point tail = it->second.spine.at(this->mLarvae[larvaIndex].getNSpinePoints()-1);
            double angleMainBodyDirectionality = Calc::calcAngleToYAxes(tail,head);
            
            double angleDiff = Calc::calcAngleDiff(movementDirection,angleMainBodyDirectionality);
            
            if (angleDiff > 50)
            {
                ++changeDirectionalityCounter;
            }
            
        }
    }
    
    if (changeDirectionalityCounter >= static_cast<int>(this->mLarvae.at(larvaIndex).parameters.size() / 2))
    {
        for (std::map<unsigned int, Larva::ValuesType>::iterator it = this->mLarvae[larvaIndex].parameters.begin(); it != this->mLarvae[larvaIndex].parameters.end(); ++it)
        {
            it->second.spine        = this->reverseVec(it->second.spine);
            it->second.spineRadii   = this->reverseVec(it->second.spineRadii);
            
            it->second.mainBodyBendingAngle = Calc::calcAngle(it->second.spine.at((this->mLarvae[larvaIndex].getNSpinePoints()-1)/2),
                                                              it->second.spine.at(0),
                                                              it->second.spine.at(this->mLarvae[larvaIndex].getNSpinePoints()-1));
            
            //it->second.movementDirection = Calc::calcCircularAngleSum(it->second.movementDirection,180);
            
            if (it->second.isLeftBended)
            {
                it->second.isLeftBended = false;
                it->second.isRightBended = true;
            }
            else if (it->second.isRightBended)
            {
                it->second.isRightBended = false;
                it->second.isLeftBended = true;
            }
        }
    }
}

void LarvaeContainer::interpolateHeadTailOverTime()
{
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        this->interpolateHeadTailOverTime(i);
    }
}

void LarvaeContainer::fillTimeSamplingGaps(uint larvaIndex)
{
    // Get Stop And Go Parameteres:
    // body bending must be +/- 30 degree deviation from 180 degree
    int angleThresh = 30;
    
    // traveled distance is calculated between one second (i.e. framesForSpeedCalc = FPS)
    int framesForSpeedCalc = static_cast<int>(CameraParameter::dFPS);
    
    // calc speed threshold
    double speedThreshTmp = static_cast<double>(GeneralParameters::iMaxLarvaeArea);
    
    // larval length is approximatly 2*sqrt(max_area)
    speedThreshTmp = 2 * std::sqrt(speedThreshTmp);
    
    // assuming a framesForSpeedCalc value equal to the frame rate: speed is calculated for one second
    // a further assumption, that a crawling larva is travaleing 10% of its body length per second (during a go) leads to:
    int speedThresh = static_cast<int>(speedThreshTmp * 0.1);
    
    if(!LarvaeExtractionParameters::StopAndGoCalculation::bUseDynamicStopAndGoParameterCalculation)
    {
        angleThresh = LarvaeExtractionParameters::StopAndGoCalculation::iAngleThreshold;
        framesForSpeedCalc = LarvaeExtractionParameters::StopAndGoCalculation::iFramesForSpeedCalculation;
        speedThresh = LarvaeExtractionParameters::StopAndGoCalculation::iSpeedThreshold;
    }
    
    // Get Movement Direction Parameters:
    int framesForMovementDirectionCalc = static_cast<int>(CameraParameter::dFPS);
    if(!LarvaeExtractionParameters::MovementDirectionParameters::bUseDynamicMovementDirectionParameterCalculation)
    {
        framesForMovementDirectionCalc = LarvaeExtractionParameters::MovementDirectionParameters::iFramesForMovementDirectionCalculation;
    }
    
    for (std::map<unsigned int, Larva::ValuesType>::iterator it = this->mLarvae[larvaIndex].parameters.begin(); it != this->mLarvae[larvaIndex].parameters.end(); ++it)
    {
        int goPhase = it->second.goPhase;
        int movementDirection = it->second.movementDirection;
        
        if(goPhase == -1)
        {
            cv::Point curMom = it->second.momentum;
            cv::Point nextMom;
            if (this->mLarvae.at(larvaIndex).getMomentumAt(it->first + framesForSpeedCalc, nextMom))
            {
                int phaseIndicator = -1;
                
                double curBending = it->second.mainBodyBendingAngle;
                double angleLowThresh = 180 - angleThresh;
                double angleHighThresh = 180 + angleThresh;
                double distance = Calc::eucledianDist(curMom,nextMom);
                
                if(distance >= speedThresh && curBending >= angleLowThresh && curBending <= angleHighThresh)
                {
                    phaseIndicator = 1;
                }
                else
                {
                    phaseIndicator = 0;
                }
                it->second.goPhase = phaseIndicator;
            }
        }
        
        if(movementDirection == -1.0)
        {
            double movementDirection = -1;
            cv::Point curMom = it->second.momentum;
            cv::Point nextMom;
            if (this->mLarvae.at(larvaIndex).getMomentumAt(it->first + framesForMovementDirectionCalc, nextMom))
            {
                // if there is no movement at all angle calculation will fail.
                // Thus, we asure that the distance between the two points is not zero!
                if((Calc::normL2(curMom - nextMom)) != 0)
                {
                    movementDirection = Calc::calcAngleToYAxes(curMom,nextMom);
                }
            }
            
            it->second.movementDirection = movementDirection;
        }
    }
}

void LarvaeContainer::fillTimeSamplingGaps()
{
    for(size_t i = 0; i < this->mLarvae.size(); ++i)
    {
        this->fillTimeSamplingGaps(i);
    }
}

double LarvaeContainer::calcMomentumDist(uint larvaIndex, uint timePoint, const cv::Point &curMomentum) const
{
    cv::Point consecutiveMomentum;
    double dist = 0;
    
    if(this->mLarvae.at(larvaIndex).getMomentumAt(timePoint-1, consecutiveMomentum))
    {
        dist = Calc::eucledianDist(curMomentum, consecutiveMomentum);
    }
    
    return dist;
}

void LarvaeContainer::changeDirection(uint larvaIndex, uint timePoint)
{
    spineType spine;
    if(this->mLarvae.at(larvaIndex).getSpineAt(timePoint, spine))
    {
        // spinesize muss be an odd number
        if(spine.size() % 2 == 0)
        {
            return;
        }
        
        cv::Point p1, p2;
        
        for(unsigned int i = 0; i < spine.size() / 2; ++i)
        {
            p1 = spine[i];
            p2 = spine[spine.size() - i - 1];
            
            spine[i] = p2;
            spine[spine.size() - i - 1] = p1;
        }
        
        this->mLarvae.at(larvaIndex).setSpineAt(timePoint, spine);
    }
}

void LarvaeContainer::eraseAt(uint larvaIndex, uint timePoint)
{
    this->mLarvae[larvaIndex].parameters.erase(timePoint);
}

double LarvaeContainer::calcDistToOrigin(uint larvaIndex, const cv::Point &curMomentum) const
{
    return Calc::eucledianDist(this->mLarvae.at(larvaIndex).getOrigin(), curMomentum);
}

int LarvaeContainer::calcGoPhaseIndicator(uint larvaIndex, 
                                          const uint timePoint, 
                                          const cv::Point &curMomentum, 
                                          const double curBending, 
                                          const uint velocityThresh, 
                                          const uint bendingThresh, 
                                          const uint timeWindow)
{
    cv::Point pastMomentum;
    int phaseIndicator = -1;
    
    if(this->mLarvae.at(larvaIndex).getMomentumAt(timePoint - timeWindow, pastMomentum))
    {
        double angleLowThresh = 180 - bendingThresh;
        double angleHighThresh = 180 + bendingThresh;
        double distance = Calc::eucledianDist(curMomentum,pastMomentum);
        
        if(distance >= velocityThresh && curBending >= angleLowThresh && curBending <= angleHighThresh)
        {
            phaseIndicator = 1;
        }
        else
        {
            phaseIndicator = 0;
        }
    }
    return phaseIndicator;
}

double LarvaeContainer::calcMovementDirection(uint larvaIndex, 
                                              const uint timePoint, 
                                              const uint timeWindow, 
                                              const cv::Point &curMomentum)
{
    double movementDirection = -1;
    cv::Point prevMomentum;
    if(this->mLarvae.at(larvaIndex).getMomentumAt(timePoint-timeWindow, prevMomentum))
    {
        // if there is no movement at all angle calculation will fail.
        // Thus, we asure that the distance between the two points is not zero!
        if((Calc::normL2(curMomentum - prevMomentum)) != 0)
        {
            movementDirection = Calc::calcAngleToYAxes(prevMomentum,curMomentum);
        }
    }
    
    return movementDirection;
}

bool LarvaeContainer::changeDirectionality(uint larvaIndex, uint timePoint, const RawLarva &rawlarva)
{
    spineType previousSpine;
    spineType curSpine = rawlarva.getDiscreteSpine();
    bool change = false;
    if(timePoint > 0 && this->mLarvae.at(larvaIndex).getSpineAt(timePoint-1,previousSpine))
    {
        cv::Point prevHead  = previousSpine.at(0);
        cv::Point prevTail  = previousSpine.at(previousSpine.size()-1);
        
        cv::Point curHead   = curSpine.at(0);
        
        double dH2H = Calc::eucledianDist(prevHead,curHead);
        double dT2H = Calc::eucledianDist(prevTail,curHead);
        if(dT2H < dH2H)
        {
            change = true;
        }
    }
    return change;
}

template<class T> std::vector<T> LarvaeContainer::reverseVec(std::vector<T> const& v)
{
    std::vector<T> retVec;
    retVec.reserve(v.size());
    std::reverse_copy(v.begin(), v.end(), std::back_inserter(retVec));
    
    return retVec;
}

void LarvaeContainer::createNewLarva(uint timePoint, const RawLarva &rawLarva, unsigned int larvaID)
{
    Larva l;
    l.setNSpinePoints(rawLarva.getDiscreteSpine().size());
    
    l.setOrigin(rawLarva.getMomentum());
    
    l.values.accDist = 0.0;
    l.setID(larvaID);
    
    this->mLarvae.push_back(l);
    
    this->insertRawLarva(l.getID(), timePoint, rawLarva);
}

uint LarvaeContainer::getLastValidLavaID() const
{
    if(this->mLarvae.empty())
        return 0;
    
    uint res = this->mLarvae.front().getID();
    
    for(size_t i = 1; i < this->mLarvae.size(); ++i)
    {
        res = qMax(res, this->mLarvae.at(i).getID());
    }
    
    return res;
}
