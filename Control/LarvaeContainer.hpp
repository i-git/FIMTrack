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

#ifndef LARVAECONTAINER_HPP
#define LARVAECONTAINER_HPP

#include "Configuration/FIMTrack.hpp"

//#include <QtCore>
//#include <QtGui>

#include "Data/Larva.hpp"
#include "InputGenerator.hpp"
#include "OutputGenerator.hpp"
#include "GUI/TrackerScene.hpp"
#include "GUI/TrackerSceneLarva.hpp"
#include "GUI/LandmarkContainer.hpp"

class LarvaeContainer : public QObject
{
    Q_OBJECT
    
private:
    std::vector < Larva >                       mLarvae;
    double                                      mMaxSpineLength;
    
    uint getLastValidLavaID() const;
    
    double calcMomentumDist(uint larvaIndex, 
                            uint timePoint, 
                            cv::Point const & curMomentum) const;
    
    void changeDirection(uint larvaIndex, uint timePoint);
    void eraseAt(uint larvaIndex, uint timePoint);
    
    double calcDistToOrigin(uint larvaIndex, 
                            cv::Point const & curMomentum) const;

    int calcGoPhaseIndicator(uint larvaIndex,
                             uint const timePoint,
                             cv::Point const & curMomentum,
                             double const curBending,
                             uint const velocityThresh,
                             uint const bendingThresh,
                             uint const timeWindow);
    double calcMovementDirection(uint larvaIndex, 
                                 uint const timePoint,
                                 uint const timeWindow,
                                 cv::Point const & curMomentum);

    bool changeDirectionality(uint larvaIndex, uint timePoint, RawLarva const & rawlarva);
    template<class T> std::vector<T> reverseVec(std::vector<T> const & v);
    
    
    bool getIndexOfLarva(uint id, size_t &index) const;
    spineType calcSpine(QPainterPath const& spinePath);
    
    void updateLarvaSpine(int index, uint time, QPainterPath const& paintSpine, const std::vector<double> &radii);
    void updateLarvaMomentum(int index, uint time, QPolygonF const& paintPolygon);
    void updateLarvaArea(int index, uint time, QPolygonF const& paintPolygon);
    void updateLarvaPerimeter(int index, uint time, QPolygonF const& paintPolygon);
    void updateLarvaDistance2Origin(int index);
    void updateLarvaAccumulatedDistance(int index);
    void updateIsCoiledIndicator(int index, uint time);
    
    void updateGoPhaseIndicator(int index, uint time);
    void updateTurnIndicator(int index, uint time);
    bool calcLeftTurnIndicator(const double curBending, const uint bendingThresh);
    bool calcRightTurnIndicator(const double curBending, const uint bendingThresh);
    void updateMovementDirection(int index, uint time);
    
    void recalculateLarvaDistanceParameter(uint larvaID);
    void recalculateLarvaDistanceToOrigin(size_t larvaIndex);
    void recalculateLarvaMomentumDistance(size_t larvaIndex);
    void recalculateLarvaVelocityAndAcceleration(size_t larvaIndex);
    void recalculateLarvaAcceleration(size_t larvaIndex);
    
    void calcLandmarkParameter(QString const& name, QPointF const& p);
    void calcLandmarkParameter(QString const& name, QLineF const& l);
    void calcLandmarkParameter(QString const& name, QRectF const& r, const bool ellipse);
    void calcBearinAngle(QString const& name, QPointF const& landmarkPoint);
    
public:
    explicit LarvaeContainer(QObject *parent = 0);
    
    Larva* createDefaultLarva(uint timeStep);
    
    /**
     * @brief createNewLarva function to initialize a larva with a given time point.
     *
     *  This function is used during tracking to initialize new larval objects
     *
     * @param timePoint specifies the first detection for this larva
     * @param rawLarva contains the uprocessed raw larva and is used to calculate all larval parameters
     * @param larvaID specifies the unique ID of this larval object
     */
    void createNewLarva(uint timePoint, RawLarva const & rawLarva, unsigned int larvaID);
    
    /**
     * @brief insert is used to add new measurements (given by rawLarva) to this larva at a given time point
     * @param timePoint the given time point (for the parameters map)
     * @param rawLarva contains all raw larval values (i.e. features; e.g. the contour).
     */
    void insertRawLarva(uint larvaID, 
                        uint timePoint, 
                        RawLarva const & rawLarva);
    
    /**
     * @brief interpolateHeadTailOverTime changes the head/tail classification if necessary
     *
     * This is a post-processing step, which is called after the whole video is processed.
     */
    void interpolateHeadTailOverTime(uint larvaIndex);
    void interpolateHeadTailOverTime();
    
    /**
     * @brief fillTimeSamplingGaps some features are impossible to calculate online for all time steps,
     * which leads to sampling gaps (see below). These gaps are filled by this function.
     *
     * This is a post-processing step, which is called after the whole video is processed.
     * Sampling gaps are caused if frame rates lower than the actual frame rate of the movie
     * are used to calculate features (e.g. for the movement direction, 10fps are to fast to
     * recognice actual movement. However 1 second is sufficient to capture movement, thus the
     * movement direction can be calculated by resampling, i.e. calculate the direction between
     * the current frame and 10 frames in the past, which leads to gaps in the beginning of the
     * movie).
     */
    void fillTimeSamplingGaps(uint larvaIndex);
    void fillTimeSamplingGaps();
    
    /**
     * @brief isAssignedAt returns true if this larva is assigned for the given time point
     * @param timePoint time point to check existence of this larva
     * @return true if the larva exists for this time point, false otherwise
     */
    bool isAssignedAt(uint larvaID, uint timePoint) const;
    
//    void assign(uint timePoint, double const minDist, const std::vector<RawLarva> &curRawLarvae);
//    std::map<uint, double> findMatches(RawLarva const & rawLarva, double const minDist, uint const curTimePoint);
//    void assignMatches2Larva(std::map<size_t ,double> matches, double const minDist, RawLarva const & rawLarva, uint timePoint);
    void interplolateLarvae();
    
    void readLarvae(QString const& ymlFileName, 
                    std::vector<std::string> &imgPaths, 
                    bool useUndist);
    
    QPair<QVector<uint>, QVector<uint> > getVisibleLarvaID(uint time);
    
    QStringList getAllTimestepsBefore(uint id, uint time);
    QStringList getAllTimestepsAfter(uint id, uint time);
    QStringList getAllContemplableLarvaeIDsForAttach(uint id);
    
    void invertLarva(uint larvaID, uint currentTime, uint toTime);
    
    void attachToLarva(uint toLarvaID, uint fromLarvaID);
    
    bool copyModel2PrevTimeStep(uint larvaID, uint currentTime);
    bool copyModel2NextTimeStep(uint larvaID, uint currentTime);
    
    bool eraseLarvaAt(uint larvaID, uint time);
    bool eraseLarva(uint larvaID);
    
    void saveResultLarvae(const std::vector<std::string> &imgPaths, 
                          QImage const& img,
                          const bool useUndist, 
                          RegionOfInterestContainer const* ROIContainer = NULL,
                          LandmarkContainer const* landmarkContainer = NULL);
    
    /// Getter
    std::vector< Larva > getAllLarvae() const {return this->mLarvae;}
    QStringList getAllLarvaeIDs() const;
    int getNumberOfLarvae() const {return this->mLarvae.size();}
    Larva* getLarvaPointer(unsigned int index)
    {
        if(index < this->mLarvae.size())
        {
            return &(this->mLarvae.at(index));
        }
        return NULL;
    }
    
    double getMaxSpineLength() const {return this->mMaxSpineLength;}
    
    bool hasLoadedLarvae() const {return !this->mLarvae.empty();}
    bool isEmpty() const  {return this->mLarvae.empty();}
    bool getLarva(const uint index, Larva& l);
    
    std::vector<uint>   getAllTimesteps(uint larvaID);
    QPair<int, int>     getStartEndTimesteps(uint larvaID);
    std::vector< int >  getAllValidLarvaeIDS(uint timePoint);
    
    bool getSpineMidPointIndex(uint larvaID, uint &index) const;
    bool getSpinePointAt(const uint larvaID, const uint timePoint, const uint index, cv::Point &spinePoint) const;
    bool getAccDistAt(const uint larvaID, const uint timePoint, double & retAccDist);
    bool getDistToOriginAt(const uint larvaID, const uint timePoint, double & retDistToOrigin);
    bool getIsCoiledIndicatorAt(const uint larvaID, uint const timePoint, bool & retIsCoiledIndicator);
    bool getMomentumAt(const uint larvaID, const uint timePoint, cv::Point & retMomentum) const;
    
    /***************** For Plotting **********************/
    QVector<cv::Point>          getAllMomentumValues(uint larvaID);
    QVector<double>             getAllAreaValues(uint larvaID);
    QVector<double>             getAllMainBodybendingAngle(uint larvaID);
    QVector<double>             getAllCoiledIndicator(uint larvaID);
    QVector<double>             getAllPerimeter(uint larvaID);
    QVector<double>             getAllDistanceToOrigin(uint larvaID);
    QVector<double>             getAllMomentumDistance(uint larvaID);
    QVector<double>             getAllAccumulatedDistance(uint larvaID);
    QVector<double>             getAllGoPhaseIndicator(uint larvaID);
    QVector<double>             getAllLeftBendingIndicator(uint larvaID);
    QVector<double>             getAllRightBendingIndicator(uint larvaID);
    QVector<double>             getAllMovementDirection(uint larvaID);
    QVector<double>             getAllDistancesToLandmark(uint larvaID, QString const& landmarkID);
    QVector<double>             getAllBearingAnglesToLandmark(uint larvaID, QString const& landmarkID);
    QVector<double>             getVelocity(uint larvaID);
    QVector<double>             getAcceleration(uint larvaID);
    
    QVector<double>             getAllTimestepsForPlotting(uint larvaID);
    QVector<QVector<double> >   getAllTimestepGaps(uint larvaID);
    
signals:
    void sendLarvaModelDeleted();
    
    void sendUpdatedResultLarvaID(uint);
    void sendRemovedResultLarvaID(uint);
    void reset();
    
public slots:
    void updateLarvaValues(TrackerSceneLarva const* tLarva);
    void removeAllLarvae();
    
    void updateLandmark(const Landmark *l);
    void removeLandmark(QString name);
    
    void removeShortTracks(uint minTrackLenght);
    
};

#endif // LARVAECONTAINER_HPP
