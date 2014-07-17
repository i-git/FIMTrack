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

#ifndef RAWLARVA_HPP
#define RAWLARVA_HPP

#include <vector>
#include <iterator>
#include <cmath>

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-register"
#include <QTime>
#pragma clang diagnostic pop

#include "Configuration/FIMTrack.hpp"
#include "Configuration/TrackerConfig.hpp"

typedef std::vector<cv::Point> contourType;
typedef std::vector<double> contourCurvatureType;
typedef std::vector<cv::Point> spineType;

/**
 * @brief The RawLarva class stores and calculates all necessary parameters from a given contour to describe a larval object.
 * The acutal (non-raw) larval object is later created from the larval object.
 */

class RawLarva
{
public:

    /**
     * @brief RawLarva constructor generates a raw larva object based on the contour.
     *
     * In this constructor, dMin is 6% of the contour size, dMax = dMin + 2 and maskSize = dMin
     * if dMin is odd, dMin+1 otherwise and distToMax = 0.25 (secondMaxCurvaturePoint must be
     * at least 25% of the contour size away from the maxCurvaturePoint).
     *
     * @param _contour the contour defining the RawLarva
     * @param img  the image to calculate the head position based on brightness information
     */
    RawLarva(contourType const _contour, cv::Mat const & img);

    /** GETTER METHODS **/

    contourCurvatureType getCurvatures(void) const {return curvatures;}

    contourType getContour(void) const {return contour;}

    cv::Point getMaxCurvaturePoint(void) const {return maxCurvaturePoint;}

    cv::Point getSecondMaxCurvaturePoint(void) const {return secondMaxCurvaturePoint;}

    unsigned int getMaxCurvatureIndex(void) const {return maxCurvatureIndex;}

    unsigned int getSecondMaxCurvatureIndex(void) const {return secondMaxCurvatureIndex;}

    spineType getSpine(void) const {return spine;}

    spineType getDiscreteSpine(void) const {return discreteSpine;}

    cv::Point getMidPoint(void) const {return discreteSpine.at((discreteSpine.size() + 1 ) / 2);}

    cv::Point getMomentum(void) const {return momentum;}

    double getArea(void) const {return area;}

    std::vector<double> getLarvalRadii(void) const {return larvalRadii;}

    double getFirstContourHalfLength(void);
    double getSecondContourHalfLength(void);
    double getSpineLength(void) const {return arcLength(discreteSpine,false);}
    double getContourPerimeter(void) const {return arcLength(contour,true);}

    bool getIsCoiledIndicator(void) const {return isCoiled;}

private:
    /**
     * @brief contour of the raw larva. The first point (with index 0) contains the contour point with highest curvature
     */
    contourType contour;
    /**
     * @brief curvatures contains the curvature values for every contour point (in the same order).
     */
    contourCurvatureType curvatures;
    /**
     * @brief spine contains the central spine points
     */
    spineType spine;

    /**
     * @brief discreteSpine contains fraction of the spine (e.g. 5 points to describe the larva)
     */
    spineType discreteSpine;

    /**
     * @brief maxCurvatureIndex stores the index to the maximal curvature point (is equal to 0).
     */
    unsigned int maxCurvatureIndex;
    /**
     * @brief maxCurvaturePoint 2D contour point with maximal curvature (eual to contour.at(0))
     */
    cv::Point maxCurvaturePoint;

    /**
     * @brief secondMaxCurvatureIndex stores the index to the second maximal curvature point
     */
    unsigned int secondMaxCurvatureIndex;
    /**
     * @brief secondMaxCurvaturePoint 2D contour point with second maximal curvature point
     */
    cv::Point secondMaxCurvaturePoint;

    /**
     * @brief momentum is the momentum of the contour
     */
    cv::Point momentum;

    /**
     * @brief area is the area of the contour
     */
    double area;

    /**
     * @brief firstDiscreteHalf stores the discrete points of the first half
     */
    spineType firstDiscreteHalf;
    /**
     * @brief reverseSecondDiscreteHalf stores the points of the second discrete half
     */
    spineType reverseSecondDiscreteHalf;

    /**
     * @brief larvalRadii stores the radii of the discrete spine points. First (head) and last (tail) radii are 0;
     *          radii inbetween are >0.
     */
    std::vector<double> larvalRadii;

    bool isCoiled;

    // curvature calculation method:
    /**
     * @brief ipanFirstPass calculates the curvatures for every contour point based onthe
     *        IPAN algorithm (Chetverikov, D. (2003). A Simple and Efficient Algorithm for
     *        Detection of High Curvature Points in Planar Curves. In Computer Analysis of
     *        Images and Patterns (Vol. 2756, pp. 746–753)).
     *        For every contour point P, a trinangle is calculated with edges P- P and P P+
     *        with length between dMin and dMax (in pixel). The smallest angle defined
     *        by the triangle (P- P P+) in range [dMin, dMax] is assigned to the respective
     *        curvature vector position.
     *
     *        This function sets the curvature vector.
     *
     * @param dMin minimal length of the line between P- P and P P+.
     * @param dMax maximal length of the line between P- P and P P+.
     */
    void ipanFirstPass(unsigned int const dMin, unsigned int const dMax);

    /**
     * @brief findMaxCurvaturePoint calculates the maximal curvature point on the contour.
     *        For every contour point, a mean curvature is calculated based on the odd
     *        maskSize: Given a maskSize of (m*2+1), m curvature values left and m
     *        curvature values right of a central point (anchor point) are summed up and
     *        devided by the maskSize.
     *        After calculating the maxCurvaturePoint, the contour and the curvature vector
     *        are reordered, to set the maxCurvatureIndex to 0 and thus the maxCurvaturePoint
     *        to the first element of the contour vector (using the method reoderParameters()).
     *
     *        This function sets the values maxCurvatureIndex and maxCurvaturePoint
     *
     * @param maskSize must be an odd size for a mask (i.e. a sliding window) to calculate
     *        the mean curvature under the anchor point (mid point of the mask).
     */
    void findMaxCurvaturePoint(int maskSize);

    /**
     * @brief findSecondMaxCurvaturePointOrdered calculates the second maximal curvature
     *        point similar to the findMaxCurvaturePoint function with a specified relative
     *        distance distanceToMax to the first contour point.
     *
     *        This function sets the values secondMaxCurvatureIndex and secondMaxCurvaturePoint
     *
     * @param maskSize must be an odd size for a mask (i.e. a sliding window) to calculate
     *        the mean curvature under the anchor point (mid point of the mask).
     * @param distanceToMax relative distance to the maxCurvatureIndex (e.g. distanceToMax=0.25
     *        forces the secondMaxCurvatureIndex to be at least 0.25*contour.size() contour
     *        points away from the maxCurvatureIndex.
     */
    void findSecondMaxCurvaturePointOrdered(int maskSize, double const distanceToMax);

    /**
     * @brief initializeSpine calculates the spine for this contour. This is done in 3 steps:
     *        1) contour is devided into firstHalf (from maxCurvaturePoint to
     *           secondMaxCurvaturePoint) and secondHalf (from secondMaxCurvaturePoint to
     *           end).
     *        2) longer contour half is reduced to a discrete set with same length than the
     *           shorter contour half. Thus, each point of the firstHalf corresponds to a
     *           specific point on the secondHalf.
     *        3) Midpoint between the corresponding contour half points is calculetes and stored
     *           in the spine vector
     *
     *        This function sets the spine vector.
     *        In addition, this funciton sets the discrete left and right contour half
     *        (i.e. firstDiscreteContourHalf and reverseSecondDiscreteContourHalf).
     */
    void initializeSpine(void);

    /**
     * @brief calcDiscreteSpine calculates a discrete spine containing nPoints points from
     *        the spine
     *
     *        This function sets the discreteSpine vector
     *
     * @param nPoints an odd number for the number of discrete points
     */
    void calcDiscreteSpine(int nPoints);


    /**
     * @brief reorderParameters is called in the findMaxCurvaturePoint method.
     *        After calculating the maxCurvaturePoint, the contour and the curvature vector
     *        are reordered, to set the maxCurvatureIndex to 0 and thus the maxCurvaturePoint
     *        to the first element of the contour vector
     */
    void reorderParameters();

    /**
     * @brief calcMomentum calculates the momentum of the contour
     *
     *        This function sets the momentum of the larva
     */
    void calcMomentum(void);

    void calcIsCoiledIndicator(double const peri2spineLengthThresh, double const midCirclePeri2PeriThresh);


    /**
     * @brief getCircularContourNeighbourIndex calculates the index of a circular contour or
     *        curvature vector (if (0 <= (pos + offset) < contourSize) return: (pos + offset);
     *        else: return the "circular index" in [0 contourSize) ).
     * @param contourSize specifies the size of the contour
     * @param pos specifies the current position
     * @param offset specifies the offest to the current position (i.e. pos + offset)
     * @return the circular index defined by (pos + offset)
     */
    unsigned int getCircularContourNeighbourIndex(size_t const contourSize,
                                                  unsigned int const pos,
                                                  int const offset) const;


    // IPAN algorithm helper methods:
    /**
     * @brief calcSuccessorPointWithDistance calculates the sucessor point on a
     *        contour for a given index curIndex with distance dist (in pixel)
     * @param curIndex current intex
     * @param dist distance to the curIndex position (in pixel)
     * @return returns the successor point to the contour point at curIndex with
     *         distance dist
     */
    cv::Point calcSuccessorPointWithDistance(int unsigned curIndex, unsigned int dist);

    /**
     * @brief calcPredecessorPointWithDistancecalculates the predecessor point on a
     *        contour for a given index curIndex with distance dist (in pixel)
     * @param curIndex current intex
     * @param dist distance to the curIndex position (in pixel)
     * @return returns the predecessor point to the contour point at curIndex with
     *         distance dist
     */
    cv::Point calcPredecessorPointWithDistance(unsigned int curIndex, unsigned int dist);


};

#endif // RAWLARVA_HPP
