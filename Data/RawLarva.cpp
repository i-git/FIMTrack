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

#include "RawLarva.hpp"

using namespace cv;
using std::vector;

typedef vector<Point> contourType;
typedef vector<double> contourCurvatureType;
typedef vector<Point> spineType;

RawLarva::RawLarva(contourType const _contour, Mat const & img)
{

    Q_UNUSED(img);
    contour=_contour;
    curvatures.reserve(contour.size());

    // number of discrete spine points
    int nPoints = 5;

    // setting the relative distToMax value for secondMaxCurvaturePoint calculation (at least
    // 25% away from the first curvature maximum
    double distToMax = 0.25;

    // dMin and dMax values for the IPAN algorithm (i.e. curvature calculation)
    int dMin = contour.size() * 0.06;
    int dMax = dMin + 2;

    // odd mask size for max curvature calculation
    int maskSize = (dMin % 2 == 0) ? dMin+1 : dMin;

    double peri2spineLengthThresh = 2.6;
    double midCirclePeri2PeriThresh = 0.5;

    if(!LarvaeExtractionParameters::bUseDefault)
    {
        nPoints = LarvaeExtractionParameters::iNumerOfSpinePoints;

        peri2spineLengthThresh = LarvaeExtractionParameters::CoiledRecognitionParameters::dPerimeterToSpinelenghtThreshold;
        midCirclePeri2PeriThresh = LarvaeExtractionParameters::CoiledRecognitionParameters::dMidcirclePerimeterToPerimeterThreshold;

    }

    if(!LarvaeExtractionParameters::IPANContourCurvatureParameters::bUseDynamicIpanParameterCalculation)
    {
        dMin = LarvaeExtractionParameters::IPANContourCurvatureParameters::iMinimalTriangelSideLenght;
        dMax = LarvaeExtractionParameters::IPANContourCurvatureParameters::iMaximalTriangelSideLenght;
        maskSize = LarvaeExtractionParameters::IPANContourCurvatureParameters::iCurvatureWindowSize;
        distToMax = LarvaeExtractionParameters::IPANContourCurvatureParameters::dMaximalCurvaturePointsDistance;
    }


    // setting the curvature vector
    ipanFirstPass(dMin, dMax);

    // setting the maxCurvatureIndex and maxCurvaturePoint and reordering the contour and curvatures vectors
    findMaxCurvaturePoint(maskSize);

    // setting the secondMaxCurvatureIndex and Point
    findSecondMaxCurvaturePointOrdered(maskSize, distToMax);

    // setting the spine
    initializeSpine();

    // setting the discrete spine
    calcDiscreteSpine(nPoints);

    // set the momentum
    calcMomentum();

    // set the contour area value
    area = contourArea(contour);

    calcIsCoiledIndicator(peri2spineLengthThresh,midCirclePeri2PeriThresh);

}

double RawLarva::getFirstContourHalfLength()
 {
    spineType firstHalf = firstDiscreteHalf;
    firstHalf.push_back(secondMaxCurvaturePoint);
    return arcLength(firstHalf,false);
 }

double RawLarva::getSecondContourHalfLength()
{
    spineType secondHalf;
    secondHalf.reserve(reverseSecondDiscreteHalf.size() +1);
    secondHalf.push_back(maxCurvaturePoint);
    for (spineType::const_iterator it=reverseSecondDiscreteHalf.begin(); it != reverseSecondDiscreteHalf.end(); ++it)
    {
        secondHalf.push_back(*it);
    }
    return arcLength(secondHalf,false);
}


void RawLarva::ipanFirstPass(unsigned int const dMin, unsigned int const dMax)
{
    // get the contour size
    size_t contourSize = contour.size();

    // for all points int hte contour
    for (int i=0; i < (int) contourSize; ++i)
    {
        // get the current point
        Point curPoint = contour.at(i);
        // initialize the angle for this point with maximum value (i.e. 360 degree)
        double angle = 360;
        // for all points in the contour with distance in [dMin, dMax]
        for (unsigned int dist = dMin; dist <= dMax; ++dist)
        {
            // get the successor point
            Point successor = calcSuccessorPointWithDistance(i, dist);
            // get the predecessor point
            Point predecessor = calcPredecessorPointWithDistance(i, dist);
            // calculate the current angle value given by the triangle (predecessor, curPoint, successor)
            double curAngle = Calc::calcAngle(curPoint, predecessor, successor);
            // is the current angle smaller?
            if (curAngle < angle)
            {
                // update the angle
                angle = curAngle;
            }
        }  // end for-loop for points within distance in [dMin, dMax]

        // push the smallest found angle into the curvatures vector
        curvatures.push_back(angle);
    }
}

void RawLarva::findMaxCurvaturePoint(int maskSize)
{
    // guarantee odd mask size
    maskSize = (maskSize % 2 == 0) ? maskSize+1 : maskSize;

    // get the amoutn of points left / right from the mask center
    int moveLeftRight = (int) maskSize / 2;

    // initialize the meanBestCurvature with maximal curvature value (i.e. 360 degree)
    double meanBestCurv = 360;
    // initialize the bestIndex (i.e. index with sharpest angle for the given mask)
    int bestIndex = 0;

    // for all contour points
    for(int contourPointNo = 0; contourPointNo < (int) contour.size(); ++contourPointNo)
    {
        // set the current mean curvature to 0
        double curMeanCurvature = 0;
        // for all points within the mask (i.e. [-moveLeftRight, moveLeftRight])
        for (int offset = -moveLeftRight; offset <= moveLeftRight; ++offset)
        {
            // get the (circular) index
            int curIndex = getCircularContourNeighbourIndex(contour.size(), contourPointNo, offset);
            // add the current curvature to the current mean curMeanCurvature
            curMeanCurvature += curvatures.at(curIndex);
        }

        // normalize the curMeanCurvature to get the real mean
        curMeanCurvature /= maskSize;

        // is the current curvature smaller than the best mean curvature?
        if (curMeanCurvature < meanBestCurv)
        {
            // update the best mean curvature
            meanBestCurv = curMeanCurvature;
            // update the bestIndex
            bestIndex = contourPointNo;
        }
    }

    // store the best index and best point (with sharpest angle within the mask)
    maxCurvatureIndex = bestIndex;
    maxCurvaturePoint = contour.at(bestIndex);

    // reorder the parameters for this RawLarva so that the first contour and curvatures element
    // is the one with the sharpest angle
    reorderParameters();
}

void RawLarva::reorderParameters(void)
{
    // get the contour size
    size_t contourSize = contour.size();

    // define parameters for the ordered contour and ordered curvatures vectors
    contourType orderedContour;
    orderedContour.reserve(contourSize);

    contourCurvatureType orderedCurvatures;
    orderedCurvatures.reserve(contourSize);

    // for all elements inside the contour and curvatures vector from maxCurvatureIndex to contourSize
    for (unsigned int i = maxCurvatureIndex; i < contourSize; ++i)
    {
        // add the respective elements into the new ordered vectors
        orderedContour.push_back(contour.at(i));
        orderedCurvatures.push_back(curvatures.at(i));
    }
    // add the remaining elements from the contour and curvatures vector to the ordered vectors
    for (unsigned int i = 0; i < maxCurvatureIndex; ++i)
    {
        orderedContour.push_back(contour.at(i));
        orderedCurvatures.push_back(curvatures.at(i));
    }

    // set the contour and curvatures vector to the ordered vectors
    contour = orderedContour;
    curvatures = orderedCurvatures;
    // reset the maxCurvatureIndex to 0
    maxCurvatureIndex = 0;
}


void RawLarva::findSecondMaxCurvaturePointOrdered(int maskSize, double const distToMax)
{
    // guarantee odd mask size
    maskSize = (maskSize % 2 == 0) ? maskSize+1 : maskSize;

    // get the distance offset ot the maxCurvaturePoint
    int distOffset = contour.size() * distToMax;

    // get the amoutn of points left / right from the mask center
    int moveLeftRight = (int) maskSize / 2;

    // initialize second best mean curvature
    double meanSecondBestCurv = 360;
    // set second best index to the opposite of the (first) maximum curvature index (wich is 0)
    int secondBestIndex = contour.size() / 2;

    // for all contour points in [distOffset, contourSize-distOffset] (to guarantee a distance to the
    // maximum curvature index)
    for(int contourPointNo = distOffset; contourPointNo <= (int) contour.size() - distOffset; ++contourPointNo)
    {
        // initialize the current mean curvature
        double meanCurvature = 0;
        // for all contour points within the mask
        for (int offset = -moveLeftRight; offset <= moveLeftRight; ++offset)
        {
            // calculate the current index
            int curIdnex = getCircularContourNeighbourIndex(contour.size(), contourPointNo, offset);
            // sum up the curvature values (within the mask)
            meanCurvature += curvatures.at(curIdnex);
        }

        // devide the summed curvature values to get the mean
        meanCurvature /= maskSize;

        // is the current mean curvature smaller than the second best mean curvature?
        if (meanCurvature < meanSecondBestCurv)
        {
            // update mean second best curvature and its index
            meanSecondBestCurv = meanCurvature;
            secondBestIndex = contourPointNo;
        }
    }

    // set the second best index and point (with second sharpest curvature with distance to the first)
    secondMaxCurvatureIndex = secondBestIndex;
    secondMaxCurvaturePoint = contour.at(secondBestIndex);
}

void RawLarva::initializeSpine(void)
{
    // reserve space for the first and the second contour half
    contourType firstHalf;
    contourType reverseSecondHalf;
    firstHalf.reserve(secondMaxCurvatureIndex);
    reverseSecondHalf.reserve(contour.size() - secondMaxCurvatureIndex);

    // push all contour points of the first half into the vector firstHalf (in [0,secondMaxCurvatureIndex) )
    for (unsigned int i = 0; i<secondMaxCurvatureIndex; ++i)
    {
        firstHalf.push_back(contour.at(i));
    }
    // push all contour points fo the second half into the vector
    // secondHalf (in [secondMaxCurvatureIndex, contourSize) )
    // and flip the order
    for (unsigned int i = contour.size() - 1; i >= secondMaxCurvatureIndex; --i)
    {
        reverseSecondHalf.push_back(contour.at(i));
    }

    // the the size of the smaller half
    int discreteSize = std::min(firstHalf.size(), reverseSecondHalf.size());

    // reserve memory for both discrete half vectors
    firstDiscreteHalf.reserve(discreteSize);
    reverseSecondDiscreteHalf.reserve(discreteSize);

    // is the first half bigger than the second half? (first half will be reduced...)
    if(firstHalf.size() > reverseSecondHalf.size())
    {
        // calculate the fraction (>1) to jump to the positions in the longer first half
        double dOffset = (double) firstHalf.size() / (double) reverseSecondHalf.size();
        // initialize the double index
        double dIndex = 0.0;

        // initialize the integer indices for the first and second contour half
        unsigned int indexFirst = 0;
        unsigned int indexSecond = 0;

        // while there are still points in the second (shorter) contour half and
        // there are still points in the first (longer) contour half
        while(indexSecond < reverseSecondHalf.size() && indexFirst < firstHalf.size())
        {
            // add the points from the first and second half to the discrete vectors
            firstDiscreteHalf.push_back(firstHalf.at(indexFirst));
            reverseSecondDiscreteHalf.push_back(reverseSecondHalf.at(indexSecond));

            // incement the double index by the double offset (>1)
            dIndex += dOffset;
            // round the double index to get the first half index
            indexFirst = cvRound(dIndex);
            // increment the second half index
            ++indexSecond;
        }
    }
    // is the second contour half bigger than the first contour half?
    else if(reverseSecondHalf.size() > firstHalf.size())
    {
        // same steps as in the if-case above...
        double dOffset = (double) reverseSecondHalf.size() / (double) firstHalf.size();
        double dIndex = 0.0;

        unsigned int indexFirst = 0;
        unsigned int indexSecond = 0;

        while(indexFirst < firstHalf.size() && indexSecond < reverseSecondHalf.size())
        {
            firstDiscreteHalf.push_back(firstHalf.at(indexFirst));
            reverseSecondDiscreteHalf.push_back(reverseSecondHalf.at(indexSecond));

            dIndex += dOffset;
            ++indexFirst;
            indexSecond = cvRound(dIndex);
        }
    }
    // both contour halfs are equally sized?
    else
    {
        // discrete contour halves are the same to the non-discrete halves
        firstDiscreteHalf = firstHalf;
        reverseSecondDiscreteHalf = reverseSecondHalf;
    }

    // reserve memory for the spine
    spine.reserve(firstDiscreteHalf.size());


    // for all discrete contour half points
    for(unsigned int i = 0; i<firstDiscreteHalf.size(); ++i)
    {
        // calculate the mid-points between the first and second half contour points
        int x = (firstDiscreteHalf.at(i).x + reverseSecondDiscreteHalf.at(i).x) / 2;
        int y = (firstDiscreteHalf.at(i).y + reverseSecondDiscreteHalf.at(i).y) / 2;

        // add the mid-points to the spine vector
        spine.push_back(Point(x,y));
    }
}

void RawLarva::calcDiscreteSpine(int nPoints)
{
    // guarantee odd number of points
    nPoints = (nPoints % 2 == 0) ? nPoints+1 : nPoints;

    // reserve memory for the discrete spine
    discreteSpine.reserve(nPoints);

    // reserve memory for the radii vector
    larvalRadii.reserve(nPoints);
    // add the first radius (which is 0)
    larvalRadii.push_back(0.0);

    // add the first point to the discrete spine
    discreteSpine.push_back(maxCurvaturePoint);

    // remaining descrete points must be ((x * (1/(nPoints-1)) * conoturSize) away from the
    // first point (with x in [1,nPoints-2]; -2 is cause by the first and last point)
    double denominator = nPoints - 1;
    double fraction = 1.0 / denominator;

    // add all discrete points form 2 to nPoints-1 to the discrete spine
    for (int i=1; i < (int) denominator; ++i)
    {
        int index = cvRound((double) spine.size() * (i*fraction));
        discreteSpine.push_back(spine.at(index));

        // calculate the larval radii:
        // a radius is defined as the smaller distance between the spine point and
        // the two corresponding contour points.
        // 1. get the spine point
        Point spinePoint = spine.at(index);
        // 2. get the corresponding contour points
        Point contourPoint1 = firstDiscreteHalf.at(index);
        Point contourPoint2 = reverseSecondDiscreteHalf.at(index);
        // 3. calculate the two eucledian distances
        double distance1 = Calc::eucledianDist(contourPoint1, spinePoint);
        double distance2 = Calc::eucledianDist(contourPoint2, spinePoint);
        // 4. select the smaller distance value
        double distance = distance1 <= distance2 ? distance1 : distance2;
        // store the distances (i.e. radii) in the vector
        larvalRadii.push_back(distance);
    }

    // add the last point to the discrete spine
    discreteSpine.push_back(secondMaxCurvaturePoint);

    // add the last radius (which is 0.0)
    larvalRadii.push_back(0.0);
}

void RawLarva::calcIsCoiledIndicator(const double peri2spineLengthThresh,
                                     const double midCirclePeri2PeriThresh)
{
    double perimeter = getContourPerimeter();
    double spineLength = getSpineLength();
    double peri2spineLengthRatio = perimeter / spineLength;

    int midPointIndex = (int) ((discreteSpine.size()-1) /2);
    double midPointRadius = larvalRadii.at(midPointIndex);
    double midCirclePeri2PeriRatio = (2*midPointRadius*CV_PI) / perimeter;

    if( peri2spineLengthRatio >= peri2spineLengthThresh &&
        midCirclePeri2PeriRatio >= midCirclePeri2PeriThresh)
    {
        isCoiled = true;
    }
    else
    {
        isCoiled = false;
    }
}

/** IPAN HELPER METHODS **/
Point RawLarva::calcSuccessorPointWithDistance(unsigned int curIndex, unsigned int dist)
{
    size_t contourSize = contour.size();
    Point curPoint = contour.at(curIndex);

    int offset = 1;
    int nextIndex = getCircularContourNeighbourIndex(contourSize, curIndex, offset);
    Point nextPoint = contour.at(nextIndex);

    double curDist = Calc::eucledianDist(curPoint, nextPoint);

    while(curDist < (double) dist && offset < (int) contourSize)
    {
        ++offset;
        nextIndex = getCircularContourNeighbourIndex(contourSize, curIndex, offset);
        nextPoint = contour.at(nextIndex);
        curDist = Calc::eucledianDist(curPoint, nextPoint);
    }

    return nextPoint;
}

Point RawLarva::calcPredecessorPointWithDistance(unsigned int curIndex, unsigned int dist)
{
    size_t contourSize = contour.size();
    Point curPoint = contour.at(curIndex);

    int offset = -1;
    int nextIndex = getCircularContourNeighbourIndex(contourSize, curIndex, offset);
    Point nextPoint = contour.at(nextIndex);

    double curDist = Calc::eucledianDist(curPoint, nextPoint);

    while(curDist < (double) dist && abs(offset) < (int) contourSize)
    {
        --offset;
        nextIndex = getCircularContourNeighbourIndex(contourSize, curIndex, offset);
        nextPoint = contour.at(nextIndex);
        curDist = Calc::eucledianDist(curPoint, nextPoint);
    }

    return nextPoint;
}

/** HELPER METHODS **/
unsigned int RawLarva::getCircularContourNeighbourIndex(size_t const contourSize, unsigned int const pos,  int const offset) const
{
    int jump = pos + offset;
    unsigned int index;
    if (jump < 0)
        index = contourSize + jump;
    else if(jump >= (int) contourSize)
        index = jump - contourSize;
    else
        index = jump;

    return index;
}

void RawLarva::calcMomentum(void)
{
    // calculate the moments of the contour
    Moments moms = moments(contour);
    // calculate center of the moments
    double cx = moms.m10 / moms.m00;
    double cy = moms.m01 / moms.m00;
    // set the momentum
    momentum = Point(cx,cy);
}
