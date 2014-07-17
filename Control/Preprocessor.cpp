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

#include "Preprocessor.hpp"

using namespace cv;
using std::vector;

Preprocessor::Preprocessor()
{
}

Mat & Preprocessor::graythresh(Mat const & src,
                                      int const thresh,
                                      Mat & dst)
{
    threshold(src,dst,thresh,255,THRESH_BINARY);
    return dst;
}

contoursType & Preprocessor::calcContours(Mat const & src, contoursType & contours)
{
    findContours(src, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE, Point());

    return contours;
}

contoursType & Preprocessor::calcPreviewContours(cv::Mat const & src, contoursType & contours)
{
    findContours(src, contours, CV_RETR_LIST, CV_CHAIN_APPROX_NONE, Point());

    return contours;
}

contoursType &Preprocessor::sizethreshold(const contoursType &contoursSrc, const int minSizeThresh, const int maxSizeThresh, contoursType &correctContoursDst, contoursType &biggerContoursDst)
{
    correctContoursDst.clear();
    biggerContoursDst.clear();

    // iterate over all contours
    for(contoursType::const_iterator it = contoursSrc.begin(); it != contoursSrc.end(); ++it)
    {
        // calculate the current size of the contour area
        double current_size = cv::contourArea((*it));
        // check the size (maxSizeThresh > current_size > minSizeThresh)
        if(current_size <= maxSizeThresh && current_size > minSizeThresh)
        {
            correctContoursDst.push_back(*it);
        }
        else if(current_size > maxSizeThresh)
        {
            biggerContoursDst.push_back(*it);
        }
    }

    // return the contours vector
    return correctContoursDst;
}

contoursType &Preprocessor::preprocessPreview2(const Mat &src,
                                               contoursType &acceptedContoursDst,
                                               contoursType &biggerContoursDst,
                                               const int gThresh,
                                               const int minSizeThresh,
                                               const int maxSizeThresh)
{
    // generate a scratch image
    Mat tmpImg = Mat::zeros(src.size().height, src.size().width, src.type());
    // generate a contours container scratch
    contoursType contours;
    // perform gray threshold
    Preprocessor::graythresh(src,gThresh,tmpImg);
    // calculate the contours
    Preprocessor::calcPreviewContours(tmpImg,contours);
    // filter the contours
    Preprocessor::sizethreshold(contours, minSizeThresh, maxSizeThresh, acceptedContoursDst, biggerContoursDst);

    // return the filtered contours
    return acceptedContoursDst;
}

contoursType &Preprocessor::preprocessTracking2(Mat const & src,
                                                        contoursType & acceptedContoursDst,
                                                        contoursType & biggerContoursDst,
                                                        int const gThresh,
                                                        int const minSizeThresh,
                                                        int const maxSizeThresh,
                                                        Backgroundsubtractor const & bs)
{
    // generate a scratch image
    Mat tmpImg = Mat::zeros(src.size().height, src.size().width, src.type());

//    bs.substract(src,tmpImg);
    bs.subtractViaThresh(src,gThresh,tmpImg);

    // generate a contours container scratch
    contoursType contours;
    // perform gray threshold
    Preprocessor::graythresh(tmpImg,gThresh,tmpImg);
    // calculate the contours
    Preprocessor::calcContours(tmpImg,contours);

    // filter the contours
    Preprocessor::sizethreshold(contours, minSizeThresh, maxSizeThresh, acceptedContoursDst, biggerContoursDst);

    // return the filtered contours
    return acceptedContoursDst;
}


//Mat & Preprocessor::medianblur(Mat const & src, Mat & dst)
//{
//    medianBlur(src,dst,3);
//    return dst;
//}

//Mat & Preprocessor::erodecircle(Mat const & src, Mat & dst)
//{
//    int elementData[] = {0,1,1,1,0,
//                         1,1,1,1,1,
//                         1,1,1,1,1,
//                         1,1,1,1,1,
//                         0,1,1,1,0 };

//    Mat element = Mat(5, 5, CV_8U, elementData).clone();

//    erode(src, dst, element, Point(-1,-1), 1, BORDER_CONSTANT, morphologyDefaultBorderValue());

//    return dst;
//}

//Mat & Preprocessor::dilatecircle(Mat const & src, Mat & dst)
//{
//    int elementData[] = {0,1,1,1,0,
//                         1,1,1,1,1,
//                         1,1,1,1,1,
//                         1,1,1,1,1,
//                         0,1,1,1,0 };

//    Mat element = Mat(5, 5, CV_8U, elementData).clone();

//    dilate(src, dst, element, Point(-1,-1), 1, BORDER_CONSTANT, morphologyDefaultBorderValue());

//    return dst;
//}
