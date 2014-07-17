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

#include "Undistorter.hpp"

Undistorter::Undistorter()
{
    this->mIsInitialised = false;
}

Undistorter::Undistorter(std::string const& path)
{
    InputGenerator::readMatrices(path, this->mCameraMatrix, this->mDistCoeffs, this->mImageSize);
    cv::initUndistortRectifyMap(this->mCameraMatrix,
                                this->mDistCoeffs,
                                cv::Mat(),
                                this->mCameraMatrix,
                                this->mImageSize,
                                CV_32FC1,
                                this->mMapX, this->mMapY);
    this->mIsInitialised = true;
    
}

Undistorter::Undistorter(cv::Mat const& cameraMatrix, 
                         cv::Mat const& distCoeffs, 
                         cv::Size const& imageSize)
{
    this->mCameraMatrix = cameraMatrix;
    this->mDistCoeffs = distCoeffs;
    this->mImageSize = imageSize;
    
    cv::initUndistortRectifyMap(this->mCameraMatrix,
                                this->mDistCoeffs,
                                cv::Mat(),
                                this->mCameraMatrix,
                                this->mImageSize,
                                CV_32FC1,
                                this->mMapX, this->mMapY);
    this->mIsInitialised = true;
}

void Undistorter::setPath(std::string const& path)
{
    InputGenerator::readMatrices(path, this->mCameraMatrix, this->mDistCoeffs, this->mImageSize);
    cv::initUndistortRectifyMap(this->mCameraMatrix,
                                this->mDistCoeffs,
                                cv::Mat(),
                                this->mCameraMatrix,
                                this->mImageSize,
                                CV_32FC1,
                                this->mMapX, this->mMapY);
    this->mIsInitialised = true;
}

void Undistorter::setParameter(cv::Mat const& cameraMatrix, 
                               cv::Mat const& distCoeffs, 
                               cv::Size const& imageSize)
{
    this->mCameraMatrix = cameraMatrix;
    this->mDistCoeffs = distCoeffs;
    this->mImageSize = imageSize;
    
    cv::initUndistortRectifyMap(this->mCameraMatrix,
                                this->mDistCoeffs,
                                cv::Mat(),
                                this->mCameraMatrix,
                                this->mImageSize,
                                CV_32FC1,
                                this->mMapX, this->mMapY);
    this->mIsInitialised = true;
}

cv::Mat& Undistorter::getUndistortImage(const cv::Mat &src, cv::Mat &dst) const
{
    if(this->mIsInitialised)
    {
        cv::remap(src, dst, this->mMapX, this->mMapY, cv::INTER_CUBIC);
    }
    
    return dst;
}

bool Undistorter::isReady() const
{
    return this->mIsInitialised;
}

void Undistorter::setReady(bool flag)
{
    this->mIsInitialised = flag;
}

void Undistorter::reset()
{
    this->mIsInitialised = false;
    this->mCameraMatrix.release();
    this->mDistCoeffs.release();
    this->mMapX.release();
    this->mMapY.release();
    this->mImageSize.height = this->mImageSize.width = 0;
}
