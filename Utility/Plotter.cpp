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

#include "Plotter.hpp"

Plotter::Plotter(QWidget *parent) 
    : QCustomPlot(parent)
{
    //    this->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
    
    this->mMarkerSize = 5;
    
    this->setLocale(QLocale(QLocale::English, QLocale::UnitedStates));
    
    // Graph for Lineplots
    this->addGraph();
    this->graph(0)->setPen(QColor(6, 48, 200, 255));
    
    // Graph to visualize current Position i.e. Timepoint as Line
    this->addGraph(this->xAxis, this->yAxis);
    this->graph(1)->setPen(QColor(255, 0, 51, 255));
    
    // Graph for Marker
    this->addGraph(this->xAxis, this->yAxis);
    this->graph(2)->setPen(QColor(6, 48, 200, 255));
    this->graph(2)->setLineStyle(QCPGraph::lsNone);
    this->graph(2)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, this->mMarkerSize));
    
    // Graph to visualize current Position in Markerplot
    this->addGraph(this->xAxis, this->yAxis);
    this->graph(3)->setPen(QColor(255, 0, 51, 255));
    this->graph(3)->setLineStyle(QCPGraph::lsNone);
    this->graph(3)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssDisc, this->mMarkerSize));
}

void Plotter::setPlotTitle(const QString &title)
{
//    qDebug() << title;
//    this->plotLayout()->insertRow(0); // inserts an empty row above the default axis rect
//    this->plotLayout()->addElement(0, 0, new QCPPlotTitle(this, title));
//    this->replot();
}

void Plotter::setXAxisLable(const QString &lable)
{
    this->xAxis->setLabel(lable);
    this->replot();
}

void Plotter::setYAxisLable(const QString &lable)
{
    this->yAxis->setLabel(lable);
    this->replot();
}

void Plotter::setXAxisRange(double from, double to)
{
    this->xAxis->setRange(from, to);
    this->replot();
}

void Plotter::setYAxisRange(double from, double to)
{
    this->yAxis->setRange(from, to);
    this->replot();
}

void Plotter::plotDataAsLine(const QVector<double> &x, const QVector<double> &y)
{
    this->mXLine    = x;
    this->mYLine    = y;
    
    this->graph(2)->clearData();
    this->graph(3)->clearData();
    
    this->graph(0)->clearData();
    this->graph(0)->setData(this->mXLine, this->mYLine);
    this->setXAxisRange(*std::min_element(this->mXLine.begin(), this->mXLine.end())-1, *std::max_element(this->mXLine.begin(), this->mXLine.end())+1);
    this->setYAxisRange(*std::min_element(this->mYLine.begin(), this->mYLine.end())-1, *std::max_element(this->mYLine.begin(), this->mYLine.end())+1);
    this->replot();
}

void Plotter::plotDataAsMarker(const QVector<double> &x, const QVector<double> &y)
{
    this->mXMarker = x;
    this->mYMarker = y;

    this->graph(0)->clearData();
    this->graph(1)->clearData();
    
    this->graph(2)->clearData();
    this->graph(2)->setData(this->mXMarker, this->mYMarker);
    //    this->setXAxisRange(*std::min_element(this->mXMarker.begin(), this->mXMarker.end())-1, *std::max_element(this->mXMarker.begin(), this->mXMarker.end())+1);
    //    this->setYAxisRange(*std::min_element(this->mYMarker.begin(), this->mYMarker.end())-1, *std::max_element(this->mYMarker.begin(), this->mYMarker.end())+1);
    this->setXAxisRange(0.0, this->mImageSize.width() + 1);
    this->setYAxisRange(0.0, this->mImageSize.height() + 1);
    this->replot();
}

void Plotter::plotCurrentTimeStepMarker(double timeStep, double firstTimeStep, bool isMarkerPlot)
{
    if(!this->graph(0)->data()->empty())
    {        
        this->graph(1)->clearData();
        this->graph(1)->setData(QVector<double>() << (timeStep+1) << (timeStep+1), QVector<double>() << *std::min_element(this->mYLine.begin(), this->mYLine.end()) << *std::max_element(this->mYLine.begin(), this->mYLine.end()));
        this->replot();
    }
    
    if(!this->graph(2)->data()->empty())
    {
        this->graph(3)->clearData();
        this->graph(3)->setData(QVector<double>() << this->mXMarker.at(timeStep - firstTimeStep + 1), QVector<double>() << this->mYMarker.at(timeStep - firstTimeStep + 1));
        this->replot();
    }
}

void Plotter::reset()
{
    this->clearFocus();
    this->clearGraphs();
    this->clearItems();
    this->clearMask();
    this->clearPlottables();
}
