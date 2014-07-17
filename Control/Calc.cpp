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

#include "Calc.hpp"

using namespace cv;
using std::vector;

namespace Calc
{

double normL2(Point const & pt)
{
    return sqrt( double(pt.x*pt.x + pt.y*pt.y) );
}

double eucledianDist(Point const & pt1, Point const & pt2)
{
    return sqrt(std::pow((double) (pt1.x - pt2.x), 2) + std::pow((double) (pt1.y - pt2.y), 2));
}

double eucledianDist(cv::Point const & pt1, QPointF const & pt2)
{
    return sqrt(std::pow((double) (pt1.x - pt2.x()), 2) + std::pow((double) (pt1.y - pt2.y()), 2));
}

double eucledianDist(QPointF const & pt1, QPointF const & pt2)
{
    double d1 = pt1.x() - pt2.x();
    double d2 = pt1.y() - pt2.y();
    return std::sqrt(d1*d1 + d2*d2);
}

double eucledianDist(QPointF const& p, QLineF const& l)
{
    /**
     * @see http://community.topcoder.com/tc?module=Static&d1=tutorials&d2=geometry1
     */
    QPointF a = l.p1();
    QPointF b = l.p2();
    
    QPointF ab = b - a;
    QPointF ba = a - b;
    QPointF ap = p - a;
    QPointF bp = p - b;
    
    double distToLine = std::abs(calcCrossProduct(ab, ap) / eucledianDist(a, b));
    double distToLineSegment;
    double dot1 = calcDotProduct(ab, bp);
    double dot2 = calcDotProduct(ba, ap);
    
    if(dot1 > 0.0)
        distToLineSegment = eucledianDist(b, p);
    else if(dot2 > 0.0)
        distToLineSegment = eucledianDist(a, p);
    else
        distToLineSegment = distToLine;
    
    //    qDebug() << "Momentum: " << p << " -> Line: " << l << " -> DistanceToLine: " << distToLine << " -> Distance to Linesegment: " << distToLineSegment;
    
    return distToLineSegment;
    
    
    /**
         * @see http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
         * for Distance to linesegment @see http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment/1501725#1501725
         */
    /*
    double l2 = eucledianDist(l.p1(), l.p2());
    l2 = std::pow(l2, 2.0);
    
    if(l2 == 0.0)
        return eucledianDist(p, l.p1());
        
    double numerator    = abs( (l.p2().x() - l.p1().x())*(l.p1().y() - p.y()) - (l.p1().x() - p.x())*(l.p2().y() - l.p1().y()) );
    double denominator  = sqrt(pow((l.p2().x() - l.p1().x()), 2) + pow((l.p2().y() - l.p1().y()), 2));
    
    double distanceToLine = (numerator / denominator);
    
    QPointF v = l.p1();
    QPointF w = l.p2();
    double distanceToLineSegment;
    
    double t = calcDotProduct(p-v, w-v) / l2;
    
    if(t < 0.0)
    {
        distanceToLineSegment = eucledianDist(p, v);
    }
    else if(t > 1.0)
    {
        distanceToLineSegment = eucledianDist(p, w);
    }
    else
    {
        QPointF projection;
        projection.setX(v.x() + t*(w.x() - v.x()));
        projection.setY(v.y() + t*(w.y() - v.y()));
        distanceToLineSegment = eucledianDist(p, projection);
    }
    
    qDebug() << "Momentum: " << p << " -> Line: " << l << " -> t: " << t << " -> DistanceToLine: " << distanceToLine << " -> Distance to Linesegment: " << distanceToLineSegment;
    
    return distanceToLineSegment;
    */
}

double eucledianDist(QPointF const& p, QRectF const& r, const bool ellipse)
{
    if(!ellipse)
    {
        double leftDist     = 0.0;
        double topDist      = 0.0;
        double rightDist    = 0.0;
        double bottomDist   = 0.0;
        
        //            leftDist    = eucledianDist(p, QLineF(r.bottomLeft(),   r.topLeft()));
        //            topDist     = eucledianDist(p, QLineF(r.topLeft(),      r.topRight()));
        //            rightDist   = eucledianDist(p, QLineF(r.topRight(),     r.bottomRight()));
        //            bottomDist  = eucledianDist(p, QLineF(r.bottomRight(),  r.bottomLeft()));
        
        /*
             * tl --------- tr
             * |            |
             * |            |
             * |            |
             * |            |
             * bl --------- br
             */
        leftDist    = eucledianDist(p, QLineF(r.topLeft(),      r.bottomLeft()));
        topDist     = eucledianDist(p, QLineF(r.topLeft(),      r.topRight()));
        bottomDist  = eucledianDist(p, QLineF(r.bottomLeft(),   r.bottomRight()));
        rightDist   = eucledianDist(p, QLineF(r.topRight(),     r.bottomRight()));
        
        return qMin(bottomDist, qMin(rightDist, qMin(leftDist, topDist)));
    }
    else
    {
        double a,b,m,c,xm,ym;
        if(r.width() > r.height())
        {
            a = r.width() / 2;
            b = r.height() / 2;
        }
        else
        {
            a = r.height() / 2;
            b = r.width() / 2;
        }
        
        xm = r.center().x();
        ym = r.center().y();
        
        m = (ym - p.y()) / (xm- p.x());
        c = ym - m*xm;
        
        double A = (a*a)*(m*m) + b*b;
        double B = 2*(a*a)*m*c - 2*(a*a)*m*ym - 2*(b*b)*xm;
        double C = (a*a)*(b*b) - (a*a)*(c*c) - (b*b)*(xm*xm) - (a*a)*(ym*ym) + 2*(a*a)*c*ym;
        
        double D = sqrt(C/A + std::pow((B/(2*A)),2));
        
        //            qDebug() << QString("xm %1; ym %2; a %3; b %4; m %5; c %6; A %7; B %8, C %9; D %10")
        //                        .arg(xm).arg(ym).arg(a).arg(b)
        //                        .arg(m).arg(c)
        //                        .arg(A).arg(B).arg(C).arg(D);
        
        QPointF p1 = QPointF(D - (B/(2*A)), m*(D - (B/(2*A))) + c);
        QPointF p2 = QPointF(-D - (B/(2*A)), m*(-D - (B/(2*A))) + c);
        
        //            qDebug() << "p1" << p1;
        //            qDebug() << "p2" << p2;
        
        return qMin(eucledianDist(p, p1), eucledianDist(p, p2));
    }
}

double calcAngle(Point const & p, Point const & p1, Point const & p2)
{
    Point l1 = p1 - p;
    Point l2 = p2 - p;
    
    double norm1 = normL2(l1);
    double norm2 = normL2(l2);
    
    double testForNan = (l1.x*l2.x + l1.y*l2.y) / (norm1 * norm2);
    
    double angle;
    
    if(testForNan >= 1 || testForNan <=-1)
        angle = 180.0;
    else
        angle = (acos( (l1.x*l2.x + l1.y*l2.y) / (norm1 * norm2) ) * 180 / CV_PI);
    
    double sgn = (l1.x * l2.y) - (l1.y * l2.x);
    
    if(sgn < 0)
        angle = 360 - angle;
    return angle;
}

double calcAngleToYAxes(Point const & p1, Point const & p2)
{
    Point yAxesPoint(p1.x,0);
    double angle = Calc::calcAngle(p1,yAxesPoint,p2);
    return angle;
}

double calcCircularAngleSum(double const angle, double const offset)
{
    double sum = angle + offset;
    double retAngle = sum;
    if(sum<0)
    {
        retAngle = 360 - sum;
    }
    else if(sum>360)
    {
        retAngle = sum - 360;
    }
    
    return retAngle;
}

double calcAngleDiff(const double angle1, const double angle2)
{
    return std::min(std::abs(angle1 - angle2), (360 - (std::abs(angle1 - angle2))));
}

double angleToRadian(double const angle)
{
    return ((angle * M_PI) / 180.0);
}

double calcPolygonArea(QPolygonF const& polygon) 
{
    double area = 0.0;
    QPointF p1;
    QPointF p2;
    
    for (int i = 0; i < polygon.size() - 1; ++i)
    {
        p1 = polygon.at(i);
        p2 = polygon.at(i+1);
        area += (p1.x()*p2.y() - p2.x()*p1.y());
    }
    
    area /= 2;
    
    return abs(area);
}

cv::Point calcPolygonCenterOfMass(QPolygonF const& polygon)
{
    cv::Point centerOfMass(0.0,0.0);
    QPointF p1;
    QPointF p2;
    double area = calcPolygonArea(polygon);
    
    for(int i = 0; i < polygon.size() - 1; ++i) 
    {
        p1 = polygon.at(i);
        p2 = polygon.at(i+1);
        
        centerOfMass.x += ( (p1.x() + p2.x()) * (p1.x()*p2.y() - p2.x()*p1.y()) );
        centerOfMass.y += ( (p1.y() + p2.y()) * (p1.x()*p2.y() - p2.x()*p1.y()) );
        
    }
    
    centerOfMass.x /= (6*area);
    centerOfMass.y /= (6*area);
    
    centerOfMass.x = abs(centerOfMass.x);
    centerOfMass.y = abs(centerOfMass.y);
    
    return centerOfMass;
}

double calcSpineLength(std::vector<Point> const& spine)
{
    return cv::arcLength(spine, false);
}

double calcPerimeter(QPolygonF const& polygon)
{
    std::vector<Point> stdPolygon;
    for(int i = 0; i < polygon.size(); ++i)
    {
        stdPolygon.push_back(cv::Point(polygon.at(i).x(), polygon.at(i).y()));
    }
    
    return arcLength(stdPolygon, true);
}

double calcDotProduct(const QPointF &p1, const QPointF &p2)
{
    return (p1.x()*p2.x() + p1.y()*p2.y());
}

double calcCrossProduct(const QPointF &p1, const QPointF &p2)
{
    return (p1.x()*p2.y() - p1.y()*p2.x());
}

double calcSmallestAngle(const Point &p, const Point &p1, const Point &p2)
{
    Point l1 = p1 - p;
    Point l2 = p2 - p;
    
    double norm1 = normL2(l1);
    double norm2 = normL2(l2);
    
    double testForNan = (l1.x*l2.x + l1.y*l2.y) / (norm1 * norm2);
    
    double angle;
    
    if(testForNan >= 1 || testForNan <=-1)
        angle = 180.0;
    else
        angle = (acos( (l1.x*l2.x + l1.y*l2.y) / (norm1 * norm2) ) * 180 / CV_PI);
    
    return angle;
}

}
