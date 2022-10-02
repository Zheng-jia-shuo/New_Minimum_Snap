#ifndef __BEZIERCURVE_H__
#define __BEZIERCURVE_H__

#include <iostream>
#include <Eigen/Dense>

class BezierCurve
{
    private:
        Eigen::Matrix<double,6,6> m_M;
    public:
        BezierCurve();
        BezierCurve(Eigen::Matrix<double,6,6> M);
        Eigen::Matrix<double,6,1> ControlPoint(Eigen::Matrix<double,6,1> &normalCoeff);
        Eigen::Matrix<double,5,1> FirstOrderControlPoint(Eigen::Matrix<double,6,1> &zeroControlPoint);
        Eigen::Matrix<double,4,1> SecondOrderControlPoint(Eigen::Matrix<double,6,1> &zeroControlPoint);
};

class BezierCurLimCon
{
    public:
        BezierCurLimCon(){};
        double UintTime(const double &no_uint_time);
        bool DisplaceBoundaryConstraints(const Eigen::Matrix<double,6,1> &displacebegincurve,const Eigen::Matrix<double,6,1> &displaceendcurve,const double &begindisplace,const double &enddisplace);
        bool VelBoundaryConstraints(const Eigen::Matrix<double,5,1> &velbegincurve,const Eigen::Matrix<double,5,1> &velendcurve,const double &beginvel,const double &endvel);
        bool AccBoundaryConstraints(const Eigen::Matrix<double,4,1> &accbegincurve,const Eigen::Matrix<double,4,1> &accendcurve,const double &beginacc,const double &endacc);
        bool ContinuityConstraints(const Eigen::Matrix<double,6,1> &lastcurvedisplace,const Eigen::Matrix<double,5,1> &lastcurvevel,const Eigen::Matrix<double,4,1> &lastcurveacc,
                                                                 const Eigen::Matrix<double,6,1> &curcurvedisplace,const Eigen::Matrix<double,5,1> &curcurvevel,const Eigen::Matrix<double,4,1> &curcurveacc);
        bool DisplaceSafetyConstraints(const Eigen::Matrix<double,6,1> &displacecurve,const double &maxdisplace,const double mindisplace);
        bool VelSafetyConstraints(const Eigen::Matrix<double,5,1> &velcurve,const double &maxvel,const double &minvel);
        bool AccSafetyConstraints(const Eigen::Matrix<double,4,1> &acccurve,const double &maxacc,const double &minacc);
};

#endif