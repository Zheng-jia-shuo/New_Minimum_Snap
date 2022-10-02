#include <iostream>
#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/bezier_curve.h"

BezierCurve::BezierCurve()
{
    m_M.setZero();
    m_M << 1,     0,      0,      0,      0,      0,
                    -5,    5,      0,      0,      0,      0,
                    10,   -20,  10,   0,      0,      0,
                    -10,  30,   -30, 10,    0,      0,
                    5,      -20,  30,  -20,  5,      0,
                    -1,     5,     -10, 10,   -5,     1;
}

BezierCurve::BezierCurve(Eigen::Matrix<double,6,6> M)
{
    m_M=M;
}

Eigen::Matrix<double,6,1> BezierCurve::ControlPoint(Eigen::Matrix<double,6,1> &normalCoeff)
{
    //P=M*C
    Eigen::Matrix<double,6,1> tmp = Eigen::Matrix<double,6,1>::Zero();
    tmp = m_M.inverse() * normalCoeff;
    return tmp;
}

Eigen::Matrix<double,5,1> BezierCurve::FirstOrderControlPoint(Eigen::Matrix<double,6,1> &zeroControlPoint)
{
    Eigen::Matrix<double,5,1> tmp = Eigen::Matrix<double,5,1>::Zero();
    double C0 = 5 * (zeroControlPoint(4,0) - zeroControlPoint(5,0));
    double C1 = 5 * (zeroControlPoint(3,0) - zeroControlPoint(4,0));
    double C2 = 5 * (zeroControlPoint(2,0) - zeroControlPoint(3,0));
    double C3 = 5 * (zeroControlPoint(1,0) - zeroControlPoint(2,0));
    double C4 = 5 * (zeroControlPoint(0,0) - zeroControlPoint(1,0));

    tmp << C4,C3,C2,C1,C0;
    return tmp;
}

Eigen::Matrix<double,4,1> BezierCurve::SecondOrderControlPoint(Eigen::Matrix<double,6,1> &zeroControlPoint)
{
    Eigen::Matrix<double,4,1> tmp = Eigen::Matrix<double,4,1>::Zero();
    double C0 = 4 * 5 * (zeroControlPoint(3,0) + zeroControlPoint(5,0) - 2 * zeroControlPoint(4,0));
    double C1 = 4 * 5 * (zeroControlPoint(2,0) + zeroControlPoint(4,0) - 2 * zeroControlPoint(3,0));
    double C2 = 4 * 5 * (zeroControlPoint(1,0) + zeroControlPoint(3,0) - 2 * zeroControlPoint(2,0));
    double C3 = 4 * 5 * (zeroControlPoint(0,0) + zeroControlPoint(2,0) - 2 * zeroControlPoint(1,0));

    tmp << C3,C2,C1,C0;
    return tmp;
}

//贝塞尔曲线的时间区间必须规一化到[0,1]中。
double BezierCurLimCon::UintTime(const double &no_uint_time)
{
    return no_uint_time * 1;
}

bool BezierCurLimCon::DisplaceBoundaryConstraints(const Eigen::Matrix<double,6,1> &displacebegincurve,const Eigen::Matrix<double,6,1> &displaceendcurve,const double &begindisplace,const double &enddisplace)
{
    double begindisp = displacebegincurve(5,0);
    double enddisp = displaceendcurve(0,0);
    if((begindisp == begindisplace) && (enddisp == enddisplace))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::VelBoundaryConstraints(const Eigen::Matrix<double,5,1> &velbegincurve,const Eigen::Matrix<double,5,1> &velendcurve,const double &beginvel,const double &endvel)
{
    double beginv = velbegincurve(4,0);
    double endv = velendcurve(0,0);
    if((beginv == beginvel) && (endv == endvel))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::AccBoundaryConstraints(const Eigen::Matrix<double,4,1> &accbegincurve,const Eigen::Matrix<double,4,1> &accendcurve,const double &beginacc,const double &endacc)
{
    double begina = accbegincurve(3,0);
    double enda = accendcurve(0,0);
    if((begina == beginacc) && (enda == endacc))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::ContinuityConstraints(const Eigen::Matrix<double,6,1> &lastcurvedisplace,const Eigen::Matrix<double,5,1> &lastcurvevel,const Eigen::Matrix<double,4,1> &lastcurveacc,
                                                                                                const Eigen::Matrix<double,6,1> &curcurvedisplace,const Eigen::Matrix<double,5,1> &curcurvevel,const Eigen::Matrix<double,4,1> &curcurveacc)
{
    double DCL5 = lastcurvedisplace(0,0);
    double DCL0 = lastcurvedisplace(5,0);
    double DCC5 = curcurvedisplace(0,0);
    double DCC0 = curcurvedisplace(5,0);

    double VCL5 = lastcurvevel(0,0);
    double VCL0 = lastcurvevel(5,0);
    double VCC5 = curcurvevel(0,0);
    double VCC0 = curcurvevel(5,0);

    double ACL5 = lastcurveacc(0,0);
    double ACL0 = lastcurveacc(5,0);
    double ACC5 = curcurveacc(0,0);
    double ACC0 = curcurveacc(5,0);

    if((DCL5 == DCC0) && (VCL5 == VCC0) && (ACL5 == ACC0))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::DisplaceSafetyConstraints(const Eigen::Matrix<double,6,1> &displacecurve,const double &maxdisplace,const double mindisplace)
{
    double C5 = displacecurve(0,0);
    double C4 = displacecurve(1,0);
    double C3 = displacecurve(2,0);
    double C2 = displacecurve(3,0);
    double C1 = displacecurve(4,0);
    double C0 = displacecurve(5,0);

    if((C0 <= maxdisplace) && (C0 >= mindisplace) && (C1 <= maxdisplace) && (C1 >= mindisplace) &&(C2 <= maxdisplace) && (C2 >= mindisplace) &&
        (C3 <= maxdisplace) && (C3 >= mindisplace) &&(C4 <= maxdisplace) && (C4 >= mindisplace) &&(C5 <= maxdisplace) && (C5 >= mindisplace))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::VelSafetyConstraints(const Eigen::Matrix<double,5,1> &velcurve,const double &maxvel,const double &minvel)
{
    double C4 = velcurve(0,0);
    double C3 = velcurve(1,0);
    double C2 = velcurve(2,0);
    double C1 = velcurve(3,0);
    double C0 = velcurve(4,0);

    if((C0 <= maxvel) && (C0 >= minvel) && (C1 <= maxvel) && (C1 >= minvel) &&(C2 <= maxvel) && (C2 >= minvel) &&
        (C3 <= maxvel) && (C3 >= minvel) &&(C4 <= maxvel) && (C4 >= minvel))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool BezierCurLimCon::AccSafetyConstraints(const Eigen::Matrix<double,4,1> &acccurve,const double &maxacc,const double &minacc)
{
    double C3 = acccurve(0,0);
    double C2 = acccurve(1,0);
    double C1 = acccurve(2,0);
    double C0 = acccurve(3,0);

    if((C0 <= maxacc) && (C0 >= minacc) && (C1 <= maxacc) && (C1 >= minacc) &&(C2 <= maxacc) && (C2 >= minacc) &&
        (C3 <= maxacc) && (C3 >= minacc))
    {
        return true;
    }
    else
    {
        return false;
    }
}
