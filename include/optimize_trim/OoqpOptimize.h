#ifndef __OOQPOPTIMIZE_H__
#define __OOQPOPTIMIZE_H__

#include "ooqp/QpGenData.h"
#include "ooqp/QpGenVars.h"
#include "ooqp/QpGenResiduals.h"
#include "ooqp/GondzioSolver.h"
#include "ooqp/QpGenSparseMa27.h"

#include <iostream>
#include <string>
#include <Eigen/Dense>

using namespace std;

class OptimizeOoqp
{
    private:
        bool EqualConstrain;
        bool UnequalConstrain;
        
    public:
        OptimizeOoqp();
        void Test_Solve();
        void ConEqualSolve(const Eigen::VectorXd& T, const int curveNums,const Eigen::Vector3d& begin,const Eigen::Vector3d& end,const Eigen::VectorXd& midlep);
        void SolveCurve(const double sumT,const int EquConNum,const int UnequConNum,const Eigen::Vector3d& begin,const Eigen::Vector3d& end,const Eigen::VectorXd& midle);
        void SingleSolveCurve(const double T,const int EqualConstrainNum,const int UnequalConNum,const double ps,const double vs,const double as,const double pe,const double ve,const double ae);
};

#endif