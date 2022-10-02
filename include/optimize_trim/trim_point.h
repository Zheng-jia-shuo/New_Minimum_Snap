#ifndef __TRIM_POINT_H__
#define __TRIM_POINT_H__

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class Trim_Point
{
    public:
    Trim_Point()
    {
        m_position=Eigen::Vector3d::Zero();
        m_velocity=Eigen::Vector3d::Zero();
        m_acceleration=Eigen::Vector3d::Zero();
        m_jerk=Eigen::Vector3d::Zero();
    }
    void SetPosition(Eigen::Vector3d position);
    void SetVelocity(Eigen::Vector3d velocity);
    void SetAcceleration(Eigen::Vector3d acceleration);
    void SetJerk(Eigen::Vector3d jerk);
    public:
    double m_pos;
    double m_vel;
    double m_acc;
    Eigen::Vector3d m_position;
    Eigen::Vector3d m_velocity;
    Eigen::Vector3d m_acceleration;
    Eigen::Vector3d m_jerk;
};

class Minimum_Snap
{
public:
    Minimum_Snap(){}
    Minimum_Snap(int segement_nums)
    {
        m_segement_nums=segement_nums;
        m_Q.resize(6*m_segement_nums,6*m_segement_nums);
        m_M.resize(6*m_segement_nums,6*m_segement_nums);
    }
    void CalSingleQ(double sum_time,Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& QMat);
    void CalSingleM(double sum_time,Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& MMat);
    void ConstructSingleC(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& CMat);
    void GainSingleKnownState(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& dFMat);
    void GainSingleUnknownState(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& dPMat);
private:
    int m_segement_nums;
    Eigen::MatrixXd m_Q;
    Eigen::MatrixXd m_M;
    Eigen::MatrixXd m_C;
    Eigen::MatrixXd m_dF;
    Eigen::MatrixXd m_dP;
    Eigen::Matrix<double,6,1> m_coeff;
};

#endif