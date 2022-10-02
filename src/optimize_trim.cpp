#include <iostream>
#include <fstream>
#include <ros/ros.h>
#include <std_msgs/Float64MultiArray.h>
#include <Eigen/Dense>
#include <pangolin/pangolin.h>
#include <cstdlib>
#include <cmath>
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/trim_point.h"
#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/bezier_curve.h"
#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/OoqpOptimize.h"
#include "/home/zjs/minimum_snap_ws/devel/include/optimize_trim/eigen2ros.h"

#define SEGEMENTNUMS 4
#define TIMEINTER 4

using namespace std;
int curvesNums = 0;
struct My_Point
{
    double x;
    double y;
    double z;
};
void WriteFile(double d1,double d2,double d3,double d4,double d5,double d6);
void DrawLine(My_Point &p1,My_Point &p2,My_Point &p3,My_Point &p4,My_Point &p5)
{
    pangolin::CreateWindowAndBind("Main",640,480);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::OpenGlRenderState s_cam(pangolin::ProjectionMatrix(640,480,420,420,320,320,0.2,100),pangolin::ModelViewLookAt(2,0,2, 0,0,0, pangolin::AxisY));

    pangolin::Handler3D handler(s_cam);
    pangolin::View& d_cam = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, 0.0, 1.0, -640.0f/480.0f)
            .SetHandler(&handler);

    while( !pangolin::ShouldQuit() )
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        d_cam.Activate(s_cam);

		//绘制直线
		glBegin(GL_LINE_STRIP);
        glLineWidth(5.0);
        glColor3f(5.0, 6.0, 7.0);
        glVertex2d(p1.x,p1.y);
        glVertex2d(p2.x,p2.y);
        glVertex2d(p3.x,p3.y);
        glVertex2d(p4.x,p4.y);
        glVertex2d(p5.x,p5.y);
        glEnd();

        pangolin::FinishFrame();
    }
}

void test01(std::vector< optimize_trim::eigen2ros > &coeff)
{
    Trim_Point startpoint;
    Trim_Point endpoint;
    std::vector<Trim_Point> v_midlepoints;
    Minimum_Snap snap;
    Eigen::MatrixXd Q;
    Eigen::MatrixXd M;
    Eigen::MatrixXd C;
    Eigen::MatrixXd dF;
    Eigen::MatrixXd dP;
    Eigen::MatrixXd dF_dP;
    Eigen::Matrix<double,6,6> Coeff_P2C = Eigen::Matrix<double,6,6>::Zero();
    double time=5;

    //P = M * C
    Coeff_P2C<<  1,     0,      0,      0,      0,      0,
                                -5,    5,      0,      0,      0,      0,
                                10,   -20,  10,   0,      0,      0,
                                -10,  30,   -30, 10,    0,      0,
                                5,      -20,  30,  -20,  5,      0,
                                -1,     5,     -10, 10,   -5,     1;

    startpoint.m_pos=0.0;
    startpoint.m_vel=0.0;
    startpoint.m_acc=0.0;
    endpoint.m_pos=10.0;
    endpoint.m_vel=0.0;
    endpoint.m_acc=0.0;
    for (int i = 0; i < 3; i++)
    {
        Trim_Point midlepoint;
        midlepoint.m_pos=2*i+1;
        std::cout<<"midleposition = "<<midlepoint.m_pos<<std::endl;
        midlepoint.m_vel=0;
        midlepoint.m_acc=0;
        v_midlepoints.push_back(midlepoint);
    }
    snap.CalSingleQ(time,startpoint,v_midlepoints,endpoint,Q);
    snap.CalSingleM(time,startpoint,v_midlepoints,endpoint,M);
    snap.ConstructSingleC(startpoint,v_midlepoints,endpoint,C);
    snap.GainSingleKnownState(startpoint,v_midlepoints,endpoint,dF);
    snap.GainSingleUnknownState(startpoint,v_midlepoints,endpoint,dP);
    dF_dP.resize((dF.rows()+dP.rows()),1);
    dF_dP.block(0,0,dF.rows(),1)=dF;
    dF_dP.block(dF.rows(),0,dP.rows(),1)=dP;
    
    std::cout<<"Q 的行列是 "<<Q.rows()<<" "<<Q.cols()<<std::endl;
    std::cout<<"Q = "<<std::endl;
    std::cout<<Q<<std::endl;
    // std::cout<<"M 的行列是 "<<M.rows()<<" "<<M.cols()<<std::endl;
    // std::cout<<"M = "<<std::endl;
    // std::cout<<M<<std::endl;
    // std::cout<<"C 的行列是 "<<C.rows()<<" "<<C.cols()<<std::endl;
    // std::cout<<"C's transpose = "<<std::endl;
    // std::cout<<C.transpose()<<std::endl;
    // std::cout<<"dF 的行列是 "<<dF.rows()<<" "<<dF.cols()<<std::endl;
    // std::cout<<"dF = "<<std::endl;
    // std::cout<<dF<<std::endl;
    // std::cout<<"dP 的行列是 "<<dP.rows()<<" "<<dP.cols()<<std::endl;
    // std::cout<<"dP = "<<std::endl;
    // std::cout<<dP<<std::endl;
    // std::cout<<"dF_dP 的行列是 "<<dF_dP.rows()<<" "<<dF_dP.cols()<<std::endl;
    // std::cout<<"dF_dP = "<<std::endl;
    // std::cout<<dF_dP<<std::endl;
    // std::cout<<"--------------------------------------"<<std::endl;
    // std::cout<<C.transpose()*dF_dP<<std::endl;
    Eigen::MatrixXd R=C * M.inverse().transpose() * Q * M.inverse() * C.transpose();
    // std::cout<<"R = "<<std::endl;
    // std::cout<<R<<std::endl;

    Eigen::MatrixXd R_pp=R.block(dF.rows(),dF.rows(),dP.rows(),dP.rows());
    Eigen::MatrixXd R_fp=R.block(0,dF.rows(),dF.rows(),dP.rows());

    Eigen::JacobiSVD<Eigen::MatrixXd> svd_pp(R_pp,Eigen::ComputeFullU|Eigen::ComputeFullV);
    int SingularMatSize=svd_pp.singularValues().size();
    Eigen::MatrixXd singularvalue=Eigen::MatrixXd::Zero(SingularMatSize,SingularMatSize);
    for (int i = 0; i < SingularMatSize; i++)
    {
        if(svd_pp.singularValues()(i) != 0)
        {
            singularvalue(i,i)=1/svd_pp.singularValues()(i);
        }
        else
        {
            singularvalue(i,i)=svd_pp.singularValues()(i);
        }
    }
    Eigen::MatrixXd S_inv=Eigen::MatrixXd::Zero(R_pp.cols(),R_pp.rows());
    S_inv.block(0,0,SingularMatSize,SingularMatSize)=singularvalue;
    Eigen::MatrixXd R_pp_inv=svd_pp.matrixV() * S_inv * svd_pp.matrixU().transpose();

    Eigen::MatrixXd dP_start=-1 * R_pp_inv * R_fp.transpose() * dF;
    dF_dP.block(dF.rows(),0,dP_start.rows(),dP_start.cols()) = dP_start;
    // std::cout<<"dP_Start = "<<std::endl;
    // std::cout<<dP_start<<std::endl;
    // std::cout<<"R_pp'sinverse = "<<std::endl;
    // std::cout<<R_pp_inv<<std::endl;
    // std::cout<<"R_fp = "<<std::endl;
    // std::cout<<R_fp<<std::endl;

    Eigen::MatrixXd coeffP=M.inverse()*C.transpose()*dF_dP;
    Eigen::Matrix<double,6,1> C_Coeff = Eigen::Matrix<double,6,1>::Zero();
    Eigen::Matrix<double,6,1> P_Coeff = Eigen::Matrix<double,6,1>::Zero();    

    int j=0;
    int k=0;
    coeff.resize((v_midlepoints.size()+1));
    curvesNums = coeff.size();
    for (int i = 0; i < (v_midlepoints.size()+1); i++)
    {
        P_Coeff << coeffP(j,0), coeffP(j+1,0), coeffP(j+2,0), coeffP(j+3), coeffP(j+4,0), coeffP(j+5);
        C_Coeff = Coeff_P2C.inverse() * P_Coeff;
        std::cout<<"第 "<<(i+1)<<" 段曲线系数是："<<coeffP(j,0)<<" "<<coeffP(j+1,0)<<" "<<coeffP(j+2,0)<<" "<<coeffP(j+3,0)<<" "<<coeffP(j+4,0)<<" "<<coeffP(j+5,0)<<std::endl;
	    // 获得贝塞尔曲线的系数
        // std::cout<<"第 "<<(i+1)<<" 段贝塞尔曲线系数是："<<C_Coeff(k,0)<<" "<<C_Coeff(k+1,0)<<" "<<C_Coeff(k+2,0)<<" "<<C_Coeff(k+3,0)<<" "<<C_Coeff(k+4,0)<<" "<<C_Coeff(k+5,0)<<std::endl;
        coeff[i].c1 = coeffP(j,0);
        coeff[i].c2 = coeffP(j+1,0);
        coeff[i].c3 = coeffP(j+2,0);
        coeff[i].c4 = coeffP(j+3,0);
        coeff[i].c5 = coeffP(j+4,0);
        coeff[i].c6 = coeffP(j+5,0);
        
        j=j+6;
    }
}

void test02()
{
    OptimizeOoqp op;
    // op.Test_Solve();
    // op.SingleSolveCurve(0.5,6,0,0,0,0,1,0,0);
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    Eigen::VectorXd midle;
    midle.resize(3);

    begin.setZero();
    end.setZero();

    begin(0) = 0;
    end(0) = 5;
    midle(0) = -1;
    midle(1) = 2;
    midle(2) = 3;

    op.SolveCurve(10,6,0,begin,end,midle);
}

void test03()
{
    OptimizeOoqp op;
    Eigen::Vector3d begin;
    Eigen::Vector3d end;
    Eigen::VectorXd midle;
    Eigen::VectorXd T;
    midle.resize(3);
    T.resize(2);

    begin.setZero();
    end.setZero();

    begin(0) = 0;
    end(0) = 5;
    midle(0) = 1;
    // midle(1) = 2;
    // midle(2) = 3;
    T(0) = 1;
    T(1) = 2;
    // T(2) = 1;
    // T(3) = 1;
    op.ConEqualSolve(T,2,begin,end,midle);
}

int main(int argc, char *argv[])
{
    ros::init(argc,argv,"pub_coeff");
    ros::NodeHandle nh;

    std::vector< optimize_trim::eigen2ros > cof;
    std_msgs::Float64MultiArray curve_arr;
    test01(cof);
    ros::Publisher pub_coeff = nh.advertise< optimize_trim::eigen2ros >("curve_coeff",100);
    ros::Publisher pub_coeff_arr = nh.advertise<std_msgs::Float64MultiArray>("curve_coef_arr",100);
    
    // for(int n = 0;n < cof.size();n++)
    // {
    //     std::cout<<cof[n].c1<<" "<<cof[n].c2<<" "<<cof[n].c3<<" "<<cof[n].c4<<" "<<cof[n].c5<<" "<<cof[n].c6<<endl;
    // }
    // std::cout<<endl;

    ros::Rate r(20);
    optimize_trim::eigen2ros tmp;
    int j = 0;
    curve_arr.data.resize(6 * cof.size());
    for(int i = 0;i < (6 * cof.size());i=i+6)
    {
        curve_arr.data[i] = cof[j].c1;
        curve_arr.data[i+1] = cof[j].c2;
        curve_arr.data[i+2] = cof[j].c3;
        curve_arr.data[i+3] = cof[j].c4;
        curve_arr.data[i+4] = cof[j].c5;
        curve_arr.data[i+5] = cof[j].c6;
        // std::cout<<cof[j].c1<<" "<<cof[j].c2<<" "<<cof[j].c3<<" "<<cof[j].c4<<" "<<cof[j].c5<<" "<<cof[j].c6<<std::endl;
        if(j < (cof.size() - 1))
        {
            j = j + 1;
        }
    }
    // std::cout<<"arr's size is "<<curve_arr.data.size()<<std::endl;

    int count = 0;

    while (ros::ok())
    {
        // std::cout<<"count = "<<count<<std::endl;
        // ROS_INFO("count = %.2f",count);
        // tmp = cof[count];
        // pub_coeff.publish(tmp);
        // count = count + 1;
        // if(count != 0 && count % curvesNums == 0)
        // {
        //     count = 0;
        // }

        pub_coeff_arr.publish(curve_arr);
        r.sleep();
        ros::spinOnce();
    }
    return 0;
}
