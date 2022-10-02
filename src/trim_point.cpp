#include "/home/zjs/minimum_snap_ws/src/optimize_trim/include/optimize_trim/trim_point.h"

void Trim_Point::SetPosition(Eigen::Vector3d position)
{
    m_position=position;
}

void Trim_Point::SetVelocity(Eigen::Vector3d velocity)
{
    m_velocity=velocity;
}

void Trim_Point::SetAcceleration(Eigen::Vector3d acceleration)
{
    m_acceleration=acceleration;
}

void Trim_Point::SetJerk(Eigen::Vector3d jerk)
{
    m_jerk=jerk;
}

void Minimum_Snap::CalSingleQ(double sum_time,Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& QMat)
{
    int num=0;
    std::vector< Eigen::Matrix<double,6,6> > Q;
    Eigen::Matrix<double,6,6> singleQ=Eigen::Matrix<double,6,6>::Zero();
    Eigen::MatrixXd MatQ = Eigen::MatrixXd::Zero(6*(1+midlep.size()),6*(1+midlep.size()));
    Eigen::Matrix2d subQ_right_bottom;
    if(std::abs((endp.m_pos - startp.m_pos))==0)
    {
        std::cout<<"请稍微调整起始点和终止点的位置！！！"<<std::endl;
        return;
    }
    for (auto it = midlep.begin() ; it != midlep.end(); it++)
    {
        if(std::abs(((*it).m_pos - startp.m_pos))==0 | (std::abs(endp.m_pos - midlep.back().m_pos))==0)
        {
            if(((*it).m_pos - startp.m_pos)==0)
            {
                std::cout<<"起始点和第二个点位置重合"<<std::endl;
            }
            if((std::abs(endp.m_pos - midlep.back().m_pos))==0)
            {
                std::cout<<"末尾点和倒数第二个点位置重合"<<std::endl;
            }
        }
        if(it == midlep.begin())
        {
            double dt;
            double portion=0;
            Eigen::Matrix<double,6,6> Q1=Eigen::Matrix<double,6,6>::Zero();
            Eigen::Matrix2d subQ1_right_bottom;
            portion = std::abs(((*it).m_pos - startp.m_pos)) / std::abs((endp.m_pos - startp.m_pos));
            //std::cout<<"time = "<<sum_time<<" "<<"portion = "<<portion<<std::endl;
            dt = sum_time * portion;
            subQ1_right_bottom << 576*(std::pow(dt,1)) , 1440 * (std::pow(dt,2)),
                                                            1440 * (std::pow(dt,2)) , 4800*(std::pow(dt,3));
            Q1.block<2,2>(4,4) = subQ1_right_bottom;
            std::cout<<"dt = "<<dt<<std::endl;
            // std::cout<<"Q1 = "<<Q1<<std::endl;
            Q.push_back(Q1);            
        }
        else
        {
            double dt;
            double portion=0;
            Eigen::Matrix<double,6,6> Q2=Eigen::Matrix<double,6,6>::Zero();
            Eigen::Matrix2d subQ2_right_bottom;
            portion = (std::abs((*it).m_pos - (*(it-1)).m_pos)) / std::abs((endp.m_pos - startp.m_pos));
            dt = sum_time * portion;
            std::cout<<"midle dt = "<<dt<<std::endl;

            subQ2_right_bottom << 576*(std::pow(dt,1)) , 1440 * (std::pow(dt,2)),
                                                            1440 * (std::pow(dt,2)) , 4800*(std::pow(dt,3));
            Q2.block<2,2>(4,4) = subQ2_right_bottom;
            Q.push_back(Q2);  
        }
        num++;
    }
    // std::cout<<"与容器大小相比的 num 为 "<<num<<std::endl;
    // std::cout<<"存放中间点的容器的大小为 "<<midlep.size()<<std::endl;
    if(num == midlep.size())
    {
        double dt;
        double portion=0;
        Eigen::Matrix<double,6,6> Q3=Eigen::Matrix<double,6,6>::Zero();
        Eigen::Matrix2d subQ3_right_bottom;
        portion = (std::abs(endp.m_pos - midlep.back().m_pos)) / std::abs((endp.m_pos - startp.m_pos));
        dt = sum_time * portion;
        std::cout<<"dt = "<<dt<<std::endl;

        subQ3_right_bottom << 576*(std::pow(dt,1)) , 1440 * (std::pow(dt,2)),
                                                            1440 * (std::pow(dt,2)) , 4800*(std::pow(dt,3));
        Q3.block<2,2>(4,4) = subQ3_right_bottom;
        // std::cout<<"Q right bottom = "<<std::endl;
        // std::cout<<Q3<<std::endl;
        Q.push_back(Q3);  
    }
    for (int i = 0; i < Q.size(); i++)
    {
        MatQ.block<6,6>(6*i,6*i)=Q[i];
    }
    QMat=MatQ;
}

void Minimum_Snap::CalSingleM(double sum_time,Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& MMat)
{
    int num=0;
    Eigen::MatrixXd MatM=Eigen::MatrixXd::Zero(6*(1+midlep.size()),6*(1+midlep.size()));
    std::vector<Eigen::Matrix<double,6,6>> M;
    for (auto it = midlep.begin() ; it != midlep.end(); it++)
    {
        if(it == midlep.begin())
        {
            double dt;
            double portion=0;
            Eigen::Matrix<double,6,6>  M1=Eigen::Matrix<double,6,6>::Zero();
            portion = std::abs(((*it).m_pos - startp.m_pos)) / std::abs((endp.m_pos - startp.m_pos));
            dt = sum_time * portion;
            M1(0,5)=1;
            M1(1,4)=1;
            M1(2,3)=2;
            M1(3,0)=std::pow(dt,5);
            M1(3,1)=std::pow(dt,4);
            M1(3,2)=std::pow(dt,3);
            M1(3,3)=std::pow(dt,2);
            M1(3,4)=std::pow(dt,1);
            M1(3,5)=1;
            M1(4,0)=5*std::pow(dt,4);
            M1(4,1)=4*std::pow(dt,3);
            M1(4,2)=3*std::pow(dt,2);
            M1(4,3)=2*std::pow(dt,1);
            M1(4,4)=1;
            M1(5,0)=20*std::pow(dt,3);
            M1(5,1)=12*std::pow(dt,2);
            M1(5,2)=6*std::pow(dt,1);
            M1(5,3)=2;
            M.push_back(M1);
        }
        else
        {
            double dt;
            double portion=0;
            Eigen::Matrix<double,6,6>  M2=Eigen::Matrix<double,6,6>::Zero();
            portion = (std::abs((*it).m_pos - (*(it-1)).m_pos)) / std::abs((endp.m_pos - startp.m_pos));
            dt = sum_time * portion;
            M2(0,5)=1;
            M2(1,4)=1;
            M2(2,3)=2;
            M2(3,0)=std::pow(dt,5);
            M2(3,1)=std::pow(dt,4);
            M2(3,2)=std::pow(dt,3);
            M2(3,3)=std::pow(dt,2);
            M2(3,4)=std::pow(dt,1);
            M2(3,5)=1;
            M2(4,0)=5*std::pow(dt,4);
            M2(4,1)=4*std::pow(dt,3);
            M2(4,2)=3*std::pow(dt,2);
            M2(4,3)=2*std::pow(dt,1);
            M2(4,4)=1;
            M2(5,0)=20*std::pow(dt,3);
            M2(5,1)=12*std::pow(dt,2);
            M2(5,2)=6*std::pow(dt,1);
            M2(5,3)=2;
            M.push_back(M2);
        }
        num++;
    }
    if(num == midlep.size())
    {
        double dt;
        double portion=0;
        Eigen::Matrix<double,6,6> M3=Eigen::Matrix<double,6,6>::Zero();
        portion = (std::abs(endp.m_pos - midlep.back().m_pos)) / std::abs((endp.m_pos - startp.m_pos));
        dt = sum_time * portion;
        M3(0,5)=1;
        M3(1,4)=1;
        M3(2,3)=2;
        M3(3,0)=std::pow(dt,5);
        M3(3,1)=std::pow(dt,4);
        M3(3,2)=std::pow(dt,3);
        M3(3,3)=std::pow(dt,2);
        M3(3,4)=std::pow(dt,1);
        M3(3,5)=1;
        M3(4,0)=5*std::pow(dt,4);
        M3(4,1)=4*std::pow(dt,3);
        M3(4,2)=3*std::pow(dt,2);
        M3(4,3)=2*std::pow(dt,1);
        M3(4,4)=1;
        M3(5,0)=20*std::pow(dt,3);
        M3(5,1)=12*std::pow(dt,2);
        M3(5,2)=6*std::pow(dt,1);
        M3(5,3)=2;
        M.push_back(M3);
    }
    for (int i = 0; i < M.size(); i++)
    {
        MatM.block<6,6>(6*i,6*i)=M[i];
    }
    MMat=MatM;
}

void Minimum_Snap::ConstructSingleC(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& CMat)
{
    int MatCRows=0;
    int MatCCols=0;
    int segements=1+midlep.size();
    Eigen::MatrixXd MatC;

    MatCRows=6*(1+midlep.size());
    MatCCols=3*(2+midlep.size());
    MatC.resize(MatCRows,MatCCols);
    MatC.setZero();

    Eigen::MatrixXd MatCStart;
    Eigen::MatrixXd MatCStartTmp;
    MatCStart.resize(6,MatCCols);
    MatCStart.setZero();
    MatCStart(0,0)=1;
    MatCStart(1,1)=1;
    MatCStart(2,2)=1;
    MatCStart(3,3)=1;
    MatCStart(4,segements+6-1)=1;
    MatCStart(5,segements+7-1)=1;
    MatC.block(0,0,6,MatCCols)=MatCStart;
    MatCStartTmp=MatCStart;

    Eigen::MatrixXd CTmp;
    Eigen::MatrixXd CLast;
    CTmp.resize(6,MatCCols);
    CLast.resize(3,MatCCols);
    CTmp.setZero();
    CLast.setZero();
    //std::cout<<"in function's segements = "<<segements<<std::endl;
    for (int i = 0; i < (segements - 2); i++)
    {
        if(i==0)
        {
            //上一个矩阵一半赋值给下一个
            CTmp.block(0,0,3,MatCCols)=MatCStartTmp.block(3,0,3,MatCCols);
            CTmp(3,5+i-1)=1;
            CTmp(4,segements+8-1+2*i)=1;
            CTmp(5,segements+9-1+2*i)=1;
            CLast=CTmp.block(3,0,3,MatCCols);
            MatC.block(6*(i+1),0,6,MatCCols)=CTmp;
        }
        else
        {
            CTmp.setZero();
            CTmp.block(0,0,3,MatCCols)=CLast;
            // std::cout<<"in function CLast = "<<std::endl;
            // std::cout<<CLast<<std::endl;
            CTmp(3,5+i-1)=1;
            CTmp(4,segements+8-1+2*i)=1;
            CTmp(5,segements+9-1+2*i)=1;
            // std::cout<<"in function CLastnext = "<<std::endl;
            // std::cout<<CTmp.block(3,0,3,MatCCols)<<std::endl;
            CLast=CTmp.block(3,0,3,MatCCols);
            MatC.block(6*(i+1),0,6,MatCCols)=CTmp;
        }
    }

    Eigen::MatrixXd MatCEnd;
    MatCEnd.resize(6,MatCCols);
    MatCEnd.setZero();
    MatCEnd.block(0,0,3,MatCCols)=CTmp.block(3,0,3,MatCCols);
    MatCEnd(3,segements+2)=1;
    MatCEnd(4,segements+3)=1;
    MatCEnd(5,segements+4)=1;
    MatC.bottomLeftCorner(6,MatCCols)=MatCEnd;

    CMat=MatC;
    // std::cout<<"CMat = "<<std::endl;
    // std::cout<<CMat<<std::endl;
    CMat=MatC.transpose();
}

void Minimum_Snap::GainSingleKnownState(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& dFMat)
{
    Eigen::MatrixXd dFTmp;
    dFTmp.resize(6+midlep.size(),1);
    
    dFTmp(0,0)=startp.m_pos;
    dFTmp(1,0)=startp.m_vel;
    dFTmp(2,0)=startp.m_acc;

    int i=0;
    for ( i ; i < midlep.size(); i++)
    {
        dFTmp(i+3,0)=midlep[i].m_pos;
    }

    dFTmp(i+3,0)=endp.m_pos;
    dFTmp(i+4,0)=endp.m_vel;
    dFTmp(i+5,0)=endp.m_acc;
    dFMat=dFTmp;
}

void Minimum_Snap::GainSingleUnknownState(Trim_Point& startp,std::vector<Trim_Point>& midlep,Trim_Point& endp,Eigen::MatrixXd& dPMat)
{
    Eigen::MatrixXd dPTmp;
    dPTmp.resize(2*midlep.size(),1);
    int j=0;
    for (int i = 0; i < midlep.size(); i++)
    {
        dPTmp(j,0)=midlep[i].m_vel;
        dPTmp(j+1,0)=midlep[i].m_acc;
        j=j+2;
    }
    dPMat=dPTmp;
}
