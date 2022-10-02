#include <ros/ros.h>
#include <std_msgs/Float64MultiArray.h>
#include <Eigen/Dense>
#include "/home/zjs/minimum_snap_ws/devel/include/optimize_trim/eigen2ros.h"

std::vector<Eigen::Matrix<double,6,1>> curve_coeficient;
std::vector< std::vector<double> > xaixs;
std::vector<double> curve_time;
std_msgs::Float64MultiArray test_points;
optimize_trim::eigen2ros curve_info;

void Show_Coef_Arr(const std_msgs::Float64MultiArray::ConstPtr& arr)
{
    int j = 0;
    int arr_size = arr->data.size();

    test_points.data.resize(51);
    xaixs.resize(arr_size);
    curve_coeficient.resize(arr_size/6);
    // curve_time.resize(arr_size/6);
    curve_time.resize(4);
    curve_time[0]=0.5;
    curve_time[1]=1;
    curve_time[2]=1;
    curve_time[3]=2.5;


    for(int i = 0;i < arr_size/6;i++)
    {
        curve_coeficient[i](0,0) = arr->data[j];
        curve_coeficient[i](1,0) = arr->data[j+1];
        curve_coeficient[i](2,0) = arr->data[j+2];
        curve_coeficient[i](3,0) = arr->data[j+3];
        curve_coeficient[i](4,0) = arr->data[j+4];
        curve_coeficient[i](5,0) = arr->data[j+5];
        j = j + 6;
    }

    for(int i = 0;i < arr_size / 6;i++)
    {
        int num = 0;
        xaixs[i].resize(100 * curve_time[i] + 1);

        // std::cout<<curve_coeficient[i](0,0)<<" "<<curve_coeficient[i](1,0)<<" "<<curve_coeficient[i](2,0)<<" "
        //                     <<curve_coeficient[i](3,0)<<" "<<curve_coeficient[i](4,0)<<" "<<curve_coeficient[i](5,0)<<" ";
        // std::cout<<std::endl;

        for(int j = 0;j <= (100 * curve_time[i]);j = j + 1)
        {
            double d_j = j / 100.f;
            curve_info.time = d_j;
            curve_info.xaixs = curve_coeficient[i](0,0) * std::pow(curve_info.time,5) + curve_coeficient[i](1,0) * std::pow(curve_info.time,4)
                                                  + curve_coeficient[i](2,0) * std::pow(curve_info.time,3) + curve_coeficient[i](3,0) * std::pow(curve_info.time,2)
                                                  + curve_coeficient[i](4,0) * std::pow(curve_info.time,1) + curve_coeficient[i](5,0);
            xaixs[i][num] = curve_info.xaixs;

            if (num == (100 * curve_time[i]))
            {
                std::cout<<"Position's section = "<<"["<<xaixs[i][0]<<","<<xaixs[i][num]<<"]"<<" num = "<<num<<std::endl;
                num = 0;
            }

            if(i == 0)
            {
                test_points.data = xaixs[i];
                // for(int cc=0;cc<test_points.data.size();cc++)
                // {
                //     std::cout<<test_points.data[cc]<<" "<<std::endl;
                // }
            }
            num = num + 1;
        }
        std::cout<<std::endl;
    }
    std::cout<<std::endl;
}

int main(int argc,char * argv[])
{
    ros::init(argc,argv,"sub_coefficient");
    ros::NodeHandle nh;
    
    ros::Subscriber sub_coef_arr = nh.subscribe<std_msgs::Float64MultiArray>("curve_coef_arr",100,Show_Coef_Arr);
    ros::Publisher pub_test_points = nh.advertise<std_msgs::Float64MultiArray>("test_points_fifty",100);
    ros::Rate r(20);

    while (ros::ok())
    {
        pub_test_points.publish(test_points);
        ros::spinOnce();
        r.sleep();
    }
    
    return 0;
}