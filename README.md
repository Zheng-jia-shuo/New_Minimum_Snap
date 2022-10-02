# Minimum-snap
# Version: First Version			Date: 2022.10.02
# Author: Zheng Jiashuo
Files in src folder need to be optimized.
In optimizing trajectory curve,trim_point.cpp that could get trajectory points is main file.You should browse it carefully.About OoqpOptimize.cpp which includes Ooqp library,if you want to read,I hope you can learn Ooqp firstly.bezier_curve.cpp include the math method that Bezier curve.subCoeff.cpp's main function is subscribing topic named "curve_coef_arr".optimize_trim.cpp's main function is calculating the coffeicient of curve and publishing topic named "curve_coef_arr".
-------------------------------------------------------------------
Optimized trajectory curve
Download source code.
If you want to use this project,just follow follow the steps below.
1.Create a workspace on your desktop
        a.  cd ~
        b.  mkdir -p minimum_snap_ws/src
        c.  cd minimum_snap_ws/src
        d.  git clone https://github.com/Zheng-jia-shuo/Minimum-snap
        e.  cd ~/minimum_snap_ws
        f.  catkin_make
2.Run executable file,in this section you can create a launch file or enter your folder that locate in minimum_snap_ws/devel/lib/your_package/
        ./your_executable_filename
In your terminal,you could find the result of curve's coefficient optimized by minimum snap.
You need to pay attention to the following points!!!
The file eigen2ros.msg in msg/ is created by yourself.You could get it by browsing my website which is https://blog.csdn.net/er_dan_love/article/details/124794614.In section 2.1,it was introduced how to create an own message type.
