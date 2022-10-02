cmake_minimum_required(VERSION 3.0.2)
project(optimize_trim)

find_package(catkin REQUIRED COMPONENTS
  roscpp
  rospy
  std_msgs
  message_generation
)

add_message_files(
  FILES
  eigen2ros.msg
)

generate_messages(
  DEPENDENCIES
  std_msgs
)

find_package(Eigen3 REQUIRED)
find_package(OpenCV 3 REQUIRED)
find_package(Pangolin REQUIRED)

include_directories(
  ${catkin_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIR}
  ${OpenCV_INCLUDE_DIRS}
  ${Pangolin_INCLUDE_DIRS}
)

include_directories("/home/zjs/OOQP_Solver/OOQP/include")
include_directories("/home/zjs/OOQP_Solver/OOQP/lib")

catkin_package(
  #  LIBRARIES optimize_trim
  CATKIN_DEPENDS roscpp rospy std_msgs message_runtime
  DEPENDS EIGEN3
  INCLUDE_DIRS 
)

add_library(trim_point src/trim_point.cpp)
target_link_libraries(trim_point ${catkin_LIBRARIES})

add_library(bezier_curve src/bezier_curve.cpp)
target_link_libraries(bezier_curve ${catkin_LIBRARIES})

add_library(OptimizeOoqp src/OoqpOptimize.cpp)
target_link_libraries(OptimizeOoqp ${catkin_LIBRARIES})

add_executable(optimize_trim_file src/optimize_trim.cpp)
target_link_libraries(optimize_trim_file trim_point  
                                                                                   bezier_curve
                                                                                   ${catkin_LIBRARIES} 
                                                                                   ${OpenCV_LIBS} 
                                                                                   ${Pangolin_LIBRARIES}
                                                                                   OptimizeOoqp
                                                                                   ooqpgensparse
                                                                                   ooqpsparse
                                                                                   ooqpgondzio
                                                                                   ooqpbase 
                                                                                   blas 
                                                                                   ma27 
                                                                                   gfortran)
add_executable(coefficient src/subCoeff.cpp)
target_link_libraries(coefficient                   trim_point  
                                                                                    bezier_curve
                                                                                    ${catkin_LIBRARIES} 
                                                                                    ${OpenCV_LIBS} 
                                                                                    ${Pangolin_LIBRARIES}
                                                                                    OptimizeOoqp
                                                                                    ooqpgensparse
                                                                                    ooqpsparse
                                                                                    ooqpgondzio
                                                                                    ooqpbase 
                                                                                    blas 
                                                                                    ma27 
                                                                                    gfortran)

add_executable(test_ooqp src/test_ooqp.cpp)
target_link_libraries(test_ooqp
                    ooqpgensparse
                    ooqpsparse
                    ooqpgondzio
                    ooqpbase 
                    blas 
                    ma27 
                    gfortran
)