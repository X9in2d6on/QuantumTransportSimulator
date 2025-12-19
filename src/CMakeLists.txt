cmake_minimum_required(VERSION 3.10)
project(NanoMOSFET_Simulator CXX)

# C++标准设置
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 设置编译类型为Debug模式
set(CMAKE_BUILD_TYPE Debug)

# 编译标志 - 添加调试信息
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# 启用OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Eigen库路径（根据实际情况修改）
include_directories(/usr/include/eigen3)

# 添加可执行文件
add_executable(device_demo devicesetupdemo.cpp)
add_executable(iv_generator ivcurvegenerator.cpp)

# 链接库
target_link_libraries(device_demo OpenMP::OpenMP_CXX)
target_link_libraries(iv_generator OpenMP::OpenMP_CXX)
