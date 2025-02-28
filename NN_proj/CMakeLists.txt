cmake_minimum_required(VERSION 3.5)

project(NN_proj LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(SRCS main.cpp matrix.cpp nn.cpp)
set(HDRS matrix.h nn.h)
add_executable(NN_proj ${SRCS} ${HDRS})
