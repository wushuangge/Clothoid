cmake_minimum_required(VERSION 3.15)
project(MXClothoid)

set(CMAKE_CXX_STANDARD 11)

add_executable(MXClothoid main.cpp
        Clothoid.cpp Clothoid.h MXClothoid.h test.h)