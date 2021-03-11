cmake_minimum_required(VERSION 3.10)
project(Rasterizer)

set(CMAKE_CXX_STANDARD 17)

include_directories(/usr/local/include)

add_executable(Rasterizer main.cpp tgaimage.cpp tgaimage.h model.h model.cpp geometry.h)
