cmake_minimum_required(VERSION 3.22)
project(251prog5)

set(CMAKE_CXX_STANDARD 14)

add_executable(251prog5 application.cpp testing.cpp dist.cpp dist.h osm.h ../../Downloads/tinyxml2.cpp ../../Downloads/tinyxml2.h ../../Downloads/tinyxml2.h graph.h)
