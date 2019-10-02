cmake_minimum_required(VERSION 3.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(ExternalProject)

project( AprilTools )

find_package( OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

add_subdirectory(apriltag)

add_executable( apriltools AprilTools.cc )
target_link_libraries( apriltools ${OpenCV_LIBS} apriltag )
