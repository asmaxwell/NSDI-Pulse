cmake_minimum_required(VERSION 2.8...3.13)

project (NSDI_Pulse
VERSION 0.0
LANGUAGES CXX)

add_library(speGrid speGrid.cpp speGrid.h)
add_library(speEqns speEqns.cpp speEqns.h)
add_library(fields fields.cpp fields.h)

add_executable(NSDI_Pulse NSDI_Pulse.cpp)

target_link_libraries(NSDI_Pulse speGrid speEqns fields)
