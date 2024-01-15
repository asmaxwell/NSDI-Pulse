cmake_minimum_required( VERSION 3.10 )

project (NSDI_Pulse
VERSION 0.0
LANGUAGES CXX)

#For testing include catch2/3
find_package(Catch2 3 REQUIRED)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

add_library(actionData actionData.cpp actionData.h)
add_library(fields fields.cpp fields.h)
add_library(multiroot multiroot.cpp multiroot.h)
add_library(speGrid speGrid.cpp speGrid.h)
add_library(speEqns speEqns.cpp speEqns.h)


add_executable(NSDI_Pulse NSDI_Pulse.cpp)
target_link_libraries(NSDI_Pulse PRIVATE actionData speGrid speEqns fields multiroot
						 Catch2::Catch2 GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)
