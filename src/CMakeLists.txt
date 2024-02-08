
add_library(actionData actionData.cpp actionData.h)
add_library(fields fields.cpp fields.h)
add_library(multiroot multiroot.cpp multiroot.h)
add_library(speGrid speGrid.cpp speGrid.h)
add_library(speEqns speEqns.cpp speEqns.h)

add_executable(NSDI_Pulse NSDI_Pulse.cpp)


target_link_libraries(NSDI_Pulse PRIVATE actionData speGrid speEqns fields multiroot
						GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX)
