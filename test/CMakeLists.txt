# cmake_minimum_required( VERSION 3.14 )

# project (NSDI_Pulse
# VERSION 0.0
# LANGUAGES CXX)

# #set c++ standard
# set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
# set(CMAKE_CXX_EXTENSIONS OFF)

# can set catch2 to download for testing
# include(FetchContent)
# FetchContent_Declare(
#   catch
#   GIT_REPOSITORY https://github.com/catchorg/Catch2.git
#   GIT_TAG v2.13.6)
# FetchContent_MakeAvailable(catch)


file(GLOB TEST_SOURCES "*.cpp")

add_executable(test_NSDI_Pulse ${TEST_SOURCES})

target_link_libraries(test_NSDI_Pulse PRIVATE actionData speGrid speEqns fields multiroot
						 Catch2::Catch2 GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX
)

include(CTest)
include(Catch)
catch_discover_tests(test_NSDI_Pulse)
enable_testing()
add_test(NAME test_NSDI_Pulse COMMAND test_NSDI_Pulse) # Command can be a target
