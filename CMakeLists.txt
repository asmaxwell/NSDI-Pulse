cmake_minimum_required( VERSION 3.18 )

project (NSDI_Pulse
VERSION 0.0
LANGUAGES CXX)

#set c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#support for folders in IDEs
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

include(FetchContent)
find_package(GSL REQUIRED)
find_package(OpenMP REQUIRED)

# The compiled library code is here
add_subdirectory(src)




# Add Testing
include(CTest)
add_subdirectory(test)

