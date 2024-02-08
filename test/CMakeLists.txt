#can set catch2 to download for testing
FetchContent_Declare(
  catch
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v2.13.6)
FetchContent_MakeAvailable(catch)


#find_package(Catch2 3 REQUIRED)
file(GLOB TEST_SOURCES "*.cpp")

add_executable(test_NSDI_Pulse ${TEST_SOURCES})

target_link_libraries(test_NSDI_Pulse PRIVATE actionData speGrid speEqns fields multiroot
						 Catch2::Catch2 GSL::gsl GSL::gslcblas OpenMP::OpenMP_CXX
)

#include(CTest)
#include(Catch)
#catch_discover_tests(test_NSDI_Pulse)
enable_testing()
add_test(NAME test_NSDI_Pulse COMMAND test_NSDI_Pulse) # Command can be a target
