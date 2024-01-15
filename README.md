# NSDI Pulse
## Introduction
This is an extension of a previous NSDI code using the strong-field approximation that was valid for monochromatic fields. This code is valid for a sin-squared laser pulse.
## Prerequisite 
In order to install this code you will need the following libraires installed:
* CMake (to compile)
* GSL (GNU Scientific Library)
* OpenMP (parallelization labrary)
* Catch2 (if you want to use the test features---if not you can simply comment out all Catch2 references from the code)

## Install and compile
First download the code from Github however you perfer. Then compile using CMake. I fyou do  this over the command line go to the directory build/default/ and run the following:
```cmake --build . --target all```
However, it may be easiet to import the project as a CMake project into eclispe, virtual studio, or another C++ friendly IDE.
