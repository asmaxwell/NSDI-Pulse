# NSDI Pulse
## Introduction
This is an extension of a previous NSDI code using the strong-field approximation that was valid for monochromatic fields. This code is valid for a sin-squared laser pulse.
## Prerequisites 
In order to install this code you will need the following libraires installed:
* CMake (to compile)
* GSL (GNU Scientific Library)
* OpenMP (parallelization labrary)
* Catch2 (if you want to use the test features---if not you can simply comment out all Catch2 references from the code)

## Install and compile
First download the code from Github however you perfer. Then compile using CMake. If you do  this over the command line go to the root directory and run the following:

```cmake -S . -B build```

This will make a build directory and configure the CMake file. Then to compile and gerneate an executable run:

```cmake --build build```

However, it may be easiet to import the project as a CMake project into eclispe, virtual studio, or another C++ friendly IDE.

## Structure of code
The ```main``` function is located in the file NSDI_Pulse.cpp, which is where the code begins it's execution, this file is what makes the final executable. Here, we have code that makes a Data/ folder for export and set some parameters. a parameter file could be used instead but due to the small size of this code I didn't implement that yet. One of the paraemters set is the field type and lasr field class, this can be swtich between monochromatic and sin2 (sin-squared).

An important class is ```speGrid```, this defines a grid of saddle point solutions. Here ```saddlePointGrid``` is delared, this is where the saddle point equations are sovled. The saddle point equations are sovled over the entire grid by first picking a sinlge point in momentum, then solving this point randomly. What this means is many (currently 2000) random complex times are generated for the single point in momentum. These times are used as guesses for a rootfinding method to solve the saddle point equations. This will exhaustively provide all solutions for a specific momentum. This is a slow way to find solutions but we only need to do it for a single point. After that the solutions for that one point are used as guess for all the adjacent points, and this can be used iteratively to propagate the grid with all the solutions. The number of solutions can vary from point-to-point, this means you should choose a starting point that has many solutions so you don't miss any.

After this the computation of the Action is straight forward, as it is just evaluating funcitions with the saddle point solutions.
