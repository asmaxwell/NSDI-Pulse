# NSDI Pulse
## Introduction
This is an extension of a previous NSDI code using the strong-field approximation that was valid for monochromatic fields. This code is valid for a sin-squared laser pulse.
## Prerequisites 
To install this code you will need the following libraries installed:
* CMake (to compile)
* GSL (GNU Scientific Library)
* OpenMP (parallelization library)
* Catch2 (if you want to use the test features---if not you can simply not build the text executable)

## Install and compile
First download the code from GitHub however you prefer. Then compile using CMake. If you do  this over the command line go to the root directory and run the following:

```cmake -S . -B build```

This will make a build directory and configure the CMake file. Then to compile and generate an executable run:

```cmake --build build```

However, it may be easiest to import the project as a CMake project into Eclipse, Visual Studio Code, or another C++-friendly IDE.

## Structure of code
The ```main``` function is located in the file ```src/NSDI_Pulse.cpp```, which is where the code begins its execution, this file is what makes the final executable. Here, we have code that makes a Data/ folder for export and sets some parameters. a parameter file could be used instead but due to the small size of this code, I didn't implement that yet. One of the parameters set is the field type and laser field class, this can be switched between monochromatic and sin2 (sin-squared).

An important class is ```speGrid```, which defines a grid of saddle point solutions. Here ```saddlePointGrid``` is declared, this is where the saddle point equations are solved. The saddle point equations are solved over the entire grid by first picking a single point in momentum, and then solving this point randomly. What this means is many (currently 2000) random complex times are generated for the single point in momentum. These times are used as guesses for a rootfinding method to solve the saddle point equations. This will exhaustively provide all solutions for a specific momentum. This is a slow way to find solutions but we only need to do it for a single point. After that the solutions for that one point are used as guesses for all the adjacent points, and this can be used iteratively to propagate the grid with all the solutions. The number of solutions can vary from point to point, this means you should choose a starting point that has many solutions so you don't miss any.

After this, the computation of the Action is straightforward, as it is just evaluating functions with the saddle point solutions.
