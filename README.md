# Multi-Objective Quantum Circuit Discovery

## Two projects - one aim

This repository contains the code for two algorithms for multi-objective quantum circuit discovery. Both use a multi-objective genetic algorithm to optimise circuits, minimising both the complexity of the circuit and the difference between circuit output and desired output. However, each takes a different approach to this task.

### Algorithm 1

If you have been reading our recently submitted paper, then the algorithm described therein is algorithm 1. The emphasis here is on making best use of the circuit evaluations performed. Routines are included for automatic circuit simplification and for full numeric optimisation of circuits with parameterised gates,  as part of the overall optimisation algorithm.

### Algorithm 2

Here the emphasis is on improving the efficiency of circuit evaluation by exploiting the results of the evaluations of previous ‘parent’ circuits. Circuit simplification routines are still used. Numeric optimisation is limited to the optimisation of single gate parameters.

## Prerequisites

The code is written in modern C++ and requires a C++17 compliant compiler. The Boost libraries are also required - specifically those providing boost::lexical_cast, and boost::to_lower. Algorithm 1 also requires Armadillo (used by QIClib) and Eigen (used by LBFGS++).

Code has been tested primarily on a Mac, though it has also been used on Linux. While the current version has not yet been tested on PC, there should be few problems - if issues do arise, please contact me.

### Other libraries

Algorithm 1 uses both QIClib (copyright Titus Chanda) for circuit simulation and LBFGS++ (copyright Yixuan Qiu) for the numerical optimization routine. Both of these libraries required some modification prior to use, so the code is provided in the repository.

## Installation

The contents of the two "Algorithm" folders are already arranged as required for successful compilation. The "Multi-objective GAs" folder is used by both algorithms - place this somewhere on your machine and make sure that the location is added compiler's header search path. Note that, in addition to numerous header files, this folder contains two files, utils.cpp and rng.cpp that will need to to manually added to your project. The "LBFGS++" code, used by algorithm 1, is header only, so just place it somewhere your compiler will be able to find it.

For algorithm 1, the compiler will also need to be told to link the binary to the libarmadillo library and the Accelerate framework.

## Command line parameters

Both algorithms are run from the command line and each takes a number of command line parameters. A brief description of these parameters can be found by running the code with command line parameter "-help". 