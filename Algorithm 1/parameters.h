// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD1.
//
// MOQCD1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD1.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <vector>
#include <iostream>
#include <limits>

class Parameters
{
public:
  ~Parameters();

  // Input
  void usage() const;
  void processCommandLine(int argc, char* argv[]);

  // Problem
  std::string problem() const;
  int gateSet() const;
  int numQbits() const;

  // GA parameters
  int populationSize() const;
  double crossoverRate() const;
  double mutationRate() const;
  double spacing() const;

  // Circuit length parameters
  double meanRandomLength() const;  // Mean length of a randomly generated circuit. (Exponential distribution.)
  int maxCircuitLength() const;

  // Evaluation parameters
  double targetError() const;      // A circuit of cost c is not evaluated if we have already found one with...
  long expectedReduction() const;  // ...error < targetError and cost <= c - expectedReduction.
  bool randomizeGateParameters() const;  // True => gate parameters are randomized before numerical optimization.
  int iterationQuota() const;  // Numerical optimization proceeds only it is expected to improve the front within a...
                               // ...quota of iterations.
  // Solution cacheing
  bool useCache() const;

  // Termination criteria and random seed
  int numGenerations() const;
  int numHalts() const;
  bool randomSeedProvided() const;
  int randomSeed() const;

  // Output
  void output(std::ostream& out) const;

private:
  void reset();
  void sanity_check();

private:
  // Program name and command line
  std::string program_name;  // Used only by usage()
  std::string command_line;  // Used in logfile output
  int num_qbits;

  // Problem
  std::string problem_;
  int gate_set = 1;  // Default value.

  // GA parameters
  int population_size = -1;
  double crossover_rate = -1.0;
  double mutation_rate = -1.0;
  double spacing_ = -1.0;

  // Problem parameters
  double mean_random_length = -1.0;
  int max_circuit_length = -1;
  double target_error = 0.00001;  // Default value.
  long expected_reduction = std::numeric_limits<long>::max();  // Default value - all circuit structures are evaluated.
  bool randomize_gate_parameters = false;  // Default value.
  int iteration_quota = 50;  // Default value.
  bool use_cache = false;  // Default value.
  

  // Termination criteria and random seed
  int num_generations = -1;
  int num_halts = 5;  // Default value.
  int random_seed = -1;  // Indicates 'no value' - set the seed randomly.
};

std::ostream& operator<<(std::ostream& out, const Parameters& parameters);

#endif
