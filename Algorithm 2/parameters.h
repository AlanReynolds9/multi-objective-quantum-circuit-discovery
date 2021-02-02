// MOQCD2 (Multi-objective quantum circuit discovery - algorithm 2.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD2.
//
// MOQCD2 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD2 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
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
  int adultPopSize() const;
  int childPopSize() const;
  double spacing() const;
  int numRandSurvivors() const;

  // Problem parameters
  int maxLength() const;  // Circuits larger thatn this will not be evaluated, but will be assigned a bad error value.
  double meanRandomLength() const;  // Mean length of a randomly generated circuit. (Exponential distribution.)
  double badErrorCutoff() const;
  double gateOptProb() const;  // Probability with which the gate that is inserted (or tweaked) is optimized, rather...
  int numInsertions() const;   // ...than randomized.
  int minInsertionLength() const;
  int maxInsertionLength() const;
  double angleTolerance() const;

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

  // Problem
  int num_qbits;
  std::string problem_;
  int gate_set = 1;  // Default value.

  // GA parameters
  int adult_pop_size = -1;
  int child_pop_size = -1;
  double spacing_ = -1.0;
  int num_rand_survivors = 0;

  // Problem parameters
  int max_length = -1;
  double mean_random_length = -1.0;
  double bad_error_cutoff = 1.1;  // Default value. Indicates that no error cutoff is used.
  double gate_opt_prob = 0.9;  // Default value.
  int num_insertions = 2;  // Default value.
  int min_insertion_length = 2;  // Default value.
  int max_insertion_length = 5;  // Default value.
  double angle_tolerance = 0.01;  // Default value.

  // Termination criteria and random seed
  int num_generations = -1;
  int num_halts = 5;  // Default value.
  int random_seed = -1;  // Indicates 'no value' - set the seed randomly.
};

std::ostream& operator<<(std::ostream& out, const Parameters& parameters);

#endif
