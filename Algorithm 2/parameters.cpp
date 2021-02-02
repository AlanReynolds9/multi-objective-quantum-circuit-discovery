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

#include <string>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include "commandLineParser.h"
#include "parameters.h"

using namespace std;
using boost::lexical_cast;      // I prefer this to std::stoi for converting from string to, for example, int, as it...
using boost::bad_lexical_cast;  // ...catches bad input like "11ddff" rather than treating it as 11.
using boost::to_lower;


Parameters::~Parameters()
{
}


void Parameters::usage() const
{
  cerr << endl;
  cerr << "Usage: " << program_name << endl;
  cerr << "\t-problem <problem>" << endl;
  cerr << "\t-gateSet <gate set (integer) (default = 1)>" << endl;
  cerr << "\t-qbits <number of qbits>" << endl;
  cerr << "\t-adults <adult population size>" << endl;
  cerr << "\t-children <child population size>" << endl;
  cerr << "\t-space <minimum spacing between adult solutions>" << endl;
  cerr << "\t-randSurvivors <number of survivors selected at random (default = 0)>" << endl;
  cerr << "\t-badError <bad error cutoff>" << endl;
  cerr << "\t-gateOptProb <gate optimization probability>" << endl;
  cerr << "\t-maxLength <maximum circuit length>" << endl;
  cerr << "\t-randLength <mean length of random circuit>" << endl;
  cerr << "\t-numInserts <number of insertions (default = 2)>" << endl;
  cerr << "\t-minInsert <minimum insertion length (default = 2)>" << endl;
  cerr << "\t-maxInsert <maximum insertion length (default = 5)>" << endl;
  cerr << "\t-angleTol <angle tolerance> (default = 0.01)>" << endl;
  cerr << "\t-generations <number of generations>" << endl;
  cerr << "\t-halts <number of halts for progress reports (default = 5)>" << endl;
  cerr << "\t-seed <seed for random number generator (optional)>" << endl;
  cerr << "\t-help" << endl << endl;

  cerr << "Option -help provides this usage summary." << endl << endl;

  cerr << "The problem to be solved is specified by the problem parameter ('fourier'," << endl;
  cerr << "'grover' or 'toffoli'), the gate set and the number of qubits. The gate set" << endl;
  cerr << "parameter is just an integer, to allow the user to choose from a small number of" << endl;
  cerr << "pre-selected gate sets." << endl << endl;

  cerr << "The child population should ideally be considerably larger than the adult" << endl;
  cerr << "population, to exploit the fact that evaluating children is cheap, but" << endl;
  cerr << "generating the data in the adults that enables this quick evaluation is not." << endl << endl;

  cerr << "The 'survival' phase of the algorithm works as in NSGA II, but with two" << endl;
  cerr << "exceptions. Firstly, during the NSGA II part of the survival routine, solutions" << endl;
  cerr << "are forbidden from surviving if a solution sufficiently close in the objective" << endl;
  cerr << "space has already been selected for survival. The cutoff distance is given by" << endl;
  cerr << "the 'space' parameter. Secondly, there is an second phase where solutions are" << endl;
  cerr << "merely selected at random. The number of solutions selected in this way is" << endl;
  cerr << "indicated by the 'randSurvivors' parameter. At present, this second phase" << endl;
  cerr << "ignores the space restriction." << endl << endl;

  cerr << "The 'bad error cutoff' indicates a value for overall error that, if exceeded," << endl;
  cerr << "results in circuit construction (from the parents) and circuit simplification" << endl;
  cerr << "being skipped." << endl << endl;

  cerr << "The 'gate optimization probability' gives the chance that the gate that has" << endl;
  cerr << "just been inserted or 'tweaked' is optimized, rather than merely assigned a" << endl;
  cerr << "random gate angle. (This clearly applies only to mutatable gates.)" << endl << endl;

  cerr << "The circuits in the initial population are produced at random, with circuit" << endl;
  cerr << "length sampled from the exponential distribution with the mean given by the" << endl;
  cerr << "-randLength parameter." << endl << endl;
  
  cerr << "Insertions are short random gate sequences, created in the Problem object" << endl;
  cerr << "once in every generation. The matrix for the sequence is pre-calculated. Then" << endl;
  cerr << "when a child solution is generated by a method that inserts a sequence of" << endl;
  cerr << "gates, one of these sequences is selected at random." << endl << endl;

  cerr << "The angle tolerance is used when reducing a parameterized gate, such as an" << endl;
  cerr << "ArbitraryPhase, to an unparameterized gate such as PiByEight, or to the" << endl;
  cerr << "identity (i.e. no gate at all). Since the gate's angle is unlikely to be" << endl;
  cerr << "exactly that required to make the gate and its replacement precisely" << endl;
  cerr << "equivalent, some tolerance is required." << endl << endl;
}


void Parameters::processCommandLine(int argc, char* argv[])
{
  reset();  // Unlikely to be necessary.

  // Get the program name and the command line string
  command_line = program_name = argv[0];
  for (int i = 1; i < argc; ++i)
  {
    command_line += " ";
    command_line += argv[i];
  }

  // Get all the parameters from the command line.
  CommandLineParser parser({"help"}, {"problem", "gateSet", "qbits", "adults", "children", "space", "randSurvivors",
                                      "maxLength", "randLength", "badError", "gateOptProb", "numInserts", "minInsert",
                                      "maxInsert", "angleTol", "generations", "halts", "seed"});
  parser.readArguments(argc, argv);

  while (!parser.finished())
  {
    auto flag = parser.flag();
    if (flag == "help")
    {
      usage();
      exit(EXIT_SUCCESS);
    }
    else if (flag == "problem")
    {
      problem_ = parser.value();
      to_lower(problem_);
    }
    else if (flag == "gateSet")
    {
      gate_set = lexical_cast<int>(parser.value());
    }
    else if (flag == "qbits")
    {
      num_qbits = lexical_cast<int>(parser.value());
    }
    else if (flag == "adults")
    {
      adult_pop_size = lexical_cast<int>(parser.value());
    }
    else if (flag == "children")
    {
      child_pop_size = lexical_cast<int>(parser.value());
    }
    else if (flag == "space")
    {
      spacing_ = lexical_cast<double>(parser.value());
    }
    else if (flag == "randSurvivors")
    {
      num_rand_survivors = lexical_cast<int>(parser.value());
    }
    else if (flag == "maxLength")
    {
      max_length = lexical_cast<int>(parser.value());
    }
    else if (flag == "randLength")
    {
      mean_random_length = lexical_cast<double>(parser.value());
    }
    else if (flag == "badError")
    {
      bad_error_cutoff = lexical_cast<double>(parser.value());
    }
    else if (flag == "gateOptProb")
    {
      gate_opt_prob = lexical_cast<double>(parser.value());
    }
    else if (flag == "numInserts")
    {
      num_insertions = lexical_cast<int>(parser.value());
    }
    else if (flag == "minInsert")
    {
      min_insertion_length = lexical_cast<int>(parser.value());
    }
    else if (flag == "maxInsert")
    {
      max_insertion_length = lexical_cast<int>(parser.value());
    }
    else if (flag == "angleTol")
    {
      angle_tolerance = lexical_cast<double>(parser.value());
    }
    else if (flag == "generations")
    {
      num_generations = lexical_cast<int>(parser.value());
    }
    else if (flag == "halts")
    {
      num_halts = lexical_cast<int>(parser.value());
    }
    else if (flag == "seed")
    {
      random_seed = lexical_cast<int>(parser.value());
    }
    else
    {
      throw CommandLineError("Unknown command line flag '" + flag + "'.");
    }

    parser.next();
  }

  sanity_check();
}


void Parameters::reset()
{
  problem_.clear();
  gate_set = 1;  // Default value
  num_qbits = -1;
  adult_pop_size = -1;
  child_pop_size = -1;
  spacing_ = -1.0;
  num_rand_survivors = 0;
  max_length = -1;
  mean_random_length = -1.0;
  bad_error_cutoff = 1.1;
  gate_opt_prob = 0.9;  // Default value.
  num_insertions = 2;  // Default value.
  min_insertion_length = 2;  // Default value.
  max_insertion_length = 5;  // Default value.
  angle_tolerance = 0.01;  // Default value.
  num_generations = -1;
  num_halts = 5;  // Default value.
  random_seed = -1;  // -1 indicates no random seed provided.
}


void Parameters::sanity_check()
{
  if (problem() != "fourier" && problem() != "grover" && problem() != "toffoli")
  {
    throw CommandLineError("The problem must be either 'fourier', 'grover' or 'toffoli'.");
  }
  if (gateSet() <= 0)
  {
    throw CommandLineError("Gate set must be a positive integer. (Upper limit depends on the problem.)");
  }
  if (numQbits() <= 0)
  {
    throw CommandLineError("Number of qubits must be a positive integer.");
  }
  if (adultPopSize() < 2)
  {
    throw CommandLineError("The adult population must contain at least two solutions.");
  }
  if (childPopSize() < 1)
  {
    throw CommandLineError("The genetic algorithm required a positive child population size.");
  }
  if (spacing() < 0.0)
  {
    throw CommandLineError("The minimum spacing between solutions must be non-negative.");
  }
  if (numRandSurvivors() < 0)
  {
    throw CommandLineError("The number of randomly selected survivors must be non-negative.");
  }
  if (numRandSurvivors() >= adultPopSize())
  {
    throw CommandLineError("The number of randomly selected survivors must be less than the adult population size.");
  }
  if (maxLength() <= 1)
  {
    throw CommandLineError("The maximum circuit length must be positive.");
  }
  if (meanRandomLength() <= 1.0)
  {
    throw CommandLineError("The mean length of a random solution must be greater than 1.0.");
  }
  if (badErrorCutoff() <= 0.0)
  {
    throw CommandLineError("The 'bad error cutoff' must be positive - perhaps the error rate of the empty circuit?");
  }
  if (gateOptProb() < 0.0 || gateOptProb() > 1.0)
  {
    throw CommandLineError("The gate optimization probability must lie between zero and one.");
  }
  if (numInsertions() < 0)
  {
    // Zero is used to indicate that sequence insertions, i.e. insertions of more than one gate, are not used.
    throw CommandLineError("Number of insertions must be non-negative.");
  }
  if (minInsertionLength() < 2)
  {
    // Now that we explicitly permit single gate insertions without optimization, we no longer permit 'insertions' of
    // length 1.
    throw CommandLineError("Minimum insertion length must be at least 2.");
  }
  if (maxInsertionLength() < minInsertionLength())
  {
    throw CommandLineError("Maximum insertion length must not be smaller than the minimum insertion length.");
  }
  if (angleTolerance() < 0)
  {
    throw CommandLineError("Angle tolerance must be non-negative.");
  }
  if (numGenerations() <= 0)
  {
    throw CommandLineError("Number of generations must be provided and be positive.");
  }
  if (numHalts() <= 0)
  {
    throw CommandLineError("Number of halts for progress reports must be provided and positive.");
  }
}


string Parameters::problem() const
{
  return problem_;
}


int Parameters::gateSet() const
{
  return gate_set;
}


int Parameters::numQbits() const
{
  return num_qbits;
}


int Parameters::adultPopSize() const
{
  return adult_pop_size;
}


int Parameters::childPopSize() const
{
  return child_pop_size;
}


double Parameters::spacing() const
{
  return spacing_;
}


int Parameters::numRandSurvivors() const
{
  return num_rand_survivors;
}


int Parameters::maxLength() const
{
  return max_length;
}


double Parameters::meanRandomLength() const
{
  return mean_random_length;
}


double Parameters::badErrorCutoff() const
{
  return bad_error_cutoff;
}


double Parameters::gateOptProb() const
{
  return gate_opt_prob;
}


int Parameters::numInsertions() const
{
  return num_insertions;
}


int Parameters::minInsertionLength() const
{
  return min_insertion_length;
}


int Parameters::maxInsertionLength() const
{
  return max_insertion_length;
}


double Parameters::angleTolerance() const
{
  return angle_tolerance;
}


int Parameters::numGenerations() const
{
  return num_generations;
}


int Parameters::numHalts() const
{
  return num_halts;
}


bool Parameters::randomSeedProvided() const
{
  return random_seed >= 0;
}


int Parameters::randomSeed() const
{
  return random_seed;
}


void Parameters::output(ostream& out) const
{
  out << "Command line:" << endl;
  out << command_line << endl << endl;

  out << "Problem: " << problem() << endl;
  out << "Gate set: " << gateSet() << endl;
  out << "Number of qubits: " << numQbits() << endl;
  out << "Adult population size: " << adultPopSize() << endl;
  out << "Child population size: " << childPopSize() << endl;
  out << "Minimum spacing between solutions: " << spacing() << endl;
  out << "Number of survivors selected at random: " << numRandSurvivors() << endl;
  out << "Maximum circuit length: " << maxLength() << endl;
  out << "Mean length of random circuits: " << meanRandomLength() << endl;
  out << "Bad error cutoff: " << badErrorCutoff() << endl;
  out << "Gate optimization probability: " << gateOptProb() << endl;
  out << "Number of insertions: " << numInsertions() << endl;
  out << "Minimum insertion size: " << minInsertionLength() << endl;
  out << "Maximum insertion size: " << maxInsertionLength() << endl;
  out << "Angle tolerance: " << angleTolerance() << endl;
  out << "Number of generations: " << numGenerations() << endl;
  out << "Number of stops for progress reports: " << numHalts() << endl;
  if (randomSeed() >= 0)
  {
    out << "Seed for random number generator: " << randomSeed() << endl;
  }
}

ostream& operator<<(ostream& out, const Parameters& parameters)
{
  parameters.output(out);
  return out;
}
