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
  cerr << "\t-pop <population size>" << endl;
  cerr << "\t-crossRate <crossover rate>" << endl;
  cerr << "\t-mutRate <mutation rate>" << endl;
  cerr << "\t-space <minimum spacing between adult solutions>" << endl;
  cerr << "\t-randLength <mean length of random circuit>" << endl;
  cerr << "\t-maxLength <maximum circuit length>" << endl;
  cerr << "\t-target <target value of overall error considered 'good enough' (default = 0.00001)>" << endl;
  cerr << "\t-expReduction <expected reduction in cost (default = 'huge')>" << endl;
  cerr << "\t-randParams" << endl;
  cerr << "\t-numOptQuota <iteration quota for numerical optimization (default = 50)>" << endl;
  cerr << "\t-useCache" << endl;
  cerr << "\t-generations <number of generations>" << endl;
  cerr << "\t-halts <number of halts for progress reports (default = 5)>" << endl;
  cerr << "\t-seed <seed for random number generator (optional)>" << endl;
  cerr << "\t-help" << endl << endl;

  cerr << "Option -help provides this usage summary." << endl << endl;

  cerr << "The problem to be solved is specified by the problem parameter ('fourier'," << endl;
  cerr << "'grover' or 'toffoli'), the gate set and the number of qubits. The gate set" << endl;
  cerr << "parameter is just an integer, to allow the user to choose from a small number of" << endl;
  cerr << "pre-selected gate sets." << endl << endl;

  cerr << "The crossover rate is the chance of crossover being performed, while the" << endl;
  cerr << "mutation rate is the probability of performing each of the 5 mutation operators." << endl << endl;

  cerr << "The randomly generations circuits in the initial population are given a length" << endl;
  cerr << "sampled from the exponential distribution with mean value given by the" << endl;
  cerr << "randLength parameter. Circuits of length greater than the maximum are not" << endl;
  cerr << "evaluated, but simply assigned maximal error values." << endl << endl;

  cerr << "The normal NSGA II survival routine is modified by requiring a minimum spacing" << endl;
  cerr << "between solutions in objective space. This is specified using the -space" << endl;
  cerr << "parameter. Note that if the algorithm finds it impossible to fill the new adult" << endl;
  cerr << "population without violating this restriction, then the value of the minimum" << endl;
  cerr << "spacing is (temporarily) halved until the population can be filled." << endl << endl;

  cerr << "The 'target value' is used for two purposes. First, it provides an extra" << endl;
  cerr << "halting condition for the numerical optimization. Second, it is used in" << endl;
  cerr << "combination with the 'expected reduction in cost' to prevent evaluation of" << endl;
  cerr << "circuits when a 'good enough' circuit of lesser cost has already been found." << endl;
  cerr << "In detail, if a circuit that meets the 'target value' for error has been" << endl;
  cerr << "found at cost c, then circuits of cost at least c + expectedReduction are" << endl;
  cerr << "not evaluated. The 'expected reduction' parameter acknowledges that, after" << endl;
  cerr << "numerical optimization, it may be possible to 'reduce' the circuit, lowering" << endl;
  cerr << "its cost. The default values ensure that numerical optimization is curtailed" << endl;
  cerr << "at a sensible point, while no solution is left unevaluated due to a 'good" << endl;
  cerr << "enough' solution having been already found." << endl << endl;
  
  cerr << "Option -randParams switches on randomization of gate parameters before" << endl;
  cerr << "numerical optimization. (If not used, gate parameters will start at the special" << endl;
  cerr << "values found to work well with the parent solution(s).)" << endl << endl;
  cerr << "The 'iteration quota' is used in numerical optimization. If it appears, from" << endl;
  cerr << "the last few iterations, that the numerical optimization will not improve" << endl;
  cerr << "the front within an extra 'quota' of iterations, then it is halted." << endl << endl;
  
  cerr << "Option -useCache is used to indicate that a solution cache is to be used to" << endl;
  cerr << "avoid re-evaluation of solutions." << endl;
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
  CommandLineParser parser({"help", "randParams", "useCache"}, {"problem", "gateSet", "qbits", "pop", "crossRate",
                            "mutRate", "space", "randLength", "maxLength", "target", "expReduction", "numOptQuota",
                            "generations", "halts", "seed"});
  parser.readArguments(argc, argv);

  while (!parser.finished())
  {
    auto flag = parser.flag();
    if (flag == "help")
    {
      usage();
      exit(EXIT_SUCCESS);
    }
    else if (flag == "randParams")
    {
      randomize_gate_parameters = true;
    }
    else if (flag == "useCache")
    {
      use_cache = true;
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
    else if (flag == "pop")  // NSGA II
    {
      population_size = lexical_cast<double>(parser.value());
    }
//    else if (flag == "adults")  // QCGA
//    {
//      adult_pop_size = lexical_cast<int>(parser.value());
//    }
//    else if (flag == "children")  // QCGA
//    {
//      child_pop_size = lexical_cast<int>(parser.value());
//    }
    else if (flag == "crossRate")
    {
      crossover_rate = lexical_cast<double>(parser.value());
    }
    else if (flag == "mutRate")
    {
      mutation_rate = lexical_cast<double>(parser.value());
    }
//    else if (flag == "elite")  // QCGA
//    {
//      num_elite = lexical_cast<int>(parser.value());
//    }
//    else if (flag == "trials")  // QCGA
//    {
//      num_trials = lexical_cast<int>(parser.value());
//    }
    else if (flag == "space")  // QCGA and crowded NSGA II
    {
      spacing_ = lexical_cast<double>(parser.value());
    }
    else if (flag == "randLength")
    {
      mean_random_length = lexical_cast<double>(parser.value());
    }
    else if (flag == "maxLength")
    {
      max_circuit_length = lexical_cast<int>(parser.value());
    }
    else if (flag == "target")
    {
      target_error = lexical_cast<double>(parser.value());
    }
    else if (flag == "expReduction")
    {
      expected_reduction = lexical_cast<long>(parser.value());
    }
    else if (flag == "numOptQuota")
    {
      iteration_quota = lexical_cast<int>(parser.value());
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
      break;
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
  population_size = -1;
  crossover_rate = -1.0;
  mutation_rate = -1.0;
  spacing_ = -1.0;
  mean_random_length = -1.0;
  max_circuit_length = -1;
  target_error = 0.00001;  // Default value
  expected_reduction = std::numeric_limits<long>::max();  // Default value - all circuit structures are evaluated.
  randomize_gate_parameters = false;  // Default value.
  use_cache = false;  // Default value.
  num_generations = -1;
  num_halts = 5;  // Default value.
  random_seed = -1;  // -1 indicates no random seed provided.
  iteration_quota = 50;  // Default value.
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
  if (populationSize() < 1)
  {
    throw CommandLineError("The genetic algorithm requires a positive population size.");
  }
  if (crossoverRate() < 0.0 || crossoverRate() > 1.0)
  {
    throw CommandLineError("The genetic algorithm requires a crossover rate between 0 and 1 inclusive.");
  }
  if (mutationRate() < 0.0 || mutationRate() > 1.0)
  {
    throw CommandLineError("The genetic algorithm requires a mutation rate between 0 and 1 inclusive.");
  }
  if (spacing() < 0.0)
  {
    throw CommandLineError("The minimum spacing between solutions must be non-negative.");
  }
  if (meanRandomLength() <= 0.0)
  {
    throw CommandLineError("The mean length of a random solution must be positive.");
  }
  if (maxCircuitLength() <= 10)
  {
    throw CommandLineError("The maximum circuit length must be at least 10.");
  }
  if (expectedReduction() < 0)
  {
    throw CommandLineError("The expected reduction must be non-negative.");
  }
  if (iterationQuota() <= 0)
  {
    throw CommandLineError("The iteration quota for the numerical optimization must be positive.");
  }

  if (numGenerations() <= 0)
  {
    throw CommandLineError("Number of generations must be provided and be positive.");
  }
  if (numHalts() <= 0)
  {
    throw CommandLineError("Number of halts for progress reports must be provided and positive.");
  }
  if (!useCache())
  {
    cerr << "Warning: no solution cache will be used." << endl;
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


int Parameters::populationSize() const
{
  return population_size;
}


double Parameters::crossoverRate() const
{
  return crossover_rate;
}


double Parameters::mutationRate() const
{
  return mutation_rate;
}


double Parameters::spacing() const
{
  return spacing_;
}


double Parameters::meanRandomLength() const
{
  return mean_random_length;
}


int Parameters::maxCircuitLength() const
{
  return max_circuit_length;
}


double Parameters::targetError() const
{
  return target_error;
}


long Parameters::expectedReduction() const
{
  return expected_reduction;
}


bool Parameters::randomizeGateParameters() const
{
  return randomize_gate_parameters;
}


int Parameters::iterationQuota() const
{
  return iteration_quota;
}


bool Parameters::useCache() const
{
  return use_cache;
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
  out << "Population size: " << populationSize() << endl;
  out << "Crossover rate: " << crossoverRate() << endl;
  out << "Mutation rate: " << mutationRate() << endl;
  out << "Minimum spacing between solutions: " << spacing() << endl;
  out << "Mean length of random circuits: " << meanRandomLength() << endl;
  out << "Maximum circuit length: " << maxCircuitLength() << endl;
  out << "Target for overall error considered 'good enough': " << targetError() << endl;
  out << "Expected reduction after numerical optimization: " << expectedReduction() << endl;

  out << "Gate parameters are ";
  if (!randomizeGateParameters())
  {
    out << "not ";
  }
  out << "randomized before numerical optimization." << endl;

  out << "Iteration quota (for numerical optimization): " << iterationQuota() << endl;

  out << "A cache for visited circuits is ";
  if (!useCache())
  {
    out << "not ";
  }
  out << "in use." << endl;

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
