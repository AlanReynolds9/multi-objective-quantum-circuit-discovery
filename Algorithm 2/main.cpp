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

#include <ctime>
#include <fstream>
#include <iomanip>
#include "checked_stream.h"
#include "QCGA2Solution.h"
#include "QCGA2Algorithm.h"
#include "fourier.h"
#include "grover.h"
#include "toffoli.h"
#include "utils.h"
#include "stores.h"
#include "parameters.h"
#include "commandLineParser.h"  // Added just to get access to CommandLineError

using namespace std;
using namespace utils::io;


// FOURIER/GROVER/TOFFOLI

template <class Algorithm>
void performExperiment(Algorithm& ga, const Parameters& parameters)
{
  // Perform an experiment with the provided genetic algorithm. In other words, run the algorithm but also manage the
  // output, sending it to the appropriate files, collect timing information etc.

  // Prepare output stream for timing information
  checked<ofstream> timeOut("time.txt");
  time_t algorithmStart;
  time(&algorithmStart);

  ga.reset(); // (Leave this out and we get a memory based run-time error. Not nice!)

  string initialFilename = "initial.sols";
  checked<ofstream> initialOut(initialFilename);
  initialOut << ga.adults() << endl << endl;
  initialOut.close();

  unsigned totalGenerations = parameters.numGenerations();
  unsigned numHalts = parameters.numHalts();
  assert(totalGenerations % numHalts == 0);

  unsigned genStepSize = totalGenerations / numHalts;
  cout << "Beginning circuit optimization: " << totalGenerations << " generations, with reports every " <<
          genStepSize << "." << endl;
  for (unsigned numGenerations = genStepSize; numGenerations <= totalGenerations; numGenerations += genStepSize)
  {
    ga.setGenerations(numGenerations);
    time_t beginTime;
    time(&beginTime);
    ga.go();
    time_t endTime;
    time(&endTime);

    // Output best solutions found
    string bestFilename = to_string(numGenerations) + "gen-best.sols";
    ga.store().sortByObjective(1);  // Sort by overall error.
    checked<ofstream> bestOut(bestFilename);
    bestOut << ga.store() << endl << endl;
    bestOut.close();

    // Output current population
    string popFilename = to_string(numGenerations) + "gen-adult.sols";
    checked<ofstream> popOut(popFilename);
    popOut << ga.adults() << endl << endl;
    popOut.close();

    // Output timings to standard output and to file.
    cout << " Generations = " << numGenerations << " -- Time taken = " << static_cast<long>(endTime - beginTime) <<
            endl;
    timeOut << "Generations = " << numGenerations << " -- Time taken = " << static_cast<long>(endTime - beginTime) <<
               endl;
  }

  time_t algorithmFinish;
  time(&algorithmFinish);
  timeOut << endl << "Total time = " << static_cast<long>(algorithmFinish - algorithmStart) << endl;
  cout << endl << "Total time = " << static_cast<long>(algorithmFinish - algorithmStart) << endl;
}


int main(int argc, char* argv[])
{
  Parameters parameters;
  parameters.processCommandLine(argc, argv);
  checked<ofstream> settingsOut("settings.txt");
  settingsOut << parameters;
  settingsOut.close();

  if (parameters.randomSeedProvided())
  {
    utils::rand::seedRand(parameters.randomSeed());
  }

  if (parameters.problem() == "fourier")
  {
    // Choose the gate set.
    double gateCostRange = 100.0;

    std::vector<GateCreator> permittedGateCreators
    {
      masterGateCreatorIndex.at("YRot"),
      masterGateCreatorIndex.at("ArbPhase"),
      masterGateCreatorIndex.at("Swap")
    };
    std::vector<PermittedControls> permittedControls
    {
      PermittedControls::none,
      PermittedControls::any,
      PermittedControls::none
    };
    std::vector<GateCost> costs{{10}, {10}, {1}};

    // Create the Problem object.
    using fourier::Problem;
    Problem problem(permittedGateCreators, permittedControls, costs, gateCostRange, parameters.numQbits(),
                    parameters.maxLength(), parameters.meanRandomLength(), parameters.badErrorCutoff(),
                    parameters.gateOptProb(), parameters.numInsertions(), parameters.minInsertionLength(),
                    parameters.maxInsertionLength(), parameters.angleTolerance());
    problem.transferQbitInputOptionsToContext();

    // Create the algorithm object and perform a set of experiments with it.
    using Solution = QCGA2Solution<fourier::Solution>;
    using Store = BasicManagedStore<Solution>;
    Store store;
    QCGA2Algorithm<Problem, Solution, Store> ga(problem, store, parameters.adultPopSize(), parameters.childPopSize(),
                                                parameters.spacing(), parameters.numRandSurvivors());
    performExperiment(ga, parameters);
  }
  else if (parameters.problem() == "toffoli")
  {
    // Choose the gate set
    double gateCostRange;
    std::vector<GateCreator> permittedGateCreators;
    std::vector<PermittedControls> permittedControls;
    std::vector<GateCost> costs;
    if (parameters.gateSet() == 1)
    {
      // Essentially any single qubit gates, expressed as rotations. I would use SU2 if I had them implemented fully.
      // Lots of gate parameters.
      gateCostRange = 1.0;  //100.0;

      permittedGateCreators =
      {
        masterGateCreatorIndex.at("XRot"),
        masterGateCreatorIndex.at("YRot"),
        masterGateCreatorIndex.at("ZRot"),
        masterGateCreatorIndex.at("X")
      };
      permittedControls =
      {
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::one
      };
      costs = {{1}, {1}, {1}, {1, 20}};
    }
    else if (parameters.gateSet() == 2)
    {
      // Non-paramaterized gates, as given in Nielsen and Chuang.
      gateCostRange = 1.0;  //100.0;

      permittedGateCreators =
      {
        masterGateCreatorIndex.at("Hadamard"),
        masterGateCreatorIndex.at("PiByEight"),
        masterGateCreatorIndex.at("PiByEightInv"),
        masterGateCreatorIndex.at("Phase"),
        masterGateCreatorIndex.at("PhaseInv"),  // Not in circuit from Nielsen and Chuang, but feels odd leaving it out.
        masterGateCreatorIndex.at("Z"),         // Aids simplification? (Apparently I used to include Y too.)
        masterGateCreatorIndex.at("X")
      };
      permittedControls =
      {
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::one
      };
      costs = {{10}, {10}, {10}, {10}, {10}, {10}, {1}/*, 20}*/};  // These costs are backward, but seem to produce...
    }                                                              // ...better results!!!
    else if (parameters.gateSet() == 3)
    {
      // Non-paramaterized gates, as given in Nielsen and Chuang.
      gateCostRange = 1.0;  //100.0;

      permittedGateCreators =
      {
        masterGateCreatorIndex.at("Hadamard"),
        masterGateCreatorIndex.at("ArbPhase"),
        masterGateCreatorIndex.at("X")
      };
      permittedControls =
      {
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::one
      };
      costs = {{10}, {10}, {1}/*, 20}*/};  // These costs are backward, but seem to produce better results!!!
    }
    else
    {
      throw CommandLineError("Unavailable gate set selected.");
    }

    // Create the Problem object.
    using toffoli::Problem;
    Problem problem(permittedGateCreators, permittedControls, costs, gateCostRange, parameters.numQbits(),
                    parameters.maxLength(), parameters.meanRandomLength(), parameters.badErrorCutoff(),
                    parameters.gateOptProb(), parameters.numInsertions(), parameters.minInsertionLength(),
                    parameters.maxInsertionLength(), parameters.angleTolerance());
    problem.transferQbitInputOptionsToContext();

    // Create the algorithm object and perform a set of experiments with it.
    using Solution = QCGA2Solution<toffoli::Solution>;
    using Store = BasicManagedStore<Solution>;
    Store store;
    QCGA2Algorithm<Problem, Solution, Store> ga(problem, store, parameters.adultPopSize(), parameters.childPopSize(),
                                                parameters.spacing(), parameters.numRandSurvivors());
    performExperiment(ga, parameters);
  }
  else if (parameters.problem() == "grover")
  {
    // Choose the gate set
    double gateCostRange;
    std::vector<GateCreator> permittedGateCreators;
    std::vector<PermittedControls> permittedControls;
    std::vector<GateCost> costs;
    if (parameters.gateSet() == 1)
    {
      // Gate set 1 is that used by Vasek, with parameterized gates.
      gateCostRange = 1000.0;
      permittedGateCreators =
      {
        masterGateCreatorIndex.at("Z"),
        masterGateCreatorIndex.at("XRot"),
        masterGateCreatorIndex.at("ArbPhase"),
        masterGateCreatorIndex.at("Oracle")
      };
      permittedControls =
      {
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::any,
        PermittedControls::none  // Anything other than none won't work with Oracle!
      };
      costs = {{9}, {10}, {10}, {500}};
    }
    else if (parameters.gateSet() == 2)
    {
      gateCostRange = 200.0;  // Was 200.0
      permittedGateCreators =
      {
        masterGateCreatorIndex.at("Z"),
        masterGateCreatorIndex.at("X"),
        masterGateCreatorIndex.at("Hadamard"),
        masterGateCreatorIndex.at("Oracle")
      };
      permittedControls =
      {
        PermittedControls::any,
        PermittedControls::none,
        PermittedControls::none,
        PermittedControls::none  // Anything other than none won't work with Oracle!
      };
      costs = {{10}, {1}, {1}, {100}};  // ZGate cost was {1, 10}.
    }
    else
    {
      throw CommandLineError("Unavailable gate set selected.");
    }

    using grover::Problem;
    Problem problem(permittedGateCreators, permittedControls, costs, gateCostRange, parameters.numQbits(),
                    parameters.maxLength(), parameters.meanRandomLength(), parameters.badErrorCutoff(),
                    parameters.gateOptProb(), parameters.numInsertions(), parameters.minInsertionLength(),
                    parameters.maxInsertionLength(), parameters.angleTolerance());
    problem.transferQbitInputOptionsToContext();

    // Create the algorithm object and perform a set of experiments with it.
    using Solution = QCGA2Solution<grover::Solution>;
    using Store = BasicManagedStore<Solution>;
    Store store;
    QCGA2Algorithm<Problem, Solution, Store> ga(problem, store, parameters.adultPopSize(), parameters.childPopSize(),
                                                parameters.spacing(), parameters.numRandSurvivors());
    performExperiment(ga, parameters);
  }
  else
  {
    throw CommandLineError("Unavailable problem type selected.");
  }

  return EXIT_SUCCESS;
}
