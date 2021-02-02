// MOMH-lib (Multi-objective metaheuristic library.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOMH-lib.
//
// MOMH-lib is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOMH-lib is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with MOQCD2.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef NSGA2GP_ALGORITHM
#define NSGA2GP_ALGORITHM

#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "survival.h"
#include "stores.h"

// Genetic programming version of NSGA II. Unlike the standard version, this one still assumes that the standard NSGA II
// survival methods will be used.
template <class Problem, class Solution, class Store>
class NSGA2GPAlgorithm
{
public:
  NSGA2GPAlgorithm(const Problem& prob, Store& store, int popSize, double crossProb);
  
  // Initialization or reset to generation 0
  void reset();

  // Changing the settings
  void setGenerations(int numGen);
  void setTimeLimit(int timeLimit);
  void setCrossoverProb(double crossProb);  // Mutation is performed whenever crossover isn't

  // Running the algorithm
  void go();     // Continue from the current generation until num_generations is reached
  void run();    // Start from a new population of adults

  // Extracting the results
  int generationsPerformed() const;
  int timeTaken() const;
  const Population<Solution>& adults() const;
  Store& store() const; // Not a const Store, as output function may remove dominated solutions 'just in time'.

private:
  void step();                // One step
  void completeStore() const; // Add the non-dominated solutions from the final population to the store.

private:
  // 'Typedefs' for the selector, breeder, evaluator and survival types
  using Selector = SortedTournamentSelector<Solution>;
  using Breeder = GPBreeder<Solution, Selector>;
  using Crowding = NormedNSGA2Crowding<Solution>;
  using Evaluator = NSGA2Evaluator<Solution, Crowding>;
  using Survival = SortedBestSurvival<Solution, Store, Evaluator>;

  // Problem
  const Problem& prob_;

  // Store
  Store& store_;

  // Termination conditions
  bool generations_limited;
  int num_generations;
  bool time_limited;
  int time_limit;

  // Algorithm settings
  const int population_size;
  double crossover_prob;

  // Populations
  Population<Solution> adults_;   // Must be before selector_ and survival_
  Population<Solution> children_; // Must be before survival_

  // Number, type and range of objectives
  Objectives objectives_;

  // Algorithm components
  Selector selector_;            // Must be before breeder_
  Crowding crowding_;            // Must be before evaluator_
  Evaluator evaluator_;
  Breeder breeder_;
  Survival survival_;

  // Algorithm status
  int generation_;
  int time_taken;
};

template <class Problem, class Solution, class Store>
NSGA2GPAlgorithm<Problem, Solution, Store>::NSGA2GPAlgorithm(const Problem& prob, Store& store, int popSize,
                                                             double crossProb) :
prob_(prob),
store_(store),
generations_limited(false),
time_limited(false),
population_size(popSize),
crossover_prob(crossProb),
adults_(popSize),
children_(popSize),
objectives_(prob),
selector_(adults_),
crowding_(objectives_),
evaluator_(objectives_, crowding_, true),
breeder_(selector_),
survival_(adults_, children_, store_, evaluator_),
generation_(0),
time_taken(0)
{
  breeder_.setCrossoverProb(crossProb);
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::reset()
{
  generation_ = 0;
  store_.reset();

  // Be careful to measure the time taken to create and evaluate the initial population. Most of the time this will just
  // be zero, but not so if we plan only a few generations on an expensive problem
  time_t beginTime;
  time(&beginTime);
  adults_.initialize(prob_);
  evaluator_.evaluate(adults_);
  time_t endTime;
  time(&endTime);
  time_taken = static_cast<int>(endTime - beginTime);
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::setGenerations(int numGen)
{
  generations_limited = true;
  num_generations = numGen;
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::setCrossoverProb(double crossProb)
{
  crossover_prob = crossProb;
  breeder_.setCrossoverProb(crossProb);
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::step()
{
  // Perform one iteration of the algorithm

  ++generation_;

  // Fill the child population
  children_.clear();
  while (!children_.full())
  {
    std::shared_ptr<Solution> child(breeder_.child());
    child->evaluateObjectives();
    children_.add(child);
  }

  // Survival of the fittest creates the new adult population
  survival_.compete();

  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations or the time limit to be increased, allowing populations to be
  // saved at various points in the search
  if (!generations_limited && !time_limited)
  {
    cerr << "Attempting to run NSGA II (GP) with no termination conditions. Aborting." << endl;
    exit(EXIT_FAILURE);
  }

  int timeLeft = time_limit - time_taken;
  time_t beginTime;
  time(&beginTime);
  time_t currentTime = beginTime;
  while ((!generations_limited || generation_ < num_generations) &&
         (!time_limited || currentTime - beginTime < timeLeft))
  {
    step();
    time(&currentTime);
  }
  time_taken += static_cast<int>(currentTime - beginTime);

  // Ensure that the non-dominated solutions from the final population are added to the store of the best found
  completeStore();
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int NSGA2GPAlgorithm<Problem, Solution, Store>::generationsPerformed() const
{
  return generation_;
}

template <class Problem, class Solution, class Store>
int NSGA2GPAlgorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& NSGA2GPAlgorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& NSGA2GPAlgorithm<Problem, Solution, Store>::store() const
{
  // Original intent was to return a constant reference, since this routine is only used to report the non-dominated
  // solutions to the user, e.g. for store output. However, for some types of store, dominated solutions are removed
  // 'just in time'. Hence dominated solutions might be eliminated during output. (Consider use of 'mutable' in such
  // stores?)
  return store_;
}

template <class Problem, class Solution, class Store>
void NSGA2GPAlgorithm<Problem, Solution, Store>::completeStore() const
{
  addWorthy(store_, adults_.cbegin(), adults_.cend());
}

#endif
