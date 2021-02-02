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

#ifndef MOGA_ALGORITHM
#define MOGA_ALGORITHM

#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "survival.h"
#include "stores.h"
#include "distances.h"

// (If the algorithm classes turn out to be pretty much the same except for the Evaluator/Survival etc. types, it might
// might sense to make just one which is templated. Currently different due to the requirement to have the niche size in
// the constructor, the need the set the OrgyBreeder size and the addition of the Crowding object.)

template <class Problem, class Solution, class Store>
class MOGAAlgorithm
{
public:
  MOGAAlgorithm(const Problem& prob, Store& store, int popSize, double crossProb, double mutationProb,
                double nicheSize);
  
  // Initialization or reset to generation 0
  void reset();

  // Changing the settings
  void setGenerations(int numGen);
  void setTimeLimit(int timeLimit);
  void setCrossoverProb(double crossProb);
  void setMutationProb(double mutationProb);

  // Running the algorithm
  void go();     // Continue from the current generation until num_generations is reached
  void run();    // Start from a new population of adults

  // Extracting the results
  int generationsPerformed() const;
  int timeTaken() const;
  const Population<Solution>& adults() const;
  Store& store() const; // Not a const Store, as output function may remove dominated solutions just in time.

private:
  void step();                // One step
  void completeStore() const; // Add the non-dominated solutions from the final population to the store.

private:
  // Typedefs for the selector, breeder, evaluator and survival types
  using Rank = std::vector<std::shared_ptr<Solution> >;
  using Selector = StochasticUniversalSelector<Solution>;
  using Breeder = OrgyBreeder<Solution, Selector>;
  using Crowding = MOGACrowding<Solution, ObjectiveNorm<Solution, 2> >;
  using Evaluator = MOGAEvaluator<Solution, Crowding>;
  using Survival = ChildSurvival<Solution, Store>;

  // Problem
  const Problem& prob_;  // (Was merely a parameter of reset. However, we may need to call reset without changing the...
                         // ...problem.)
  // Store
  Store& store_;

  // Number, type and range of objectives
  Objectives objectives_;

  // Termination conditions
  bool generations_limited;
  int num_generations;
  bool time_limited;
  int time_limit;

  // Algorithm settings
  const int population_size;
  double crossover_prob;
  double mutation_prob;

  // Populations
  Population<Solution> adults_;   // Must be before selector_ and survival_
  Population<Solution> children_; // Must be before survival_

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
MOGAAlgorithm<Problem, Solution, Store>::MOGAAlgorithm(const Problem& prob, Store& store, int popSize, double crossProb,
                                                       double mutationProb, double nicheSize) :
prob_(prob),
store_(store),
generations_limited(false),
time_limited(false),
population_size(popSize),
crossover_prob(crossProb),
mutation_prob(mutationProb),
adults_(popSize),
children_(popSize),
objectives_(prob),
selector_(adults_),
crowding_(objectives_, nicheSize),
evaluator_(crowding_),
breeder_(selector_),
survival_(adults_, children_, store_),
generation_(0),
time_taken(0)
{
  breeder_.setSize(popSize);  // A bit that is different to NSGA II!
  breeder_.setCrossoverProb(crossProb);
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::reset()
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
void MOGAAlgorithm<Problem, Solution, Store>::setGenerations(int numGen)
{
  generations_limited = true;
  num_generations = numGen;
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::setCrossoverProb(double crossProb)
{
  crossover_prob = crossProb;
  breeder_.setCrossoverProb(crossProb);
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::setMutationProb(double mutationProb)
{
  mutation_prob = mutationProb;
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::step()   // CHECK: Is the child population ever 'evaluated'??!!??
{                                         // (Not necessary for survival, but necessary for sensible parent selection!)
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

  // Calculate a single value of fitness for each solution. (Unlike NSGA II and SPEA 2, this is not handled by the
  // Survival class)
  evaluator_.evaluate(adults_);

  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations or time limit to be increased, allowing populations to be
  // saved at various points in the search
  if (!generations_limited && !time_limited)
  {
    std::cerr << "Attempting to run MOGA with no termination conditions. Aborting." << std::endl;
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
void MOGAAlgorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int MOGAAlgorithm<Problem, Solution, Store>::generationsPerformed() const
{
  return generation_;
}

template <class Problem, class Solution, class Store>
int MOGAAlgorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& MOGAAlgorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& MOGAAlgorithm<Problem, Solution, Store>::store() const
{
  // Function to return a constant reference to the store of best solutions, allowing these to be output
  // Not a constant member, since printing the store may require dominated solutions to be eliminated
  return store_;
}

template <class Problem, class Solution, class Store>
void MOGAAlgorithm<Problem, Solution, Store>::completeStore() const
{
  addWorthy(store_, adults_.cbegin(), adults_.cend());
}

#endif
