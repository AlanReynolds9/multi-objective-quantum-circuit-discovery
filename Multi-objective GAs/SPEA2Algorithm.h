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

#ifndef SPEA2_ALGORITHM
#define SPEA2_ALGORITHM

#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "survival.h"
#include "stores.h"
#include "distances.h"

template <class Problem, class Solution, class Store>
class SPEA2Algorithm
{
public:
  SPEA2Algorithm(const Problem& prob, Store& store, int k, int adultSize, int childSize, double crossProb,
                 double mutationProb);
  
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
  using Selector = TournamentSelector<Solution>;  // NOTE: SortedTournamentSelector should also work: check
  using Breeder = StandardBreeder<Solution, Selector>;
  using Crowding = SPEA2Crowding<Solution, ObjectiveNorm<Solution, 2> >;
  using Evaluator = SPEA2Evaluator<Solution, Crowding>;
  using Survival = SPEA2Survival<Solution, Store, Evaluator>;

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
  const int k_;
  const int adult_population_size;
  const int child_population_size;
  double crossover_prob;
  double mutation_prob;

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
SPEA2Algorithm<Problem, Solution, Store>::SPEA2Algorithm(const Problem& prob, Store& store, int k, int adultSize,
                                                         int childSize, double crossProb, double mutationProb) :
prob_(prob),
store_(store),
k_(k),
generations_limited(false),
time_limited(false),
adult_population_size(adultSize),
child_population_size(childSize),
crossover_prob(crossProb),
mutation_prob(mutationProb),
adults_(adultSize),
children_(childSize),
objectives_(prob),
selector_(adults_),
crowding_(objectives_, k),
evaluator_(crowding_),
breeder_(selector_),
survival_(adults_, children_, store_, evaluator_, objectives_),
generation_(0),
time_taken(0)
{
  breeder_.setCrossoverProb(crossProb);
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::reset()
{
  // Unlike NSGA II or MOGA, SPEA2 starts with a child population which is evaluated and which suffers competition to
  // get the first adults.
  generation_ = 0;
  store_.reset();

  // Be careful to measure the time taken to create and evaluate the initial population. Most of the time this will just
  // be zero, but not so if we plan only a few generations on an expensive problem
  time_t beginTime;
  time(&beginTime);
  children_.initialize(prob_);
  survival_.compete();
  time_t endTime;
  time(&endTime);
  time_taken = static_cast<int>(endTime - beginTime);
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::setGenerations(int numGen)
{
  generations_limited = true;
  num_generations = numGen;
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::setCrossoverProb(double crossProb)
{
  crossover_prob = crossProb;
  breeder_.setCrossoverProb(crossProb);
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::setMutationProb(double mutationProb)
{
  mutation_prob = mutationProb;
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::step()
{
  // Perform one iteration of the algorithm

  ++generation_;

  // Fill the child population
  children_.clear();
  while (!children_.full())
  {
    std::shared_ptr<Solution> child(breeder_.child());
    children_.add(child);
  }
  children_.evaluateObjectives();

  // Survival of the fittest creates the new adult population
  survival_.compete();

  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations or the time limit to be increased, allowing populations to be
  // saved at various points in the search.
  if (!generations_limited && !time_limited)
  {
    std::cerr << "Attempting to run SPEA 2 with no termination conditions. Aborting." << std::endl;
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
void SPEA2Algorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int SPEA2Algorithm<Problem, Solution, Store>::generationsPerformed() const
{
  return generation_;
}

template <class Problem, class Solution, class Store>
int SPEA2Algorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& SPEA2Algorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& SPEA2Algorithm<Problem, Solution, Store>::store() const
{
  // Function to return a constant reference to the store of best solutions, allowing these to be output
  // Not a constant member, since printing the store may require dominated solutions to be eliminated
  return store_;
}

template <class Problem, class Solution, class Store>
void SPEA2Algorithm<Problem, Solution, Store>::completeStore() const
{
  addWorthy(store_, adults_.cbegin(), adults_.cend());
}

#endif
