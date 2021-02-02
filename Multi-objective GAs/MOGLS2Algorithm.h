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

#ifndef MOGLS2_ALGORITHM
#define MOGLS2_ALGORITHM

#include <limits>
#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "survival.h"
#include "stores.h"
#include "localSearch.h"

template <class Problem, class Solution, class Store>
class MOGLS2Algorithm
{
public:
  MOGLS2Algorithm(const Problem& prob, Store& store, int popSize, double crossProb, double mutationProb, int numNeighs,
                  int numElite);
  
  // Initialization or reset to generation 0
  void reset();

  // Changing the settings
  void setEvaluations(int numEval);  // Rough limit, to allow for output midway through a run. A generation will not...
                                     // ...be interrupted due to this limit
  void setTotalEvaluations(int numEval);  // Hard limit, to allow comparison with other algorithms. A generation may...
                                          // ...be cut short, to prevent unfair comparisons.
  void setTimeLimit(int timeLimit);  // Will not cut a generation short. (Time based comparisons are more flaky anyway.)

  // Running the algorithm
  void go();     // Continue from the current generation until num_generations is reached
  void run();    // Start from a new population of adults

  // Extracting the results
  int evaluationsPerformed() const;
  int timeTaken() const;
  const Population<Solution>& adults() const;
  Store& store() const; // Not a const Store, as output function may remove dominated solutions just in time.

private:
  void step();                // One step
  void completeStore() const; // Add the non-dominated solutions from the final population to the store.

private:
  // Typedefs for the selector, breeder, evaluator and survival types
  using Selector = ScaledRouletteSelector<Solution>;
  using Breeder = SingleChildBreeder<Solution, Selector>;
  using Evaluator = RandomWeightedEvaluator<Solution>;
  using Survival = ChildSurvival<Solution, Store>;

  // Problem
  const Problem& prob_;  // (Was merely a parameter of reset. However, we may need to call reset without changing the...
                         // ...problem.)
  // Store
  Store& store_;

  // Termination conditions
  bool evaluations_limited;
  int num_evaluations;    // Number of evaluations that must occur before the algorithm pauses. LS is not interrupted.
  int total_evaluations;  // Hard limit which will result in the algorithm being halted, even within local search.
  bool time_limited;
  int time_limit;         // Local search not interrupted.

  // Algorithm settings
  const int population_size;
  double crossover_prob;
  double mutation_prob;
  int num_neighbours;  // Number of neighbours to examine during local search
  int num_elite;       // Number of elite solutions taken from the store when creating the new population

  // Populations
  Population<Solution> adults_;   // Must be before selector_ and survival_
  Population<Solution> children_; // Must be before survival_

  // Number, type and range of objectives
  Objectives objectives_;        // Must be before evaluator_

  // Algorithm components
  Selector selector_;            // Must be before breeder_
  Evaluator evaluator_;
  Breeder breeder_;
  Survival survival_;

  // Algorithm status
  int evaluations_performed;
  int time_taken;
};

template <class Problem, class Solution, class Store>
MOGLS2Algorithm<Problem, Solution, Store>::MOGLS2Algorithm(const Problem& prob, Store& store, int popSize,
                                                           double crossProb, double mutationProb, int numNeighs,
                                                           int numElite) :
prob_(prob),
store_(store),
evaluations_limited(false),
num_evaluations(std::numeric_limits<int>::max()),   // Not great, but adding another bool (evaluations_hard_limited)...
total_evaluations(std::numeric_limits<int>::max()), // ...is not appealing at present
time_limited(false),
population_size(popSize),
crossover_prob(crossProb),
mutation_prob(mutationProb),
num_neighbours(numNeighs),
num_elite(numElite),
adults_(popSize),
children_(popSize),
objectives_(prob),
selector_(adults_),
evaluator_(objectives_),
breeder_(selector_),
survival_(adults_, children_, store_),
evaluations_performed(0),
time_taken(0)
{
  breeder_.setCrossoverProb(crossProb);
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::reset()
{
  store_.reset();

  // Be careful to measure the time taken to create and evaluate the initial population. Most of the time this will just
  // be zero, but not so if we plan only a few generations on an expensive problem
  time_t beginTime;
  time(&beginTime);

  adults_.initialize(prob_);
  evaluations_performed = population_size;  // It is fairly safe to assume that the population size does not exceed...
                                            // ...the evaluation limit!
  // Add the newly created adults to the store - must be done now, as we need the 'elite' solutions
  for (int i = 0; i < adults_.size(); ++i)
  {
    store_.add(adults_[i]);
  }

  time_t endTime;
  time(&endTime);
  time_taken = static_cast<int>(endTime - beginTime);
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::setEvaluations(int numEval)
{
  evaluations_limited = true;
  num_evaluations = numEval;
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::setTotalEvaluations(int numEval)
{
  evaluations_limited = true;
  total_evaluations = numEval;
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::step()
{
  // Perform one generation of the algorithm, which will involve an unknown number of evaluations

  // Fill the child population
  children_.clear();

  // ...first copy some elite solutions from the store
  assert(!store_.empty());
  for (int i = 0; i < num_elite; ++i)
  {
    auto child = make_shared<Solution>(*store_.randomElement());  // Looks right not to clone here.
    children_.add(child);
  }
  while (!children_.full() && evaluations_performed < total_evaluations)
  {
    // ...then create the rest of the children by breeding from the adult population
    evaluator_.evaluate(adults_); // Before each child is created, re-evaluate with different random weighted objective.
    std::shared_ptr<Solution> child(breeder_.child());
    child->evaluateObjectives();
    ++evaluations_performed;
    children_.add(child);

    // Add child to the store. (Not in algorithm as described by Ishibuchi and Murata - they only add the solution after
    // local search.)
    store_.add(child);
  }

  // Perform first found local search on each of the solutions, using a equally spaced linear combination of the scaled
  // objectives, accepting only better solutions. The local search object also updates the number of evaluations
  // performed.
  AngledLocalSearch<Solution, Store> localSearch(children_, objectives_, store_, num_neighbours, evaluations_performed,
                                                 total_evaluations);
  localSearch.optimize();

  // Survival of the fittest creates the new adult population
  survival_.compete();

  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations to be increased, allowing populations to be saved at various
  // points in the search
  if (!evaluations_limited && !time_limited)
  {
    cerr << "Attempting to run MOGLS2 with no termination conditions. Aborting." << endl;
    exit(EXIT_FAILURE);
  }

  int timeLeft = time_limit - time_taken;
  time_t beginTime;
  time(&beginTime);
  time_t currentTime = beginTime;
  while ((!evaluations_limited || evaluations_performed < std::min(num_evaluations, total_evaluations)) &&
         (!time_limited || currentTime - beginTime < timeLeft))  // Playing safe - num_evaluations shouldn't be...
  {                                                              // ...greater than total_evaluations.
    step();
    time(&currentTime);
  }
  time_taken += static_cast<int>(currentTime - beginTime);

  // Ensure that the non-dominated solutions from the final population are added to the store of the best found
  completeStore();
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int MOGLS2Algorithm<Problem, Solution, Store>::evaluationsPerformed() const
{
  return evaluations_performed;
}

template <class Problem, class Solution, class Store>
int MOGLS2Algorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& MOGLS2Algorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& MOGLS2Algorithm<Problem, Solution, Store>::store() const
{
  // Function to return a constant reference to the store of best solutions, allowing these to be output
  // Not a constant member, since printing the store may require dominated solutions to be eliminated
  return store_;
}

template <class Problem, class Solution, class Store>
void MOGLS2Algorithm<Problem, Solution, Store>::completeStore() const
{
  // Store is already complete, since solutions have been added upon creation. No need to call addWorthy, though it wouldn't add anything
  // anyway, since 'worthy()' always returns false.
}

#endif
