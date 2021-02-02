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

#ifndef JMOGLS_ALGORITHM
#define JMOGLS_ALGORITHM

#include <limits>
#include <rng.h>
#include <memory>
#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "survival.h"
#include "stores.h"
#include "localSearch.h"

template <class Problem, class Solution, class Store>
class JMOGLSAlgorithm
{
public:
  JMOGLSAlgorithm(const Problem& prob, Store& store, int tempSize, int numStart, double crossProb, double mutationProb,
                  int numNeighs);
  
  // Initialization or reset to generation 0
  void reset();

  // Changing the settings
  void setEvaluations(int numEval);       // Set a rough evaluation limit, to allow for output midway through a run.
                                          // A generation is not interrupted due to this limit
  void setTotalEvaluations(int numEval);  // Set a hard evaluation limit, to allow comparison with other algorithms.
                                          // A generation may be cut short, to prevent unfair comparisons.
  void setTimeLimit(int timeLimit);       // Will not cut a generation short. (Time based comparisons are somewhat...
                                          // ...more flaky anyway.)
  // Running the algorithm
  void go();     // Continue from the current generation until num_generations is reached
  void run();    // Start from a new population of adults

  // Extracting the results
  int evaluationsPerformed() const;
  int numLSPerformed() const;
  int timeTaken() const;
  const Population<Solution>& adults() const;
  Store& store() const; // Not a const Store, as output function may remove dominated solutions just in time.

private:
  void step();                // One step
  void completeStore() const; // Add the non-dominated solutions from the final population to the store.

private:
  // 'Typedefs' for the selector, breeder, evaluator and survival types
  using TempSelector = SortedBestSelector<Solution>;         // Selector used to create the temporary population
  using ParentSelector = UniformDifferentManySelector<Solution>; // Selector used to choose parents from the...
  using Breeder = SingleChildBreeder<Solution, ParentSelector>;  // ...temporary population.
  using Evaluator = WeightedEvaluator<Solution>;
  using Survival = JMOGLSSurvival<Solution, Store>;
  using LocalSearch = WeightedLocalSearch<Solution, Store>;

  // Problem
  const Problem& prob_;  // (Was merely a parameter of reset. However, we may need to call reset without changing the...
                         // ...problem
  // Store
  Store& store_;

  // Termination conditions
  bool evaluations_limited;
  int num_evaluations;    // Number of evaluations that must occur before the algorithm pauses. LS is not interrupted
  int total_evaluations;  // Hard limit on the number of evaluations - halts the algorithm even within local search.
  bool time_limited;
  int time_limit;

  // Algorithm settings
  const int temp_size;
  const int start_size;
  const int pop_size;
  double crossover_prob;
  double mutation_prob;
  int num_neighbours;  // Number of neighbours to examine during local search

  // Populations
  Population<Solution> adults_;         // Must be before selector_ and survival_
  Population<Solution> temp_population; // Must be before parent_selector
  std::shared_ptr<Solution> child_;

  // Number, type and range of objectives
  Objectives objectives_;        // Must be before evaluator_

  // Algorithm components
  TempSelector temp_selector;
  ParentSelector parent_selector;      // Must be before breeder_
  Evaluator evaluator_;
  Breeder breeder_;
  Survival survival_;

  // Local search weights
  std::vector<double> weights_;

  // Algorithm status
  int evaluations_performed;
  int num_LS_performed;
  int time_taken;
};

template <class Problem, class Solution, class Store>
JMOGLSAlgorithm<Problem, Solution, Store>::JMOGLSAlgorithm(const Problem& prob, Store& store, int tempSize,
                                                           int numStart, double crossProb, double mutationProb,
                                                           int numNeighs) :
prob_(prob),
store_(store),
evaluations_limited(false),
num_evaluations(std::numeric_limits<int>::max()),
total_evaluations(std::numeric_limits<int>::max()),
time_limited(false),
temp_size(tempSize),
start_size(numStart),
pop_size(tempSize * numStart),
crossover_prob(crossProb),
mutation_prob(mutationProb),
num_neighbours(numNeighs),
adults_(pop_size),
temp_population(tempSize),
child_(new Solution(prob_)),                   // Need to initialize the child before the survival object
objectives_(prob),
temp_selector(adults_),
parent_selector(temp_population),
evaluator_(objectives_),
breeder_(parent_selector),
survival_(adults_, temp_population,/* child_,*/ store_),
weights_(objectives_.numObj(), 0),
evaluations_performed(0),
num_LS_performed(0),
time_taken(0)
{
  temp_selector.setNumSelections(tempSize);
  parent_selector.setNumSelections(2);
  breeder_.setCrossoverProb(crossProb);
  breeder_.setMutationProb(mutationProb);
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::reset()
{
  store_.reset();

  // Be careful to measure the time taken to create and evaluate the initial population. Most of the time this will just
  // be zero, but not so if we plan only a few generations on an expensive problem
  time_t beginTime;
  time(&beginTime);

  adults_.initialize(prob_, start_size);
  evaluations_performed = start_size;  // It is fairly safe to assume that the population size does not exceed the...
                                       // ...evaluation limit!
  // Add the newly created adults to the store (for this algorithm, all solutions are added to the store immediately...
  // ...after creation and evaluation.
  for (int i = 0; i < adults_.size(); ++i)
  {
    store_.add(adults_[i]);
  }

  // Apply local search to each of the solutions.
  for (int i = 0; i < adults_.size(); ++i)
  {
    // For now, we only handle two objectives, though the algorithm can handle more.
    assert(objectives_.numObj() == 2);
    weights_[0] = utils::rand::rand01();
    weights_[1] = 1 - weights_[0];

    LocalSearch localSearch(adults_[i], objectives_, store_, weights_, num_neighbours, evaluations_performed,
                            total_evaluations);
    localSearch.optimize();
  }
  num_LS_performed = start_size;

  time_t endTime;
  time(&endTime);
  time_taken = static_cast<int>(endTime - beginTime);
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::setEvaluations(int numEval)
{
  evaluations_limited = true;
  num_evaluations = numEval;
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::setTotalEvaluations(int numEval)
{
  evaluations_limited = true;
  total_evaluations = numEval;
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::step()
{
  // Perform one generation of the algorithm, which will involve an unknown number of evaluations

  // Randomize the weights
  weights_[0] = utils::rand::rand01();
  weights_[1] = 1 - weights_[0];

  // Evaluate the population
  evaluator_.setWeights(weights_);
  evaluator_.evaluate(adults_);

  // Fill the temporary population
  temp_population.clear();
  while (!temp_population.full())
  {
    temp_population.add(new Solution(*temp_selector.select()));  // I would like to be able to just pass what was...
  }                                                              // ...selected, but problems with 'const' occur

  // Create the child using the breeder which should already be referencing the temporary population
  child_ = breeder_.child();
  child_->evaluateObjectives();
  ++evaluations_performed;

  // Add child to the store. (Not in algorithm as described by Jaszkiewicz.)
  store_.add(child_);

  // Perform first found local search on the child, using the same weights as used to select the parents, accepting only
  // better solutions. The local search object also updates the number of evaluations performed.
  LocalSearch localSearch(child_, objectives_, store_, weights_, num_neighbours, evaluations_performed, total_evaluations);
  localSearch.optimize();
  ++num_LS_performed;

  // Survival of the fittest creates the new adult population - not a major change in this case, since the most that can
  // happen is one solution replacement.
  survival_.compete(child_);

  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations to be increased, allowing populations to be saved at various
  // points in the search
  if (!evaluations_limited && !time_limited)
  {
    cerr << "Attempting to run JMOGLS with no termination conditions. Aborting." << endl;
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
void JMOGLSAlgorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int JMOGLSAlgorithm<Problem, Solution, Store>::evaluationsPerformed() const
{
  return evaluations_performed;
}

template <class Problem, class Solution, class Store>
int JMOGLSAlgorithm<Problem, Solution, Store>::numLSPerformed() const
{
  return num_LS_performed;
}

template <class Problem, class Solution, class Store>
int JMOGLSAlgorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& JMOGLSAlgorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& JMOGLSAlgorithm<Problem, Solution, Store>::store() const
{
  // Function to return a constant reference to the store of best solutions, allowing these to be output
  // Not a constant member, since printing the store may require dominated solutions to be eliminated
  return store_;
}

template <class Problem, class Solution, class Store>
void JMOGLSAlgorithm<Problem, Solution, Store>::completeStore() const
{
  // Store is already complete, since solutions have been added upon creation. No need to call addWorthy, though it
  // wouldn't add anything anyway, since 'worthy()' always returns false.
}

#endif
