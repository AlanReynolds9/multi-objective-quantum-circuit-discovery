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

#ifndef QCGA_ALGORITHM
#define QCGA_ALGORITHM

#include "population.h"
#include "selectors.h"
#include "breeders.h"
#include "evaluators.h"
#include "QCGA2Survival.h"
#include "stores.h"

// This second attempt is basically just the crowded NSGA2 of the parameterized gate code, but with different adult and
// child population sizes and 'always on' mutation and crossover.

template <class Problem, class Solution, class Store>
class QCGA2Algorithm
{
public:
  QCGA2Algorithm(Problem& prob, Store& store, int adultPopSize, int childPopSize, double spacing, int numRand);

  // Initialization or reset to generation 0
  void reset();

  // Changing the settings
  void setGenerations(int numGen);
  void setTimeLimit(int timeLimit);

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
  // 'Typedefs' for the selector, breeder and survival types
  using Selector = UniformDifferentFewSelector<Solution>;
  using Breeder = SingleChildBreeder<Solution, Selector>;
  using Crowding = NSGA2Crowding<Solution>;
  using Evaluator = NSGA2Evaluator<Solution, Crowding>;
  using Survival = QCGA2Survival<Solution, Store, Evaluator>;

  // Problem
  Problem& prob_;  // Was const, but we now need to use prob_.prepareForStep();  (Create a Workspace class?)

  // Store
  Store& store_;

  // Termination conditions
  bool generations_limited;
  int num_generations;
  bool time_limited;
  int time_limit;

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
QCGA2Algorithm<Problem, Solution, Store>::QCGA2Algorithm(Problem& prob, Store& store, int adultPopSize,
                                                         int childPopSize, double spacing, int numRand) :
prob_(prob),
store_(store),
generations_limited(false),
time_limited(false),
adults_(adultPopSize),
children_(childPopSize),
objectives_(prob),
selector_(adults_),
crowding_(objectives_),
evaluator_(objectives_, crowding_, true),  // In standard NSGA II we assume that standard dominance is used. Hence...
breeder_(selector_),                       // ...the 'true'.
survival_(adults_, children_, store_, evaluator_, objectives_, spacing, numRand),
generation_(0),
time_taken(0)
{
  selector_.setNumSelections(2);  // For use with UniformDifferentFewSelector.
  breeder_.setCrossoverProb(1.1);  // Always perform crossover, which merely tells the child who its second parent is...
  breeder_.setMutationProb(1.1);   // ...and mutation, which determines the details of the mutation or crossover...
}                                  // ...operation to be performed.

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::reset()
{
  generation_ = 0;
  store_.reset();

  // Be careful to measure the time taken to create and evaluate the initial population. Most of the time this will just
  // be zero, but not so if we plan only a few generations on an expensive problem.
  time_t beginTime;
  time(&beginTime);
  adults_.initialize(prob_);
  evaluator_.evaluate(adults_);
  time_t endTime;
  time(&endTime);
  time_taken = static_cast<int>(endTime - beginTime);
}

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::setGenerations(int numGen)
{
  generations_limited = true;
  num_generations = numGen;
}

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::setTimeLimit(int timeLimit)
{
  time_limited = true;
  time_limit = timeLimit;
}

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::step()
{
  // Perform one iteration of the algorithm

  // Increment generation count.
  ++generation_;

  // Perform whatever problem specific tasks must be performed each step.
  prob_.prepareForStep();  // (I will want to be able to pass a gate sequence, solution or set of solutions to this...
                           // ...function.)
  // Fill the child population
  children_.clear();
  while (!children_.full())
  {
    std::shared_ptr<Solution> child(breeder_.child());
    children_.add(child);
  }
  children_.evaluateObjectives();  // These should be 'quick' evaluations. (Solution::evaluateObjectives() should be...
                                   // ...able to work this out.)
  // Survival of the fittest creates the new adult population
  survival_.compete();

  adults_.evaluateObjectives();  // These should be 'full' evaluations, if the solution has not already been fully...
                                 // ...evaluated.
  // Ensure that the breeder does not have any children left inside it
  breeder_.reset();
}

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::go()
{
  // Start without reset. Allows the number of generations or the time limit to be increased, allowing populations to be
  // saved at various points in the search
  if (!generations_limited && !time_limited)
  {
    std::cerr << "Attempting to run QCGA with no termination conditions. Aborting." << std::endl;
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
void QCGA2Algorithm<Problem, Solution, Store>::run()
{
  // Restart
  reset();
  go();
}

template <class Problem, class Solution, class Store>
int QCGA2Algorithm<Problem, Solution, Store>::generationsPerformed() const
{
  return generation_;
}

template <class Problem, class Solution, class Store>
int QCGA2Algorithm<Problem, Solution, Store>::timeTaken() const
{
  return time_taken;
}

template <class Problem, class Solution, class Store>
const Population<Solution>& QCGA2Algorithm<Problem, Solution, Store>::adults() const
{
  // Function to return a constant reference to the adult population, allowing for the output of the final population of
  // solutions.
  return adults_;
}

template <class Problem, class Solution, class Store>
Store& QCGA2Algorithm<Problem, Solution, Store>::store() const
{
  // Function to return a constant reference to the store of best solutions, allowing these to be output
  // Not a constant member, since printing the store may require dominated solutions to be eliminated
  return store_;
}

template <class Problem, class Solution, class Store>
void QCGA2Algorithm<Problem, Solution, Store>::completeStore() const
{
  addWorthy(store_, adults_.cbegin(), adults_.cend());
}

#endif
