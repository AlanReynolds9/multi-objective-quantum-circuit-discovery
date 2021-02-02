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

#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <cassert>
#include <vector>
#include <algorithm>
#include <memory>
#include "objectives.h"

const double PI = 3.14159265358979323846;

class LesserAngleP
{
  // Function object to determine whether the scaled objectives for one solution, when plotted, make a lesser angle to
  // the first objective axis than those for a second solution. This should only be used with two objectives. Works with
  // pointers to solutions to allow sorting of populations.
public:
  template <class Pointer>
  bool operator()(Pointer lhs, Pointer rhs) const
  {
    return lhs->angle() < rhs->angle();
  }
};

class WeightedObjectiveImproves
{
  // Function object to determine whether one solution is better than another, according to a weighted objective.
  // Objectives are scaled before combination into the weighted objective, so that the resulting problem is a
  // maximization problem, with objectives from 0 to 1.
public:
  WeightedObjectiveImproves(const Objectives& objectives, std::vector<double> weights) :
  objectives_(objectives),
  weights_(weights)
  {
  }

  template <class Solution>
  bool operator()(const Solution& lhs, const Solution& rhs) const
  {
    double lhsValue = 0;
    double rhsValue = 0;
    for (int obj = 0; obj < static_cast<int>(objectives_.numObj()); ++obj)
    {
      lhsValue += objectives_.scaleObjective(obj, lhs.objective(obj)) * weights_[obj];
      rhsValue += objectives_.scaleObjective(obj, rhs.objective(obj)) * weights_[obj];
    }

    return lhsValue > rhsValue;
  }

private:
  const Objectives& objectives_;
  std::vector<double> weights_;
};

template <class Solution, class Store>
class AngledLocalSearch
{
  // Class that handles the second form of local search employed by Ishibuchi and Murata. Only works with two
  // objectives. Plotting the first objective along the x-axis and the second along the y-axis, each solution subtends a
  // certain angle to the x-axis. The population on which the search is performed is sorted according to this angle.
  // Then equally spaced angles are used to give the direction of the local search for each solution. It is assumed that
  // the objectives are scaled, so that they are all maximized and take values from 0 to 1.
public:
  AngledLocalSearch(Population<Solution>& population, const Objectives& objectives, Store& store, int numNeighbours,
                    int& evaluationsPerformed, int evaluationLimit);

  void optimize();

private:
  void calc_directions();
  bool single_step(int sol, const WeightedObjectiveImproves& better);
  void single_run(int sol);

private:
  Population<Solution>& population_;
  const Objectives& objectives_;
  Store& store_;
  int& evaluations_performed;
  int evaluation_limit;

  int num_neighbours;
  std::vector<std::vector<double> > weights_;  // Indexed on solution number, then objective
  std::vector<std::shared_ptr<Solution> > tried_;
};

template <class Solution, class Store>
class WeightedLocalSearch
{
  // Class that handles the local search required by Jaszkiewicz's version of genetic local search, though with the
  // modification that only a certain number of neighbours are examined. This is applied to a single solution and
  // objective weights are provided
public:
  WeightedLocalSearch(std::shared_ptr<Solution>& solution, const Objectives& objectives, Store& store,
                      const std::vector<double>& weights, int numNeighbours, int& evaluationsPerformed,
                      int evaluationLimit);

  void optimize();

private:
  bool single_step(const WeightedObjectiveImproves& better);

private:
  std::shared_ptr<Solution>& solution_;  // A reference means that we can copy the pointer to the new solution, rather
  // than the whole solution
  const Objectives& objectives_;
  Store& store_;
  std::vector<double> weights_;
  int& evaluations_performed;
  int evaluation_limit;

  int num_neighbours;
  std::vector<std::shared_ptr<Solution> > tried_;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store>
void AngledLocalSearch<Solution, Store>::calc_directions()
{
  for (int i = 0; i < population_.size(); ++i)
  {
    double angle = PI / 2 / (population_.size() - 1) * i;
    weights_[i][0] = cos(angle);
    weights_[i][1] = sin(angle);
  }
}

template <class Solution, class Store>
AngledLocalSearch<Solution, Store>::AngledLocalSearch(Population<Solution>& population, const Objectives& objectives,
                                                      Store& store, int numNeighbours, int& evaluationsPerformed,
                                                      int evaluationLimit) :
population_(population),
objectives_(objectives),
store_(store),
evaluations_performed(evaluationsPerformed),
evaluation_limit(evaluationLimit),
num_neighbours(numNeighbours),
weights_(population.size(), std::vector<double>(objectives.numObj()))
{
  assert(objectives.numObj() == 2);

  tried_.reserve(num_neighbours);
  calc_directions();
}

template <class Solution, class Store>
bool AngledLocalSearch<Solution, Store>::single_step(int sol, const WeightedObjectiveImproves& better)
{
  tried_.resize(0);  // Visual C++: clear deallocates memory

  // Need to ensure that the same move is not selected more than once. Also need to keep track of the number of...
  // ...evaluations. (This is SO hacky! Sort this out if we ever use this again.)
  for (int i = 0; i < num_neighbours && evaluations_performed < evaluation_limit; ++i)
  {
    // Get a new neighbouring solution and evaluate
    std::shared_ptr<Solution> neighbour;
    int attempts = 0;
    do
    {
      neighbour = make_shared<Solution>(*population_[sol]);
      neighbour->randomMove();
      ++attempts;
    }
    while (std::find_if(tried_.begin(), tried_.end(),
                        PointerEqual<std::shared_ptr<Solution> >(neighbour)) != tried_.end() &&
           attempts <= 5000);

    neighbour->evaluateObjectives();
    ++evaluations_performed;

    // Make a note of the neighbour to avoid retrying.
    tried_.push_back(neighbour);

    // Add neighbour to the store. (Not in algorithm as described by Ishibuchi and Murata - they only add the solution
    // after completion of LS.)
    store_.add(neighbour);  // Must be after evaluation, or chaos in the store ensues!

    if (better(*neighbour, *population_[sol]))
    {
      // Found a better neighbour.
      population_[sol] = neighbour;
      return true;
    }
  }

  // Failed to find a better neighbour in the permitted number of evaluations.
  return false;
}


template <class Solution, class Store>
void AngledLocalSearch<Solution, Store>::single_run(int sol)
{
  // Perform single steps of the search, until either an attempt to make a step fails or we reach the evaluation limit.
  WeightedObjectiveImproves better(objectives_, weights_[sol]);
  while (evaluations_performed < evaluation_limit && single_step(sol, better))
  {
  }
}

template <class Solution, class Store>
void AngledLocalSearch<Solution, Store>::optimize()
{
  // Sort the solutions in order of angle with the x-axis, where x is the first objective.

  // First determine the angles
  // (Couldn't we get the MOGLS solution to calculate these when evaluating the objectives? (Would require access to the
  // Objectives object.))
  // (Use of atan2 is not efficient - we only use the angles for the purpose of sorting the solutions.)
  // (Used to be done by using a function object that compared tangents - however, it is difficult to see how to deal
  // with the case where both objectives are zero.)
  for (int i = 0; i < population_.size(); ++i)
  {
    double x = objectives_.scaleObjective(0, population_[i]->objective(0));
    double y = objectives_.scaleObjective(1, population_[i]->objective(1));
    if (x == 0 && y == 0)
    {
      population_[i]->setAngle(utils::rand::randDouble(0, PI / 2));
    }
    else
    {
      population_[i]->setAngle(atan2(y, x));
    }
  }

  sort(population_.begin(), population_.end(), LesserAngleP());

  for (int sol = 0; sol < population_.size() && evaluations_performed < evaluation_limit; ++sol)
  {
    single_run(sol);
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store>
WeightedLocalSearch<Solution, Store>::WeightedLocalSearch(std::shared_ptr<Solution>& solution,
                                                          const Objectives& objectives, Store& store,
                                                          const std::vector<double>& weights, int numNeighbours,
                                                          int& evaluationsPerformed, int evaluationLimit) :
solution_(solution),
objectives_(objectives),
store_(store),
weights_(weights),
evaluations_performed(evaluationsPerformed),
evaluation_limit(evaluationLimit),
num_neighbours(numNeighbours)
{
  tried_.reserve(num_neighbours);
}

template <class Solution, class Store>
bool WeightedLocalSearch<Solution, Store>::single_step(const WeightedObjectiveImproves& better)
{
  tried_.resize(0);  // Visual C++: clear deallocates memory

  // (This is SO hacky! Sort this out if we ever use this again.)
  for (int i = 0; i < num_neighbours && evaluations_performed < evaluation_limit; ++i)
  {
    // Get a new neighbouring solution and evaluate
    std::shared_ptr<Solution> neighbour;
    int attempts = 0;
    do
    {
      neighbour = make_shared<Solution>(*solution_);
      neighbour->randomMove();
      ++attempts;
    }
    while (std::find_if(tried_.begin(), tried_.end(),
                        PointerEqual<std::shared_ptr<Solution> >(neighbour)) != tried_.end() &&
           attempts <= 10000);

    neighbour->evaluateObjectives();
    ++evaluations_performed;

    // Make a note of the neighbour to avoid retrying.
    tried_.push_back(neighbour);

    // Add neighbour to the store. (Not in algorithm as described by Jaszkiewicz - they only add the solution after...
    // ...completion of LS.)
    store_.add(neighbour);  // Must be after evaluation, or chaos in the store ensues!

    if (better(*neighbour, *solution_))
    {
      // Found a better neighbour.
      solution_ = neighbour;
      return true;
    }
  }

  // Failed to find a better neighbour in the permitted number of evaluations.
  return false;
}


template <class Solution, class Store>
void WeightedLocalSearch<Solution, Store>::optimize()
{
  // Perform single steps of the search, until either an attempt to make a step fails or we reach the evaluation limit.
  WeightedObjectiveImproves better(objectives_, weights_);
  while (evaluations_performed < evaluation_limit && single_step(better))
  {
  }
}

#endif
