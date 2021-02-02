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

// Selectors are used to select members from a Population for breeding.
// Typically, fitter members are selected in preference to the less fit, though this need not be the case.

#ifndef SELECTORS_H
#define SELECTORS_H

#include <cassert>
#include <vector>
#include <algorithm>
#include <memory>
#include <rng.h>
#include "utils.h"
#include "population.h"

// (Both tournament and rank should be parameterized. (Tournament size and selection probability.))

template <class Solution>
class SortedTournamentSelector
{
  // Class to perform binary tournament selection on a sorted Population.
  // Currently has no checks to see if the Population is indeed sorted.
public:
  explicit SortedTournamentSelector(const Population<Solution>& population);
  void setNumEntrants(int newSize);
  std::shared_ptr<const Solution> select() const;

private:
  const Population<Solution>& population_;
  int num_entrants;
};


template <class Solution>
class TournamentSelector
{
  // Class to perform binary tournament selection on an unsorted Population.
  // A `fitterThan(lhs, rhs)' method must exist for comparing Solutions. (Note that this is NOT domination)
public:
  explicit TournamentSelector(const Population<Solution>& population);
  void setNumEntrants(int newSize);
  std::shared_ptr<const Solution> select() const;

private:
  const Population<Solution>& population_;
  int num_entrants;
};


template <class Solution>
class SortedProbabilisticSelector
{
  // Choose the best solution with a certain probability. If not selected, repeat with the next best, and then the next
  // and so on.
public:
  explicit SortedProbabilisticSelector(const Population<Solution>& population);
  void setProb(double newProb);
  std::shared_ptr<const Solution> select() const;

private:
  const Population<Solution>& population_;
  double prob_;
};


template <class Solution>
class StochasticUniversalSelector
{
  // A form of proportional representation used by MOGA. Whenever a solution is selected, a ball is rolled on a roulette
  // wheel. The sectors of the wheel correspond to the solutions, with widths proportional to the fitness. However, the
  // ball rolls only 1/N times the circumference of the wheel at each iteration, where N is the number of children to be
  // selected. The starting point of the ball is selected at random, but then the selector deterministically selects
  // equally spaced solutions.

  // Note: Fitness must be maximized, and appropriately scaled.
public:
  explicit StochasticUniversalSelector(const Population<Solution>& population);
  void setNumChildren(int newNumChildren);
  std::shared_ptr<const Solution> select();

private:
  const Population<Solution>& population_;
  int num_children;
  double increment_;
  double roulette_pos;
  double total_fitness;
  double fitness_so_far;
  int solution_num;
};


template <class Solution>
class ScaledRouletteSelector
{
  // Roulette selection, but with the fitness scaled so that the minimum fitness achieved is zero.
public:
  explicit ScaledRouletteSelector(const Population<Solution>& population);
  std::shared_ptr<const Solution> select();

private:
  void calc_selection_probs();
  int play_roulette();

private:
  const Population<Solution>& population_;

  std::vector<double> selection_probs;
};


template <class Solution>
class SortedBestSelector
{
  // Class that selects the best solutions according to the fitness value, one after another
public:
  explicit SortedBestSelector(const Population<Solution>& population);
  void setNumSelections(int numSelections);
  std::shared_ptr<const Solution> select();

private:
  const Population<Solution>& population_;
  int num_selections;
  int selection_num;
};


template <class Solution>
class UniformDifferentFewSelector
{
  // Class that selects solutions with equal probability, but does not select the same solution twice. After a set
  // number of different solutions has been selected, the selector resets ready for selecting a new set of different
  // solutions. Used for when few but different solutions are required, e.g. if we wish to ensure that the same solution
  // isn't selected twice for crossover.
public:
  explicit UniformDifferentFewSelector(const Population<Solution>& population);
  void setNumSelections(int numSelections);
  std::shared_ptr<const Solution> select();

private:
  const Population<Solution>& population_;
  int num_selections;
  int selection_num;
  std::vector<int> selected_solutions;
};


template <class Solution>
class UniformDifferentManySelector
{
  // Class that selects solutions with equal probability, but does not select the same solution twice. Used for when
  // many different solutions are required.
public:
  explicit UniformDifferentManySelector(const Population<Solution>& population);
  void setNumSelections(int numSelections);
  std::shared_ptr<const Solution> select();

private:
  const Population<Solution>& population_;
  int num_selections;
  int selection_num;
  std::vector<bool> selected_;   // Warning: vector<bool>
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
SortedTournamentSelector<Solution>::SortedTournamentSelector(const Population<Solution>& population) :
population_(population),
num_entrants(2)
{
}

template <class Solution>
void SortedTournamentSelector<Solution>::setNumEntrants(int newSize)
{
  num_entrants = newSize;
}

template <class Solution>
std::shared_ptr<const Solution> SortedTournamentSelector<Solution>::select() const
{
  // Hold a tournament of 'num_entrants' competitors and select the winner
  // Likely to be inefficient for large tournaments
  // 2018: The same solution might be picked twice - is this a problem?
  int best = utils::rand::randInt(0, population_.size());
  for (int i = 1; i < num_entrants; ++i)
  {
    best = std::min(best, utils::rand::randInt(0, population_.size()));
  }

  return population_[best];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
TournamentSelector<Solution>::TournamentSelector(const Population<Solution>& population) :
population_(population),
num_entrants(2)
{
}

template <class Solution>
void TournamentSelector<Solution>::setNumEntrants(int newSize)
{
  num_entrants = newSize;
}

template <class Solution>
std::shared_ptr<const Solution> TournamentSelector<Solution>::select() const
{
  // Hold a tournament of 'num_entrants' competitors and select the winner
  // Likely to be inefficient for large tournaments
  // 2018: The same solution might be picked twice - is this a problem?
  int best = utils::rand::randInt(0, population_.size());
  for (int i = 1; i < num_entrants; ++i)
  {
    int contender = utils::rand::randInt(0, population_.size());
    if (fitterThanP(population_[contender], population_[best]))  // Was just fitterThan
    {
      best = contender;
    }
  }

  return population_[best];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
SortedProbabilisticSelector<Solution>::SortedProbabilisticSelector(const Population<Solution>& population) :
population_(population),
prob_(0.5)
{
  // The repeated constructors seem to suggest that inheritance might be more suitable.
}

template <class Solution>
void SortedProbabilisticSelector<Solution>::setProb(double newProb)
{
  prob_ = newProb;
}

template <class Solution>
std::shared_ptr<const Solution> SortedProbabilisticSelector<Solution>::select() const
{
  // Select the best with probability, 0.5.
  // If not selected, select the second best with prob 0.5. (Gives 0.25 overall).
  // Repeat until a solution is selected or only the last one is left.
  int i;
  for (i = 0; i < population_.size() - 1; ++i)
  {
    if (utils::rand::rand01() < prob_)
    {
      return population_[i];
    }
  }

  return population_[i];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
StochasticUniversalSelector<Solution>::StochasticUniversalSelector(const Population<Solution>& population) :
population_(population),
num_children(0)
{
  // Calculate the total fitness of the population. This will be used to scale the fitnesses so that they sum to one.
  total_fitness = 0;
  for (auto i = population_.cbegin(); i != population_.cend(); ++i)  // (Ranged loop?)
  {
    const Solution& sol = **i;  // Iterator into a population of shared pointers
    total_fitness += sol.fitness();
  }
}

template <class Solution>
void StochasticUniversalSelector<Solution>::setNumChildren(int newNumChildren)
{
  // Not only set the number of children, but also find the size of the steps taken around
  // the roulette wheel and the starting position.
  num_children = newNumChildren;
  increment_ = total_fitness / num_children;  // Set wheel circumference to be the sum of the fitnesses of the solutions

  // Set the wheel ready to go, by setting the position of the ball and the solution pointer
  roulette_pos = utils::rand::randDouble(0, increment_) - increment_;  // Subtraction so that we can roll to the...
  solution_num = 0;                                                    // ...first position, then select.
  fitness_so_far = population_[solution_num]->fitness();
}

template <class Solution>
std::shared_ptr<const Solution> StochasticUniversalSelector<Solution>::select()
{
  assert(num_children > 0);

  // Roll the ball
  roulette_pos += increment_;

  // Find the selected solution. If the roulette position has gone beyond the current solution, increment the solution
  // iterator until the correct solution is found
  while (roulette_pos > fitness_so_far && solution_num < population_.size() - 1)  // Be careful with floating point...
  {                                                                               // ...arithmetic!
    ++solution_num;
    fitness_so_far += population_[solution_num]->fitness();
  }

  if (roulette_pos > fitness_so_far)
  {
    std::cerr << "Roulette ball has gone flying!" << std::endl;
  }

  return population_[solution_num];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
ScaledRouletteSelector<Solution>::ScaledRouletteSelector(const Population<Solution>& population):
population_(population),
selection_probs(population.capacity())  // Has to be capacity, because current size is zero!
{
}

template <class Solution>
void ScaledRouletteSelector<Solution>::calc_selection_probs()
{
  // Get the minimum weighted objective
  double minFitness = population_[0]->fitness();
  for (int i = 1; i < population_.size(); ++i)
  {
    minFitness = std::min(minFitness, population_[i]->fitness());
  }

  // Subtract the minimum from all the weighted_objectives then normalize to get the selection probabilities
  for (int i = 0; i < population_.size(); ++i)
  {
    selection_probs[i] = population_[i]->fitness() - minFitness;
  }
  normalize(selection_probs.begin(), selection_probs.end());
}

template <class Solution>
int ScaledRouletteSelector<Solution>::play_roulette()
{
  // Spin the wheel
  double roulettePos = utils::rand::rand01();
  int solutionNum = -1;
  while (roulettePos >= 0 && solutionNum < population_.size() - 1)
  {
    ++solutionNum;
    roulettePos -= selection_probs[solutionNum];  // Be careful with floating point arithmetic!
  }

  if (roulettePos >= 0)
  {
    std::cerr << "Roulette ball has gone flying!!" << std::endl;
  }

  return solutionNum;
}


template <class Solution>
std::shared_ptr<const Solution> ScaledRouletteSelector<Solution>::select()
{
  calc_selection_probs();
  int solutionNum = play_roulette();

  return population_[solutionNum];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
SortedBestSelector<Solution>::SortedBestSelector(const Population<Solution>& population):
population_(population),
num_selections(0),
selection_num(0)
{
}

template <class Solution>
void SortedBestSelector<Solution>::setNumSelections(int numSelections)
{
  num_selections = numSelections;
  selection_num = 0;
}

template <class Solution>
std::shared_ptr<const Solution> SortedBestSelector<Solution>::select()
{
  if (selection_num == num_selections)
  {
    selection_num = 0;
  }

  return population_[selection_num++];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
UniformDifferentFewSelector<Solution>::UniformDifferentFewSelector(const Population<Solution>& population):
population_(population),
num_selections(0),
selection_num(0)
{
}

template <class Solution>
void UniformDifferentFewSelector<Solution>::setNumSelections(int numSelections)
{
  num_selections = numSelections;
  selection_num = 0;
  selected_solutions.clear();
  selected_solutions.reserve(num_selections);
}

template <class Solution>
std::shared_ptr<const Solution> UniformDifferentFewSelector<Solution>::select()
{
  // Check to see if we are making a new set of selections - if so, reset the selection number and the array of selected
  // items
  if (num_selections < 1)
  {
    throw std::runtime_error("Set number of selection to a positive value before using UniformDifferentFewSelector.");
  }
  if (selection_num >= num_selections)
  {
    selection_num = 0;
    selected_solutions.resize(0);  // Visual C++: clear() deallocates memory, or at least used to.
  }

  // Find a random unselected item
  int i = utils::rand::randInt(0, population_.size());
  while (std::find(begin(selected_solutions), end(selected_solutions), i) != end(selected_solutions))
  {
    i = utils::rand::randInt(0, population_.size());
  }

  // Add to the list of selected solutions and return the solution
  selected_solutions.push_back(i);
  ++selection_num;

  return population_[i];
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
UniformDifferentManySelector<Solution>::UniformDifferentManySelector(const Population<Solution>& population):
population_(population),
num_selections(0),
selection_num(0),
selected_(population.capacity(), false)   // Using capacity just in case the population grows during the lifetime of...
{                                         // ...the selector.
}

template <class Solution>
void UniformDifferentManySelector<Solution>::setNumSelections(int numSelections)
{
  num_selections = numSelections;
  selection_num = 0;
}

template <class Solution>
std::shared_ptr<const Solution> UniformDifferentManySelector<Solution>::select()
{
  // Check to see if we are making a new set of selections - if so, reset the selection number and the array of selected
  // items
  if (selection_num == num_selections)
  {
    selection_num = 0;
    selected_.assign(population_.capacity(), false);
  }

  // Find a random unselected item
  int i = utils::rand::randInt(0, population_.size());
  while (selected_[i])
  {
    i = utils::rand::randInt(0, population_.size());
  }

  // Mark as selected and return the solution
  selected_[i] = true;
  ++selection_num;

  return population_[i];
}

#endif
