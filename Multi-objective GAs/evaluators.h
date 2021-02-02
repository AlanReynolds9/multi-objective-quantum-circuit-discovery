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

// The Evaluator classes for each of the algorithms. Note that these classes do NOT perform evaluateObjectives() on the
// solutions. Rather, given the objective values, these Evaluators peform the ranking, crowding, etc., required to
// complete evaluation.

#ifndef EVALUATORS_H
#define EVALUATORS_H

#include <vector>
#include <iostream>
#include <functional>
#include <algorithm>
#include <memory>
#include <rng.h>
#include "utils.h"
#include "objectives.h"
#include "population.h"
#include "crowding.h"

template <class Solution>
class ParetoRankEvaluator
{
  // Class to calculate the pareto ranks
public:
  void evaluate(Population<Solution>& population);
};

template <class Solution, class Crowding>
class NSGA2Evaluator
{
  // Class to calculate the fronts that each solution belongs to and the crowding distance
  using Front = std::vector<std::shared_ptr<Solution> >;

public:
  // (Considering the following constructor, can't we have the evaluator create the crowding object, rather than have us
  // pass it in as a constructor parameter?)
  NSGA2Evaluator(const Objectives& objectives, Crowding& crowding, bool standardDominance);

  void evaluate(Population<Solution>& population);

private:
  void initialize_dominance();
  void remove_solution(int solution);
  void identify_fronts();
  void two_objective_identify_fronts();
  void general_identify_fronts();
  void identify_best_front(Available& available, int frontNum);
  void calculate_crowding();

private:
  const Objectives& objectives_;
  Crowding& crowding_;
  bool standard_dominance{true};

  Population<Solution>*          population_{nullptr};
  std::vector<int>               num_superiors;  // Number of solutions that dominate solution i.
  std::vector<std::vector<int> > inferiors_;     // List of solutions that solution i dominates.
  std::vector<Front>             fronts_;        // Could we just number each solution instead?
};

template <class Solution, class Crowding>
class MOGAEvaluator
{
  // The MOGA evaluator does not extract fronts by ranks. The rank of each solution begin equal to the number of
  // dominating solutions
  using Rank = std::vector<std::shared_ptr<Solution> >;

public:
  explicit MOGAEvaluator(Crowding& crowding);

  void evaluate(Population<Solution>& population);

private:
  void calculate_ranks();
  void assign_basic_fitness();
  void calculate_crowding();

private:
  Population<Solution>* population_{nullptr};
  Crowding& crowding_;
  int max_rank{0};  // Not really necessary to initialize here, but added initialization for 'good form'.
  std::vector<Rank> ranks_;  // Indexed on rank, from 0 being the best.
};

template <class Solution, class Crowding>
class SPEA2Evaluator
{
  // Class to calculate the raw fitness of each solution and the crowding distance.
  // Survival and mating selection are different to each other, so sorting the resulting population might not be best.
public:
  explicit SPEA2Evaluator(Crowding& crowding);

  void evaluate(Population<Solution>& population);

private:
  void calculate_strengths();
  void calculate_unadjusted_fitness();
  void calculate_crowding();

private:
  Population<Solution>* population_{nullptr};
  Crowding& crowding_;
  std::vector<int> strength_;
};

template <class Solution>
class RandomWeightedEvaluator
{
  // Class to calculate the fitness of solutions based on a single weighted objective.
public:
  explicit RandomWeightedEvaluator(const Objectives& objectives);  // No crowding required, but we do need to know...
                                                                 // ...whether each objective is minimized or maximized.
  void evaluate(Population<Solution>& population);

  void outputWeights(std::ostream& out) const;

private:
  void randomize_weights();

private:
  const Objectives& objectives_;
  std::vector<double> weights_;
};

template <class Solution>
class WeightedEvaluator  // (Should we have both forms of weighted evaluator?)
{
  // Class to calculate the fitness of solutions based on a single weighted objective.
public:
  explicit WeightedEvaluator(const Objectives& objectives);  // No crowding, but we do need to know whether each...
                                                             // ...objective is minimized or maximized.
  WeightedEvaluator(const Objectives& objectives, const std::vector<double>& weights);

  void setWeights(const std::vector<double>& weights);

  void evaluate(Population<Solution>& population);

private:
  const Objectives& objectives_;
  std::vector<double> weights_;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
void ParetoRankEvaluator<Solution>::evaluate(Population<Solution>& population)
{
  // Reset the ranks
  for (auto i = population.begin(); i != population.end(); ++i)
  {
    auto sol = *i;  // (CHECK: Should be a std::shared_ptr<Solution>?)
    sol->resetRank(); // (Requires the solution class to know about ranks. Doesn't strike me as the right place for it.)
  }

  // Calculate the Pareto ranks of the solutions in the population
  for (auto i = population.begin(); i != population.end(); ++i)
  {
    auto sol1 = *i;
    sol1->incrementRank();
    for (auto j = population.begin(); j != population.end(); ++j)
    {
      auto sol2 = *j;
      if (i != j && sol2->beats(*sol1))
      {
        sol1->incrementRank();
      }
    }
    sol1->setRankEvaluated();
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Crowding>
NSGA2Evaluator<Solution, Crowding>::NSGA2Evaluator(const Objectives& objectives, Crowding& crowding,
                                                   bool standardDominance) :
objectives_(objectives),
crowding_(crowding),
standard_dominance(standardDominance)
{
  // (The standard NSGA II algorithm classes (NSGA2Algorithm and NSGA2GPAlgorithm) assume that standard dominance is
  // used. Otherwise the NSGA II crowding algorithm is likely to lead to rather odd results. However, using alternative
  // crowding methods allows us to use modified dominance relations. In this case, the NSGA2Evaluator class needs to
  // know that non-standard dominance is used in order to prevent it using more efficient non-dominated sorting
  // algorithms that only work with standard dominance.)
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::evaluate(Population<Solution>& population)
{
  population_ = &population;
  fronts_.clear();

  // Extract fronts of non-dominated solutions
  identify_fronts();

  // Determine how crowded each solution is in its front
  calculate_crowding();

  // Sort each front so that the least crowded solutions come first
  for (auto i = begin(fronts_); i != end(fronts_); ++i)
  {
    // We don't want any residual sorting of solutions according to the last objective to survive through the sort on
    // crowding, on to the survival stage, as this may create bias towards one particular part of the front.
    shuffle(i->begin(), i->end(), utils::rand::rng);
    sort(i->begin(), i->end(), lessCrowdedP<std::shared_ptr<Solution> >);
  }

  // Rearrange the population so that the fittest come first
  auto destIt = population.begin();
  for (auto i = begin(fronts_); i != end(fronts_); ++i)
  {
    destIt = copy(i->begin(), i->end(), destIt);
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::initialize_dominance()
{
  // Function that calculates, for each solution, the number that dominate it and a list of those it dominates.
  // (For large populations, this takes a lot of time! Consider both algorithmic improvements and parallelization.)
  
  // Start by resetting the values in each of the vectors
  num_superiors.assign(population_->size(), 0);
  inferiors_.assign(population_->size(), std::vector<int>());
  
  // For each combination of i and j, see whether solution i dominates solution j. If so, increase 'num_superiors' for j
  // and add j to i's dominance list.
  for (int i = 0; i < population_->size(); ++i)
  {
    auto& sol1 = (*population_)[i];
    for (int j = i + 1; j < population_->size(); ++j)
    {
      auto& sol2 = (*population_)[j];
      if (sol1->beats(*sol2))
      {
        ++num_superiors[j];
        inferiors_[i].push_back(j);
      }
      else if (sol2->beats(*sol1))
      {
        ++num_superiors[i];
        inferiors_[j].push_back(i);
      }
    }
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::remove_solution(int solution)
{
  for (auto i = inferiors_[solution].begin(); i != inferiors_[solution].end(); ++i)
  {
    --num_superiors[*i];
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::identify_fronts()
{
  if (standard_dominance && objectives_.numObj() == 2)  // (What about only one objective?)
  {
    two_objective_identify_fronts();  // Only works with two objectives and 'standard' dominance
  }
  else
  {
    general_identify_fronts();  // Should work with any sensible dominance relation
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::two_objective_identify_fronts()
{
  // First sort the solutions according to the first objective, breaking ties according to the second.
  LexBetterByObjective<Solution> betterByObjective(objectives_, 0);
  sort(population_->begin(), population_->end(), betterByObjective);

  // Add the first solution, which must be non-dominated, to the first front. 
  Front newFront;
  fronts_.push_back(newFront);
  auto solution = (*population_)[0];
  solution->setRank(0);  // To allow use of TournamentSelector, rather than SortedTournamentSelector.
  fronts_[0].push_back(solution);

  // We will also keep track of the best value of objective(1) in each front.
  std::vector<double> bests;
  bests.push_back(solution->objective(1));

  // Also keep a note of the last solution and the front to which it was added, in case the next solution is identical
  // in both objectives.
  auto lastSolution = solution;
  int frontNum = 0;

  // Scan through the rest of the solutions, adding to the appropriate front.
  for (int sol = 1; sol < population_->size(); ++sol)
  {
    solution = (*population_)[sol];

    // If there is a tie on both objectives, we add the solution to the same front as last time. Otherwise, we find the
    // front to which the solution belongs by finding the first front who's 'best' solution, according to objective(1),
    // is worse than this solution.
    // (Should we have a 'ties' function?)
    if (solution->objective(0) != lastSolution->objective(0) || solution->objective(1) != lastSolution->objective(1))
    {
      // Find the first front with 'best' value strictly worse than objective(1) for this solution.
      std::vector<double>::const_iterator i;
      if (objectives_.maximized(1))
      {
        i = std::upper_bound(bests.cbegin(), bests.cend(), solution->objective(1), std::greater<double>());
      }
      else
      {
        i = std::upper_bound(bests.cbegin(), bests.cend(), solution->objective(1), std::less<double>());
      }
      frontNum = static_cast<int>(i - bests.begin());
    }

    // If necessary, create a new front.
    if (frontNum == static_cast<int>(fronts_.size()))
    {
      Front newFront;
      fronts_.push_back(newFront);
      bests.push_back(0);
    }

    // Add solution to the front.
    solution->setRank(frontNum);  // NEW: To allow use of TournamentSelector (rather than sorted version).
    fronts_[frontNum].push_back(solution);
    bests[frontNum] = solution->objective(1);
    lastSolution = solution;
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::general_identify_fronts()
{
  // Calculate, for each solution, the number that dominate it. Also create a list of the solutions it dominates.
  initialize_dominance();

  // Keep a record of the solutions that are left
  Available available(population_->size());

  // While there are some, extract fronts
  // (Added use of frontNum so that we can label Solutions with their rank. This, in turn, is so that we can use
  // SortedTournamentSelector, rather than assuming that the population is sorted. Finally, this is so that we can use a
  // different Survival class - specifically SortedBestUniqueSurvival.)
  int frontNum = 0;
  while (!available.empty())
  {
    identify_best_front(available, frontNum++);
  }
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::identify_best_front(Available& available, int frontNum)
{
  // Create an array of ints to hold the indices of those solutions that will be in the front
  std::vector<int> newFrontIndices;

  auto i = available.begin();
  while (i != available.end())
  {
    // If the solution is non-dominated, place it in the front, removing it from the available list.
    int solution = *i;
    if (num_superiors[solution] == 0)
    {
      newFrontIndices.push_back(solution);
      i = available.erase(i);
    }
    else
    {
      ++i;
    }
  }

  if (newFrontIndices.empty())
  {
    std::cerr << "  Failed to find front. Aborting!" << std::endl;
    exit(EXIT_FAILURE);
  }

  // We have found the front, storing the indices of the solutions in a vector. What we really want is a Front of shared
  // pointers to solutions
  Front newFront;
  for (auto i = newFrontIndices.begin(); i != newFrontIndices.end(); ++i)
  {
    auto solutionPtr = (*population_)[*i];
    solutionPtr->setRank(frontNum);  // (To allow use of unsorted Selectors.)
    newFront.push_back(solutionPtr);
    
    // Remove the solution in the new front from the population, updating the dominance relations.
    remove_solution(*i);
  }

  // Finally, add the new front to the front list  // (We need to give each solution the front number!!!)
  fronts_.push_back(newFront);
}

template <class Solution, class Crowding>
void NSGA2Evaluator<Solution, Crowding>::calculate_crowding()
{
  // Calculate the crowding in each front
  // HERE: Must eventually rewrite so that this is only performed as necessary
  for (auto& front : fronts_)
  {
    crowding_.calculateCrowding(front.begin(), front.end());
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Crowding>
MOGAEvaluator<Solution, Crowding>::MOGAEvaluator(Crowding& crowding) :
crowding_(crowding)
{
}

template <class Solution, class Crowding>
void MOGAEvaluator<Solution, Crowding>::evaluate(Population<Solution>& population)
{
  population_ = &population;

  calculate_ranks();
  assign_basic_fitness();
  calculate_crowding();  // In the case of MOGA, this adjusts the fitness values
}

template <class Solution, class Crowding>
void MOGAEvaluator<Solution, Crowding>::calculate_ranks()
{
  auto& pop = *population_;
  auto popSize = pop.size();

  // Create a temporary vector for storing the ranks of each solution
  std::vector<int> rank(popSize, 0);

  // For each i that dominates j, increase the rank of j
  for (int i = 0; i < popSize; ++i)
  {
    for (int j = 0; j < popSize; ++j)
    {
      if (i != j && pop[i]->beats(*pop[j]))
      {
        ++rank[j];
      }
    }
  }

  // Find the maximum rank
  max_rank = 0;
  for (int i = 0; i < popSize; ++i)
  {
    max_rank = std::max(max_rank, rank[i]);
  }

  // Create a vector of solutions of each rank
  ranks_.assign(max_rank + 1, Rank());
  for (int i = 0; i < popSize; ++i)
  {
    ranks_[rank[i]].push_back(pop[i]);
  }

  // Solutions should now all be grouped according to rank
}

template <class Solution, class Crowding>
void MOGAEvaluator<Solution, Crowding>::assign_basic_fitness()
{
  auto& pop = *population_;
  auto popSize = pop.size();

  int numBetter = 0;
  for (int rank = 0; rank <= max_rank; ++rank)
  {
    int numOfRank = static_cast<int>(ranks_[rank].size());
    double fitnessOfRank = popSize - numBetter - (numOfRank - 1.0) / 2.0;
    for (auto i = ranks_[rank].begin(); i != ranks_[rank].end(); ++i)
    {
      Solution& sol = **i;         // i is an iterator to a shared pointer
      sol.setFitness(fitnessOfRank);
    }

    numBetter += numOfRank;
  }
}

template <class Solution, class Crowding>
void MOGAEvaluator<Solution, Crowding>::calculate_crowding()
{
  // Calculate the crowding for each rank
  // HERE: Must eventually rewrite so that this is only performed as necessary
  for (auto& rank : ranks_)
  {
    crowding_.calculateCrowding(rank.begin(), rank.end());
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Crowding>
SPEA2Evaluator<Solution, Crowding>::SPEA2Evaluator(Crowding& crowding) :
crowding_(crowding)
{
}

template <class Solution, class Crowding>
void SPEA2Evaluator<Solution, Crowding>::evaluate(Population<Solution>& population)
{
  population_ = &population;

  calculate_strengths();
  calculate_unadjusted_fitness();
  calculate_crowding();           // This is probably the messy one, as crowding is used in two different ways
}

template <class Solution, class Crowding>
void SPEA2Evaluator<Solution, Crowding>::calculate_strengths()
{
  // For each solution, count the number it dominates
  strength_.assign(population_->size(), 0);

  for (int i = 0; i < population_->size(); ++i)
  {
    auto sol1 = (*population_)[i];
    for (int j = 0; j < population_->size(); ++j)
    {
      if (i != j)
      {
        auto sol2 = (*population_)[j];
        if (sol1->beats(*sol2))
        {
          ++strength_[i];
        }
      }
    }
  }
}

template <class Solution, class Crowding>
void SPEA2Evaluator<Solution, Crowding>::calculate_unadjusted_fitness()
{
  // For each solution, sum the strengths of those solutions that dominate it
  // The results are stored in the solutions, which does not seem ideal. However, the SPEA2 survival class needs access
  // to the unadjusted fitnesses, as well as the distances to each other solution. (Yick!)
  std::vector<int> unadjustedFitness(population_->size(), 0);

  for (int i = 0; i < population_->size(); ++i)
  {
    auto sol1 = (*population_)[i];
    for (int j = 0; j < population_->size(); ++j)
    {
      if (i != j)
      {
        auto sol2 = (*population_)[j];
        if (sol1->beats(*sol2))
        {
          unadjustedFitness[j] += strength_[i];
        }
      }
    }
  }

  // Copy unadjusted fitnesses into the solutions
  for (int i = 0; i < population_->size(); ++i)
  {
    auto sol = (*population_)[i];
    sol->setFitness(unadjustedFitness[i]);
  }
}

template <class Solution, class Crowding>
void SPEA2Evaluator<Solution, Crowding>::calculate_crowding()
{
  // Calculate the crowding over the entire population
  // (We must eventually rewrite so that this is only performed as necessary.)
  crowding_.calculateCrowding(population_->begin(), population_->end());
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
RandomWeightedEvaluator<Solution>::RandomWeightedEvaluator(const Objectives& objectives) :
objectives_(objectives),
weights_(objectives.numObj())
{
}

template <class Solution>
void RandomWeightedEvaluator<Solution>::randomize_weights()
{
  // Create a random set of weights, normalized so that they add to one.
  for (int obj = 0; obj < static_cast<int>(objectives_.numObj()); ++obj)
  {
    weights_[obj] = utils::rand::rand01();
  }
  normalize(weights_.begin(), weights_.end());
}

template <class Solution>
void RandomWeightedEvaluator<Solution>::evaluate(Population<Solution>& population)
{
  // Given a random set of objective weights, calculate the weighted objective values for each solution
  randomize_weights();

  for (int i = 0; i < population.size(); ++i)
  {
    double fitness = 0;
    for (int obj = 0; obj < static_cast<int>(objectives_.numObj()); ++obj)
    {
      fitness += weights_[obj] * objectives_.scaleObjective(obj, population[i]->objective(obj));  // (Might prefer to...
    }                                                                                    // ...scale wrt the population.

    population[i]->setFitness(fitness);
  }
}

template <class Solution>
void RandomWeightedEvaluator<Solution>::outputWeights(std::ostream& out) const
{
  out << "Weights:" << std::endl;
  for (int obj = 0; obj < static_cast<int>(objectives_.numObj()); ++obj)   // Use copy and a stream iterator instead
  {
    out << weights_[obj] << std::endl;
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
WeightedEvaluator<Solution>::WeightedEvaluator(const Objectives& objectives) :
objectives_(objectives),
weights_(objectives.numObj(), 1/objectives.numObj())  // Set all weights to be equal by default
{
}

template <class Solution>
WeightedEvaluator<Solution>::WeightedEvaluator(const Objectives& objectives, const std::vector<double>& weights) :
objectives_(objectives),
weights_(weights)
{
}

template <class Solution>
void WeightedEvaluator<Solution>::setWeights(const std::vector<double>& weights)
{
  weights_ = weights;
}

template <class Solution>
void WeightedEvaluator<Solution>::evaluate(Population<Solution>& population)
{
  // Given the providedset of objective weights, calculate the weighted objective values for each solution
  for (int i = 0; i < population.size(); ++i)
  {
    double fitness = 0;
    for (int obj = 0; obj < static_cast<int>(objectives_.numObj()); ++obj)
    {
      fitness += weights_[obj] * objectives_.scaleObjective(obj, population[i]->objective(obj));  // NOTE: Might...
    }                                                                          // ...prefer to scale wrt the population.

    population[i]->setFitness(fitness);
  }
}

#endif
