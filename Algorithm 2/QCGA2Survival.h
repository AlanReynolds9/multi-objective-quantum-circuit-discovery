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

#ifndef QCGA2SURVIVAL_H
#define QCGA2SURVIVAL_H

#include <algorithm>
#include "utils.h"
#include <rng.h>
#include "population.h"
#include "distances.h"
#include "stores.h"
#include "nonDominatedStrimmer.h"


template <class Solution, class Store, class Evaluator>
class QCGA2Survival
{
  // As for BestSurvival, but uses an Evaluator that sorts the population in order of quality. It is therefore
  // unnecessary to resort. (Of course, the Evaluator selected must be one that sorts the population.) Used to be called
  // SortedBestCrowdedSurvival.
public:
  using Problem = typename Solution::Problem;

  QCGA2Survival(Population<Solution>& adults, Population<Solution>& children, Store& store, Evaluator& evaluator,
                const Objectives& objectives, double spacing, int numRand);

  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
  ObjectiveNorm<Solution, 1> distance_;
  double preferred_spacing;
  int num_rand;  // Number of solutions to be selected at random, rather than via the usual Survival methods.
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store, class Evaluator>
QCGA2Survival<Solution, Store, Evaluator>::QCGA2Survival(Population<Solution>& adults, Population<Solution>& children,
                                                         Store& store, Evaluator& evaluator,
                                                         const Objectives& objectives, double spacing, int numRand):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator),
distance_(objectives),
preferred_spacing(spacing),
num_rand(numRand)
{
}

template <class Solution, class Store, class Evaluator>
void QCGA2Survival<Solution, Store, Evaluator>::compete()
{
  // Create a merged population, let the evaluator sort it, then take solutions from the front, ignoring solutions that
  // are too close, in objective space, to already selected solutions.

  // Create a merged population, let the Evaluator evaluate the solutions and place them in quality order.
  Population<Solution> pool(adults_, children_);
  evaluator_.evaluate(pool);

  // Set the current spacing limit to the preferred value set by the user.
  auto currentSpacing = preferred_spacing;

  // Keep a count of the number of solutions transferred to the new adult population, and a list (in quality order) of
  // those solutions discarded due to crowding.
  int numTransferred{0};
  auto numNonRand = adults_.size() - num_rand;  // Number of solutions to be selected by standard NSGA2-like methods.
  Population<Solution> discards(adults_.size() + children_.size());  // Over-reserving (but safe).

  // Main loop
  auto i = pool.begin();
  while (numTransferred < numNonRand)
  {
    for (i = pool.begin(); numTransferred < numNonRand && i != pool.end(); ++i)
    {
      auto& candidate = *i;
      bool crowded{false};
      for (auto j = 0; j < numTransferred && !crowded; ++j)
      {
        if (distance_(*candidate, *adults_[j]) < currentSpacing)
        {
          crowded = true;
        }
      }

      if (crowded)
      {
        discards.add(candidate);
      }
      else if (!candidate->bad())
      {
        adults_[numTransferred++] = candidate;
      }
    }

    // Unable to find sufficiently many solutions using the current spacing restriction. Reduce crowding distance by a
    // factor of 2. (The extra bit that subtracts 0.000001 ensures that the crowding distance eventually becomes non-
    // positive, meaning that it no longer applies. Hence this eventually guarantees that the population can be filled.)
    if (numTransferred < numNonRand)
    {
      currentSpacing = currentSpacing / 2 - 0.000001;
      std::cout << "Failed to fill adult population. Reducing spacing threshold. New spacing = " << currentSpacing <<
                   "." << std::endl;
      pool = discards;  // (Adding a move constructor to the Population class might simplify this and make it more...
      discards.clear(); // ...efficient.)
    }
  }

  // We have transferred the good solutions to the adult population. We now transfer a random selection of the remaining
  // solutions. At present, this ignores spacing. We may wish to consider spacing in future.
  // First merge the remaining solutions into a single vector
  std::vector<std::shared_ptr<Solution> > remaining(discards.begin(), discards.end());
  remaining.insert(remaining.end(), i, pool.end());

  // Then randomly select solutions, removing them from the 'remaining' vector and adding them to the adult population.
  while (numTransferred < adults_.size())
  {
    auto r = utils::rand::randInt(0, static_cast<int>(remaining.size()));
    auto& candidate = remaining[r];
    if (!candidate->bad())
    {
      adults_[numTransferred++] = remaining[r];
    }
    remaining.erase(remaining.begin() + r);  // (Highly inefficient!!)
  }

  // Finally save any worthy solution that still remains on to the store.
  addWorthy(store_, remaining.begin(), remaining.end());
}


#endif  // QCGA2SURVIVAL
