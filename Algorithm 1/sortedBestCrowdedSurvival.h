// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
//
// This file is part of MOQCD1.
//
// MOQCD1 is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
//
// MOQCD1 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with MOQCD1.  If not, see
// <http://www.gnu.org/licenses/>.

//----------------------------------------------------------------------------------------------------------------------

#ifndef SORTEDBESTCROWDEDSURVIVAL_H
#define SORTEDBESTCROWDEDSURVIVAL_H

#include <algorithm>
#include "utils.h"
#include "population.h"
#include "distances.h"
#include "stores.h"
#include "nonDominatedStrimmer.h"


template <class Solution, class Store, class Evaluator>
class SortedBestCrowdedSurvival
{
  // As for BestSurvival, but uses an Evaluator that sorts the population in order of quality. It is therefore
  // unnecessary to resort. (Of course, the Evaluator selected must be one that sorts the population.)
public:
  using Problem = typename Solution::Problem;

  SortedBestCrowdedSurvival(Population<Solution>& adults, Population<Solution>& children, Store& store,
                            Evaluator& evaluator, const Objectives& objectives, double spacing);

  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
  ObjectiveNorm<Solution, 1> distance_;
  double preferred_spacing;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store, class Evaluator>
SortedBestCrowdedSurvival<Solution, Store, Evaluator>::SortedBestCrowdedSurvival(Population<Solution>& adults,
                                                                                 Population<Solution>& children,
                                                                                 Store& store, Evaluator& evaluator,
                                                                                 const Objectives& objectives,
                                                                                 double spacing):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator),
distance_(objectives),
preferred_spacing(spacing)
{
}

template <class Solution, class Store, class Evaluator>
void SortedBestCrowdedSurvival<Solution, Store, Evaluator>::compete()
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
  Population<Solution> discards(adults_.size() + children_.size());  // Over-reserving (but safe).

  // Main loop
  auto i = pool.begin();
  while (numTransferred < adults_.size())
  {
    for (i = pool.begin(); numTransferred < adults_.size() && i != pool.end(); ++i)
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
      else
      {
        adults_[numTransferred++] = candidate;
      }
    }

    // Unable to find sufficiently many solutions using the current spacing restriction. Temporarily reduce crowding
    // distance by a factor of 2.
    if (numTransferred < adults_.size())
    {
      currentSpacing = currentSpacing / 2 - 0.000001;   // The -0.000001 ensures that, eventually, this becomes...
                                                        // non-positive.
      pool = discards;  // (Adding a move constructor to the Population class might simplify this and make it more...
      discards.clear(); // ...efficient.)
    }
  }

  // Save any of those solutions that are worthy, but marked as being duplicates, on to the store. (This might seem
  // silly, since copies still exist in the population. However, it is possible that the Solution class's equality
  // operator allows for some tolerance, i.e. merely indicates that the solutions are similar. We would like to ensure
  // that the Store contains the best of these.)
  addWorthy(store_, discards.begin(), discards.end());

  // Also ensure that any worthy solutions remaining in mergedPop but about to be eliminated are saved to store.
  addWorthy(store_, i, pool.end());
}



#endif
