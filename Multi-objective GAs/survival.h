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

#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <algorithm>
#include "utils.h"
#include "population.h"
#include "distances.h"
#include "stores.h"
#include "nonDominatedStrimmer.h"

// May wish to create a Survival base class

template <class Solution, class Store>
class ChildSurvival
{
  // Children simply replace the adults.
public:
  ChildSurvival(Population<Solution>& adults, Population<Solution>& children, Store& store);

  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;
};

template <class Solution, class Store, class Evaluator>
class BestSurvival
{
  // The adult and child populations are combined and evaluated, after which only the best survive.
public:
  BestSurvival(Population<Solution>& adults, Population<Solution>& children, Store& store, Evaluator& evaluator);

  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
};

template <class Solution, class Store, class Evaluator>
class SortedBestSurvival
{
  // As for BestSurvival, but uses an Evaluator that sorts the population in order of quality. It is therefore
  // unnecessary to resort. (Of course, the Evaluator selected must be one that sorts the population.)
public:
  SortedBestSurvival(Population<Solution>& adults, Population<Solution>& children, Store& store, Evaluator& evaluator);

  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
};

template <class Solution, class Store, class Evaluator>
class SortedBestUniqueSurvival
{
  // As for SortedBestSurvival, but duplicates are not allowed to survive unless absolutely necessary to get the correct
  // population size.
public:
  SortedBestUniqueSurvival(Population<Solution>& adults, Population<Solution>& children, Store& store,
                           Evaluator& evaluator);
  void compete();

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
};

// A horrible mess. Must sort this lot out!
template <class Solution, class Store, class Evaluator>
class SPEA2Survival
{
  // If the number of non-dominated solutions does not exceed the adult population size, this works in the same way as
  // BestSurvival. However, when there are too many non-dominated solutions, different crowding is used to determine
  // which survive.
public:
  SPEA2Survival(Population<Solution>& adults, Population<Solution>& children, Store& store, Evaluator& evaluator,
                const Objectives& objectives);

  void compete();

private:
  void truncate_non_dominated(Population<Solution>& mergedPop); // Remove enough non-dominated solutions so that...
                                                   // ...those that remain will fit in adults_. Then add them to adults_

private:
  Population<Solution>& adults_;
  Population<Solution>& children_;
  Store& store_;

  Evaluator& evaluator_;
  const Objectives& objectives_;
};

template <class Solution, class Store>
class JMOGLSSurvival
{
  // A rather odd survival class for the MOGLS algorithm developed by Jaszkiewicz. Solutions are assumed to be already
  // evaluated with respect to the correct weighted objective. Only one child is available for addition and it is added
  // if it is better than the worst solution in the temporary population. If the adult population is full, then the
  // oldest solution is removed.
public:
  JMOGLSSurvival(Population<Solution>& adults, Population<Solution>& tempPopulation, Store& store);

  void compete(const std::shared_ptr<const Solution>& child);  // Made this a parameter, because it wasn't being...
                       // ...initialized correctly by the constructor - at least not when it was a boost::intrusive_ptr.
private:
  Population<Solution>& adults_;
  Population<Solution>& temp_population;
  Store& store_;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store>
ChildSurvival<Solution, Store>::ChildSurvival(Population<Solution>& adults, Population<Solution>& children,
                                              Store& store):
adults_(adults),
children_(children),
store_(store)
{
}

template <class Solution, class Store>
void ChildSurvival<Solution, Store>::compete()
{
  // Ensure that any worthy solutions in the adult population are saved to the store, before elimination
  addWorthy(store_, adults_.begin(), adults_.end()); 

  // Replace adults with children
  adults_ = children_;
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store, class Evaluator>
BestSurvival<Solution, Store, Evaluator>::BestSurvival(Population<Solution>& adults, Population<Solution>& children,
                                                       Store& store, Evaluator& evaluator):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator)
{
}

template <class Solution, class Store, class Evaluator>
void BestSurvival<Solution, Store, Evaluator>::compete()
{
  // Create a merged pool of solutions and evaluate
  Population<Solution> mergedPop(adults_, children_);
  evaluator_.evaluate(mergedPop);

  // Sort solutions
  // (Erm, the STL has a algorithm for doing a partial sort!!!)
  std::sort(mergedPop.begin(), mergedPop.end(), fitterThanP<std::shared_ptr<Solution> >);

  // Ensure that any worthy solutions that are about to be eliminated are saved to store
  addWorthy(store_, mergedPop.begin() + adults_.size(), mergedPop.end());

  // Copy the best solutions into the adult pool, which eliminates the worst.
  copy(mergedPop.begin(), mergedPop.begin() + adults_.size(), adults_.begin());
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store, class Evaluator>
SortedBestSurvival<Solution, Store, Evaluator>::SortedBestSurvival(Population<Solution>& adults,
                                                                   Population<Solution>& children, Store& store,
                                                                   Evaluator& evaluator):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator)
{
}

template <class Solution, class Store, class Evaluator>
void SortedBestSurvival<Solution, Store, Evaluator>::compete()
{
  // Create a merged population, let the evaluator order it and take the first half
  Population<Solution> mergedPop(adults_, children_);
  evaluator_.evaluate(mergedPop);

  // Ensure that any worthy solutions that are about to be eliminated are saved to store
  addWorthy(store_, mergedPop.begin() + adults_.size(), mergedPop.end());

  // Copy the best solutions into the adult pool, which eliminates the worst.
  copy(mergedPop.begin(), mergedPop.begin() + adults_.size(), adults_.begin());
}

//------------------------------------------------------

template <class Solution, class Store, class Evaluator>
SortedBestUniqueSurvival<Solution, Store, Evaluator>::SortedBestUniqueSurvival(Population<Solution>& adults,
                                                                               Population<Solution>& children,
                                                                               Store& store, Evaluator& evaluator):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator)
{
}

template <class Solution, class Store, class Evaluator>
void SortedBestUniqueSurvival<Solution, Store, Evaluator>::compete()
{
  // Create a merged population, let the evaluator sort it, then take solutions from the front, ignoring duplicates.
  // (Currently inefficient - O(n^2). However, it is difficult to see how to create an efficient version that is
  // sufficiently general. For example, we might provide a Solution class where equality is tested up to some tolerance
  // (in the floating point parts). I cannot see how to sort such solutions such that those considered 'equal' are
  // guaranteed to occur next to each other in the new sequence.)
  // (We also assume, at present, that there are sufficient unique solutions to fill the adult population. If the
  // previous adult population contained only unique solutions, then there should be enough. However, it is just about
  // possible that the initial population contains duplicates, so this is NOT SAFE.)

  // Create a merged population and let the Evaluator evaluate the solutions and sort them.
  Population<Solution> mergedPop(adults_, children_);
  evaluator_.evaluate(mergedPop);

  int numTransferred{0};  // Number of solutions transferred to the new adult population.
  std::vector<std::shared_ptr<Solution> > discards;
  discards.reserve(children_.size());    // Over-reserving (but safe).
  auto i = mergedPop.begin();
  for (; numTransferred < adults_.size(); ++i)  // Assumes that the new adult population is filled before we reach...
  {                                             // ...the end of the candidates.
    auto& candidate = *i;
    bool alreadyPresent{false};
    for (auto j = 0; j < numTransferred && !alreadyPresent; ++j)
    {
      if (*candidate == *adults_[j])  // If the solutions (not the pointers!) are the same
      {
        alreadyPresent = true;
      }
    }

    if (alreadyPresent)
    {
      discards.push_back(candidate);
    }
    else
    {
      adults_[numTransferred++] = candidate;
    }
  }

  // Save any of those solutions that are worthy, but marked as being duplicates, on to the store. (This might seem
  // silly, since copies still exist in the population. However, it is possible that the Solution class's equality
  // operator allows for some tolerance, i.e. merely indicates that the solutions are similar. We would like to ensure
  // that the Store contains the best of these.)
  addWorthy(store_, discards.begin(), discards.end());

  // Also ensure that any worthy solutions remaining in mergedPop but about to be eliminated are saved to store.
  addWorthy(store_, i, mergedPop.end());
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store, class Evaluator>
SPEA2Survival<Solution, Store, Evaluator>::SPEA2Survival(Population<Solution>& adults, Population<Solution>& children,
                                                         Store& store, Evaluator& evaluator,
                                                         const Objectives& objectives):
adults_(adults),
children_(children),
store_(store),
evaluator_(evaluator),
objectives_(objectives)    // I don't like passing the objectives data in here
{
}

template <class Solution, class Store, class Evaluator>
void SPEA2Survival<Solution, Store, Evaluator>::compete()
{
  // Create a merged pool of solutions and evaluate
  Population<Solution> mergedPop(adults_, children_);
  evaluator_.evaluate(mergedPop);

  // Sort solutions
  // (Erm, the STL has a algorithm for doing a partial sort!!!)
  sort(mergedPop.begin(), mergedPop.end(), fitterThanP<std::shared_ptr<Solution> >);  // Uses the generic fitterThan...
                                                                                      // ...defined on smart pointers.
  // If the first half is assumed to be the survivors, check that the first solution in the second half is dominated. If
  // not, then a more complicated survival method is required
  if (mergedPop.size() <= adults_.capacity() || mergedPop[adults_.capacity()]->fitness() > 0)  // First bit for the...
  {                                                                                  // ...first generation (No adults.)
    // Either the merged population is smaller than the capacity for adults, e.g. in the first generation when
    // childCapacity <= adultCapacity, or there are insufficient non-dominated solutions to fill the adult population.
    // In this case, just take the best solutions.

    // copy(mergedPop.begin(), mergedPop.begin() + adults_.capacity(), adults_.begin());
    // I would like to use the above call to copy, but cannot since at the first iteration, adults has size of zero.
    // Also cannot use resize unless we add it to the population class
    adults_.clear();
    int i = 0;
    while (i < mergedPop.size() && !adults_.full())  // Added the first bit in case total number of solutions is less...
    {                                                // ...than adult capacity.
      adults_.add(mergedPop[i++]);
    }
  }
  else
  {
    // Too many non-dominated solutions. Perform the routine that iteratively removes the most crowded solution,
    // according to lexicographical distance comparisons.
    truncate_non_dominated(mergedPop);
  }
}

template <class Solution, class Store, class Evaluator>
void SPEA2Survival<Solution, Store, Evaluator>::truncate_non_dominated(/*const */Population<Solution>& mergedPop)
{
  // There are too many non-dominated solutions. Iteratively remove the most crowded solutions until there are few
  // enough to fit into the new adult population. The most crowded solution is found by using lexicographical
  // comparisons of the sorted distances to the other solutions.

  // Find the position of the first dominated solution in the merged population
  auto nonDominatedEnd = mergedPop.begin();  // (Was briefly cbegin(), but cbegin doesn't exist yet for populations!)
  for (; nonDominatedEnd != mergedPop.end(); ++nonDominatedEnd)  // (Was briefly cend().)
  {
    if ((*nonDominatedEnd)->fitness() > 0)
    {
      break;
    }
  }
  int numNonDominated = static_cast<int>(nonDominatedEnd - mergedPop.begin());
  int numSurvivors = numNonDominated;

  // Create the object for handling the distance matrix and strimming out the most crowded solutions.
  NonDominatedStrimmer strimmer(numNonDominated);
  std::vector<bool> survives(numNonDominated, true);
  for (int i = 0; i < numNonDominated; ++i)
  {
    const Solution& sol1 = *mergedPop[i];
    for (int j = i + 1; j < numNonDominated; ++j)
    {
      const Solution& sol2 = *mergedPop[j];

      strimmer.addDistance(i, j, ObjectiveNorm<Solution, 2>(objectives_)(sol1, sol2));  // Adds both i->j and j->i
    }
  }

  // Need to initialize, otherwise this will not work!
  strimmer.initializeIndices();

  // Find the most crowded item, remove and repeat until the number of remaining solutions equals the adult population
  // size.
  while (numSurvivors > adults_.capacity())
  {
    // Find the most crowded solution
    int mostCrowded = strimmer.mostCrowded();

    // Remove it from the distance matrix handler and mark it as a non-survivor
    strimmer.remove(mostCrowded);
    survives[mostCrowded] = false;
    --numSurvivors;

    // Ensure that any worthy solutions are added to the store
    if (mergedPop[mostCrowded]->worthy())
    {
      store_.add(mergedPop[mostCrowded]);  // (This one adds a shared pointer to a constant solution, and is ok.)
    }
  }

  // Copy the remaining solutions to the adult population
  adults_.clear();
  for (int i = 0; i < numNonDominated; ++i)
  {
    if (survives[i])
    {
      adults_.add(mergedPop[i]);   // (This line is the reason the merged population is not const.)
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Store>
JMOGLSSurvival<Solution, Store>::JMOGLSSurvival(Population<Solution>& adults, Population<Solution>& tempPopulation,
                                                Store& store) :
adults_(adults),
temp_population(tempPopulation),
store_(store)
{
}

template <class Solution, class Store>
void JMOGLSSurvival<Solution, Store>::compete(const std::shared_ptr<const Solution>& child)
{
  // Is the child better than at least one of the solutions, and different from all the solutions, in the temporary
  // population?
  bool better = false;
  bool same = false;
  for (int i = 0; i < temp_population.size() && !same; ++i)
  {
    if (child->fitness() > temp_population[i]->fitness())
    {
      better = true;
    }
    if (*child == *temp_population[i])
    {
      same = true;
    }
  }

  if (!same && better)
  {
    // Child should be added to the population.
    // Is the population currently full - if so, find the oldest solution and replace with the new one, otherwise simply
    // add the new one.
    if (adults_.full())
    {
      // Find the oldest
      auto oldest = min_element(adults_.begin(), adults_.end(), olderThanP<std::shared_ptr<Solution> >);
      *oldest = new Solution(*child);  // (Would like to just do *oldest = child_, but problems with const occur.)
    }
    else
    {
      adults_.add(new Solution(*child));  // (Would like to just do *oldest = child_, but problems with const occur.)
    }
  }
}

#endif
