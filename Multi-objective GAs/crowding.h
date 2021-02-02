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

#ifndef CROWDING_H
#define CROWDING_H

#include <cassert>
#include <exception>
#include <limits>
#include <functional>
#include <vector>
#include <algorithm>
#include <memory>
#include <rng.h>
#include "objectives.h"
#include "distances.h"

// (The way the Solution and the Crowding classes interact seems, to me, to be less than ideal, since it requires the
// solution to provide both an objective value function and functions for setting and getting the crowding distances.
// One possibility is to use the private inheritance trick to create an additional interface, but this would require me
// to upcast all the items pointed to by the Iterators. How do I know what to upcast to?)

const double MAXDOUBLE = std::numeric_limits<double>::max();

//----------------------------------------------------------------------------------------------------------------------

// Crowding classes

template <class Solution>
class NSGA2Crowding
{
  // Class to perform crowding of solutions in a front/population.
  // (Can we not make the crowding class a friend of the solution class (e.g. NSGA2Solution), allowing us to hide the
  // objective and crowding distances within the solution.)
  using Iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;

public:
  NSGA2Crowding(const Objectives& objectives);
  void calculateCrowding(Iterator first, Iterator last) const;

private:
  const Objectives* objectives_;
};

template <class Solution>
class NormedNSGA2Crowding
{
  // Version that scales objectives by the difference between the greatest and smallest values in the front.
  using Iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;

public:
  NormedNSGA2Crowding(const Objectives& objectives);
  void calculateCrowding(Iterator first, Iterator last) const;

private:
  const Objectives* objectives_;
};

template <class Solution, class Distance>
class NearestObjCrowding
{
  // Solutions are randomly ordered. Crowding distance for the first is the distance of the closest other solution, by
  // objective (deltaX + deltaY). This solution is removed and then the next considered. (Avoids 'mutual crowding' and
  // the elimination of sections of the front
  using Iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;

public:
  NearestObjCrowding(const Objectives& objectives);
  void calculateCrowding(Iterator first, Iterator last) const;

private:
  const Objectives* objectives_;
};

template <class Solution, class Distance>
class MOGACrowding
{
  // Class to perform crowding of solutions in a rank/population. Rather that calculating a crowding value, this class
  // adjusts the fitness of solutions
  using Iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;

public:
  MOGACrowding(const Objectives& objectives, double nicheSize);
  void calculateCrowding(Iterator first, Iterator last) const;

private:
  const Objectives* objectives_;
  double niche_size;
};

template <class Solution, class Distance>
class SPEA2Crowding
{
  // Class to perform crowding of solutions using the distance to the kth nearest neighbour, as in SPEA2. Different
  // crowding techniques are used in the survival routine if the number of non-dominated solutions exceeds the
  // population size.
  using Iterator = typename std::vector<std::shared_ptr<Solution> >::iterator;

public:
  SPEA2Crowding(const Objectives& objectives, int k);
  void calculateCrowding(Iterator first, Iterator last) const;

private:
  const Objectives* objectives_;
  const int k_;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
NSGA2Crowding<Solution>::NSGA2Crowding(const Objectives& objectives) :
objectives_(&objectives)
{
  // (Would prefer some way of using a compile time assertion. Alternatively, some way to automatically use an
  // alternative.)
  for (auto obj = 0; obj < objectives_->numObj(); ++obj)
  {
    if (objectives_->upperBound(obj) - objectives_->lowerBound(obj) == 0)
    {
      throw std::domain_error("Attempting to use standard NSGA II crowding when the objective ranges are unknown."
                              " Suggest using normalized NSGA II crowding.");
    }
  }
}

template <class Solution>
void NSGA2Crowding<Solution>::calculateCrowding(Iterator first, Iterator last) const
{
  // Function to calculate the crowding of solutions in the range provided by the iterators.
  // (I suspect that it might be better to assume that the iterators passed are proper STL iterators, inheriting from
  // some random access iterator base class.)
  
  // First, reset the crowding distances of the solutions 
  for (auto i = first; i != last; ++i)
  {
    (*i)->setCrowdingDistance(0);
  }

  // For each objective
  for (auto objNum = 0; objNum < objectives_->numObj(); ++objNum)
  {
    // Create a comparator and sort the solutions by objective i
    CompareByObjective<Solution, std::less> lessByObjective(objNum);
    sort(first, last, lessByObjective);

    // Solutions with extreme objective values are given the maximum crowding distance
    (*first)->setCrowdingDistance(MAXDOUBLE);
    (*(last - 1))->setCrowdingDistance(MAXDOUBLE);
    
    // Each solution in the middle of the front...
    if (last - first > 2)
    {
      for (Iterator i = first + 1; i != last - 1; ++i)
      {
        Solution& prevSol = **(i - 1);
        Solution& thisSol = **i;
        Solution& nextSol = **(i + 1);
        // ...has a crowding distance obtained by summing the difference in objective values of those solutions...
        // ...surrounding it in the front
        if (thisSol.crowdingDistance() != MAXDOUBLE)
        {
          // (Note that the difference in objective values is normalized using the objective range information from the
          // ...Solution class.)
          double objDiff = static_cast<double>(nextSol.objective(objNum) - prevSol.objective(objNum));
          thisSol.setCrowdingDistance(thisSol.crowdingDistance() +
                                      objDiff / (objectives_->upperBound(objNum) - objectives_->lowerBound(objNum)));
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution>
NormedNSGA2Crowding<Solution>::NormedNSGA2Crowding(const Objectives& objectives) :
objectives_(&objectives)
{
}

template <class Solution>
void NormedNSGA2Crowding<Solution>::calculateCrowding(Iterator first, Iterator last) const
{
  // Function to calculate the crowding of solutions in the range provided by the iterators.
  // (I suspect that it might be better to assume that the iterators passed are proper STL iterators, inheriting from
  // some random access iterator base class.)
  
  // First, reset the crowding distances of the solutions 
  for (auto i = first; i != last; ++i)
  {
    (*i)->setCrowdingDistance(0);
  }

  // Fer each objective
  for (auto objNum = 0; objNum < objectives_->numObj(); ++objNum)
  {
    // Create a comparator and sort the solutions by objective i
    CompareByObjective<Solution, std::less> lessByObjective(objNum);
    sort(first, last, lessByObjective);

    // Objectives will be scaled according to the distance between the smallest and largest objective values
    auto minObj = static_cast<double>((*first)->objective(objNum));
    auto maxObj = static_cast<double>((*(last - 1))->objective(objNum));

    // If the maximum and the minimum are the same, set all crowding distances to the max
    if (maxObj == minObj)
    {
      for (auto i = first; i != last; ++i)
      {
        Solution& thisSol = **i;
        thisSol.setCrowdingDistance(MAXDOUBLE);
      }
    }
    else
    {
      // Solutions with extreme objective values are given the maximum crowding distance
      (*first)->setCrowdingDistance(MAXDOUBLE);
      (*(last - 1))->setCrowdingDistance(MAXDOUBLE);

      // Each solution in the middle of the front...
      if (last - first > 2)
      {
        for (Iterator i = first + 1; i != last - 1; ++i)
        {
          Solution& prevSol = **(i - 1);
          Solution& thisSol = **i;
          Solution& nextSol = **(i + 1);
          // ...has a crowding distance obtained by summing the difference in objective values of those solutions...
          // ...surrounding it in the front
          if (thisSol.crowdingDistance() != MAXDOUBLE)
          {
            // (Note that the difference in objective values is normalized using the range of objective values in the...
            // ...front.)
            double objDiff = static_cast<double>(nextSol.objective(objNum) - prevSol.objective(objNum));
            thisSol.setCrowdingDistance(thisSol.crowdingDistance() + objDiff / (maxObj - minObj));
          }
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Distance>
NearestObjCrowding<Solution, Distance>::NearestObjCrowding(const Objectives& objectives) :
objectives_(&objectives)
{
}

template <class Solution, class Distance>
void NearestObjCrowding<Solution, Distance>::calculateCrowding(Iterator first, Iterator last) const
{
  // Create a random ordering for the solutions
  int numSolutions = static_cast<int>(last - first);
  std::vector<int> randSequence(numSolutions);
  for (int i = 0; i < numSolutions; ++i)
  {
    randSequence[i] = i;
  }
  shuffle(randSequence.begin(), randSequence.end(), utils::rand::rng);

  // Calculate the crowding distance for each solution in turn. Once a solution's crowding distance is calculated, it is
  // not considered in the crowding of other solutions. This means that if two solutions form a tight but isolated pair,
  // only one is crowded out.
  for (int i = 0; i < numSolutions; ++i)
  {
    Solution& crowdedSolution = **(first + randSequence[i]);

    // Initially set the distance to be huge
    double minDistance = MAXDOUBLE;

    // Find the distance to the closest remaining solution
    for (int j = i + 1; j < numSolutions; ++j)
    {
      Solution& crowdingSolution = **(first + randSequence[j]);
      minDistance = std::min(minDistance, Distance(*objectives_)(crowdedSolution, crowdingSolution));
    }
    crowdedSolution.setCrowdingDistance(minDistance);
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Distance>
MOGACrowding<Solution, Distance>::MOGACrowding(const Objectives& objectives, double nicheSize) :
objectives_(&objectives),
niche_size(nicheSize)
{
}

template <class Solution, class Distance>
void MOGACrowding<Solution, Distance>::calculateCrowding(Iterator first, Iterator last) const
{
  int numSolutions = static_cast<int>(last - first);

  // Calculate the niche count for each solution, and the sum of the fitness of each solution
  // (In the standard MOGA algorithm, crowding is performed on a rank by rank basis, with all solutions in a rank having
  // the same fitness, so we could just multiply the fitness of the first by the number of solutions.)
  std::vector<double> nicheCount(numSolutions, 0.0); 
  double originalFitnessSum = 0.0;
  for (auto i = first; i != last; ++i)
  {
    const Solution& sol1 = **i;   // iterator to a collection of shared_ptrs
    originalFitnessSum += sol1.fitness();
    for (auto j = first; j != last; ++j)
    {
      const Solution& sol2 = **j;

      double distance = Distance(*objectives_)(sol1, sol2);
      if (distance < niche_size)
      {
        nicheCount[i - first] += 1 - distance / niche_size;
      }
    }
  }
  
  // Adjust the fitness of each solution be dividing by the niche count. Also note the new sum of fitness
  double adjustedFitnessSum = 0.0;
  for (auto i = first; i != last; ++i)
  {
    Solution& sol = **i;

    // Adjust...
    double fitness = sol.fitness();
    fitness /= nicheCount[i - first];
    sol.setFitness(fitness);

    // Sum...
    adjustedFitnessSum += fitness;
  }

  // Normalise so that the sum of the solution fitness within the rank is unchanged
  double scaleFactor = originalFitnessSum / adjustedFitnessSum;
  for (auto i = first; i != last; ++i)
  {
    Solution& sol = **i;

    double fitness = sol.fitness();
    fitness *= scaleFactor;
    sol.setFitness(fitness);
  }
}

//----------------------------------------------------------------------------------------------------------------------

template <class Solution, class Distance>
SPEA2Crowding<Solution, Distance>::SPEA2Crowding(const Objectives& objectives, int k) :
objectives_(&objectives),
k_(k)
{
  assert(k > 0);
}

template <class Solution, class Distance>
void SPEA2Crowding<Solution, Distance>::calculateCrowding(Iterator first, Iterator last) const
{
  // Function to get the distance of the kth nearest solution to each solution in the range
  int numSolutions = static_cast<int>(last - first);
  assert(numSolutions > k_); // e.g need at least 3 solutions to be able to find the second closest

  // First create a distance matrix 
  std::vector<std::vector<double> > dist(numSolutions, std::vector<double>(numSolutions, 0.0));

  // Fill with the distances
  for (int i = 0; i < numSolutions; ++i)
  {
    Solution& sol1 = **(first + i);
    for (int j = i + 1; j < numSolutions; ++j)
    {
      Solution& sol2 = **(first + j);
      dist[i][j] = dist[j][i] = Distance(*objectives_)(sol1, sol2);
    }
  }

  // Find the distance to the kth nearest neighbour 
  for (int i = 0; i < numSolutions; ++i)
  {
    Solution& crowdedSolution = **(first + i);
    nth_element(dist[i].begin(), dist[i].begin() + k_, dist[i].end());

    // Rather than create an adjusted fitness, store the crowding distance, which is the distance to the kth neighbour.
    crowdedSolution.setCrowdingDistance(dist[i][k_]);
  }
}

#endif
