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

#include <cassert>
#include <algorithm>
#include <iostream>
#include "nonDominatedStrimmer.h"

using namespace std;

class IndexLess
{
public:
  explicit IndexLess(const vector<double>& values) :
  values_(values)
  {
  }

  bool operator()(int lhs, int rhs) const
  {
    return values_[lhs] < values_[rhs];
  }

private:
  const vector<double>& values_;
};


NonDominatedStrimmer::NonDominatedStrimmer(int numNonDominated) :
num_items(numNonDominated),
distances_(numNonDominated, vector<double>(numNonDominated, -1)),
index_by_distance(numNonDominated),
available_(numNonDominated, true)
{
  // Ensure that diagonal elements are zero
  for (int item = 0; item < numNonDominated; ++item)
  {
    distances_[item][item] = 0;
  }

  // Initialize the indices so that they hold all the other objects, in whatever order!
  for (int item = 0; item < numNonDominated; ++item)
  {
    for (int index = 0; index < numNonDominated; ++index)
    {
      if (index != item)
      {
        index_by_distance[item].push_back(index);
      }
    }
  }
}

void NonDominatedStrimmer::addDistance(int i, int j, double dist)
{
  distances_[i][j] = distances_[j][i] = dist;
}

// For now, we completely sort the arrays of indices. We will later consider lazy sorting
void NonDominatedStrimmer::initializeIndices()
{
  for (int item = 0; item < num_items; ++item)
  {
    sort(index_by_distance[item].begin(), index_by_distance[item].end(), IndexLess(distances_[item]));
  }
}

int NonDominatedStrimmer::mostCrowded() const
{
  int best = -1;  // Added a value to avoid 'best might be uninitialized' below.
  int i;
  for (i = 0; i < num_items; ++i)
  {
    if (available_[i])
    {
      best = i;
      break;
    }
  }

  assert(best >= 0);  // Added, just in case some bug results in all solutions being unavailable.

  for (++i; i < num_items; ++i)
  {
    if (available_[i] && moreCrowded(i, best))  // If best is not initialized above, we get a warning here.
    {
      best = i;
    }
  }

  return best;
}

void NonDominatedStrimmer::remove(int toRemove)
{
  available_[toRemove] = false;
}

bool NonDominatedStrimmer::moreCrowded(int lhs, int rhs) const
{
  // Lexicographical compare does not work, as we need to use two reference distance vectors. Therefore we do this
  // manually
  // (Is there a better way?)
  DistanceIterator lhsIt = begin(lhs);
  DistanceIterator rhsIt = begin(rhs);
  while (lhsIt != end(lhs) && rhsIt != end(rhs))
  {
    if (distances_[lhs][*lhsIt] < distances_[rhs][*rhsIt])
    {
      return true;
    }
    if (distances_[lhs][*lhsIt] > distances_[rhs][*rhsIt])
    {
      return false;
    }

    ++lhsIt;
    ++rhsIt;
  }

  // The arrays being compared ought to be the same length.
  if (lhsIt == end(lhs) && rhsIt != end(rhs))
  {
    cerr << "Left hand array is shorter than the right hand array. Odd." << endl;
    return true;
  }
  if (lhsIt != end(lhs) && rhsIt == end(rhs))
  {
    cerr << "Left hand array is longer than the right hand array. Odd." << endl;
    return false;
  }

  return false;
}

NonDominatedStrimmer::DistanceIterator::DistanceIterator(const NonDominatedStrimmer& comp, int base, int i) :
comp_(&comp),
base_(base),
i_(i)
{
  // Ensure that i points to an available item
  while (i_ < comp_->num_items - 1 && !comp_->available_[comp_->index_by_distance[base_][i_]])
  {
    ++i_;
  }
}

NonDominatedStrimmer::DistanceIterator::DistanceIterator(const NonDominatedStrimmer::DistanceIterator& rhs) :
comp_(rhs.comp_),
base_(rhs.base_),
i_(rhs.i_)
{
}

bool NonDominatedStrimmer::DistanceIterator::operator==(const NonDominatedStrimmer::DistanceIterator& rhs) const
{
  return comp_ == rhs.comp_ && base_ == rhs.base_ && i_ == rhs.i_;
}

bool NonDominatedStrimmer::DistanceIterator::operator!=(const NonDominatedStrimmer::DistanceIterator& rhs) const
{
  return !(*this == rhs);
}

int NonDominatedStrimmer::DistanceIterator::operator*() const
{
  return comp_->index_by_distance[base_][i_];
}

NonDominatedStrimmer::DistanceIterator& NonDominatedStrimmer::DistanceIterator::operator++()
{
  ++i_;
  while (i_ < comp_->num_items - 1 && !comp_->available_[comp_->index_by_distance[base_][i_]])
  {
    ++i_;
  }

  return *this;
}

NonDominatedStrimmer::DistanceIterator NonDominatedStrimmer::DistanceIterator::operator++(int)
{
  NonDominatedStrimmer::DistanceIterator old(*this);
  ++*this;

  return old;
}

NonDominatedStrimmer::DistanceIterator NonDominatedStrimmer::begin(int base) const
{
  return DistanceIterator(*this, base, 0);
}

NonDominatedStrimmer::DistanceIterator NonDominatedStrimmer::end(int base) const
{
  return DistanceIterator(*this, base, num_items - 1);
}
