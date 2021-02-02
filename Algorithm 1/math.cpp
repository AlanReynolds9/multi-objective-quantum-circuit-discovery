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

#include <algorithm>
#include "math.h"

using std::vector;

size_t binCoeff(size_t n, size_t r)
{
  // Calculates the binomial coefficient nCr. Only useful for small n!
  // Only used in calculations associated with the hash function for circuits, so use size_t rather than int.

  // If the number of items we are attempting to select exceeds the set size, then there are zero ways of doing so.
  if (r > n)
  {
    return 0;
  }

  // If r is large (i.e. more than half of n), it will be quicker to calculate nC(n - r)
  if (n < 2 * r)
  {
    return binCoeff(n, n - r);
  }

  // Calculate using a series of multiplications and divisions, ordered so as to reduce the chance of overflow.
  size_t result{1};
  for (size_t i = 1; i <= r; ++i)
  {
    result *= n - r + i;
    result /= i;
  }

  return result;
}


size_t combinationRank(const vector<int>& combination)
{
  // Assign a rank to a combination of natural numbers. Each k-combination is assigned a unique rank. This is a
  // bijection from k-combinations to the natural numbers, including zero. (Note that combinations of different sizes
  // may be mapped to the same number.) The 'combination' should be pre-sorted in increasing order.
  size_t rank{0};
  for (auto i = 0; i < combination.size(); ++i)
  {
    rank += binCoeff(combination[i], i + 1);
  }

  return rank;
}


vector<int> combinationUnrank(int rank, int n, int k)
{
  // Find the k-combination that, on application of combinationRank(), produces the rank provided. The number n is
  // provided to help, and gives the number of qbits the combination was selected from. (We could implement the function
  // without n, but it would require more thought.)
  vector<int> combination;
  combination.reserve(k);

  for (int candidate = n - 1; candidate >= 0; --candidate)
  {
    size_t bc = binCoeff(candidate, k);
    if (bc <= rank)
    {
      combination.push_back(candidate);
      rank -= bc;
      --k;
    }
  }

  // We found the qbits starting with that with the highest index. We therefore need to reverse the vector to get the
  // bits in increasing order.
  std::reverse(combination.begin(), combination.end());

  return combination;
}
