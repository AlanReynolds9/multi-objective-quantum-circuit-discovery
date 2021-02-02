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
#include <stdexcept>
#include <algorithm>
#include <rng.h>

namespace utils
{
  namespace rand
  {
    Rng rng(static_cast<unsigned>(std::random_device()()));


    void seedRand(uint_fast32_t seed)
    {
      rng.seed(seed);
    }


    double rand01()
    {
      // Function to produce a random double in [0, 1).
      // (Unnecessary 'static', but might be slightly faster constructing the distribution just once.)
      static std::uniform_real_distribution<> dist {0, 1};
      return dist(rng);
    }


    int randInt(int rangeBegin, int rangeEnd)  // Parameters were lowerBound and upperBound
    {
      // Function to produce a random integer in [rangeBegin, rangeEnd).
      // Consider changing to [lowerBound, upperBound], i.e. allow upperBound? (Seems the standard style for random
      // number generators nowadays, despite conflicting with the usual C++ idiom for ranges.)
      assert(rangeEnd > rangeBegin);
      std::uniform_int_distribution<> dist {rangeBegin, rangeEnd - 1};
      return dist(rng);
    }


    std::pair<int, int> randIntPair(int rangeBegin, int rangeEnd)
    {
      // Function to select, at random, two different integers from the range.
      assert(rangeEnd > rangeBegin + 1);
      auto i = randInt(rangeBegin, rangeEnd);  // Select random element.
      auto j = randInt(rangeBegin, rangeEnd - 1);  // Select random remaining element.
      if (j >= i)
      {
        ++j;  // We have selected the jth remaining element, which is the j + 1th original element.
        return std::make_pair(i, j);
      }
      else
      {
        return std::make_pair(j, i);
      }
    }


    std::vector<int> randIntSet(int rangeBegin, int rangeEnd)
    {
      // Function to produce a random subset of the integers in the range [rangeBegin, rangeEnd).
      assert(rangeEnd > rangeBegin);
      std::vector<int> out;
      out.reserve(rangeEnd - rangeBegin);
      for (auto i = rangeBegin; i < rangeEnd; ++i)  // Would like to use a range-based for loop, but have problems...
      {                                             // ...with vector<bool>.
        if (randBool())
        {
          out.push_back(i);
        }
      }
      return out;
    }


    std::vector<int> randIntSet(int rangeBegin, int rangeEnd, int count)
    {
      // Function to produce a random sorted vector of 'count' integers, selected from [rangeBegin, rangeEnd).
      // EFFICIENCY: An alternative version would be preferable if the range is large, but 'count' is small.
      auto rangeSize = rangeEnd - rangeBegin;
      if (count > rangeSize)
      {
        // (Should we use assert or throw instead?)
        throw std::invalid_argument("Attempting to use randIntSet(int rangeBegin, int rangeEnd, int count) to get more"
                                    " (count) unique integers than is available in the range provided.");
      }

      // Create a permutation of the permissible values.
      std::vector<int> permutation(rangeSize);
      std::iota(permutation.begin(), permutation.end(), rangeBegin);

      // Shuffle
      std::shuffle(permutation.begin(), permutation.end(), rng);  // (Rather excessive if 'count' is small compared...
                                                                  // ...to 'rangeSize'.
      // Return the first 'count' elements
      permutation.erase(permutation.begin() + count, permutation.end());
      std::sort(permutation.begin(), permutation.end());
      return permutation;
    }


    int randGeometric(double mean)
    {
      // Produces a positive integer. Probability of n is given by p(1 - p)^(n - 1), where p = 1 / mean.
      // Note that std::geometric_distribution produces a non-negative integer, so we must add 1.
      std::geometric_distribution<> dist {1 / mean};
      return dist(rng) + 1;
    }


    double randDouble(double rangeBegin, double rangeEnd)
    {
      // Function to produce a random double in [rangeBegin, rangeEnd).
      std::uniform_real_distribution<> dist {rangeBegin, rangeEnd};
      return dist(rng);
    }


    double randNormal(double mean, double stdev)
    {
      // Function to produce a random normally distributed double.
      std::normal_distribution<> dist {mean, stdev};
      return dist(rng);
    }


    bool randBool()
    {
      // Function to produce a random bool
      static std::bernoulli_distribution dist;  // Default probability is 0.5.
      return dist(rng);
    }
  }
}
