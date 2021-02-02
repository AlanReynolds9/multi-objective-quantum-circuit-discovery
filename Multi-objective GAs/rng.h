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

#ifndef RNG_H
#define RNG_H

// Used to include <random> in rng.cpp, but now need to make the base random number generator available for use by
// shuffle.
#include <random>
#include <utility>
#include <vector>

namespace utils
{
  namespace rand
  {
    // (Should we put all this random number generation stuff into a class?)
    
    // The base random number generator. Rng and rng used to be hidden away in an unnamed namespace in rng.cpp, but
    // shuffle needs rng.
    using Rng = std::mt19937;
    extern Rng rng;
    
    void seedRand(uint_fast32_t seed);
    double rand01();
    int randInt(int rangeBegin, int rangeEnd);  // Was lowerBound and upperBound.
    std::pair<int, int> randIntPair(int rangeBegin, int rangeEnd);  // Each different. First < second.
    std::vector<int> randIntSet(int rangeBegin, int rangeEnd);  // Selects unique sorted elements.
    std::vector<int> randIntSet(int rangeBegin, int rangeEnd, int count);
    int randGeometric(double mean);  // Values start at one.
    double randDouble(double rangeBegin, double rangeEnd);
    double randNormal(double mean, double stdev);
    bool randBool();
  }
}

#endif
