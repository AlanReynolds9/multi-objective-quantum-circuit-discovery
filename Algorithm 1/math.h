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

#ifndef MATH_H
#define MATH_H

#include <vector>

// Binomial coefficient calculation. For small n only!
size_t binCoeff(size_t n, size_t r);

// Function to assign a rank to a combination of natural numbers. Each k-combination is assigned a unique rank. This is
// a bijection (assuming we include 0). Note that 'combination' should be sorted in increasing order.
size_t combinationRank(const std::vector<int>& combination);
std::vector<int> combinationUnrank(int rank, int n, int k);

#endif // MATH_H
