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

#ifndef MATH_H
#define MATH_H

#include <vector>
#include "constants.h"  // NEW: To get cmplx.

// Binomial coefficient calculation. For small n only!
size_t binCoeff(size_t n, size_t r);  // Was used in hashing - hence 'size_t'. Feel free to change to int or long.

// Function to assign a rank to a combination of natural numbers. Each k-combination is assigned a unique rank. This is
// a bijection (assuming we include 0). Note that 'combination' should be sorted in increasing order.
size_t combinationRank(const std::vector<int>& combination);  // Currently unused. Feel free to change size_t to int...
std::vector<int> combinationUnrank(int rank, int n, int k);   // ...or long.

// Matrix sparsifier
std::vector<std::vector<std::pair<int, cmplx> > > sparse(const std::vector<std::vector<cmplx> >& matrix);

#endif // MATH_H
