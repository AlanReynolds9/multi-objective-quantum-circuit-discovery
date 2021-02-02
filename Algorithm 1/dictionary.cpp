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

#include <cassert>
#include "dictionary.h"

arma::uword qicLibQbitIndex(int ourQbitIndex, int numQbits)
{
  // QIClib numbers its qbits from 1. Moreover, qbit 1 seems to correspond with the MOST significant qbit, while I'd
  // prefer it to be the least. Hence our qbit indices 0, 1, 2, 3 become 4, 3, 2, 1.
  assert(ourQbitIndex >= 0);
  assert(ourQbitIndex < numQbits);
  return static_cast<arma::uword>(numQbits - ourQbitIndex);
}
