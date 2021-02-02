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

// At present, this just contains a function for converting our qbit indices into those that QIClib uses.
// Alas, this doesn't quite hide all of the qbit indexing peculiarities from the rest of the code - there is also
// swapPermutation() (in state.cpp) that uses the fact that QIClib qbit indices start from 1, but must be placed in a
// uvec with indexing starting from 0 when using qic::sysperm().

#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <armadillo>

// We could change the return type to auto (I suspect), but only by moving the function definition here too. (Code that
// uses this function needs to know the return type, which can only be deduced from the definition.
arma::uword qicLibQbitIndex(int ourQbitIndex, int numQbits);

#endif // DICTIONARY_H
