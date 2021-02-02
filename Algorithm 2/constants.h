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

#ifndef CONSTANTS_H
#define CONSTANTS_H

// Used to be just constants. Now contains some algorithm parameters that are 'constant once set'

#include <complex>

using cmplx = std::complex<double>;

namespace constants
{
  constexpr std::complex<double> i{0, 1};
  constexpr double pi = 3.141592653589793238462643383279502884;  // Excessive? :-)
  constexpr double tolerance = 0.0000000001; // Used in the dominance relation only. (Thus far.)
  constexpr int totalNumGateTypes = 14;
}


#endif // CONSTANTS_H
