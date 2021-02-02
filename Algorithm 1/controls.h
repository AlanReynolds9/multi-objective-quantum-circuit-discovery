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

// Policy classes that handle the creation of Controls on Gates. This means that a single Gate implementation can handle
// not only, e.g., an uncontrolled x-rotation gate type, but also singly controlled x-rotations and multiply controlled
// x-rotations as separate types. This allows a finer degree of control over the set of permitted gates.

#ifndef CONTROLS_H
#define CONTROLS_H

#include <iostream>
#include <armadillo>

enum class PermittedControls
{
  invalid,
  none,
  one,
  notMany,  // i.e. no more than one
  many,     // i.e. at least two
  notOne,   // i.e. none or at least two
  atLeastOne,
  any
};

class Controls
{
public:
  Controls(int numQbits);
  Controls(int numQbits, const std::vector<int> qbits);  // Added for use by simplification code.

  bool operator==(const Controls& rhs) const;
  bool operator!=(const Controls& rhs) const;
  bool operator<(const Controls& rhs) const;  // Only needed for calculations of population entropy.
                                              // (Rename to sortBefore()?)
  void swapBits(int i, int j);
  void remove(int i);
  void add(int i);

  int numControls() const;
  std::vector<int> qbits() const;  // Output controls as a sorted vector of ints. (Inefficient and best avoided.)
  bool isControl(int qbit) const;
  int first() const;  // First control qbit. If there are none, returns num_qbits.

  // Output the controls in a format suitable for the Transformation object.
  arma::uvec transControls() const;

  void output(std::ostream& out) const;

private:
  int num_qbits;
  std::vector<char> is_control;  // vector<bool> breaks swap(), and we wish to be able to swap bits straightforwardly.
};

std::ostream& operator<<(std::ostream& out, const Controls& controls);

#endif // CONTROLS_H
