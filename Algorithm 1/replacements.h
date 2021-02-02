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

#ifndef REPLACEMENTS_H
#define REPLACEMENTS_H

//#include "circuit.h"

class Replacement
{
public:
  Replacement();
  Replacement(int leftPos, int rightPos, int meetPos);

  bool valid() const;

  int leftPos() const;
  int rightPos() const;
  int meetPos() const;

private:
  int left_pos;
  int right_pos;
  int meet_pos;
};



#endif // REPLACEMENTS_H
