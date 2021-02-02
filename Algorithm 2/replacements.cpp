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

#include "replacements.h"


Replacement::Replacement() :
left_pos(-2),  // -1 is used to indicate that a gate can be moved to the beginning and then eliminated.
right_pos(-2),
meet_pos(-2)
{
}


Replacement::Replacement(int leftPos, int rightPos, int meetPos) :
left_pos(leftPos),
right_pos(rightPos),
meet_pos(meetPos)
{
}


bool Replacement::valid() const
{
  return right_pos >= 0;
}


int Replacement::leftPos() const
{
  return left_pos;
}


int Replacement::rightPos() const
{
  return right_pos;
}


int Replacement::meetPos() const
{
  return meet_pos;
}
