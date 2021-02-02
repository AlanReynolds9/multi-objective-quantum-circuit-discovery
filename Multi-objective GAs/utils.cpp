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

#include "utils.h"

Available::Available(int size)
{
  for (int i = 0; i < size; ++i)
  {
    available_.push_back(i);
  }
}

Available::iterator Available::begin()
{
  return available_.begin();
}

Available::iterator Available::end()
{
  return available_.end();
}

Available::iterator Available::erase(iterator i)
{
  return available_.erase(i);
}

bool Available::empty() const
{
  return available_.empty();
}
