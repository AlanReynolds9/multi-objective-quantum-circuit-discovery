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

#include <cassert>
#include "rng.h"
#include "controls.h"

using utils::rand::randInt;
using utils::rand::randBool;


Controls::Controls(int numQbits) :
num_qbits(numQbits),
is_control(numQbits, false)
{
}


Controls::Controls(int numQbits, const std::vector<int> qbits) :
num_qbits(numQbits),
is_control(numQbits, false)
{
  for (auto qbit : qbits)
  {
    assert(0 <= qbit && qbit < numQbits);
    is_control[qbit] = true;
  }
}


bool Controls::operator==(const Controls& rhs) const
{
  // We're unlikely to ever be comparing Controls for different values of num_qbits, and besides, if num_qbits differs
  // then length of is_control should also differ. Hence the first part is not really necessary.
  return num_qbits == rhs.num_qbits && is_control == rhs.is_control;
}


bool Controls::operator!=(const Controls& rhs) const
{
  return !operator==(rhs);
}


bool Controls::operator<(const Controls &rhs) const
{
  // We're unlikely to ever be comparing Controls for different values of num_qbits, and besides, if num_qbits differs
  // then length of is_control should also differ. Hence the first part is not really necessary.
  if (num_qbits != rhs.num_qbits)
  {
    return num_qbits < rhs.num_qbits;
  }
  return is_control < rhs.is_control;
}


void Controls::swapBits(int i, int j)
{
  std::swap(is_control[i], is_control[j]);
}


void Controls::remove(int i)
{
  // If i is in the list of controls, remove it.
  is_control[i] = false;
}


void Controls::add(int i)
{
  // If i is not in the list of controls, add it.
  is_control[i] = true;
}


int Controls::numControls() const
{
  // We don't attempt to use an algorithm here (or a ranged for loop) because we are using vector<bool>, which has
  // issues.
  int count {0};
  for (auto i = 0; i < num_qbits; ++i)
  {
    if (is_control[i])
    {
      ++count;
    }
  }

  return count;
}
    


std::vector<int> Controls::qbits() const
{
  // Output the control qbits in a sorted vector.
  std::vector<int> qbitSet;
  qbitSet.reserve(numControls());
  for (auto i = 0; i < num_qbits; ++i)
  {
    if (is_control[i])
    {
      qbitSet.push_back(i);
    }
  }
  return qbitSet;
}


bool Controls::isControl(int qbit) const
{
  return is_control[qbit];
}


int Controls::first() const
{
  for (auto qbit = 0; qbit < num_qbits; ++qbit)
  {
    if (is_control[qbit])
    {
      return qbit;
    }
  }
  return num_qbits;
}


std::vector<int> Controls::asList() const
{
  // Translate the controls stored here into a vector of control indices.
  std::vector<int> controls;
  for (auto i = 0; i < num_qbits; ++i)
  {
    if (is_control[i])
    {
      controls.push_back(i);
    }
  }
  return controls;
}


void Controls::output(std::ostream& out) const
{
  std::vector<int> controlBits;
  controlBits.reserve(num_qbits);
  for (auto i = 0; i < num_qbits; ++i)
  {
    if (is_control[i])
    {
      controlBits.push_back(i);
    }
  }

  if (controlBits.size() > 1)
  {
    out << ", with control qbits " << controlBits.front();
    for (auto i = controlBits.begin() + 1; i < controlBits.end() - 1; ++i)
    {
      out << ", " << *i;
    }
    out << " and " << controlBits.back();
  }
  else if (controlBits.size() == 1)
  {
    out << ", with control qbit " << controlBits.front();
  }
}


std::ostream& operator<<(std::ostream& out, const Controls& controls)
{
  controls.output(out);
  return out;
}
