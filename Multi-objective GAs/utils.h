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

#ifndef UTILS_H
#define UTILS_H

#include <list>
#include <algorithm>
#include <functional>
#include <numeric>

class Available
{
  // Class to indicate which items in an array, vector or population are still available for use.
  // Could definitely be done more generically, I should think. 
  // (Rather than allowing the user to determine whether an element is available, we could allow access only to those
  // available elements. Consider using a linked list of references into the container. Marking an element as
  // unavailable could be awkward though.)
public:
  explicit Available(int size);

  using iterator = std::list<int>::iterator;
  iterator begin();
  iterator end();

  iterator erase(iterator i);
  bool empty() const;

private:
  std::list<int> available_;
};

// Can we do all these pointer based comparison functions with one function templated on the comparison? Yes, but not
// without wrapping the overloaded fitterThan function into a templated FitterThan class, and the same for lessCrowded,
// etc. This would not save us a great deal!
template <class Pointer>
bool fitterThanP(Pointer lhs, Pointer rhs)  // (Changed to fitterThanP from fitterThan, as there were ambiguity...
{                                           // ...issues with overloading fitterThan.)
  return fitterThan(*lhs, *rhs);
}

template <class Pointer>
bool olderThanP(Pointer lhs, Pointer rhs)
{
  return olderThan(*lhs, *rhs);
}

template <class Pointer>
bool lessCrowdedP(Pointer lhs, Pointer rhs)  // (Changed to lessCrowdedP from lessCrowded - ambiguity issues.)
{
  return lessCrowded(*lhs, *rhs);
}

template <class Pointer>
bool sortBeforeP(Pointer lhs, Pointer rhs)  // (Changed to sortBeforeP from sortBefore - ambiguity issues.)
{
  return sortBefore(*lhs, *rhs);
}

template <class Pointer>
class PointerEqual
{
public:
  explicit PointerEqual(Pointer rhs) :
  rhs_(rhs)
  {
  }

  bool operator()(Pointer lhs)
  {
    return *lhs == *rhs_;
  }

private:
  Pointer rhs_;
};

template <class Pointer>
class PointerNotEqual
{
public:
  explicit PointerNotEqual(Pointer rhs) :
  rhs_(rhs)
  {
  }

  bool operator()(Pointer lhs)
  {
    return *lhs != *rhs_;
  }

private:
  Pointer rhs_;
};

//----------------------------------------------------------------------------------------------------------------------

template <class Iterator>
void normalize(const Iterator& first, const Iterator& last)
{
  // Normalize should not be called on an array with negative values.
  assert(*std::min_element(first, last) >= 0);

  double sum = std::accumulate(first, last, 0.0);
  if (sum == 0.0)
  {
    // All elements are (hopefully) zero, which would result in a division by zero.
    // when calculating selection probabilities, this means that each solution is equally fit. They should therefore be
    // selected with equal probability, so each element is replaced with 1/n, where n is the number of items.
    double newValue = 1 / static_cast<double>(last - first);
    std::fill(first, last, newValue);
  }
  else
  {
    std::transform(first, last, first, [sum](int value){return value / sum;}); // The container MUST contain doubles.
  }
}

#endif
