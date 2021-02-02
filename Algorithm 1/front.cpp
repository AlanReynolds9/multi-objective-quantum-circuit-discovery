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

#include <limits>
#include <algorithm>
#include "front.h"
#include <iostream>  // TEMP


void Front::clear()
{
  front_.clear();
}


double Front::best_unthreaded(long cost) const
{
  // Get the best error value found for a circuit with cost equal or less than 'cost'. Not thread-safe and therefore
  // private, but called by the thread-safe functions.

  // Find first element with cost greater than 'cost'.
  auto it = front_.upper_bound(cost);

  // If this is the first element, then no solutions have yet been found at 'cost' or less.
  if (it == front_.begin())
  {
    return std::numeric_limits<double>::infinity();
  }

  // If not, then the previous element contains the value we are looking for.
  return (--it)->second;
}


double Front::best(long cost) const
{
  // Get the best error value found for a circuit with cost equal or less than 'cost'. Thread-safe.
  std::shared_lock<std::shared_timed_mutex> lck{mtx_};  // Lock for reading.
  return best_unthreaded(cost);
}


void Front::add(long cost, double error)
{
  // If the new entry improves the front, add it.
  // (There are two searches in the map, one in best_unthreaded() and one in front_.lower_bound(cost). We could
  // duplicate code from best(), which would allow us to use insert() with a hint.)
  std::unique_lock<std::shared_timed_mutex> lck{mtx_};  // Lock for writing
  if (error < best_unthreaded(cost))
  {
    // Insert new element, getting the iterator
    auto i = front_.lower_bound(cost);
    if (i != front_.end() && i->first == cost)
    {
      i->second = error;
    }
    else
    {
      i = front_.insert(i, std::make_pair(cost, error));
    }

    // Delete any elements that follow the new one and have error at least as large
    ++i;  // Make iterator i point to the element after the new one.
    auto j = std::find_if(i, front_.end(), [error](auto val){return val.second < error;});  // Find next element that...
                                                                                            // ...has better error.
    front_.erase(i, j);  // Erase the elements that are worse in error and cost.
  }
}

