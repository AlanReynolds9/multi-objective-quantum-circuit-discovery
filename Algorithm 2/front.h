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

#ifndef FRONT_H
#define FRONT_H

#include <map>
#include <shared_mutex>
#include <limits>
#include <algorithm>

// Class for storing a Front, i.e. the quality of the solutions in a non-dominated set. I have templated it just in case
// we meet a scenario where, for example, it might be beneficial to use 'long' rather than 'double' as an objective
// type. (Numerical accuracy?)


template <typename Cost, typename Error>  // Cost and Error can, of course, be any pair of objective types.
class Front
{
  // Class that stores the best overall error at each possible value of cost. This is, of course, a non-dominated front,
  // but often a different one from that considered by the GA, which often includes a measure of worst case error too.
  // This is used to allow this information to be passed to the numerical algorithm code, which can then determine
  // whether the numerical search looks promising or not.
  // (Rather problem specific.)
public:
  void clear();
  double best(Cost cost) const;  // Best overall error found with circuit cost no greater than 'cost'.
  void add(Cost cost, Error error);
  void output(std::ostream& out) const;

private:
  double best_unthreaded(Cost cost) const;

private:
  std::map<Cost, Error> front_;
  mutable std::shared_timed_mutex mtx_;  // C++: All I want is a shared_mutex, but the Mac compiler can't find it!
};

template <typename Cost, typename Error>
std::ostream& operator<<(std::ostream& out, const Front<Cost, Error>& front);

//----------------------------------------------------------------------------------------------------------------------

template <typename Cost, typename Error>
void Front<Cost, Error>::clear()
{
  front_.clear();
}


template <typename Cost, typename Error>
double Front<Cost, Error>::best_unthreaded(Cost cost) const
{
  // Get the best error value found for a circuit with cost equal or less than 'cost'. Not thread-safe and therefore
  // private, but called by the thread-safe functions.

  // Find first element with cost greater than 'cost'.
  auto it = front_.upper_bound(cost);

  // If this is the first element, then no solutions have yet been found at 'cost' or less.
  if (it == front_.begin())
  {
    return std::numeric_limits<Error>::max();
  }

  // If not, then the previous element contains the value we are looking for.
  return (--it)->second;
}


template <typename Cost, typename Error>
double Front<Cost, Error>::best(Cost cost) const
{
  // Get the best error value found for a circuit with cost equal or less than 'cost'. Thread-safe.
  std::shared_lock<std::shared_timed_mutex> lck{mtx_};  // Lock for reading.
  return best_unthreaded(cost);
}


template <typename Cost, typename Error>
void Front<Cost, Error>::add(Cost cost, Error error)
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
    auto j = std::find_if(i, front_.end(),
                          [error](auto val){return val.second < error;});  // Find next element that has better error.
    front_.erase(i, j);  // Erase the elements that are worse in error and cost.
  }
}


template <typename Cost, typename Error>
void Front<Cost, Error>::output(std::ostream& out) const
{
  for (const auto& point : front_)
  {
    out << "(" << point.first << ", " << point.second << ")" << std::endl;
  }
}


template <typename Cost, typename Error>
std::ostream& operator<<(std::ostream& out, const Front<Cost, Error>& front)
{
  front.output(out);
  return out;
}

#endif  // FRONT_H
