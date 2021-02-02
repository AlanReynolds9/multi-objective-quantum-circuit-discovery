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

#ifndef FRONT_H
#define FRONT_H

#include <map>
#include <shared_mutex>

class Front
{
  // Class that stores the best overall error at each possible value of cost. This is, of course, a non-dominated front,
  // but often a different one from that considered by the GA, which often includes a measure of worst case error too.
  // This is used to allow this information to be passed to the numerical algorithm code, which can then determine
  // whether the numerical search looks promising or not.
public:
  void clear();
  double best(long cost) const;  // Best overall error found with circuit cost no greater than 'cost'.
  void add(long cost, double error);

private:
  double best_unthreaded(long cost) const;

private:
  std::map<long, double> front_;
  mutable std::shared_timed_mutex mtx_;  // C++: All I want is a shared_mutex, but the Mac compiler can't find it!
};

#endif  // FRONT_H
