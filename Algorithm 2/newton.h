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

#ifndef NEWTON_H
#define NEWTON_H

#include <functional>

class Newton
{
public:
  Newton(const std::function<double(double)>& f, const std::function<double(double)>& df,
         const std::function<double(double)>& ddf);

  void minimize(double firstGuess, double epsilon);  // Performs iterations until size of the gradient is less than...
                                                     // ...epsilon.
  double minimum() const;
  double minimumF() const;

private:
  std::function<double(double)> f_;
  std::function<double(double)> df_;
  std::function<double(double)> ddf_;

  double min_;
};

#endif // NEWTON_H
