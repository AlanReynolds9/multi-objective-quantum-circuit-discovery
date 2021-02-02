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

#include <stdexcept>
#include <cmath>
#include "newton.h"


Newton::Newton(const std::function<double(double)>& f, const std::function<double(double)>& df,
               const std::function<double(double)>& ddf) :
f_(f),
df_(df),
ddf_(ddf),
min_(0)
{
}


void Newton::minimize(double firstGuess, double epsilon)
{
  double guess = firstGuess;
  while (abs(df_(guess)) > epsilon)
  {
    double ddf = ddf_(guess);
    if (ddf > 0.0)
    {
      throw std::runtime_error("Second derivative is not positive.");
    }
    guess -= df_(guess) / ddf;
  }

  min_ = guess;
}


double Newton::minimum() const
{
  return min_;
}


double Newton::minimumF() const
{
  return f_(min_);
}
