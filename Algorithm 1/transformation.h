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

#ifndef TRANSFORMATION_H
#define TRANSFORMATION_H

#include <iostream>
#include <armadillo>
#include "constants.h"

class Transformation
{
  // A Transformation is created by a Gate and then passed to the State so that it can apply it to itself.
  // At present, this handles only simple gate transformations, i.e. multiple control qbits but a single target.
  // (Data members can all be constant. Given that this is such a lightweight class, one wonders whether it should
  // simply be a struct or whether we should just make it a friend of State.)
public:
  Transformation(cmplx u11, cmplx u12, cmplx u21, cmplx u22, arma::uword target, const arma::uvec& controlBits);

  //  State apply(const State& state) const;

  const arma::cx_mat22& matrix() const;  // Inline?
  const arma::uword& target() const;      // Perhaps should simply return by value, but may wish to have multiple...
  const arma::uvec& controlBits() const;  // ...target bits in future.

  void output(std::ostream& out) const;

private:
  const arma::cx_mat22 matrix_;
  const arma::uword target_;
  const arma::uvec control_bits;
};

std::ostream& operator<<(std::ostream& out, const Transformation& transformation);


#endif // TRANSFORMATION_H
