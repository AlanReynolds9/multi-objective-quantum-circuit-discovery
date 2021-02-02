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

#include "transformation.h"

Transformation::Transformation(cmplx u11, cmplx u12, cmplx u21, cmplx u22, arma::uword target,
                               const arma::uvec& controlBits) :
matrix_{u11, u21, u12, u22},  // Note the reordering - QIClib stores matrices in column-major order
target_(target),
control_bits(controlBits)
{
}


const arma::cx_mat22& Transformation::matrix() const
{
  return matrix_;
}


const arma::uword& Transformation::target() const
{
  return target_;
}


const arma::uvec& Transformation::controlBits() const
{
  return control_bits;
}


void Transformation::output(std::ostream& out) const
{
  // Assumes that we start at the beginning of a line
  out << "Matrix" << std::endl;
  matrix_.print(out);
  out << "applied to qbit " << target_;
  if (control_bits.n_elem == 1)
  {
    out << " with control qbit " << control_bits[0] << "." << std::endl;
  }
  else if (control_bits.n_elem > 1)
  {
    out << " with control qbits " << control_bits[0];
    for (auto i = 1; i != control_bits.n_elem - 1; ++i)
    {
      out << ", " << control_bits[i];
    }
    if (control_bits.n_elem > 1)
    {
      out << " and " << control_bits[control_bits.size() - 1] << "." << std::endl;
    }
  }
  else  // No controls
  {
    out << "." << std::endl;
  }
}


std::ostream& operator<<(std::ostream& out, const Transformation& gate)
{
  gate.output(out);
  return out;
}

