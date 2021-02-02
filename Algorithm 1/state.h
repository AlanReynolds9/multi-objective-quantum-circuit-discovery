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

#ifndef STATE_H
#define STATE_H

#include <armadillo>  // We don't need all of QIClib yet.
#include "transformation.h"

extern arma::uword dim(int numQbits);

class State
{
public:
  explicit State(int numQbits, int init = 0);  // Here 'init' indicates the selected basis state.
  State(int numQbits, arma::cx_vec&& stateVector);  // Added to allow a Gate to apply itself to a State. Used by...
                                                    // ...transform() and swapQbits().
  // Inner product of two states.
  cmplx overlap(const State& rhs) const;

  // Apply a unitary transformation, specified by a matrix, to a specified set of qbits, with a specified set of
  // controls. Doesn't actually change the state, but creates a new one.
  State transform(const Transformation& trans) const;

  // Apply a gradient transformation, specified by a matrix, to a specified set of qbits, with a specified set of
  // controls. Here, we are wishing to find the gradient of some other transformation. Hence, given a basis state with
  // unset control bits, the gradient is ZERO, i.e. rather than simply leaving the state alone when control bits are not
  // set, we annihilate it instead.
  State gradTransform(const Transformation& trans) const;

  // Transform the state by swapping the contents of two qbits. (Note that QIClib has functionality for more general
  // permutations.) Doesn't actually change the state, but creates a new one.
  State swapQbits(int bit1, int bit2) const;

  // Apply a mark to a component of the state vector, i.e. multiply that component by -1.
  State mark(int toMark) const;

  // Apply the fast fourier transform to the state.
  // (Consider allowing const access to the state_vector - otherwise we could end up adding a member function to the
  // class whenever we work on a new problem. (Alternatively, we could create single member function that takes an
  // armadillo function to apply to state_vector.)
  State fourier() const;

  void output(std::ostream& out) const;

private:
  // While it might be nice to use a friendly std::vector here, this would necessitate converting back and forth
  // whenever we apply a gate to the state. Given the number of simulations required, this seems inefficient.
  int num_qbits;
  arma::cx_vec state_vector;
};

cmplx stateOverlap(const State& lhs, const State& rhs);  // (Make a static member of State?)
std::ostream& operator<<(std::ostream& out, const State& state);

#endif  // STATE_H
