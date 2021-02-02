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

#define QICLIB_DONT_USE_NLOPT  // We haven't installed NLopt and don't want to use NLopt dependent features
#define ARMA_DONT_USE_WRAPPER

#include <cassert>
#include <iostream>
#include "dictionary.h"
#include "state.h"
#include "qiclib/include/QIClib"
#include "QIClib-additions/apply_grad_ctrl.hpp"


arma::uword dim(int numQbits)
{
  // Get the dimensionality of the Hilbert space from the number of qbits.
  return arma::uword(1) << numQbits;
}


arma::uvec swapPermutation(int numQbits, int bit1, int bit2)
{
  // Return, as an arma::uvec for QIClib, the permutation associated with swapping qbits bit1 and bit2. So, with 6
  // qbits, swapping qbits 1 and 3 transforms 012345 into 032145. Each of these qbit indices then must be converted to
  // an index for QIClib (i.e. i -> 6 - i.) This means that the transformation, from QIClib's perspective, is
  // 654321 -> 634521. Of course, this is 123456 -> 125436. It is this last permutation, i.e. 125436, that we want. Of
  // course, just to make things even worse, 125436 must be stored in a uvec that starts at index 0! (QIClib's qbit
  // indexing system is a pain!)

  // (To avoid the '+ 1' and '- 1's that appear below, it might make sense to create a vector of size numQbits + 1 and
  // ignore the zeroth element. Then, at the end, copy the elements (apart from element 0) one by one into the
  // arma::uvec. This might be clearer.)

  // Create the 'non-permutation', i.e. a uvec storing 123456
  arma::uvec permutation(numQbits);
  for (auto i = 0; i < numQbits; ++i)
  {
    permutation(i) = i + 1;
  }

  // Adjust the two qbits that are swapped.
  permutation(qicLibQbitIndex(bit1, numQbits) - 1) = qicLibQbitIndex(bit2, numQbits);
  permutation(qicLibQbitIndex(bit2, numQbits) - 1) = qicLibQbitIndex(bit1, numQbits);

  return permutation;
}


State::State(int numQbits, int init) :
num_qbits{numQbits},
state_vector(dim(numQbits), arma::fill::zeros)
{
  state_vector(init) = 1;  // Playing safe for now. (Armadillo uses () for safe element access and [] for quick/unsafe.
}


State::State(int numQbits, arma::cx_vec&& stateVector) :
num_qbits(numQbits),
state_vector(std::move(stateVector))
{
}


cmplx State::overlap(const State& rhs) const
{
  return arma::cdot(state_vector, rhs.state_vector);
}


State State::transform(const Transformation& trans) const
{
  // Note that apply_ctrl makes no change to state_vector. It creates a new state and returns that.
  return {num_qbits, qic::apply_ctrl(state_vector, trans.matrix(), trans.controlBits(), {trans.target()})};
}


State State::gradTransform(const Transformation& trans) const
{
  // Note that apply_grad_ctrl makes no change to state_vector. It creates a new state and returns that.
  return {num_qbits, qic::apply_grad_ctrl(state_vector, trans.matrix(), trans.controlBits(), {trans.target()})};
}


State State::swapQbits(int bit1, int bit2) const
{
  // Note that sysperm makes no change to state_vector. It creates a new state and returns that.
  return {num_qbits, qic::sysperm(state_vector, swapPermutation(num_qbits, bit1, bit2))};
}


State State::mark(int toMark) const
{
  assert(0 <= toMark && toMark < dim(num_qbits));
  arma::cx_vec newStateVector = state_vector;
  newStateVector(toMark) = -newStateVector(toMark);  // Playing safe for now. (Armadillo uses () for safe element...
  return {num_qbits, std::move(newStateVector)};     // ...access, [] for unsafe.)
}


State State::fourier() const
{
  // I am following the conventions of N. David Mermin. As a result the fft in armadillo gives the complex conjugate of
  // the result I want. I need ifft.
  return {num_qbits, arma::ifft(state_vector) * std::sqrt(dim(num_qbits))};
}


void State::output(std::ostream& out) const
{
  out << state_vector(0);  // Playing safe for now. (Armadillo uses () for safe element access.)
  for (auto i = 1; i < dim(num_qbits); ++i)
  {
    out << " " << state_vector(i);  // Playing safe for now.
  }
}


cmplx stateOverlap(const State& lhs, const State& rhs)
{
  return lhs.overlap(rhs);
}


std::ostream& operator<<(std::ostream& out, const State& state)
{
  state.output(out);
  return out;
}
