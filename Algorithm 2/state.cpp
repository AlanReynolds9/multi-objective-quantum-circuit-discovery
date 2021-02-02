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

#include <cassert>
#include <iostream>
#include <armadillo>
#include "state.h"
#include "simulator.h"

using std::vector;


int dim(int numQbits)  // Should this be 'unsigned'? Maybe, maybe not - getting a negative value for dim is a sure...
{                      // ...sign of overflow!
  // Get the dimensionality of the Hilbert space from the number of qbits.
  return 1 << numQbits;
}


State::State(int numQbits, int init) :
num_qbits{numQbits},
state_vector(dim(numQbits), 0.0)
{
  state_vector[init] = 1;
}


State::State(int numQbits, vector<cmplx>&& stateVector) :
num_qbits(numQbits),
state_vector(std::move(stateVector))
{
}


State::State(int numQbits, const vector<cmplx>& stateVector) :
num_qbits(numQbits),
state_vector(stateVector)
{
}


cmplx State::operator[](int i) const
{
  return state_vector[i];
}


cmplx State::overlap(const State& rhs) const
{
  return dotProduct(state_vector, rhs.state_vector);
}


void State::transform(const GateSimulator& simulator, const SimulatorIndexManager& indexManager)
{
  state_vector = simulator.apply(state_vector, indexManager);
}


void State::swapQbits(int bit1, int bit2)
{
  SwapSimulator simulator(num_qbits, bit1, bit2);
  state_vector = simulator.apply(state_vector);
}


void State::mark(int toMark)
{
  assert(0 <= toMark && toMark < dim(num_qbits));
  state_vector[toMark] = -state_vector[toMark];
}


State State::fourier() const  // For use when the State is represented by a std::vector.
{
  // I am following the conventions of N. David Mermin. As a result the fft in armadillo gives the complex conjugate of
  // the result I want. We will need to convert from a std::vector to an arma::cx_vec before applying arma::ifft and
  // then transform the result back again.

  arma::cx_vec armaState(state_vector);
  armaState = arma::ifft(armaState) * std::sqrt(dim(num_qbits));

  vector<cmplx> result(dim(num_qbits));
  for (int i = 0; i < dim(num_qbits); ++i)
  {
    result[i] = armaState(i);  // Playing safe with the Armadillo vector index for now.
  }

  return {num_qbits, result};
}


void State::output(std::ostream& out) const
{
  out << state_vector[0];
  for (auto i = 1; i < dim(num_qbits); ++i)
  {
    out << " " << state_vector[i];
  }
}


cmplx stateOverlap(const State& lhs, const State& rhs)
{
  return lhs.overlap(rhs);
}


vector<vector<cmplx> > State::bigMatrix(const vector<State>& outputStates)  // Used?
{
  // (All this does is take the transpose.)
  vector<vector<cmplx> > matrix(outputStates.size(), vector<cmplx>(outputStates.size()));
  for (int col = 0; col < outputStates.size(); ++col)
  {
    for (int row = 0; row < outputStates.size(); ++row)
    {
      matrix[row][col] = outputStates[col][row];
    }
  }

  return matrix;
}


std::ostream& operator<<(std::ostream& out, const State& state)
{
  state.output(out);
  return out;
}
