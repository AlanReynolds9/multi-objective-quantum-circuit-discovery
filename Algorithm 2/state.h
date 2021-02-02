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

#ifndef STATE_H
#define STATE_H

#include <vector>
#include "constants.h"

class GateSimulator;
class SimulatorIndexManager;


int dim(int numQbits);


class State
{
public:
  explicit State(int numQbits, int init = 0);  // Here 'init' indicates the selected basis state.
  State(int numQbits, std::vector<cmplx>&& stateVector);
  State(int numQbits, const std::vector<cmplx>& stateVector);

  cmplx operator[](int i) const;  // Do we really need this?

  // Inner product of two states.
  cmplx overlap(const State& rhs) const;

  // Apply a GateSimulator to the state. Doesn't actually change the state, but creates a new one.
  void transform(const GateSimulator& simulator, const SimulatorIndexManager& indexManager);

  // Transform the state by swapping the contents of two qbits.
  // Doesn't actually change the state, but creates a new one.
  void swapQbits(int bit1, int bit2);

  // Apply a mark to a component of the state vector, i.e. multiply that component by -1.
  void mark(int toMark);

  // Apply the fast fourier transform to the state.
  // (Consider allowing const access to the state_vector - otherwise we could end up adding a member function to the
  // class whenever we work on a new problem. Alternatively, we could create single member function that takes a
  // function to apply to state_vector.)
  State fourier() const;

  void output(std::ostream& out) const;

  // Function for creating the big matrix that transforms basis states into the provided output states.
  static std::vector<std::vector<cmplx> > bigMatrix(const std::vector<State>& outputStates);

private:
  int num_qbits;
  std::vector<cmplx> state_vector;
};

cmplx stateOverlap(const State& lhs, const State& rhs);  // (Make a static member of State?)
std::ostream& operator<<(std::ostream& out, const State& state);

#endif  // STATE_H
