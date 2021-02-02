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

#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <climits>
#include <vector>
#include "constants.h"

namespace Simulator
{
  constexpr int maxBits{sizeof(unsigned) * CHAR_BIT};
};


cmplx dotProduct(const std::vector<cmplx>& lhs, const std::vector<cmplx>& rhs);


class SimulatorIndexManager
{
  // A helper class for the GateSimulators. Uses the stored identities of the target and control bits to determine the
  // row of the matrix to be used in calculating a specified output coefficient index, and the index of the input
  // coefficient that corresponds with a specified output coefficient - column combination. Also provides access to the
  // control bits.
  // This part of the Simulator does not need to be updated after gate mutation, as it does not depend on gate
  // parameters.
public:
  SimulatorIndexManager(const std::vector<int>& targets, const std::vector<int>& controls);

  int matrixSize() const;
  const std::vector<int>& controls() const;
  unsigned rowIndex(const std::bitset<Simulator::maxBits>& outBits) const;
  unsigned inputIndex(const std::bitset<Simulator::maxBits>& outBits, unsigned column) const;

  void output(std::ostream& out) const;

private:
  void set_masks();

private:
  const int matrix_size;
  const std::vector<int> targets_;
  const std::vector<int> controls_;

  std::bitset<Simulator::maxBits> in_mask;  // Ones give input bits whose values are determined by the choice of...
                                            // ...column of the matrix.
  std::bitset<Simulator::maxBits> out_mask;  // Ones give input bits whose values are the same as the output bits.
  std::vector<std::bitset<Simulator::maxBits> > in_bits_set_by_column;  // Gives values of input bits covered by the...
                                                                        // ...in_mask.
  // (As far as I can tell, in_mask might be unecessary, except that it is used to initialize out_mask.)
};

std::ostream& operator<<(std::ostream& out, const SimulatorIndexManager& simulator);


class GateSimulator
{
  // Base class for the two types of GateSimulator.
  // (There is probably a better way to do this. Templates, following a policy class approach??)
public:
  // (I should probably add a virtual destructor?)
  virtual std::vector<cmplx> apply(const std::vector<cmplx>& state,
                                   const SimulatorIndexManager& indexManager) const = 0;
  virtual void output(std::ostream& out) const = 0;
};

std::ostream& operator<<(std::ostream& out, const GateSimulator& simulator);


class DenseGateSimulator : public GateSimulator
{
public:
  DenseGateSimulator(std::vector<std::vector<cmplx> >&& matrix, bool copyIfControlUnset = true);
  DenseGateSimulator(const std::vector<std::vector<cmplx> >& matrix, bool copyIfControlUnset = true);

  virtual std::vector<cmplx> apply(const std::vector<cmplx>& state,
                                   const SimulatorIndexManager& indexManager) const override;

  void output(std::ostream& out) const override;

private:
  std::vector<std::vector<cmplx> > matrix_;  // (Efficiency: copying.)
  bool copy_if_control_unset;  // If false, output state is zero when control is unset. Used, for example, when...
};           // ...considering just the 'cos' part of a Rotation matrix. (Could also be used for gradient calculations.)


class SparseGateSimulator : public GateSimulator
{
public:
  SparseGateSimulator(std::vector<std::vector<std::pair<int, cmplx> > >&& matrix, bool copyIfControlUnset = true);
  SparseGateSimulator(const std::vector<std::vector<std::pair<int, cmplx> > >& matrix, bool copyIfControlUnset = true);

  virtual std::vector<cmplx> apply(const std::vector<cmplx>& state,
                                   const SimulatorIndexManager& indexManager) const override;

  void output(std::ostream& out) const override;

private:
  std::vector<std::vector<std::pair<int, cmplx> > > matrix_;  // Each row is sparse, rather than entire matrix.
  bool copy_if_control_unset;  // If false, output state is zero when control is unset. Used, for example, when...
};           // ...considering just the 'cos' part of a Rotation matrix. (Could also be used for gradient calculations.)


class SwapSimulator
{
public:
  SwapSimulator(int numQbits, int swapBit1, int swapBit2);

  std::vector<cmplx> apply(const std::vector<cmplx>& state) const;

private:
  const int num_qbits;
  const int swap_bit1;
  const int swap_bit2;
};


#endif  // SIMULATOR_H
