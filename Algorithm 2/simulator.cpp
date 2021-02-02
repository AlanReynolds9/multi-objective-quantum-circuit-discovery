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

#include <climits>
#include <bitset>
#include <numeric>  // For std::inner_product.
#include "simulator.h"
#include "state.h"  // Just for use of dim(), at present.

using std::pair;
using std::vector;
using std::bitset;
using std::endl;

using Simulator::maxBits;


cmplx dotProduct(const vector<cmplx>& lhs, const vector<cmplx>& rhs)
{
  return std::inner_product(lhs.begin(), lhs.end(), rhs.begin(), cmplx(), std::plus<>(),
                            [](cmplx lhs, cmplx rhs){return conj(lhs) * rhs;});
}

//----------------------------------------------------------------------------------------------------------------------

SimulatorIndexManager::SimulatorIndexManager(const vector<int>& targets, const vector<int>& controls) :
matrix_size(dim(static_cast<int>(targets.size()))),
targets_(targets),
controls_(controls),
in_mask(0),
in_bits_set_by_column(matrix_size)
{
  set_masks();
}


void SimulatorIndexManager::set_masks()
{
  // Applying a gate involves calculating the coefficient, in the output state, of all the basis states. Each output
  // basis state is considered in turn. For 5 qbits, basis states are indexed from 0 to 31, i.e. |0> = |00000>, |1> =
  // |00001>, |2> = |00010>, ..., |31> = |11111>. The coefficient of more than one input basis state may contribute to
  // that for the output basis state under consideration. The set of input basis state indices that contribute to an
  // output coefficient is constructed by inserting all possible combinations of the states of the target qbits into the
  // output basis state being calculated. To do this, we need a mask that covers the target bits and another that covers
  // the rest.
  for (auto target : targets_)
  {
    in_mask.set(target);
  }
  out_mask = ~in_mask;

  for (unsigned column = 0; column < matrix_size; ++column)
  {
    bitset<maxBits> columnBits(column);  // Just the integer, 'column', in binary.
    for (auto i = 0; i < targets_.size(); ++i)
    {
      if (columnBits[i])
      {
        in_bits_set_by_column[column].set(targets_[i]);
      }
    }
  }
}


int SimulatorIndexManager::matrixSize() const
{
  return matrix_size;
}


const vector<int>& SimulatorIndexManager::controls() const
{
  return controls_;
}


unsigned SimulatorIndexManager::rowIndex(const bitset<maxBits>& outBits) const
{
  // Get the row of the transformation matrix that must be used to calculate contributions to the coefficient of the
  // basis state represented by 'outBits'.
  //
  // Example:
  // Suppose that are 5 qbits and the targets are qbits 3, 0, 2, in that order. The matrix is 8x8. Row zero contributes
  // to the output coefficients of basis states where qbits 3, 0 and 2 are all zero. Row one maps to basis states where
  // qbits 3, 0 and 2 are 1, 0, 0 respectively, and so on. To find the required row index, we simply find the binary
  // number that matches the settings of qbits 3, 0 and 2 in 'outBits'. So if outBits is 11011, then qbits 3, 0 and 2
  // take the values 1, 1 and 0, so the row index is 011 = 3.
  bitset<maxBits> rowBits;
  for (auto i = 0; i < targets_.size(); ++i)
  {
    if (outBits[targets_[i]])
    {
      rowBits.set(i);
    }
  }

  return static_cast<unsigned>(rowBits.to_ulong());
}


unsigned SimulatorIndexManager::inputIndex(const bitset<maxBits>& outBits, unsigned column) const
{
  // Get the input index that, when combined with the specified column of the transformation matrix, produces a
  // contribution to the coefficient of the basis state represented by 'outBits'.
  bitset<maxBits> inBits = in_bits_set_by_column[column];
  inBits = (inBits & in_mask) | (outBits & out_mask);  // '& in_mask' might be unecessary - no bit unset in in_mask...
  return static_cast<unsigned>(inBits.to_ulong());     // ...will be set in inBits.
}


void SimulatorIndexManager::output(std::ostream& out) const
{
  out << "Targets:";
  for (auto target : targets_)
  {
    out << " " << target;
  }
  out << endl;

  out << "Controls:";
  if (controls_.empty())
  {
    out << " none";
  }
  for (auto control : controls_)
  {
    out << " " << control;
  }
  out << endl;
}


std::ostream& operator<<(std::ostream& out, const SimulatorIndexManager& indexManager)
{
  indexManager.output(out);
  return out;
}

//----------------------------------------------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const GateSimulator& simulator)
{
  simulator.output(out);
  return out;
}

//----------------------------------------------------------------------------------------------------------------------

DenseGateSimulator::DenseGateSimulator(vector<vector<cmplx> >&& matrix, bool copyIfControlUnset) :
matrix_(std::move(matrix)),
copy_if_control_unset(copyIfControlUnset)
{
}


DenseGateSimulator::DenseGateSimulator(const vector<vector<cmplx> >& matrix, bool copyIfControlUnset) :
matrix_(matrix),
copy_if_control_unset(copyIfControlUnset)
{
}


vector<cmplx> DenseGateSimulator::apply(const vector<cmplx>& state, const SimulatorIndexManager& indexManager) const
{
  // Check that the index manager is compatible.
  assert(matrix_.size() == indexManager.matrixSize());

  // Create output state.
  // (A note about efficiency. In serial, we could knock about 8% off total run time (for 4 qbit Grover) by giving the
  // Simulator a 'new_state' member, with memory allocation only happening once, rather than every time the simulator is
  // applied. However, this is not thread-safe - a simulator is likely shared between multiple threads, so providing
  // write access to such a member is dangerous. While we could try locking and unlocking new_state, this seems
  // excessive.)
  int dimensions = static_cast<int>(state.size());
  vector<cmplx> newState(dimensions, cmplx(0.0));

  // Calculate.
  for (unsigned outIndex = 0; outIndex < dimensions; ++outIndex)  // outIndex identifies which basis state's coeff is...
  {                                                               // ...being calculated.
    bitset<maxBits> outBits{outIndex};  // (Separate this block (before 'if (allControlsSet())') into a function?)
    bool allControlsSet = true;         // (Perhaps put into SimulatorIndexManager.)
    const auto& controls = indexManager.controls();
    for (auto control : controls)
    {
      if (!outBits[control])
      {
        allControlsSet = false;
        break;
      }
    }

    if (allControlsSet)
    {
      auto row = indexManager.rowIndex(outIndex);
      for (unsigned col = 0; col < matrix_.size(); ++col)
      {
        auto inIndex = indexManager.inputIndex(outIndex, col);

        // We can now combine the coefficient of the input basis state with the correct member of the matrix to get the
        // contribution to the coefficient of the output basis state.
        newState[outIndex] += matrix_[row][col] * state[inIndex];
      }
    }
    else if (copy_if_control_unset)
    {
      // A control bit is not set. Basis state is unaffected by the gate, so it's coefficient is copied to the output.
      newState[outIndex] = state[outIndex];
    }
  }

  return newState;
}


void DenseGateSimulator::output(std::ostream& out) const
{
  for (const auto& row : matrix_)
  {
    for (auto element : row)
    {
      out << element << " ";
    }
    out << endl;
  }
  out << endl;
}

//----------------------------------------------------------------------------------------------------------------------

SparseGateSimulator::SparseGateSimulator(vector<vector<pair<int, cmplx> > >&& matrix, bool copyIfControlUnset) :
matrix_(std::move(matrix)),
copy_if_control_unset(copyIfControlUnset)
{
}


SparseGateSimulator::SparseGateSimulator(const vector<vector<pair<int, cmplx> > >& matrix, bool copyIfControlUnset) :
matrix_(matrix),
copy_if_control_unset(copyIfControlUnset)
{
}


vector<cmplx> SparseGateSimulator::apply(const vector<cmplx>& state, const SimulatorIndexManager& indexManager) const
{
  // Check that the index manager is compatible.
  assert(matrix_.size() == indexManager.matrixSize());

  // Create output state.
  // (A note about efficiency. In serial, we could knock about 8% off total run time (for 4 qbit Grover) by giving the
  // Simulator a 'new_state' member, with memory allocation only happening once, rather than every time the simulator is
  // applied. However, this is not thread-safe - a simulator is likely shared between multiple threads, so providing
  // write access to such a member is dangerous. While we could try locking and unlocking new_state, this seems
  // excessive.)
  int dimensions = static_cast<int>(state.size());
  vector<cmplx> newState(dimensions, cmplx(0.0));

  // Calculate.
  for (unsigned outIndex = 0; outIndex < dimensions; ++outIndex)  // outIndex identifies which basis state's coeff is...
  {                                                               // ...being calculated.
    bitset<maxBits> outBits{outIndex};  // (Separate block into function. Also repeats code from DenseGateSimulator.)
    bool allControlsSet = true;
    const auto& controls = indexManager.controls();
    for (auto control : controls)
    {
      if (!outBits[control])
      {
        allControlsSet = false;
        break;
      }
    }

    if (allControlsSet)
    {
      auto row = indexManager.rowIndex(outIndex);
      for (const auto& [col, val] : matrix_[row])  // (The 'const &' makes a noticeable difference to code speed when...
      {                                            // using sparse simulators for normal gates. (About 8.5% overall.)
        auto inIndex = indexManager.inputIndex(outIndex, col);

        // We can now combine the coefficient of the input basis state with the correct member of the matrix to get the
        // contribution to the coefficient of the output basis state.
        newState[outIndex] += val * state[inIndex];
      }
    }
    else if (copy_if_control_unset)
    {
      // A control bit is not set. Basis state is unaffected by the gate, so it's coefficient is copied to the output.
      newState[outIndex] = state[outIndex];
    }
  }

  return newState;
}


void SparseGateSimulator::output(std::ostream& out) const
{
  for (const auto& row : matrix_)
  {
    int col = 0;
    for (auto [nonZeroCol, val] : row)
    {
      for (; col < nonZeroCol; ++col)
      {
        out << "0 ";
      }
      out << val << " ";
      ++col;
    }
    for (; col < static_cast<int>(matrix_.size()); ++col)
    {
      out << "0 ";
    }
    out << endl;
  }
}

//----------------------------------------------------------------------------------------------------------------------

SwapSimulator::SwapSimulator(int numQbits, int swapBit1, int swapBit2) :
num_qbits(numQbits),
swap_bit1(swapBit1),
swap_bit2(swapBit2)
{
}


vector<cmplx> SwapSimulator::apply(const vector<cmplx>& state) const
{
  // (Much like the case with vector<bool>, swap WILL NOT WORK on the bits in a bitset.)
  vector<cmplx> newState(dim(num_qbits));
  for (unsigned outIndex = 0; outIndex < dim(num_qbits); ++outIndex)  // outIndex identifies which basis state's...
  {                                                                   // ...coefficient is being calculated.
    bitset<maxBits> outBits{outIndex};
    if (outBits[swap_bit1] == outBits[swap_bit2])
    {
      newState[outIndex] = state[outIndex];
    }
    else
    {
      auto inBits = outBits;

      // The next two lines do the swap, since we know that the swaps bits are taking opposite values for this basis
      // state.
      inBits.flip(swap_bit1);
      inBits.flip(swap_bit2);

      unsigned inIndex = static_cast<unsigned>(inBits.to_ulong());
      newState[outIndex] = state[inIndex];
    }
  }

  return newState;
}
