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

#ifndef FOURIER_H
#define FOURIER_H

#include "circuit.h"

namespace fourier
{
  class Problem : public circuit::Problem
  {
  public:
    // (At present, gate costs are stored in the CircuitContext. However, the gateCostRange is sent, via the Problem
    // class, into the circuit::Problem object. This is a little odd.)
    Problem(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
            const std::vector<GateCost>& gateCosts, int numQbits, double meanRandomLength, double gateCostRange,
            int maxCircuitLength, double targetError, long expectedReduction, bool randomizeParameters,
            int iterationQuota, bool useCache);

    QbitOptions qbitInputOptions(int qbit) const override;

    int numObj() const override;
    bool objMaximized(int objNum) const override;
    double objLowerBound(int objNum) const override;
    double objUpperBound(int objNum) const override;
  };
  

  class Solution : public circuit::Solution
  {
  public:
    using Problem = fourier::Problem;  // Required by the algorithm solution classes such as NSGA2Solution.

    explicit Solution(const Problem& problem);
    // Default destructor and copy constructor are fine. No assignment operator - what would it do with problem_!?
    
    // Create a 'clone' solution. Not necessarily an exact copy - the parent solution may have collected some
    // non-genetic baggage that the child does not need, or it may be useful for the child to 'know' its parent.
    Solution clone() const;

    double overallError() const;    // For accessing solution quality information by name, rather than by using...
    double worstCaseError() const;  // ...Solution::objective(x).

    // Two versions of the dominance relation, required by the MOMH code. Function 'beats' is used within the algorithm.
    // Function 'dominates' is used only by the Store object that stores the best solutions seen. A solution may be
    // eliminated from consideration for storage if it is beaten during the operation of the algorithm. Hence it is
    // important that solution that is not dominated by another (and that we might wish to be stored) should not be
    // beaten either. I.e. we need A beats B => A dominates B. These two functions will usually be the same, but we may
    // wish to use a modified dominance relation within the algorithm to encourage population diversity.
    inline bool dominates(const Solution& rhs) const;
    inline bool beats(const Solution& rhs) const;

    void output(std::ostream& out) const;
    void outputQuality(std::ostream& out) const;

  private:
    void assign_worst_error() override;
    void evaluate_error() override;  // The hard bit - requires multiple circuit simulations.
    void evaluate_error_and_gradient() override;  // More efficient than calculating the two separately.
  };


  // bool sortBefore(const Solution& lhs, const Solution& rhs);  // Use base class version
  std::ostream& operator<<(std::ostream& out, const Solution& solution);

  
  bool Solution::dominates(const Solution& rhs) const  // Inlined, as this gives a significant speed boost. (Clang)
  {
    // The simplest dominance relation would simply return true if *this is at least as good as rhs on all objectives,
    // and better in at least one. However, the use of floating point objectives complicates matters. It would be
    // possible for solutions to survive because one of the objective appears to be better by some tiny amount, when in
    // fact this is just numerical error (less likely?). Similarly, it is possible that solutions may be killed due to
    // one of their objectives having some small rounding error added to it (more likely?). Hence we require *this to be
    // better, on at least one objective, by at least 'tolerance' (to avoid the second problem), and at least as good,
    // up to a possible error of 'tolerance' (to avoid the first problem) in all other objectives.
    // NOTE: This does NOT give a partial order! We could choose to use my 'weighted objectives' approach from previous
    // research, which would deal with the problem and produce a partial order.
    // (If using NSGA II and large population sizes, it may be important that this function be efficient, as it is
    // called O(n^2) times when sorting the population, where n is the number of solutions being sorted. I tried moving
    // the cost comparison first, seeing as that is only comparing integers, but it made little difference.)

    assert(evaluated_);
    assert(rhs.evaluated_);
    using constants::tolerance;

    // First a special case - the circuit beats the target error.
    if (overallError() < problem_.targetError())
    {
      // A solution that meets the 'target' error is considered to have unbeatable error values. Hence if rhs has larger
      // cost, or has equal cost but does not beat the target, then this solution dominates.
      // ALTERNATIVE: We could, instead, have objective functions that correct the error to zero (or even -1.0, given
      // that I've seen negative error values crop up). However, at present, I would also want the objective functions
      // to return the 'calculated' error values too, just for interest's sake.
      return total_cost < rhs.total_cost || total_cost == rhs.total_cost &&
             rhs.overallError() >= problem_.targetError();
    }

    if (total_cost > rhs.total_cost || overallError() > rhs.overallError() + tolerance ||
        worstCaseError() > rhs.worstCaseError() + tolerance)
    {
      return false;
    }
    return total_cost < rhs.total_cost || overallError() < rhs.overallError() - tolerance ||
    worstCaseError() < rhs.worstCaseError() - tolerance;
  }


  bool Solution::beats(const Solution& rhs) const  // Inlined, as this seems to give a significant speed boost. (Clang)
  {
    return dominates(rhs);
  }

}  // namespace fourier


#endif  // FOURIER_H
