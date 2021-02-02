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

#include <cassert>
#include <algorithm>
#include "state.h"
#include "grover.h"

using std::endl;


namespace grover
{
  Problem::Problem(const std::vector<GateCreator>& gateCreator, const std::vector<PermittedControls>& permittedControls,
                   const std::vector<GateCost>& gateCosts, int numQbits, double meanRandomLength, double gateCostRange,
                   int maxCircuitLength, double targetError, long expectedReduction, bool randomizeParameters,
                   int iterationQuota, bool useCache) :
  circuit::Problem(gateCreator, permittedControls, gateCosts, numQbits, meanRandomLength, gateCostRange,
                   maxCircuitLength, targetError, expectedReduction, randomizeParameters, iterationQuota, useCache)
  {
  }


  QbitOptions Problem::qbitInputOptions(int qbit) const
  {
    return QbitOptions::alwaysZero;
  }


  int Problem::numObj() const
  {
    return 3;
  }


  bool Problem::objMaximized(int objNum) const
  {
    // All objectives are minimized at present.
    return false;
  }


  double Problem::objLowerBound(int objNum) const
  {
    // All objectives have a minimum value of 0.0.
    return 0.0;
  }


  double Problem::objUpperBound(int objNum) const
  {
    assert(0 <= objNum && objNum < numObj());
    if (objNum == 0 || objNum == 1)
    {
      return 1.0;  // Maximum value for both error measures.
    }
    else  // objNum == 2
    {
      return gate_cost_range;  // A fake value for the maximum network cost. Provided by the user.
    }
  }

  //--------------------------------------------------------------------------------------------------------------------

  Solution::Solution(const Problem& problem) :
  circuit::Solution(problem)
  {
  }

  
  Solution Solution::clone() const
  {
    return Solution(*this);
  }

  double Solution::overallError() const
  {
    return error_[0];
  }


  double Solution::worstCaseError() const
  {
    return error_[1];
  }


  void Solution::output(std::ostream& out) const
  {
    out << "Circuit:" << endl << circuit_;

    if (evaluated_)
    {
      out << "Overall error:    " << overallError() << endl;
      out << "Worst case error: " << worstCaseError() << endl;
      out << "Total cost:       " << cost() << endl;
    }
    else
    {
      out << "Solution unevaluated." << endl;
    }
  }


  void Solution::outputQuality(std::ostream& out) const
  {
    if (evaluated_)
    {
      out << overallError() << ", " << worstCaseError() << ", " << cost();
    }
    else
    {
      out << "Solution unevaluated." << endl;
    }
  }


  void Solution::assign_worst_error()
  {
    // Assign really bad error values! (Used when there is little point in evaluating the solution.)
    error_ = {1.0, 1.0};
  }


  void Solution::evaluate_error()
  {
    // The difficult bit of evaluating a circuit - performing simulations to determine error.

    // Initialize
    double& overall_error = error_[0];
    double& worst_case_error = error_[1];
    overall_error = worst_case_error = -1.0;
    double errorTotal{0};

    // For each basis state...
    auto numQbits = problem_.numQbits();
    for (auto i = 0; i < dim(numQbits); ++i)
    {
      // ...set the oracle gates to mark the selected basis state |i> as special
      circuit_.setMarkedState(i);

      // ...set the state state to |0> and the target state to |i>
      State startState{numQbits, 0};
      State targetState{numQbits, i};

      // ...calculate the overlap of the target state with the results of simulating with the indicated startState.
      auto overlap = circuit_.simulatedOverlap(startState, targetState);

      // ...check if this is the worst so far
      double thisError = 1 - std::norm(overlap);  // The probability that a measurement will produce the wrong value.
      if (thisError > worst_case_error)
      {
        worst_case_error = thisError;
      }

      // ...create an 'errorTotal', to be used in determining the 'overall_error'.
      errorTotal += thisError;
    }

    overall_error = errorTotal / dim(numQbits);
  }


  void Solution::evaluate_error_and_gradient()
  {
    // Calculates both errors (overall and worst case) and the gradient of the overall error. This should be more
    // efficient than calculating the error and the gradient separately, at least if there is more than one
    // parameterized gate present. (It shouldn't be much slower even if there is only one parameterized gate.)
    // NOTE: Worst case error is unnecessary at present.

    // Initialize
    double& overall_error = error_[0];
    double& worst_case_error = error_[1];
    overall_error = worst_case_error = -1.0;
    double errorTotal{0};
    gradient_.assign(numParameters(), 0.0);
    std::vector<double> gradTotalError(numParameters(), 0.0);

    // For each basis state...
    auto numQbits = problem_.numQbits();
    for (auto i = 0; i < dim(numQbits); ++i)
    {
      // ...set the oracle gates to mark the selected basis state |i> as special
      circuit_.setMarkedState(i);

      // ...set the state state to |0> and the target state to |i>
      State startState{numQbits, 0};
      State targetState{numQbits, i};

      // ...calculate the overlap of the target state with the results of simulating with the indicated startState,
      // and the gradient of this overlap
      auto [overlap, gradOverlap] = circuit_.simulatedOverlapAndGrad(startState, targetState);

      // ...check if the overlap is the worst so far
      double thisError = 1 - std::abs(overlap);  // The probability that a measurement will produce the wrong value.
      if (thisError > worst_case_error)
      {
        worst_case_error = thisError;
      }

      // ...add to the error total, to be used in determining the 'overall_error'
      errorTotal += thisError;

      // ...add the gradient of this error to the gradTotalError
      if (abs(overlap) != 0)  // Prevents division by zero. Zero overlap => zero gradient contribution.
      {                       // Do we need a tolerance??
        for (int paramNum = 0; paramNum < circuit_.numParameters(); ++paramNum)
        {
          gradTotalError[paramNum] -= real(conj(overlap) * gradOverlap[paramNum]) / abs(overlap);
          if (std::isnan(gradTotalError[paramNum]) || std::isinf(gradTotalError[paramNum]))
          {
            throw std::runtime_error("Last contribution to gradTotalError[paramNum] has resulted in 'infinity' or"
                                     "'not a number'!");
          }
        }
      }
    }

    overall_error = errorTotal / dim(numQbits);
    for (int paramNum = 0; paramNum < circuit_.numParameters(); ++paramNum)
    {
      gradient_[paramNum] = gradTotalError[paramNum] / dim(numQbits);
    }

    evaluated_ = true;
  }


  std::ostream& operator<<(std::ostream& out, const Solution& solution)
  {
    solution.output(out);
    return out;
  }

}  // namespace grover
