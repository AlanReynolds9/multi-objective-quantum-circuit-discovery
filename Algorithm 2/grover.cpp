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
#include <algorithm>
#include "state.h"
#include "grover.h"

using std::vector;


namespace grover
{
  Problem::Problem(const vector<GateCreator>& gateCreator, const vector<PermittedControls>& permittedControls,
                   const vector<GateCost>& gateCosts, double gateCostRange, int numQbits, int maxLength,
                   double meanRandomLength, double badErrorCutoff, double gateOptProb, int numInsertions,
                   int minInsertionLength, int maxInsertionLength, double angleTolerance) :
  circuit::Problem(gateCreator, permittedControls, gateCosts, gateCostRange, numQbits, maxLength, meanRandomLength,
                   badErrorCutoff, gateOptProb, numInsertions, minInsertionLength, maxInsertionLength, angleTolerance)
  {
  }


  QbitOptions Problem::qbitInputOptions(int qbit) const
  {
    return QbitOptions::alwaysZero;
  }


  int Problem::numObj() const
  {
    return 2;
//    return 3;
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
    if (objNum == 0)
    {
      return gate_cost_range;  // A fake value for the maximum network cost. Provided by the user.
    }
    else  // objNum == 1
    {
      return 1.0;  // Maximum value for both error measures.
    }
  }

  //--------------------------------------------------------------------------------------------------------------------

  Solution::Solution(const Problem& problem) :
  circuit::Solution(problem)
  {
  }


  Solution Solution::clone() const
  {
    Solution child{static_cast<const Problem&>(problem_)};  // Yick!
    child.clone_from(*this);

    return child;
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
    using std::endl;

    if (circuit_)
    {
      out << "Circuit:" << endl << *circuit_;
    }
    else
    {
      out << "Circuit not yet constructed." << endl;
    }

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
    using std::endl;

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
    error_ = {1.1, 1.1};
  }


  int Solution::num_simulations() const
  {
    // Perform a simulation for each basis state.
    return dim(problem_.numQbits());
  }


  void Solution::prepare_(int sim)
  {
    // To prepare for simulation number 'sim', simply set the Oracles to mark basis state |sim>.
    *marked_state = sim;
  }


  State Solution::start_state(int sim) const
  {
    // The start state is always |0>.
    return State(problem_.numQbits(), 0);
  }


  State Solution::target_state(int sim) const
  {
    // The target state is |sim>.
    return State(problem_.numQbits(), sim);
  }


  void Solution::calculate_error()
  {
    // Overlaps for the various inputs have already been calculated. We just need to use these to construct the correct
    // error values.

    // Initialize
    double& overall_error = error_[0];
    double& worst_case_error = error_[1];
    overall_error = worst_case_error = -1.0;
    double errorTotal{0};

    // For each basis state...
    for (auto i = 0; i < num_simulations(); ++i)  // num_simulations() is the number of basis states.
    {
      // ...check if the overlap is the worst so far
      double thisError = 1 - std::norm(overlaps_[i]);  // The probability that measurement will produce the wrong value.
      if (thisError > worst_case_error)
      {
        worst_case_error = thisError;
      }

      // ...create an 'errorTotal', to be used in determining the 'overall_error'.
      errorTotal += thisError;
    }

    overall_error = errorTotal / num_simulations();
  }


  std::ostream& operator<<(std::ostream& out, const Solution& solution)
  {
    solution.output(out);
    return out;
  }


  void Solution::calculate_symbolic()
  {
    // Coefficients for the various coefficients in the (symbolic) overlap have already been calculated. These are
    // coefficients for the constant term, cosine term, sine term, and possibly more in future. The cosines and sines
    // are of an angle, i.e. a gate parameter, that is currently unfixed and that we wish to optimize. This function
    // converts these into coefficients in a function of the gate angle that is related to (and in the case of the
    // Grovers problem, actually is) the overall error, to be minimized. Starting with the 3 coefficients (constant,
    // cosine, sine), this error surrogate will have 6 coefficients (const, cosine, sine, cos^2, sin^2, cosSin). These
    // are stored in primary_error_coeffs.

    using std::norm, std::real;

    // Initialize
    primary_error_coeffs.assign(3, vector<double>(3, 0.0));  // We use 6 of the 9 slots here. (Triangular matrix)

    // Sum the errors for each input, ignoring the leading 1 in each case. In other words, starting from zero, we
    // subtract the norm of the overlap for each input.
    for (auto i = 0; i < num_simulations(); ++i)
    {
      const auto& overlap = symbolic_overlaps[i];
      primary_error_coeffs[0][0] -= norm(overlap[0]);  // Constant part, -a*a. (* is conjugation, not multiplication!)
      primary_error_coeffs[0][1] -= 2 * real(conj(overlap[0]) * overlap[1]);  // Cosine part, -(a*b + b*a).
      primary_error_coeffs[0][2] -= 2 * real(conj(overlap[0]) * overlap[2]);  // Sine part, -(a*c + c*a).
      primary_error_coeffs[1][1] -= norm(overlap[1]);  // Cos^2 part, -b*b.
      primary_error_coeffs[1][2] -= 2 * real(conj(overlap[1]) * overlap[2]);  // SinCos part, -(b*c + c*b).
      primary_error_coeffs[2][2] -= norm(overlap[2]);  // Sin^2 part, -c*c.
    }

    // Average by dividing by the number of simulations. Unnecessary, but we'll leave this for now.
    for (int i = 0; i < 3; ++i)
    {
      for (int j = i; j < 3; ++j)
      {
        primary_error_coeffs[i][j] /= num_simulations();
      }
    }
  }

}  // namespace grover
