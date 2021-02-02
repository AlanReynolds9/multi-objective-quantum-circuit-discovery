// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Modifications to LBFGS++ code copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
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

// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name> & Dirk Toewe <DirkToewe@GoogleMail.com>
// Under MIT license

#ifndef LINE_SEARCH_BRACKETING_H
#define LINE_SEARCH_BRACKETING_H

#include <Eigen/Core>
#include <stdexcept>  // std::runtime_error

namespace LBFGSpp
{
  ///
  /// The bracketing line search algorithm for LBFGS. Mainly for internal use.
  ///
  template <typename Scalar>
  class LineSearchBracketing
  {
  private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

  public:
    ///
    /// Line search by bracketing. Similar to the backtracking line search
    /// except that it actively maintains an upper and lower bound of the
    /// current search range.
    ///
    /// \param f      A function object such that `f(x, grad)` returns the
    ///               objective function value at `x`, and overwrites `grad` with
    ///               the gradient.
    /// \param fx     In: The objective function value at the current point.
    ///               Out: The function value at the new point.
    /// \param x      Out: The new point moved to.
    /// \param grad   In: The current gradient vector. Out: The gradient at the
    ///               new point.
    /// \param step   In: The initial step length. Out: The calculated step length.
    /// \param drt    The current moving direction.
    /// \param xp     The current point.
    /// \param param  Parameters for the LBFGS algorithm
    ///
    template <typename Foo>
    static void LineSearch(Foo& f, Scalar& fx, Vector& x, Vector& grad,
                           Scalar& step,
                           const Vector& drt, const Vector& xp,
                           const LBFGSParam<Scalar>& param)
    {
      // Check the value of step
      if (step <= Scalar(0))
        std::invalid_argument("'step' must be positive");

      // Save the function value at the current x
      const Scalar fx_init = fx;
      // Projection of gradient on the search direction
      const Scalar dg_init = grad.dot(drt);
      // Make sure d points to a descent direction
      if(dg_init > 0)
        std::logic_error("the moving direction increases the objective function value");

      const Scalar dg_test = param.ftol * dg_init;

      // Upper and lower end of the current line search range
      Scalar step_lo = 0,
      step_hi = std::numeric_limits<Scalar>::infinity();

      int iter;
      for (iter = 0; iter < param.max_linesearch; iter++)
      {
        // x_{k+1} = x_k + step * d_k
        x.noalias() = xp + step * drt;
        // Evaluate this candidate
        fx = f(x, grad);

        if (fx > fx_init + step * dg_test)
        {
          step_hi = step;
        }
        else
        {
          // Armijo condition is met
          if(param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
            break;

          const Scalar dg = grad.dot(drt);
          if (dg < param.wolfe * dg_init)
          {
            step_lo = step;
          }
          else
          {
            // Regular Wolfe condition is met
            if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE)
              break;

            if (dg > -param.wolfe * dg_init)
            {
              step_hi = step;
            }
            else
            {
              // Strong Wolfe condition is met
              break;
            }
          }
        }

        assert(step_lo < step_hi);

//        if (iter >= param.max_linesearch)  // This CANNOT happen here. It can happen AFTER exiting the for loop.
//          throw std::runtime_error("the line search routine reached the maximum number of iterations");

        if (step < param.min_step)
          throw std::runtime_error("the line search step became smaller than the minimum value allowed");

        if (step > param.max_step)
          throw std::runtime_error("the line search step became larger than the maximum value allowed");

        // continue search in mid of current search range
        step = std::isinf(step_hi) ? 2*step : step_lo/2 + step_hi/2;
      }

      if (iter >= param.max_linesearch)
        std::cerr << "Line search has reached the maximum number of iterations!" << std::endl;
//        throw std::runtime_error("the line search routine reached the maximum number of iterations");
    }
  };


} // namespace LBFGSpp

#endif // LINE_SEARCH_BRACKETING_H

