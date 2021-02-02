// MOQCD1 (Multi-objective quantum circuit discovery - algorithm 1.)
//
// Modification to LBFGS++ code copyright (c) 2020  Alan Reynolds (alanreynolds9@gmail.com)
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

// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef LINE_SEARCH_BACKTRACKING_H
#define LINE_SEARCH_BACKTRACKING_H

#include <Eigen/Core>
#include <stdexcept>  // std::runtime_error
#include <iostream>


namespace LBFGSpp
{
  //
  // The backtracking line search algorithm for LBFGS. Mainly for internal use.
  //
  template <typename Scalar>
  class LineSearchBacktracking
  {
  private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;

  public:
    //
    // Line search by backtracking.
    //
    // f      A function object such that `f(x, grad)` returns the objective function value at `x`, and overwrites
    //        `grad` with the gradient.
    // fx     In: The objective function value at the current point. Out: The function value at the new point.
    // x      Out: The new point moved to.
    // grad   In: The current gradient vector. Out: The gradient at the new point.
    // step   In: The initial step length. Out: The calculated step length.
    // drt    The current moving direction.
    // xp     The current point.
    // param  Parameters for the LBFGS algorithm
    //
    template <typename Foo>
    static void LineSearch(Foo& f, Scalar& fx, Vector& x, Vector& grad, Scalar& step, const Vector& drt,
                           const Vector& xp, const LBFGSParam<Scalar>& param)
    {
      // Decreasing and increasing factors
      const Scalar dec = 0.5;
      const Scalar inc = 2.1;

      // Check the value of step
      if (step <= Scalar(0))
        throw std::invalid_argument("'step' must be positive");

      // Save the function value at the current x
      const Scalar fx_init = fx;
      // Projection of gradient on the search direction
      const Scalar dg_init = grad.dot(drt);
      // Make sure d points to a descent direction
      if (dg_init > 0)
        throw std::logic_error("the moving direction increases the objective function value");

      const Scalar dg_test = param.ftol * dg_init;
      Scalar width;

      int iter;
      for (iter = 0; iter < param.max_linesearch; iter++)
      {
        // x_{k+1} = x_k + step * d_k
        x.noalias() = xp + step * drt;
        // Evaluate this candidate
        fx = f(x, grad);

        if (fx > fx_init + step * dg_test)
        {
          width = dec;
        }
        else
        {
          // Armijo condition is met
          if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_ARMIJO)
            break;

          const Scalar dg = grad.dot(drt);
          if (dg < param.wolfe * dg_init)
          {
            width = inc;
          }
          else
          {
            // Regular Wolfe condition is met
            if (param.linesearch == LBFGS_LINESEARCH_BACKTRACKING_WOLFE)
              break;

            if (dg > -param.wolfe * dg_init)
            {
              width = dec;
            }
            else
            {
              // Strong Wolfe condition is met
              break;
            }
          }
        }

//        if (iter >= param.max_linesearch)  // This CANNOT happen here. It can happen AFTER exiting the loop.
//          throw std::runtime_error("the line search routine reached the maximum number of iterations");

        if (step < param.min_step)
          throw std::runtime_error("the line search step became smaller than the minimum value allowed");

        if (step > param.max_step)
          throw std::runtime_error("the line search step became larger than the maximum value allowed");

        step *= width;
      }
      
      if(iter >= param.max_linesearch)
        std::cerr << "Line search has reached the maximum number of iterations!" << std::endl;
//        throw std::runtime_error("the line search routine reached the maximum number of iterations");
    }
  };


} // namespace LBFGSpp

#endif // LINE_SEARCH_BACKTRACKING_H
