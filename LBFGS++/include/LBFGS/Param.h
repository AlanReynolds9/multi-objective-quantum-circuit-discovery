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

// Copyright (C) 2016-2019 Yixuan Qiu <yixuan.qiu@cos.name>
// Under MIT license

#ifndef PARAM_H
#define PARAM_H

#include <Eigen/Core>
#include <stdexcept>  // std::invalid_argument


namespace LBFGSpp
{
  enum LINE_SEARCH_ALGORITHM
  {
    // Backtracking method with the Armijo condition, which attempts to find a step length that satisfies the sufficient
    // decrease (Armijo) condition, f(x + a d) <= f(x) + ftol a g(x).d, where x is the current point, d is the current
    // search direction, a is the step length, and ftol is specified by the user. f and g are the function and gradient
    // values respectively.
    LBFGS_LINESEARCH_BACKTRACKING_ARMIJO = 1,

    // The backtracking method with the regular Wolfe condition. An alias of `LBFGS_LINESEARCH_BACKTRACKING_WOLFE`.
    LBFGS_LINESEARCH_BACKTRACKING = 2,

    // Backtracking method with regular Wolfe condition, which attempts to find a step length that satisfies both the
    // Armijo condition (see above) and the curvature condition, g(x + a d).d >= wolfe g(x).d, where wolfe is a value
    // specified by the user.
    LBFGS_LINESEARCH_BACKTRACKING_WOLFE = 2,

    // Backtracking method with strong Wolfe condition. Attempts to find a step length tha satisfies both the Armijo
    // condition (see above) and the strong curvature condition, |g(x + a d).d| <= |wolfe g(x).d|, where wolfe is a
    // value specified by the user.
    LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE = 3
  };


  //
  // Parameters to control the LBFGS algorithm.
  //
  template <typename Scalar = double>
  class LBFGSParam
  {
  public:
    // The number of corrections to approximate the inverse hessian matrix. The L-BFGS routine stores the computation
    // results of previous m iterations to approximate the inverse hessian matrix of the current iteration. This
    // parameter controls the size of the limited memories (corrections). The default value is 6. Values less than 3 are
    // not recommended. Large values will result in excessive computing time.
    int m;

    // Tolerance for the gradient based convergence test. A minimization terminates when ||g|| < \epsilon, where ||.||
    // denotes the Euclidean (L2) norm. The default value is 1e-5.
    Scalar epsilon;

    // Number of past iterations, d, kept in store for both delta-based and target based convergence tests. The default
    // value is 0, indicating that none of these convergence tests are applied.
    int    past;

    // Convergence test based on absolute improvement in objective function. The algorithm stops when the following
    // condition is met: f(x_{k-d})-f(x_k) < delta, where k is the current iteration and x_i is the position evaluated
    // in iteration i. Default values are 'false' and zero.
    bool absoluteDeltaTest;
    Scalar absoluteDelta;

    // Convergence test based on relative improvement in objective function. This is the same as the above, but the
    // improvement in objective is compared against delta * f(x_k).
    bool relativeDeltaTest;
    Scalar relativeDelta;

    // Halting condition based on a target and the number of iterations expected to reach it at the current rate. The
    // algorithm stops if quota * (f(x_{k-d})-f(x_k)) / d < f(x_k) - target, i.e. if the current rate of solution
    // improvement is insufficient to reach the target within the quota of iterations allowed. (This does not test for
    // convergence, but rather halts unpromising searches.)
    bool targetTest;
    Scalar target;
    int iterationQuota;

    // Halting condition based on a objective value that is considered 'good enough'.
    bool goodEnoughTest;
    Scalar goodEnough;

    // The maximum number of iterations. Zero (the default value) indicates that there is no maximum number and the
    // algorithm will continue until stopped by one of the convergence tests.
    int    max_iterations;

    // The line search algorithm used by the LBFGS routine. The default value is `LBFGS_LINESEARCH_BACKTRACKING_ARMIJO`.
    int    linesearch;

    // The maximum number of trials for the line search. The default value is 20.
    int    max_linesearch;

    // The minimum step length allowed in the line search. The default value is 1e-20. Usually this value does not need
    // to be modified.
    Scalar min_step;

    // The maximum step length allowed in the line search. The default value is \c 1e+20. Usually this value does not
    // need to be modified.
    Scalar max_step;

    // A parameter that determines what the line search routine considers to be a sufficient decrease to satisfy the
    // Armijo condition. The default value is 1e-4. This parameter should be greaterthan zero and smaller than 0.5.
    Scalar ftol;

    // The coefficient for the Wolfe condition, if used. The default value is 0.9. This parameter should be greater
    // the ftol parameter and smaller than 1.0.
    Scalar wolfe;

  public:
    ///
    /// Constructor for LBFGS parameters.
    /// Default values for parameters will be set when the object is created.
    ///
    LBFGSParam()
    {
      m                 = 6;
      epsilon           = Scalar(1e-5);
      past              = 0;
      absoluteDeltaTest = false;
      absoluteDelta     = Scalar(0);
      relativeDeltaTest = false;
      relativeDelta     = Scalar(0);
      targetTest        = false;
      target            = Scalar(0);
      goodEnoughTest    = false;
      goodEnough        = Scalar(0);
      iterationQuota    = 0;
      max_iterations    = 0;
      linesearch        = LBFGS_LINESEARCH_BACKTRACKING_ARMIJO;
      max_linesearch    = 20;
      min_step          = Scalar(1e-20);
      max_step          = Scalar(1e+20);
      ftol              = Scalar(1e-4);
      wolfe             = Scalar(0.9);
    }

    ///
    /// Checking the validity of LBFGS parameters.
    /// An `std::invalid_argument` exception will be thrown if some parameter
    /// is invalid.
    ///
    inline void check_param()
    {
      if (m <= 0)
        throw std::invalid_argument("'m' must be positive.");
      if (epsilon <= 0)
        throw std::invalid_argument("'epsilon' must be positive.");

      if (absoluteDeltaTest && past <= 0)
        throw std::invalid_argument("For the absolute 'delta' test, 'past' must be positive.");
      if (relativeDeltaTest && past <= 0)
        throw std::invalid_argument("For the relative 'delta' test, 'past' must be positive.");
      if (targetTest && past <= 0)
        throw std::invalid_argument("For the 'target' test, 'past' must be positive.");

      if (absoluteDeltaTest && absoluteDelta <= 0)
        throw std::invalid_argument("For the absolute 'delta' test, 'absoluteDelta' must be positive");
      if (relativeDeltaTest && relativeDelta <= 0)
        throw std::invalid_argument("For the relative 'delta' test, 'relativeDelta' must be positive");
      if (targetTest && iterationQuota <= 0)
        throw std::invalid_argument("For the target test, 'iterationQuota must be positive");

      if (!(absoluteDeltaTest || relativeDeltaTest || targetTest))
        past = 0;  // No convergence test in use requires the storage of past objective values.

      if (max_iterations < 0)
        throw std::invalid_argument("'max_iterations' must be non-negative");

      if (linesearch < LBFGS_LINESEARCH_BACKTRACKING_ARMIJO || linesearch > LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE)
        throw std::invalid_argument("unsupported line search algorithm");
      if (max_linesearch <= 0)
        throw std::invalid_argument("'max_linesearch' must be positive");
      if (min_step < 0)
        throw std::invalid_argument("'min_step' must be positive");
      if (max_step < min_step )
        throw std::invalid_argument("'max_step' must be greater than 'min_step'");
      if (ftol <= 0 || ftol >= 0.5)
        throw std::invalid_argument("'ftol' must satisfy 0 < ftol < 0.5");
      if (wolfe <= ftol || wolfe >= 1)
        throw std::invalid_argument("'wolfe' must satisfy ftol < wolfe < 1");
    }
  };


} // namespace LBFGSpp

#endif // PARAM_H
