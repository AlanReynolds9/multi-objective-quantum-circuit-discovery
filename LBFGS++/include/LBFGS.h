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

#ifndef LBFGS_H
#define LBFGS_H

#include <cmath>  // Added to get isnan and isinf working on linux. (Might not be needed?)
#include <Eigen/Core>
#include "LBFGS/Param.h"
#include "LBFGS/LineSearchBacktracking.h"
#include "LBFGS/LineSearchBracketing.h"
#include "LBFGS/LineSearchNocedalWright.h"


namespace LBFGSpp
{
  ///
  /// LBFGS solver for unconstrained numerical optimization
  ///
  template < typename Scalar,
  template<class> class LineSearch = LineSearchBracketing >  // LineSearchBacktracking >
  class LBFGSSolver
  {
  private:
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    typedef Eigen::Map<Vector> MapVec;

    LBFGSParam<Scalar>& m_param;  // Parameters to control the LBFGS algorithm
    Matrix              m_s;      // History of the s vectors
    Matrix              m_y;      // History of the y vectors
    Vector              m_ys;     // History of the s'y values
    Vector              m_alpha;  // History of the step lengths
    Vector              m_fx;     // History of the objective function values
    Vector              m_xp;     // Old x
    Vector              m_grad;   // New gradient
    Vector              m_gradp;  // Old gradient
    Vector              m_drt;    // Moving direction

    inline void reset(int n)
    {
      const int m = m_param.m;
      m_s.resize(n, m);
      m_y.resize(n, m);
      m_ys.resize(m);
      m_alpha.resize(m);
      m_xp.resize(n);
      m_grad.resize(n);
      m_gradp.resize(n);
      m_drt.resize(n);
      if(m_param.past > 0)
        m_fx.resize(m_param.past);
    }

  public:
    ///
    /// Constructor for LBFGS solver.
    ///
    /// \param param An object of \ref LBFGSParam to store parameters for the
    ///        algorithm
    ///
    LBFGSSolver(LBFGSParam<Scalar>& param) :
    m_param(param)
    {
      m_param.check_param();
    }

    ///
    /// Minimizing a multivariate function using LBFGS algorithm.
    /// Exceptions will be thrown if error occurs.
    ///
    /// \param f  A function object such that `f(x, grad)` returns the
    ///           objective function value at `x`, and overwrites `grad` with
    ///           the gradient.
    /// \param x  In: An initial guess of the optimal point. Out: The best point
    ///           found.
    /// \param fx Out: The objective function value at `x`.
    ///
    /// \return Number of iterations used.
    ///
    template <typename Foo>
    inline int minimize(Foo& f, Vector& x, Scalar& fx)
    {
      const int n = x.size();
      const int fpast = m_param.past;
      reset(n);
      int numStored = 0;  // Number of y, s and ys stored.

      // Evaluate function and compute gradient
      fx = f(x, m_grad);
//      std::cout << fx;   //********
      Scalar xnorm = x.norm();
      Scalar gnorm = m_grad.norm();
      if (fpast > 0)
        m_fx[0] = fx;

      // Early exit if the initial x is already a minimizer
      if (gnorm <= m_param.epsilon)// * std::max(xnorm, Scalar(1.0)))  // I don't need this scaling. (If I did, I'd...
      {                                                                // ...do it differently anyway.)
//        std::cout << std::endl;  //********
        return 0;  // This was 1, which didn't make sense to me.
      }

      // Initial direction
      m_drt.noalias() = -m_grad;
      // Initial step
      Scalar step = Scalar(1.0) / m_drt.norm();
      if (std::isnan(step) || std::isinf(step))  // Added std:: to get these working on linux.
      {
        std::cerr << "step is either infinite or not a number." << std::endl;
        std::cerr << "m_drt.norm() = " << m_drt.norm() << std::endl;
        std::cerr << "step = " << step << std::endl;
        exit(EXIT_FAILURE);
      }

      int k = 1;
      int end = 0;
      for ( ; ; )
      {
        // Save the curent x and gradient
        m_xp.noalias() = x;
        m_gradp.noalias() = m_grad;

        // Line search to update x, fx and gradient
        LineSearch<Scalar>::LineSearch(f, fx, x, m_grad, step, m_drt, m_xp, m_param);
//        std::cout << ", " << fx;  //********

        // New x norm and gradient norm
        xnorm = x.norm();
        gnorm = m_grad.norm();

        // Convergence test -- gradient
        if (gnorm <= m_param.epsilon)// * std::max(xnorm, Scalar(1.0)))  // I don't need this scaling, and would do...
        {                                                                // ...it differently if I did.
//          std::cout << std::endl;  //********
          return k;
        }
        // Convergence test -- absolute objective function improvement
        if (m_param.absoluteDeltaTest)
        {
          if (k >= fpast && m_fx[k % fpast] - fx <= m_param.absoluteDelta)
          {
//            std::cout << std::endl;  //********
            return k;
          }

          m_fx[k % fpast] = fx;
        }
        // Convergence test -- relative objective function improvement
        if (m_param.relativeDeltaTest)
        {
          if (k >= fpast && m_fx[k % fpast] - fx <= m_param.relativeDelta * std::abs(fx))
          {
//            std::cout << std::endl;  //********
            return k;
          }

          m_fx[k % fpast] = fx;
        }
        // Halting test -- is search likely to meet the target within the iteration quota. (If target is already met,
        // just continue!)
        if (m_param.targetTest)
        {
          if (k >= fpast &&
              fx > m_param.target && m_param.iterationQuota * (m_fx[k % fpast] - fx) / fpast < fx - m_param.target)
          {
//            std::cout << std::endl;  //********
            return k;
          }

          m_fx[k % fpast] = fx;
        }
        // Convergence test -- reached a target value considered 'good enough'
        if (m_param.goodEnoughTest)
        {
          if (fx < m_param.goodEnough)
          {
//            std::cout << std::endl;  //********
            return k;
          }
        }
        // Maximum number of iterations
        if (m_param.max_iterations != 0 && k >= m_param.max_iterations)
        {
//          std::cout << std::endl;  //********
          return k;
        }

        // Update s and y
        // s_{k+1} = x_{k+1} - x_k
        // y_{k+1} = g_{k+1} - g_k
        Vector svec = x - m_xp;
        Vector yvec = m_grad - m_gradp;

        // ys = y's = 1/rho
        // yy = y'y
        Scalar ys = yvec.dot(svec);
        Scalar yy = yvec.squaredNorm();

        // If ys > 0 (which should be the case unless the line search failed), update the stored versions in m_s, m_y
        // and m_ys
        if (ys > 0)
        {
          MapVec(&m_s(0, end), n).noalias() = svec;
          MapVec(&m_y(0, end), n).noalias() = yvec;
          m_ys[end] = ys;
          end = (++end) % m_param.m;
          ++numStored;
        }
        else
        {
          std::cerr << "Warning: ys <= 0 (" << ys << "). In this step we do not store y, s or ys." << std::endl;
          std::cerr << "svec = " << svec << std::endl;
          std::cerr << "yvec = " << yvec << std::endl;

          if (numStored == 0)
          {
            std::cerr << "...and this is happening before any values for y, s and ys have been stored." << std::endl;
            ys = yy = 1;
          }
          else
          {
            int old = (end + m_param.m - 1) % m_param.m;
            MapVec oldY(&m_y(0, old), n);
            ys = m_ys[old];
            yy = oldY.dot(oldY);
          }
        }

        // Recursive formula to compute d = -H * g
        m_drt.noalias() = -m_grad;
        int bound = std::min(m_param.m, numStored);
        int j = end;
        for(int i = 0; i < bound; i++)
        {
          j = (j + m_param.m - 1) % m_param.m;
          MapVec sj(&m_s(0, j), n);
          MapVec yj(&m_y(0, j), n);
          m_alpha[j] = sj.dot(m_drt) / m_ys[j];
          if (std::isnan(m_alpha[j]) || std::isinf(m_alpha[j]))  // Added std:: to get these working on linux.
          {
            std::cerr << "m_alpha[j] is either infinite or not a number." << std::endl;
            std::cerr << "sj.dot(m_drt) = " << sj.dot(m_drt) << std::endl;
            std::cerr << "m_ys[j] = " << m_ys[j] << std::endl;
            std::cerr << "m_alpha[j] = " << m_alpha[j] << std::endl;
            std::cerr << "gnorm = " << gnorm << std::endl;
            std::cerr << "numStored = " << numStored << std::endl;
            std::cerr << "j = " << j << std::endl;
            exit(EXIT_FAILURE);
          }
          m_drt.noalias() -= m_alpha[j] * yj;
        }

        if (std::isnan(ys / yy) || std::isinf(ys / yy))  // Added std:: to get these working on linux.
        {
          std::cerr << "ys / yy is either infinite or not a number." << std::endl;
          std::cerr << "ys = " << ys << std::endl;
          std::cerr << "yy = " << yy << std::endl;
          exit(EXIT_FAILURE);
        }
        m_drt *= (ys / yy);

        for (int i = 0; i < bound; i++)
        {
          MapVec sj(&m_s(0, j), n);
          MapVec yj(&m_y(0, j), n);
          Scalar beta = yj.dot(m_drt) / m_ys[j];
          if (std::isnan(beta) || std::isinf(beta))
          {
            std::cerr << "beta is either infinite or not a number." << std::endl;
            std::cerr << "yj.dot(m_drt) = " << yj.dot(m_drt) << std::endl;
            std::cerr << "m_ys[j] = " << m_ys[j] << std::endl;
            exit(EXIT_FAILURE);
          }
          m_drt.noalias() += (m_alpha[j] - beta) * sj;
          j = (j + 1) % m_param.m;
        }

        // step = 1.0 as initial guess
        step = Scalar(1.0);
        k++;
      }

//      std::cout << std::endl;  //********
      return k;
    }
  };


} // namespace LBFGSpp

#endif // LBFGS_H
