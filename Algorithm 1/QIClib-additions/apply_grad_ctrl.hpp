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

namespace qic {

//******************************************************************************

  template <typename T1, typename T2,
            typename TR = typename std::enable_if<is_floating_point_var<trait::pT<T1>, trait::pT<T2> >::value &&
                                                  is_same_pT_var<T1, T2>::value,
                                                  arma::Mat<typename eT_promoter_var<T1, T2>::type> >::type>
  inline TR apply_grad_ctrl(const T1& rho1, const T2& A, arma::uvec ctrl, arma::uvec sys, arma::uvec dim)
  {
    // APR: My own addition to the QIClib library - a function that allows one to apply the deriviated of a controlled
    // gate, given the controls, target and the matrix of the derivative of the uncontrolled gate. This is just an
    // adaptation of the (somewhat icky) code in apply_ctrl().
    // (The temptation to rewrite much of this is strong! Variable names are not friendly, the code could do with some
    // comments, and the repeated requirement to add or subtract one to correct indices from 1, 2, 3 to 0, 1, 2 and vice
    // versa gives one a headache. (There are points where I still haven't convinced myself that these corrections are,
    // in fact, correct. This conversion of indices should really be done in a single location, wrapped in a function.)
    using eTR = typename eT_promoter_var<T1, T2>::type;

    const auto& p = as_Mat(rho1);
    const auto& A1 = as_Mat(A);

    bool checkV = true;
    if (p.n_cols == 1)
      checkV = false;

    arma::uword d = ctrl.n_elem > 0 ? dim.at(ctrl.at(0) - 1) : 1;  // If controls exist, this had better be equal to 2.

#ifndef QICLIB_NO_DEBUG
    if (p.n_elem == 0)
      throw Exception("qic::apply_grad_ctrl", Exception::type::ZERO_SIZE);

    if (A1.n_elem == 0)
      throw Exception("qic::apply_grad_ctrl", Exception::type::ZERO_SIZE);

    if (checkV)
//      if (p.n_rows != p.n_cols)  // I haven't attempted to adapt the density matrix code - we assume checkV == false.
      throw Exception("qic::apply_grad_ctrl", Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    if (A1.n_rows != A1.n_cols)
      throw Exception("qic::apply_grad_ctrl", Exception::type::MATRIX_NOT_SQUARE);

    for (arma::uword i = 1; i < ctrl.n_elem; ++i)
      if (dim.at(ctrl.at(i) - 1) != d)
        throw Exception("qic::apply_grad_ctrl", Exception::type::DIMS_NOT_EQUAL);

    if (d > 2)  // My addition - I don't know how to handle ternaty control bits!
      throw Exception("qic::apply_grad_ctrl", Exception::type::INVALID_DIMS);

    if (dim.n_elem == 0 || arma::any(dim == 0))
      throw Exception("qic::apply_grad_ctrl", Exception::type::INVALID_DIMS);

    if (arma::prod(dim) != p.n_rows)
      throw Exception("qic::apply_grad_ctrl", Exception::type::DIMS_MISMATCH_MATRIX);

    if (arma::prod(dim(sys - 1)) != A1.n_rows)
      throw Exception("qic::apply_grad_ctrl", Exception::type::DIMS_MISMATCH_MATRIX);

    const arma::uvec ctrlsys = arma::join_cols(sys, ctrl);

    if (ctrlsys.n_elem > dim.n_elem || arma::unique(ctrlsys).eval().n_elem != ctrlsys.n_elem ||
        arma::any(ctrlsys > dim.n_elem) || arma::any(ctrlsys == 0))
      throw Exception("qic::apply_grad_ctrl", Exception::type::INVALID_SUBSYS);
#endif

    _internal::dim_collapse_sys_ctrl(dim, sys, ctrl);  // Removes any 'bit' with dim = 1 and combines consecutive
                                                       // 'uninvolved' qbits into a single 'qdit'.
    const arma::uword n = dim.n_elem;
    const arma::uword m = sys.n_elem;
    const arma::uword o = ctrl.n_elem;

    // 'keep' stores the uninvolved qdits.
    arma::uvec keep(n - m);
    arma::uword keep_count(0);
    for (arma::uword run = 0; run < n; ++run) {
      if (!arma::any(sys == run + 1)) {
        keep.at(keep_count) = run + 1;
        ++keep_count;
      }
    }

    // The 'product' array allow for simple calculation of the number/index of the basis state (i.e. 0, 1, 2, 3, etc)
    // from the status of the individual qdits.
    arma::uword product[_internal::MAXQDIT];
    product[n - 1] = 1;
    for (arma::sword i = n - 2; i >= 0; --i)
      product[i] = product[i + 1] * dim.at(i + 1);

    // The 'productr' array does similar, but considers the qdits of the system affected by the gate only. This is used
    // to calculate indices into the transformation matrix.
    arma::uword productr[_internal::MAXQDIT];
    productr[m - 1] = 1;
    for (arma::sword i = m - 2; i >= 0; --i)
      productr[i] = productr[i + 1] * dim.at(sys(i) - 1);

    arma::uword p_num = std::max(static_cast<arma::uword>(1), d - 1);

    // 'Ap' seems only to be necessary if one is using more general quantum digits that qbits. (It appears that, when
    // using qtrits, if the control trits are set to 2, the square of the matrix is applied. This does not seem to be
    // standard.)
    arma::field<arma::Mat<trait::eT<T2> > > Ap(p_num + 1);
    for (arma::uword i = 0; i <= p_num; ++i)
      Ap.at(i) = powm_gen(A1, i);

    if (!checkV)
    {
      arma::Col<eTR> rho(p.n_rows, arma::fill::zeros);  // Where we will put the output state.

      // The first n digits of 'loop_counter' will indicate which basis state we are calculating/adjusting the
      // coefficient of. The next n digits indicate the basis state whose coefficient we are using in calculating this
      // adjustment. The last digit merely indicates when we have finished. (Terrible name, and not the most transparent
      // data structure either.)
      const arma::uword loop_no = 2 * n;
      constexpr auto loop_no_buffer = 2 * _internal::MAXQDIT + 1;
      arma::uword loop_counter[loop_no_buffer] = {0};
      arma::uword MAX[loop_no_buffer];

      for (arma::uword i = 0; i < n; ++i)
      {
        MAX[i] = dim.at(i);
        if (arma::any(keep == i + 1))
          MAX[i + n] = 1;
        else
          MAX[i + n] = dim.at(i);
      }
      MAX[loop_no] = 2;

      arma::uword p1 = 0;

      while (loop_counter[loop_no] == 0)
      {
        arma::uword count1(0), count2(0);

        for (arma::uword i = 0; i < n; ++i)
        {
          count1 += (arma::any(ctrl == i + 1) && loop_counter[i] != 0) ? 1 : 0;  // The number of qbits set in the...
                                                                                 // ...state given by first n digits.
          count2 += loop_counter[i + n] == 0 ? 1 : 0;  // Only used to determine whether we have started setting any...
        }                                              // ...of the next n digits.

        if ((count1 != o) && (count2 == n))  // If not all control qbits are set (and we are still only changing the...
        {                                    // ...first n digits).
          arma::uword I(0);
          for (arma::uword i = 0; i < n; ++i)
            I += product[i] * loop_counter[i];
          rho.at(I) = static_cast<eTR>(0);  // The important change from apply_ctrl. If control bits are not set,
                                            // gradient is ZERO.
        }
        else if (count1 == o)  // All control bits are set.
        {
          arma::uword I(0), J(0), K(0), L(0);
          arma::uword power = o == 0 ? 1 : 0;

          for (arma::uword i = 0; i < n; ++i)
          {
            if (arma::any(keep == i + 1))
            {
              I += product[i] * loop_counter[i];  // 'I' indicates which basis state is having its coefficient changed.
              J += product[i] * loop_counter[i];  // 'J' indicates which basis state's coefficient is influencing...
            }                                     // ...this change.
            else
            {
              I += product[i] * loop_counter[i];
              J += product[i] * loop_counter[i + n];
            }

            // Calculate the power of the matrix to be applied. In my code, which only uses qbits, this will always be
            // 1. Note that this is the bit that is probably broken when using qtrits - the derivative of the square of
            // a matrix is not the square of the derivative. Moreover, it looks like QIClib only applys the
            // transformation when all qtrits have the same value, in which case different values should result in zero
            // gradient - the code below applies the identity matrix instead.
            if (o != 0)
            {
              arma::uword counter_1(1);
              for (arma::uword j = 1; j < o; ++j)
                counter_1 += loop_counter[ctrl.at(0) - 1] == loop_counter[ctrl.at(j) - 1] ? 1 : 0;
              power = counter_1 == o ? loop_counter[ctrl.at(0) - 1] : 0;
            }

            arma::uword counter(0);
            while (any(sys == i + 1))  // YUCK: Use 'if' here, 'while' below to get counter to the correct value.
            {                          //       Then there is no need for 'break'.
              if (sys.at(counter) != i + 1)
              {
                ++counter;
              }
              else
              {
                K += productr[counter] * loop_counter[i];      // 'K' and 'L' gives indices into the matrix.
                L += productr[counter] * loop_counter[i + n];
                break;
              }
            }
          }
          rho.at(I) += Ap.at(power).at(K, L) * p.at(J);
        }

        // The last bit is essentially basisState++, i.e. move to next setting. (000->001->010->011->020->021->100 etc,
        // if dim = {2, 3, 2}.)
        ++loop_counter[0];
        while (loop_counter[p1] == MAX[p1]) {
          loop_counter[p1] = 0;
          loop_counter[++p1]++;
          if (loop_counter[p1] != MAX[p1])
            p1 = 0;
        }
      }
      return rho;

    }
    else
    {
      // I haven't attempted to adapt the density matrix code - we assume checkV == false.
      throw Exception("qic::apply_grad_ctrl", Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);
    }
  }

//******************************************************************************

  template <typename T1, typename T2,
            typename TR = typename std::enable_if<is_floating_point_var<trait::pT<T1>, trait::pT<T2> >::value &&
                                                  is_same_pT_var<T1, T2>::value,
                                                  arma::Mat<typename eT_promoter_var<T1, T2>::type> >::type>
  inline TR apply_grad_ctrl(const T1& rho1, const T2& A, arma::uvec ctrl, arma::uvec sys, arma::uword dim = 2)
  {
    // APR: My own addition to the QIClib library - a function that allows one to apply the deriviated of a controlled
    // gate, given the controls, target and the matrix of the derivative of the uncontrolled gate. This is just an
    // adaptation of the (somewhat icky) code in apply_ctrl().
    const auto& rho = as_Mat(rho1);

#ifndef QICLIB_NO_DEBUG
    bool checkV = true;
    if (rho.n_cols == 1)
      checkV = false;

    if (rho.n_elem == 0)
      throw Exception("qic::apply_grad_ctrl", Exception::type::ZERO_SIZE);

    if (checkV)
      if (rho.n_rows != rho.n_cols)
        throw Exception("qic::apply_grad_ctrl",
                        Exception::type::MATRIX_NOT_SQUARE_OR_CVECTOR);

    if (dim == 0)  // I am not convinced that this code will work correctly with dim != 2. (How does QIClib deal...
      throw Exception("qic::apply_grad_ctrl", Exception::type::INVALID_DIMS);  // ...with ternary control bits??)
#endif

    arma::uword n = static_cast<arma::uword>(std::llround(std::log(rho.n_rows) / std::log(dim)));

    arma::uvec dim2(n);
    dim2.fill(dim);
    return apply_grad_ctrl(rho, A, std::move(ctrl), std::move(sys), std::move(dim2));
  }

//******************************************************************************

}  // namespace qic
