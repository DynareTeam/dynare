/*
 * Oct-file for bringing the complex QZ decomposition to Octave.
 * Simple wrapper around LAPACK's zgges.
 *
 * Written by Sébastien Villemot <sebastien@dynare.org>.
 */

/*
 * Copyright (C) 2010-2016 Dynare Team
 *
 * This file is part of Dynare.
 *
 * Dynare is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Dynare is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __MINGW32__
#define __CROSS_COMPILATION__
#endif

#ifdef __MINGW64__
#define __CROSS_COMPILATION__
#endif

#ifdef __CROSS_COMPILATION__
#define M_PI 3.14159265358979323846
#endif

#include <octave/oct.h>
#include <octave/f77-fcn.h>

extern "C"
{
  F77_RET_T
  F77_FUNC(zgges, ZGGES) (F77_CONST_CHAR_ARG_DECL, F77_CONST_CHAR_ARG_DECL, F77_CONST_CHAR_ARG_DECL,
                          octave_idx_type (*)(Complex *, Complex *), const octave_idx_type &,
                          Complex *, const octave_idx_type &, Complex *, const octave_idx_type &,
                          octave_idx_type &, Complex *, Complex *, Complex *, const octave_idx_type &,
                          Complex *, const octave_idx_type &, Complex *, const octave_idx_type &,
                          double *, octave_idx_type *, octave_idx_type &);
}

DEFUN_DLD(qzcomplex, args, nargout, "-*- texinfo -*-\n\
@deftypefn {Loadable Function} [ @var{aa}, @var{bb}, @var{q}, @var{z} ] = qzcomplex (@var{a}, @var{b})\n\
\n\
Computes the complex QZ decomposition of @math{(A, B)}, satisfying:\n\
@example\n\
    AA = Q'*A*Z, BB = Q'*B*Z\n\
@end example\n\
@end deftypefn\n\
")
{
  int nargin = args.length();
  octave_value_list retval;

  if (nargin != 2 || nargout != 4)
    {
      print_usage();
      return retval;
    }

  ComplexMatrix A = args(0).complex_matrix_value();
  ComplexMatrix B = args(1).complex_matrix_value();
  if (error_state)
    return retval;

  dim_vector dimA = A.dims();
  dim_vector dimB = B.dims();
  octave_idx_type n = dimA(0);
  if (n != dimA(1) || n != dimB(0) || n != dimB(1))
    {
      error("qzcomplex: input matrices must be square and of same size");
      return retval;
    }

  octave_idx_type lwork = 2*n;
  OCTAVE_LOCAL_BUFFER(Complex, alpha, n);
  OCTAVE_LOCAL_BUFFER(Complex, beta, n);
  OCTAVE_LOCAL_BUFFER(Complex, work, lwork);
  OCTAVE_LOCAL_BUFFER(double, rwork, 8*n);
  ComplexMatrix vsl(n, n), vsr(n, n);
  octave_idx_type sdim, info;

  F77_XFCN(zgges, ZGGES, (F77_CONST_CHAR_ARG("V"), F77_CONST_CHAR_ARG("V"),
                          F77_CONST_CHAR_ARG("N"), NULL,
                          n, A.fortran_vec(), n, B.fortran_vec(), n, sdim,
                          alpha, beta, vsl.fortran_vec(), n, vsr.fortran_vec(), n,
                          work, lwork, rwork, NULL, info));

  if (info != 0)
    {
      error("qzcomplex: zgges failed");
      return retval;
    }

  retval(0) = octave_value(A);
  retval(1) = octave_value(B);
  retval(2) = octave_value(vsl);
  retval(3) = octave_value(vsr);
  return retval;
}

/*

  %!test
  %! A = [ 1 2 3+1i; 4 5-1i 0; -1 -5 3];
  %! B = [ -2 -8i 4; 1 5+3i 5; 7 -10 -2];
  %! [AA,BB,Q,Z] = qzcomplex(A,B);
  %! assert(Q'*A*Z, AA, sqrt(eps))
  %! assert(Q'*B*Z, BB, sqrt(eps))

*/
