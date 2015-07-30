dnl Process this file with autoconf to produce a configure script.

dnl Copyright (C) 2009-2011 Dynare Team
dnl
dnl This file is part of Dynare.
dnl
dnl Dynare is free software: you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl Dynare is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.
dnl
dnl You should have received a copy of the GNU General Public License
dnl along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.61])
AC_INIT([dynare], [4.3.3])
AC_CONFIG_SRCDIR([configure.ac])
AM_INIT_AUTOMAKE([-Wall -Wno-portability -Werror foreign subdir-objects])

AC_CHECK_PROG([MKOCTFILE], [mkoctfile], [mkoctfile])

if test "x$MKOCTFILE" != "x"; then
  CC=`$MKOCTFILE -p CC`
  CXX=`$MKOCTFILE -p CXX`
  CFLAGS=`$MKOCTFILE -p CFLAGS`
  FFLAGS=`$MKOCTFILE -p FFLAGS`
  CXXFLAGS=`$MKOCTFILE -p CXXFLAGS`
fi

CFLAGS="$CFLAGS -Wall -Wno-parentheses"
FFLAGS="$FFLAGS -Wall"
CXXFLAGS="$CXXFLAGS -Wall -Wno-parentheses"

AC_PROG_F77([gfortran g77 f77])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AX_PROG_LN_S
AC_PROG_MKDIR_P
AM_PROG_AR
AX_PTHREAD
# Check for boost libraries used in estimation DLL
AC_LANG_PUSH([C++])
AX_BOOST_BASE([1.36])
CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADER([boost/random/linear_congruential.hpp], [], [AC_MSG_ERROR([Can't find Boost Random Library])])
AC_CHECK_HEADER([boost/math/distributions/beta.hpp], [], [AC_MSG_ERROR([Can't find Boost Math Library])])
CPPFLAGS="$CPPFLAGS_SAVED"
AC_LANG_POP([C++])
# Check for dlopen(), needed by k_order_perturbation DLL
AC_CHECK_LIB([dl], [dlopen], [LIBADD_DLOPEN="-ldl"], [])
AC_SUBST([LIBADD_DLOPEN])
# Check for GSL, needed by MS-SBVAR
AX_GSL
AM_CONDITIONAL([HAVE_GSL], [test "x$has_gsl" = "xyes"])

# Check for libmatio, needed by estimation DLL (only under Octave, as an alternative to MATLAB's libmat)
AX_MATIO
AM_CONDITIONAL([HAVE_MATIO], [test "x$has_matio" = "xyes"])

# Check for libslicot, needed by kalman_steady_state
AX_SLICOT([octave])
AM_CONDITIONAL([HAVE_SLICOT], [test "x$has_slicot" = "xyes"])

AM_CONDITIONAL([DO_SOMETHING], [test "x$MKOCTFILE" != "x"])

if test "x$MKOCTFILE" != "x"; then
  if test x"$ax_pthread_ok" = "xyes"; then
     BUILD_MEX_OCTAVE="yes"
  else
     BUILD_MEX_OCTAVE="yes (without POSIX threads)"
  fi
else
  BUILD_MEX_OCTAVE="no (missing mkoctfile)"
fi

if  test "x$MKOCTFILE" != "x" -a "x$has_slicot" = "xyes"; then
   BUILD_KALMAN_STEADY_STATE_OCTAVE="yes"
else
   BUILD_KALMAN_STEADY_STATE_OCTAVE="no (missing SLICOT)"
fi

if test "x$MKOCTFILE" != "x" -a "x$has_gsl" = "xyes" -a "x$has_matio" = "xyes"; then
   BUILD_MS_SBVAR_MEX_OCTAVE="yes"
else
   BUILD_MS_SBVAR_MEX_OCTAVE="no (missing GSL or MatIO library)"
fi

AC_ARG_ENABLE([openmp], AS_HELP_STRING([--enable-openmp], [use OpenMP for parallelization of some MEX files]), [
  if test "x$enable_openmp" = "xyes"; then
    CPPFLAGS="$CPPFLAGS -DUSE_OMP"
    CFLAGS="$CFLAGS -fopenmp"
    CXXFLAGS="$CXXFLAGS -fopenmp"
  fi
])

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 MEX files for Octave (except MS-SBVAR and Kalman Steady State):     $BUILD_MEX_OCTAVE
 MS-SBVAR MEX files for Octave:                                      $BUILD_MS_SBVAR_MEX_OCTAVE
 Kalman Steady State MEX file for Octave:                            $BUILD_KALMAN_STEADY_STATE_OCTAVE

])

AC_CONFIG_FILES([Makefile
                 mjdgges/Makefile
                 kronecker/Makefile
                 bytecode/Makefile
                 libdynare++/Makefile
                 gensylv/Makefile
                 k_order_perturbation/Makefile
                 dynare_simul_/Makefile
                 qzcomplex/Makefile
                 ordschur/Makefile
                 kalman_steady_state/Makefile
                 ms_sbvar/Makefile
                 block_kalman_filter/Makefile
		 sobol/Makefile
		 local_state_space_iterations/Makefile])

AC_OUTPUT
