dnl Process this file with autoconf to produce a configure script.

dnl Copyright (C) 2009-2012 Dynare Team
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

AC_PREREQ([2.62])
AC_INIT([dynare], [4.3.3])
AC_CONFIG_SRCDIR([preprocessor/DynareMain.cc])
AM_INIT_AUTOMAKE([1.11 -Wall -Wno-portability -Werror foreign no-dist-gzip dist-xz tar-pax subdir-objects])

AC_PROG_CC
AC_PROG_CXX

AC_CANONICAL_HOST
case ${host_os} in
  *mingw32*)
    # On mingw32, we don't want dynamic libgcc
    # Note that static-libstdc++ is only supported since GCC 4.5 (but generates no error on older versions)
    LDFLAGS="$LDFLAGS -static-libgcc -static-libstdc++ -static-libgfortran"
    have_windows="yes"
    ;;
  *cygwin*)
    AC_MSG_WARN([You are compiling for the Cygwin target. This means that the preprocessor will])
    AC_MSG_WARN([not run from MATLAB unless you add the Cygwin DLL to the path.])
    AC_MSG_WARN([This is probably not what you want. Consider using a MinGW cross-compiler.])
    if test "x$F77" = "x"; then
      # On Cygwin 1.7, g77 comes has version 3, and is not compatible with default gcc/g++ which has version 4
      # And by default, the AC_PROG_F77 will pick up g77 if it is present (even if gfortran is also here)
      F77=gfortran
    fi
    have_windows="yes"
    ;;
esac

# Use C++ for testing headers
AC_LANG([C++])

CXXFLAGS="$CXXFLAGS -Wall -Wno-parentheses"

# If default 'ar' is not available, try to find one with a host prefix (see ticket #145)
AC_CHECK_PROGS([AR], [ar ${host_alias}-ar])

AC_PROG_RANLIB
AM_PROG_AR

AX_PROG_LN_S

AC_PROG_MKDIR_P

AM_PROG_LEX

AC_CHECK_PROG([YACC], [bison], [bison])
if test "x$YACC" = "x"; then
  unset YACC # AM_MISSING_PROG needs an unset variable: an empty variable won't do
  AM_MISSING_PROG([YACC], [bison])
fi

# We need 1.36 because of unordered_{set,hash} used by Dynare++
AX_BOOST_BASE([1.36])

CPPFLAGS_SAVED="$CPPFLAGS"
CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
AC_CHECK_HEADERS([boost/graph/adjacency_list.hpp], [], [AC_MSG_ERROR([Can't find Boost Graph Library])])
AC_CHECK_HEADERS([boost/algorithm/string/trim.hpp], [], [AC_MSG_ERROR([Can't find Boost String Library])])
AC_CHECK_HEADERS([boost/algorithm/string/split.hpp], [], [AC_MSG_ERROR([Can't find Boost String Library])])
AC_CHECK_HEADERS([boost/lexical_cast.hpp], [], [AC_MSG_ERROR([Can't find Boost Lexical Cast Header])])
CPPFLAGS="$CPPFLAGS_SAVED"

# Don't use deprecated hash structures
AC_DEFINE([BOOST_NO_HASH], [], [Don't use deprecated STL hash structures])

AC_CHECK_PROG([MAKEINFO], [makeinfo], [makeinfo])

AC_CHECK_PROG([PDFTEX], [pdftex], [pdftex])
AM_CONDITIONAL([HAVE_PDFTEX], [test "x$PDFTEX" != "x"])

AC_CHECK_PROG([PDFETEX], [pdfetex], [pdfetex])
AM_CONDITIONAL([HAVE_PDFETEX], [test "x$PDFETEX" != "x"])

if test "x$PDFTEX" != "x"; then
  AC_MSG_CHECKING([for eplain])
  AX_TEX_TEST([\input eplain
\end
], [ax_tex_have_eplain])
  AC_MSG_RESULT([$ax_tex_have_eplain])
fi
AM_CONDITIONAL([HAVE_EPLAIN], [test "x$ax_tex_have_eplain" = "xyes"])

AC_CHECK_PROG([PDFLATEX], [pdflatex], [pdflatex])
AM_CONDITIONAL([HAVE_PDFLATEX], [test "x$PDFLATEX" != "x"])

AC_CHECK_PROG([BIBTEX], [bibtex], [bibtex])
AM_CONDITIONAL([HAVE_BIBTEX], [test "x$BIBTEX" != "x"])

AC_CHECK_PROG([TEXI2HTML], [texi2html], [texi2html])
AM_CONDITIONAL([HAVE_TEXI2HTML], [test "x$TEXI2HTML" != "x"])

AC_CHECK_PROG([LATEX2HTML], [latex2html], [latex2html])
AM_CONDITIONAL([HAVE_LATEX2HTML], [test "x$LATEX2HTML" != "x"])

if test "x$PDFLATEX" != "x" -a "x$BIBTEX" != "x"; then
  AC_MSG_CHECKING([for econometrica BST])
  AX_LATEX_BIBTEX_TEST([\documentclass{article}
\begin{document}
\cite{Juillard1996}
\bibliography{conftest}
\bibliographystyle{econometrica}
\end{document}
], [@techreport{Juillard1996,
	Author = {Juillard, Michel},
	Institution = {CEPREMAP},
	Number = {9602},
	Title = {Dynare : a program for the resolution and simulation of dynamic models with forward variables through the use of a relaxation algorithm},
	Type = {CEPREMAP working papers},
	Year = {1996}}
], [ax_latex_have_econometrica])
  AC_MSG_RESULT([$ax_latex_have_econometrica])

  AX_LATEX_CLASS([beamer], [ax_latex_have_beamer])
fi
AM_CONDITIONAL([HAVE_ECONOMETRICA], [test "x$ax_latex_have_econometrica" = "xyes"])
AM_CONDITIONAL([HAVE_BEAMER], [test "x$ax_latex_have_beamer" = "xyes"])

AC_CHECK_PROG([DOXYGEN], [doxygen], [doxygen])
AM_CONDITIONAL([HAVE_DOXYGEN], [test "x$DOXYGEN" != "x"])

AC_CHECK_PROG([OCTAVE], [octave], [octave])
AM_CONDITIONAL([HAVE_OCTAVE], [test "x$OCTAVE" != "x"])

AC_CHECK_PROG([CTANGLE], [ctangle], [ctangle])
AM_CONDITIONAL([HAVE_CTANGLE], [test "x$CTANGLE" != "x"])
if test "x$CTANGLE" = "x"; then
  unset CTANGLE # AM_MISSING_PROG needs an unset variable: an empty variable won't do
  AM_MISSING_PROG([CTANGLE], [ctangle])
fi

AC_CHECK_PROG([CWEAVE], [cweave], [cweave])
AM_CONDITIONAL([HAVE_CWEAVE], [test "x$CWEAVE" != "x"])

AC_PROG_F77
AC_F77_LIBRARY_LDFLAGS
case ${host_os} in
  *mingw32*)
    # Hack to get static linking of libgfortran on MinGW
    # (-static-libgfortran doesn't act on gcc/g++)
    FLIBS=`echo $FLIBS | sed 's/-lgfortran/-Wl,-Bstatic -lgfortran -Wl,-Bdynamic/'`
    ;;
  *darwin*)
    # * OS X doesn't need to link against the Fortran libraries because
    #   everything is contained within the vecLib framework.
    FLIBS=''
    ;;
esac

if test "x$F77" != "x"; then
   AX_BLAS
   AX_LAPACK
fi
AM_CONDITIONAL([HAVE_BLAS], [test x"$ax_blas_ok" = "xyes"])
AM_CONDITIONAL([HAVE_LAPACK], [test x"$ax_lapack_ok" = "xyes"])

AX_PTHREAD

AC_CONFIG_FILES([Makefile
                 preprocessor/macro/Makefile
                 preprocessor/Makefile
                 doc/Makefile
                 doc/preprocessor/Makefile
                 doc/macroprocessor/Makefile
                 doc/userguide/Makefile
                 doc/parallel/Makefile
                 doc/internals/Makefile
                 doc/gsa/Makefile
                 tests/Makefile
                 matlab/dynare_version.m
                 windows/dynare-version.nsi
                 dynare++/Makefile
                 dynare++/parser/cc/Makefile
                 dynare++/sylv/Makefile
                 dynare++/sylv/cc/Makefile
                 dynare++/sylv/testing/Makefile
                 dynare++/tl/Makefile
                 dynare++/tl/cc/Makefile
                 dynare++/tl/testing/Makefile
                 dynare++/doc/Makefile
                 dynare++/utils/cc/Makefile
                 dynare++/integ/Makefile
                 dynare++/integ/cc/Makefile
                 dynare++/integ/src/Makefile
                 dynare++/integ/testing/Makefile
                 dynare++/kord/Makefile
                 dynare++/src/Makefile
                 mex/sources/Makefile
])

AC_ARG_ENABLE([matlab], AS_HELP_STRING([--disable-matlab], [disable compilation of MEX files for MATLAB]), [], [enable_matlab=yes])
if test "x$enable_matlab" = "xyes"; then
  AC_CONFIG_SUBDIRS([mex/build/matlab])
  AX_MATLAB
fi
AM_CONDITIONAL([ENABLE_MATLAB], [test "x$enable_matlab" = "xyes"])
AM_CONDITIONAL([HAVE_CMD_LINE_MATLAB], [test "x$ax_enable_matlab" = "xyes" -a "x$have_windows" = "x"])

AC_ARG_ENABLE([octave], AS_HELP_STRING([--disable-octave], [disable compilation of MEX files for Octave]), [], [enable_octave=yes])
if test "x$enable_octave" = "xyes"; then
  AC_CONFIG_SUBDIRS([mex/build/octave])
fi
AM_CONDITIONAL([ENABLE_OCTAVE], [test "x$enable_octave" = "xyes"])

# Enable exporting of Org files
# The clean way would be to test for Emacs, Org-mode, latex, dvipng...
AC_ARG_ENABLE([org-export], AS_HELP_STRING([--enable-org-export], [enable exporting of Org files (requires Emacs, org-mode and other external programs)]))
AM_CONDITIONAL([ENABLE_ORG_EXPORT], [test "x$enable_org_export" != "x"])

# Construct final output message

BUILD_PREPROCESSOR="yes"

if test "x$ax_blas_ok" = "xyes" -a "x$ax_lapack_ok" = "xyes"; then
  if test x"$ax_pthread_ok" = "xyes"; then
     BUILD_DYNAREPLUSPLUS="yes"
  else
     BUILD_DYNAREPLUSPLUS="yes (without POSIX threads)"
  fi
else
  BUILD_DYNAREPLUSPLUS="no (missing one of: BLAS, LAPACK)"
fi

if test "x$CWEAVE" != "x" -a x"$PDFTEX" != "x" -a "x$ax_tex_have_eplain" = "xyes"; then
  BUILD_DYNAREPLUSPLUS_SRCDOC="yes"
else
  BUILD_DYNAREPLUSPLUS_SRCDOC="no (missing one of: cweave, pdftex, eplain)"
fi

if test "x$PDFLATEX" != "x" -a "x$BIBTEX" != "x" -a "x$ax_latex_have_econometrica" = "xyes"; then
  BUILD_DYNARE_USERGUIDE="yes"
else
  BUILD_DYNARE_USERGUIDE="no (missing one of: pdflatex, bibtex, econometrica BST)"
fi

if test "x$PDFLATEX" != "x" -a "x$ax_latex_have_beamer" = "xyes"; then
  BUILD_BEAMER_DOC="yes"
else
  BUILD_BEAMER_DOC="no (missing one of: pdflatex, beamer)"
fi

if test "x$PDFLATEX" != "x"; then
  BUILD_OTHER_PDF_DOC="yes"
else
  BUILD_OTHER_PDF_DOC="no (missing pdflatex)"
fi

if test "x$DOXYGEN" != "x"; then
  BUILD_DYNARE_PREPROC_DOC="yes"
else
  BUILD_DYNARE_PREPROC_DOC="no (missing doxygen)"
fi

if test "x$enable_org_export" != "x"; then
  BUILD_DYNARE_INTERNAL_DOC="yes"
else
  BUILD_DYNARE_INTERNAL_DOC="no (Org export not enabled)"
fi

if test "x$MAKEINFO" != "x"; then
  BUILD_DYNARE_INFO="yes"
  if test "x$TEXI2HTML" != "x" -a "x$LATEX2HTML" != "x"; then
    BUILD_DYNARE_HTML_MANUAL="yes"
  else
    BUILD_DYNARE_HTML_MANUAL="yes (but with ugly math formulas, missing texi2html or latex2html)"
  fi
  BUILD_DYNARE_PDF_MANUAL="yes"
else
  BUILD_DYNARE_INFO="no (missing makeinfo)"
  BUILD_DYNARE_HTML_MANUAL="no (missing makeinfo)"
  BUILD_DYNARE_PDF_MANUAL="no (missing makeinfo)"
fi

if test "x$OCTAVE" != "x"; then
  TESTSUITE_OCTAVE="yes"
else
  TESTSUITE_OCTAVE="no"
fi

if test "x$ax_blas_ok" = "xyes" -a "x$ax_lapack_ok" = "xyes"; then
  TESTSUITE_DYNAREPLUSPLUS="yes"
else
  TESTSUITE_DYNAREPLUSPLUS="no"
fi

AC_MSG_NOTICE([

Dynare is now configured for building the following components...

Binaries (with "make"):
 Dynare preprocessor:                   $BUILD_PREPROCESSOR
 Dynare++:                              $BUILD_DYNAREPLUSPLUS

PDF documentation (with "make pdf"):
 Dynare reference manual:               $BUILD_DYNARE_PDF_MANUAL
 Dynare user guide:                     $BUILD_DYNARE_USERGUIDE
 Dynare++ developer documentation:      $BUILD_DYNAREPLUSPLUS_SRCDOC
 Beamer presentations:                  $BUILD_BEAMER_DOC
 Various other documents:               $BUILD_OTHER_PDF_DOC

HTML documentation (with "make html"):
 Dynare reference manual:               $BUILD_DYNARE_HTML_MANUAL
 Dynare preprocessor developer doc:     $BUILD_DYNARE_PREPROC_DOC
 Dynare internal doc:                   $BUILD_DYNARE_INTERNAL_DOC

Info documentation (with "make info"):  $BUILD_DYNARE_INFO

Testsuites (run with "make check"):
 Dynare for Octave:                     $TESTSUITE_OCTAVE
 Dynare++:                              $TESTSUITE_DYNAREPLUSPLUS
])

AC_OUTPUT
