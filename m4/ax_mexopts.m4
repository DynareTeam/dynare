dnl Copyright (C) 2009-2016 Dynare Team
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

AC_DEFUN([AX_MEXOPTS],
[dnl
AC_REQUIRE([AX_MEXEXT])
AC_REQUIRE([AX_MATLAB_ARCH])
AC_REQUIRE([AX_MATLAB_VERSION])
AC_REQUIRE([AC_PROG_SED])

AX_COMPARE_VERSION([$MATLAB_VERSION], [lt], [7.5], [AC_MSG_ERROR([Your MATLAB is too old, please upgrade to 7.5 (R2007b) at least.])])

AC_MSG_CHECKING([for options to compile MEX for MATLAB])

MATLAB_CPPFLAGS="-I$MATLAB/extern/include"

case ${MATLAB_ARCH} in
  glnx86 | glnxa64)
    MATLAB_DEFS="$MATLAB_DEFS -D_GNU_SOURCE -DNDEBUG"
    MATLAB_CFLAGS="-fexceptions -fPIC -pthread -g -O2"
    MATLAB_CXXFLAGS="-fPIC -pthread -g -O2"
    MATLAB_FFLAGS="-fPIC -g -O2 -fexceptions"
    MATLAB_LDFLAGS_NOMAP="-shared -Wl,--no-undefined -Wl,-rpath-link,$MATLAB/bin/${MATLAB_ARCH} -L$MATLAB/bin/${MATLAB_ARCH}"
    MATLAB_LDFLAGS="$MATLAB_LDFLAGS_NOMAP -Wl,--version-script,$MATLAB/extern/lib/${MATLAB_ARCH}/mexFunction.map"
    MATLAB_LIBS="-lmx -lmex -lmat -lm -lstdc++ -lmwlapack -lmwblas"
    if test "${MATLAB_ARCH}" = "glnx86"; then
      MATLAB_DEFS="$MATLAB_DEFS -D_FILE_OFFSET_BITS=64"
      MATLAB_CFLAGS="$MATLAB_CFLAGS -m32"
      MATLAB_CXXFLAGS="$MATLAB_CXXFLAGS -m32"
    else # glnxa64
      MATLAB_CFLAGS="$MATLAB_CFLAGS -fno-omit-frame-pointer"
      MATLAB_CXXFLAGS="$MATLAB_CXXFLAGS -fno-omit-frame-pointer"
    fi
    ax_mexopts_ok="yes"
    ;;
  win32 | win64)
    MATLAB_CFLAGS="-fexceptions -g -O2"
    MATLAB_CXXFLAGS="-g -O2"
    MATLAB_FFLAGS="-fexceptions -g -O2 -fno-underscoring"
    MATLAB_DEFS="$MATLAB_DEFS -DNDEBUG"
    # Note that static-libstdc++ is only supported since GCC 4.5 (but generates no error on older versions)
    MATLAB_LDFLAGS_NOMAP="-static-libgcc -static-libstdc++ -shared -L$MATLAB/bin/${MATLAB_ARCH}"
    MATLAB_LDFLAGS="$MATLAB_LDFLAGS_NOMAP $(pwd)/$srcdir/mex.def"
    MATLAB_LIBS="-lmex -lmx -lmat -lmwlapack -lmwblas"
    ax_mexopts_ok="yes"
    ;;
  maci | maci64)
    MACOSX_DEPLOYMENT_TARGET='10.7'
    if test "${MATLAB_ARCH}" = "maci"; then
        ARCHS='i386'
    else
        ARCHS='x86_64'
    fi
    MATLAB_DEFS="$MATLAB_DEFS -DNDEBUG"
    MATLAB_CFLAGS="-fno-common -arch $ARCHS -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -fexceptions"
    MATLAB_CXXFLAGS="-fno-common -fexceptions -arch $ARCHS -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET"
    MATLAB_FFLAGS="-fexceptions -fbackslash -arch $ARCHS"
    MATLAB_LDFLAGS_NOMAP="-L$MATLAB/bin/${MATLAB_ARCH} -Wl,-twolevel_namespace -undefined error -arch $ARCHS -mmacosx-version-min=$MACOSX_DEPLOYMENT_TARGET -bundle"
    MATLAB_LDFLAGS="$MATLAB_LDFLAGS_NOMAP -Wl,-exported_symbols_list,$(pwd)/$srcdir/mexFunction-MacOSX.map"
    MATLAB_LIBS="-lmx -lmex -lmat -lstdc++ -lmwlapack -lmwblas"
    ax_mexopts_ok="yes"
    ;;
  *)
    ax_mexopts_ok="no"
    ;;
esac

# Starting from MATLAB 7.8, on 64-bit platforms, BLAS and LAPACK expect 64-bit integers, so make it the default for integers in Fortran code
if test "${MATLAB_ARCH}" = "glnxa64" -o "${MATLAB_ARCH}" = "win64" -o "${MATLAB_ARCH}" = "maci64"; then
  AX_COMPARE_VERSION([$MATLAB_VERSION], [ge], [7.8], [MATLAB_FFLAGS="$MATLAB_FFLAGS -fdefault-integer-8"])
fi

# Converts the MATLAB version number into comparable integers with only major and minor version numbers
# For example, 7.4.2 will become 0704
ax_matlab_ver=`echo "$MATLAB_VERSION" | $SED -e 's/\([[0-9]]*\)\.\([[0-9]]*\).*/Z\1ZZ\2Z/' \
                                             -e 's/Z\([[0-9]]\)Z/Z0\1Z/g' \
                                             -e 's/[[^0-9]]//g'`

MATLAB_DEFS="$MATLAB_DEFS -DMATLAB_VERSION=0x${ax_matlab_ver}"

if test "$ax_mexopts_ok" = "yes"; then
  AC_MSG_RESULT([ok])
else
  AC_MSG_RESULT([unknown])
fi

# Allow user to override default Matlab compilation flags
# Github ticket #1121
if test "x$MATLAB_MEX_CPPFLAGS" != "x"; then
  MATLAB_CPPFLAGS="$MATLAB_CPPFLAGS $MATLAB_MEX_CPPFLAGS"
fi

if test "x$MATLAB_MEX_DEFS" != "x"; then
  MATLAB_DEFS="$MATLAB_DEFS $MATLAB_MEX_DEFS"
fi

if test "x$MATLAB_MEX_CFLAGS" != "x"; then
  MATLAB_CFLAGS="$MATLAB_CFLAGS $MATLAB_MEX_CFLAGS"
fi

if test "x$MATLAB_MEX_CXXFLAGS" != "x"; then
  MATLAB_CXXFLAGS="$MATLAB_CXXFLAGS $MATLAB_MEX_CXXFLAGS"
fi

if test "x$MATLAB_MEX_LDFLAGS" != "x"; then
  MATLAB_LDFLAGS="$MATLAB_LDFLAGS $MATLAB_MEX_LDFLAGS"
fi

if test "x$MATLAB_MEX_LIBS" != "x"; then
  MATLAB_LIBS="$MATLAB_LIBS $MATLAB_MEX_LIBS"
fi

AC_SUBST([MATLAB_CPPFLAGS])
AC_SUBST([MATLAB_DEFS])
AC_SUBST([MATLAB_CFLAGS])
AC_SUBST([MATLAB_CXXFLAGS])
AC_SUBST([MATLAB_LDFLAGS])
AC_SUBST([MATLAB_LIBS])
])
