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

AC_DEFUN([AX_MATLAB_VERSION],
[dnl
AC_REQUIRE([AX_MATLAB])

AC_MSG_CHECKING([for MATLAB version])
if test "x$MATLAB_VERSION" != "x"; then
  case $MATLAB_VERSION in
    *2016a | *2016A)
      MATLAB_VERSION="9.0"
      ;;
    *2015b | *2015B)
      MATLAB_VERSION="8.6"
      ;;
    *2015a | *2015A)
      MATLAB_VERSION="8.5"
      ;;
    *2014b | *2014B)
      MATLAB_VERSION="8.4"
      ;;
    *2014a | *2014A)
      MATLAB_VERSION="8.3"
      ;;
    *2013b | *2013B)
      MATLAB_VERSION="8.2"
      ;;
    *2013a | *2013A)
      MATLAB_VERSION="8.1"
      ;;
    *2012b | *2012B)
      MATLAB_VERSION="8.0"
      ;;
    *2012a | *2012A)
      MATLAB_VERSION="7.14"
      ;;
    *2011b | *2011B)
      MATLAB_VERSION="7.13"
      ;;
    *2011a | *2011A)
      MATLAB_VERSION="7.12"
      ;;
    *2010b | *2010B)
      MATLAB_VERSION="7.11"
      ;;
    *2010a | *2010A)
      MATLAB_VERSION="7.10"
      ;;
    *2009b | *2009B)
      MATLAB_VERSION="7.9"
      ;;
    *2009a | *2009A)
      MATLAB_VERSION="7.8"
      ;;
    *2008b | *2008B)
      MATLAB_VERSION="7.7"
      ;;
    *2008a | *2008A)
      MATLAB_VERSION="7.6"
      ;;
    *2007b | *2007B)
      MATLAB_VERSION="7.5"
      ;;
    *2007a | *2007A)
      MATLAB_VERSION="7.4"
      ;;
    *2006b | *2006B)
      MATLAB_VERSION="7.3"
      ;;
    *2006a | *2006A)
      MATLAB_VERSION="7.2"
      ;;
    *14SP3 | *14sp3)
      MATLAB_VERSION="7.1"
      ;;
    *14SP2 | *14sp2)
      MATLAB_VERSION="7.0.4"
      ;;
    *14SP1 | *14sp1)
      MATLAB_VERSION="7.0.1"
      ;;
    R14 | r14)
      MATLAB_VERSION="7.0.0"
      ;;
  esac
  AC_MSG_RESULT([${MATLAB_VERSION}])
  ax_matlab_version_ok="yes"
else
  AC_MSG_RESULT([unknown])
  ax_matlab_version_ok="no"
fi

AC_ARG_VAR([MATLAB_VERSION], [MATLAB version])
])
