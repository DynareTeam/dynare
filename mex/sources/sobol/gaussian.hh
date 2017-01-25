/* Generates gaussian random deviates from uniform random deviates.
** 
** Pseudo code of the algorithm is given at http://home.online.no/~pjacklam/notes/invnorm  
**
** Copyright (C) 2010-2016 Dynare Team
**
** This file is part of Dynare.
**
** Dynare is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation, either version 3 of the License, or
** (at your option) any later version.
**
** Dynare is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
** 
** AUTHOR(S): stephane DOT adjemian AT univ DASH lemans DOT fr  
*/

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include <dynblas.h>

using namespace std;

#define lb  .02425
#define ub  .97575

#ifdef USE_OMP
# include <omp.h>
#endif
#define DEBUG_OMP 0

#ifdef __MINGW32__
#define __CROSS_COMPILATION__
#endif

#ifdef __MINGW64__
#define __CROSS_COMPILATION__
#endif

#ifdef __CROSS_COMPILATION__
#define M_PI 3.14159265358979323846
#endif

template<typename T> T icdf( const T uniform )
/*
**  This function invert the gaussian cumulative distribution function.
**
*/
{ 
  static T A[6] = 
    {
      -3.969683028665376e+01,
       2.209460984245205e+02,
      -2.759285104469687e+02,
       1.383577518672690e+02,
      -3.066479806614716e+01,
       2.506628277459239e+00
    };
  static T B[5] = 
    {
      -5.447609879822406e+01,
       1.615858368580409e+02,
      -1.556989798598866e+02,
       6.680131188771972e+01,
      -1.328068155288572e+01
    };
  static T C[6] = 
    {
      -7.784894002430293e-03,
      -3.223964580411365e-01,
      -2.400758277161838e+00,
      -2.549732539343734e+00,
       4.374664141464968e+00,
       2.938163982698783e+00
    };
  static T D[4] = 
    {
      7.784695709041462e-03,
      3.224671290700398e-01,
      2.445134137142996e+00,
      3.754408661907416e+00
    };
  T gaussian = (T)0.0;
  if ( (0<uniform)  && (uniform<lb) )
    {
      T tmp;
      tmp = sqrt(-2*log(uniform));
      gaussian = (((((C[0]*tmp+C[1])*tmp+C[2])*tmp+C[3])*tmp+C[4])*tmp+C[5])/((((D[0]*tmp+D[1])*tmp+D[2])*tmp+D[3])*tmp+1);
    }
  else
    {
      if ((lb <= uniform) && (uniform <= ub))
        {
          T tmp, TMP;
          tmp = uniform - .5;
          TMP = tmp*tmp;
          gaussian = (((((A[0]*TMP+A[1])*TMP+A[2])*TMP+A[3])*TMP+A[4])*TMP+A[5])*tmp/(((((B[0]*TMP+B[1])*TMP+B[2])*TMP+B[3])*TMP+B[4])*TMP+1);
        }
      else
        {
          if ((ub < uniform) && (uniform < 1))
            {
              T tmp;
              tmp = sqrt(-2*log(1-uniform));
              gaussian = -(((((C[0]*tmp+C[1])*tmp+C[2])*tmp+C[3])*tmp+C[4])*tmp+C[5])/((((D[0]*tmp+D[1])*tmp+D[2])*tmp+D[3])*tmp+1);
            }
        }
    }
  if ( (0<uniform) && (uniform<1) )
    {
      T tmp, tmp_;
      tmp = .5*erfc(-gaussian/sqrt(2.0))-uniform;
      tmp_ = tmp*sqrt(2*M_PI)*exp(.5*gaussian*gaussian);
      gaussian = gaussian - tmp_/(1+.5*gaussian*tmp_);
    }
  if ( uniform==0)
    {
      gaussian = -INFINITY;
    }
  if ( uniform==1)
    {
      gaussian = INFINITY;
    }
 return(gaussian);
}

template<typename T> void icdfm( const int n, T *U)
{ 
  #if USE_OMP
  #pragma omp parallel for num_threads(omp_get_num_threads())
  #endif
  for(int i=0; i<n; i++)
    {
      U[i] = icdf(U[i]);
    }
  return;
}

template<typename T> void icdfmSigma( const int d, const int n, T *U, const double *LowerCholSigma)
{ 
  double one = 1.0;
  double zero = 0.0;
  blas_int dd(d);
  blas_int nn(n);
  icdfm(n*d, U);
  double tmp[n*d];
  dgemm("N","N",&dd,&nn,&dd,&one,LowerCholSigma,&dd,U,&dd,&zero,tmp,&dd);
  memcpy(U,tmp,d*n*sizeof(double));
  return;
}

template<typename T> void usphere( const int d, const int n, T *U)
{ 
  icdfm(n*d, U);
  #if USE_OMP
  #pragma omp parallel for num_threads(omp_get_num_threads())
  #endif
  for (int j=0; j<n; j++)// sequence index.
    {
      int k = j*d;
      double norm = 0.0;
      for(int i=0; i<d; i++)// dimension index.
        {
          norm = norm + U[k+i]*U[k+i];
        }
      norm = sqrt(norm);
      for(int i=0; i<d; i++)// dimension index.
        {
          U[k+i] = U[k+i]/norm;
        }
    }
  return;
}

template<typename T> void usphereRadius( const int d, const int n, double radius, T *U)
{ 
  icdfm(n*d, U);
  #if USE_OMP
  #pragma omp parallel for num_threads(omp_get_num_threads())
  #endif
  for (int j=0; j<n; j++)// sequence index.
    {
      int k = j*d;
      double norm = 0.0;
      for(int i=0; i<d; i++)// dimension index.
        {
          norm = norm + U[k+i]*U[k+i];
        }
      norm = sqrt(norm);
      for(int i=0; i<d; i++)// dimension index.
        {
          U[k+i] = radius*U[k+i]/norm;
        }
    }
  return;
}
