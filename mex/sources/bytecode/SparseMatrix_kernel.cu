/*
 * Copyright (C) 2007-2012 Dynare Team
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

#ifndef SPARMATRIX_KERNEL
#define SPARMATRIX_KERNEL

// Kernel definition of vector division
__global__ void
VecDiv(double* A, double* B, double* C, int n)
{
  int i = blockIdx.x * 1024 + threadIdx.x;
  if (i < n)
    C[i] = (B[i] != 0.0 ? A[i] / B[i] : A[i]);
}

__global__ void
VecAdd(double* res, double* r, double alpha, double* x, int n)
{
  int i = blockIdx.x * 1024 + threadIdx.x;
  if (i < n)
    res[i] = r[i] + alpha * x[i];
}

__global__ void
VecInc(double* res, double alpha, double* x, int n)
{
  int i = blockIdx.x * 1024 + threadIdx.x;
  if (i < n)
    res[i] += alpha * x[i];
}


__global__ void
update_x(double* x, double alpha, double* y, double omega, double *z)
{
  int i = threadIdx.x;
  x[i] += alpha * y[i] + omega * z[i];
}

__global__ void
Get_LU_dim(int *n, int* A_tild_i, int *A_tild_p, int *nnz_l, int *nnz_u)
{
  nnz_u[0] = 0;
  nnz_l[0] = 0;
  for (int i = 0; i < n[0]; i++)
    {
      for (int j = A_tild_p[i]; j < A_tild_p[i+1]; j++)
        {
          if (A_tild_i[j] < i)
            nnz_l[0]++;
          else if (A_tild_i[j] == i)
            {
              nnz_u[0]++;
              //nnz_l[0]++;
            }
          else
            nnz_u[0]++;
        }
    }
}

__global__ void
Get_LU1_dim(int* n, int *nnz_l, int *nnz_u)
{
  nnz_u[0] = 3+n[0];
  nnz_l[0] = 1+n[0];
}


__global__ void
Get_L_and_U(int *n, double* A_tild_x, int* A_tild_i, int *A_tild_p, double* Lx, int* Li, int *Lp, double* Ux, int* Ui, int* Up)
{
  int nnz_u = 0, nnz_l = 0;
  Lp[0] = 0;
  Up[0] = 0;
  for (int i = 0; i < n[0]; i++)
    {
      for (int j = A_tild_p[i]; j < A_tild_p[i+1]; j++)
        {
          if (A_tild_i[j] < i)
            {
              Lx[nnz_l] = A_tild_x[j];
              Li[nnz_l] = A_tild_i[j];
              nnz_l++;
            }
          else if (A_tild_i[j] == i)
            {
              Ux[nnz_u] = A_tild_x[j];
              Lx[nnz_l] = 1.0;
              Li[nnz_l] = Ui[nnz_u] = A_tild_i[j];
              nnz_u++;
              //nnz_l++;
            }
          else
            {
              Ux[nnz_u] = A_tild_x[j];
              Ui[nnz_u] = A_tild_i[j];
              nnz_u++;
            }
        }
      Lp[i+1] = nnz_l;
      Up[i+1] = nnz_u;
    }
}
#endif
