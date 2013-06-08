  /*
 * Copyright (C) 2007-2013 Dynare Team
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

//define _GLIBCXX_USE_C99_FENV_TR1 1
//include <cfenv>

#include <cstring>
#include <ctime>
#include <sstream>
//#include <gsl/gsl_min.h>
//#include <minimize.h>
#include "SparseMatrix.hh"

#ifdef CUDA
#include "SparseMatrix_kernel.cu"
#endif

using namespace std;
#ifdef _MSC_VER
#include <windows.h>
HINSTANCE hinstLib;

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
/* used in all UMFPACK_report_* routines: */
#define UMFPACK_PRL 0			/* print level */
/* returned by all routines that use Info: */
#define UMFPACK_OK (0)
#define UMFPACK_STATUS 0	/* UMFPACK_OK, or other result */

typedef void (*t_umfpack_dl_free_numeric)(void **Numeric);
t_umfpack_dl_free_numeric umfpack_dl_free_numeric;
typedef void (*t_umfpack_dl_free_symbolic)(void **Symbolic);
t_umfpack_dl_free_symbolic umfpack_dl_free_symbolic;
typedef int64_t (*t_umfpack_dl_solve)(int64_t sys,
                                      const int64_t Ap [ ],
                                      const int64_t Ai [ ],
                                      const double Ax [ ],
                                      double X [ ],
                                      const double B [ ],
                                      void *Numeric,
                                      const double Control [UMFPACK_CONTROL],
                                      double Info [UMFPACK_INFO]);
t_umfpack_dl_solve umfpack_dl_solve;
typedef int64_t (*t_umfpack_dl_numeric)(const int64_t Ap [ ],
                                        const int64_t Ai [ ],
                                        const double Ax [ ],
                                        void *Symbolic,
                                        void **Numeric,
                                        const double Control [UMFPACK_CONTROL],
                                        double Info [UMFPACK_INFO]);
t_umfpack_dl_numeric umfpack_dl_numeric;
typedef int64_t (*t_umfpack_dl_symbolic)(int64_t n_row,
    int64_t n_col,
    const int64_t Ap [ ],
    const int64_t Ai [ ],
    const double Ax [ ],
    void **Symbolic,
    const double Control [UMFPACK_CONTROL],
    double Info [UMFPACK_INFO]);
t_umfpack_dl_symbolic umfpack_dl_symbolic;
typedef void (*t_umfpack_dl_report_info)(const double Control [UMFPACK_CONTROL],
    const double Info [UMFPACK_INFO]);
t_umfpack_dl_report_info umfpack_dl_report_info;
typedef void (*t_umfpack_dl_report_status)(const double Control [UMFPACK_CONTROL],
    int64_t status);
t_umfpack_dl_report_status umfpack_dl_report_status;
typedef void (*t_umfpack_dl_defaults)(double Control [UMFPACK_CONTROL]);
t_umfpack_dl_defaults umfpack_dl_defaults;

#endif




dynSparseMatrix::dynSparseMatrix()
{
  pivotva = NULL;
  g_save_op = NULL;
  g_nop_all = 0;
  mem_mngr.init_Mem();
  symbolic = true;
  alt_symbolic = false;
  alt_symbolic_count = 0;
  max_u = 0;
  min_u = 0x7FFFFFFF;
  res1a = 9.0e60;
  tbreak_g = 0;
  start_compare = 0;
  restart = 0;
  IM_i.clear();
  lu_inc_tol = 1e-10;
  Symbolic = NULL;
  Numeric = NULL;
#ifdef _MSC_VER
  // Get a handle to the DLL module.
  hinstLib = LoadLibrary(TEXT("libmwumfpack.dll"));
  // If the handle is valid, try to get the function address.
  if (hinstLib)
    {
      umfpack_dl_free_numeric = (t_umfpack_dl_free_numeric) GetProcAddress(hinstLib, "umfpack_dl_free_numeric");
      if (!umfpack_dl_free_numeric)
        {
          mexPrintf("umfpack_dl_free_numeric not found\n");
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_free_numeric is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_free_symbolic = (t_umfpack_dl_free_symbolic) GetProcAddress(hinstLib, "umfpack_dl_free_symbolic");
      if (!umfpack_dl_free_symbolic)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_free_symbolic is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_solve = (t_umfpack_dl_solve) GetProcAddress(hinstLib, "umfpack_dl_free_solve");
      if (!umfpack_dl_solve)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_solve is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_numeric = (t_umfpack_dl_numeric) GetProcAddress(hinstLib, "umfpack_dl_numeric");
      if (!umfpack_dl_numeric)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_numeric is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_symbolic = (t_umfpack_dl_symbolic) GetProcAddress(hinstLib, "umfpack_dl_symbolic");
      if (!umfpack_dl_symbolic)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_symbolic is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_report_info = (t_umfpack_dl_report_info) GetProcAddress(hinstLib, "umfpack_dl_report_info");
      if (!umfpack_dl_report_info)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_report_info is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_report_status = (t_umfpack_dl_report_status) GetProcAddress(hinstLib, "umfpack_dl_report_status");
      if (!umfpack_dl_report_status)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_report_status is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_defaults = (t_umfpack_dl_defaults) GetProcAddress(hinstLib, "umfpack_dl_defaults");
      if (!umfpack_dl_defaults)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_defaults is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  else
    {
      mexPrintf("library loading error\n");
      ostringstream tmp;
      tmp << " in main, libmwumfpack.dll not found. \n Check that \\Program files\\MATLAB\\RXXXXX\\bin\\win64 is in the current path.";
      throw FatalExceptionHandling(tmp.str());
    }
#endif
}

dynSparseMatrix::dynSparseMatrix(const int y_size_arg, const int y_kmin_arg, const int y_kmax_arg, const bool print_it_arg, const bool steady_state_arg, const int periods_arg,
                           const int minimal_solving_periods_arg, const double slowc_arg
#ifdef CUDA
                           , const int CUDA_device_arg, cublasHandle_t cublas_handle_arg, cusparseHandle_t cusparse_handle_arg, cusparseMatDescr_t descr_arg
#endif
                           ):
  Evaluate(y_size_arg, y_kmin_arg, y_kmax_arg, print_it_arg, steady_state_arg, periods_arg, minimal_solving_periods_arg, slowc_arg)
{
  pivotva = NULL;
  g_save_op = NULL;
  g_nop_all = 0;
  mem_mngr.init_Mem();
  symbolic = true;
  alt_symbolic = false;
  alt_symbolic_count = 0;
  max_u = 0;
  min_u = 0x7FFFFFFF;
  res1a = 9.0e60;
  tbreak_g = 0;
  start_compare = 0;
  restart = 0;
  IM_i.clear();
  lu_inc_tol = 1e-10;
  Symbolic = NULL;
  Numeric = NULL;
#ifdef CUDA
  CUDA_device = CUDA_device_arg;
  cublas_handle = cublas_handle_arg;
  cusparse_handle = cusparse_handle_arg;
  CUDA_descr = descr_arg;
#endif
#ifdef _MSC_VER
  // Get a handle to the DLL module.
  hinstLib = LoadLibrary(TEXT("libmwumfpack.dll"));
  // If the handle is valid, try to get the function address.
  if (hinstLib != NULL)
    {
      umfpack_dl_free_numeric = (t_umfpack_dl_free_numeric) GetProcAddress(hinstLib, "umfpack_dl_free_numeric");
      if (!umfpack_dl_free_numeric)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_free_numeric is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_free_symbolic = (t_umfpack_dl_free_symbolic) GetProcAddress(hinstLib, "umfpack_dl_free_symbolic");
      if (!umfpack_dl_free_symbolic)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_free_symbolic is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_report_info = (t_umfpack_dl_report_info) GetProcAddress(hinstLib, "umfpack_dl_report_info");
      if (!umfpack_dl_report_info)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_report_info is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_solve = (t_umfpack_dl_solve) GetProcAddress(hinstLib, "umfpack_dl_solve");
      if (!umfpack_dl_solve)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_solve is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_numeric = (t_umfpack_dl_numeric) GetProcAddress(hinstLib, "umfpack_dl_numeric");
      if (!umfpack_dl_numeric)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_numeric is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_symbolic = (t_umfpack_dl_symbolic) GetProcAddress(hinstLib, "umfpack_dl_symbolic");
      if (!umfpack_dl_symbolic)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_symbolic is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_report_status = (t_umfpack_dl_report_status) GetProcAddress(hinstLib, "umfpack_dl_report_status");
      if (!umfpack_dl_report_status)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_report_status is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
      umfpack_dl_defaults = (t_umfpack_dl_defaults) GetProcAddress(hinstLib, "umfpack_dl_defaults");
      if (!umfpack_dl_defaults)
        {
          ostringstream tmp;
          tmp << " in libmwumfpack.dll, the function umfpack_dl_defaults is not found.";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  else
    {
      mexPrintf("library loading error\n");
      ostringstream tmp;
      tmp << " in main, libmwumfpack.dll not found. \n Check that \\Program files\\MATLAB\\RXXXXX\\bin\\win64 in the current path.";
      throw FatalExceptionHandling(tmp.str());
    }
#endif
}


int
dynSparseMatrix::NRow(int r)
{
  return NbNZRow[r];
}

int
dynSparseMatrix::NCol(int c)
{
  return NbNZCol[c];
}

int
dynSparseMatrix::At_Row(int r, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  return NbNZRow[r];
}

int
dynSparseMatrix::Union_Row(int row1, int row2)
{
  NonZeroElem *first1, *first2;
  int n1 = At_Row(row1, &first1);
  int n2 = At_Row(row2, &first2);
  int i1 = 0, i2 = 0, nb_elem = 0;
  while (i1 < n1 && i2 < n2)
    {
      if (first1->c_index == first2->c_index)
        {
          nb_elem++;
          i1++;
          i2++;
          first1 = first1->NZE_R_N;
          first2 = first2->NZE_R_N;
        }
      else if (first1->c_index < first2->c_index)
        {
          nb_elem++;
          i1++;
          first1 = first1->NZE_R_N;
        }
      else
        {
          nb_elem++;
          i2++;
          first2 = first2->NZE_R_N;
        }
    }
  return nb_elem;
}

int
dynSparseMatrix::At_Pos(int r, int c, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  while ((*first)->c_index != c)
    (*first) = (*first)->NZE_R_N;
  return NbNZRow[r];
}

int
dynSparseMatrix::At_Col(int c, NonZeroElem **first)
{
  (*first) = FNZE_C[c];
  return NbNZCol[c];
}

int
dynSparseMatrix::At_Col(int c, int lag, NonZeroElem **first)
{
  (*first) = FNZE_C[c];
  int i = 0;
  while ((*first)->lag_index != lag && (*first))
    (*first) = (*first)->NZE_C_N;
  if ((*first))
    {
      NonZeroElem *firsta = (*first);
      if (!firsta->NZE_C_N)
        i++;
      else
        {
          while (firsta->lag_index == lag && firsta->NZE_C_N)
            {
              firsta = firsta->NZE_C_N;
              i++;
            }
          if (firsta->lag_index == lag)
            i++;
        }
    }
  return i;
}

void
dynSparseMatrix::Delete(const int r, const int c)
{
  NonZeroElem *first = FNZE_R[r], *firsta = NULL;

  while (first->c_index != c)
    {
      firsta = first;
      first = first->NZE_R_N;
    }
  if (firsta != NULL)
    firsta->NZE_R_N = first->NZE_R_N;
  if (first == FNZE_R[r])
    FNZE_R[r] = first->NZE_R_N;
  NbNZRow[r]--;

  first = FNZE_C[c];
  firsta = NULL;
  while (first->r_index != r)
    {
      firsta = first;
      first = first->NZE_C_N;
    }

  if (firsta != NULL)
    firsta->NZE_C_N = first->NZE_C_N;
  if (first == FNZE_C[c])
    FNZE_C[c] = first->NZE_C_N;

  u_liste.push_back(first->u_index);
  mem_mngr.mxFree_NZE(first);
  NbNZCol[c]--;
}

void
dynSparseMatrix::Print(int Size, int *b)
{
  int a, i, j, k, l;
  mexPrintf("   ");
  for (k = 0; k < Size*periods; k++)
    mexPrintf("%-2d ", k);
  mexPrintf("    |    ");
  for (k = 0; k < Size*periods; k++)
    mexPrintf("%8d", k);
  mexPrintf("\n");
  for (i = 0; i < Size*periods; i++)
    {
      NonZeroElem *first = FNZE_R[i];
      j = NbNZRow[i];
      mexPrintf("%-2d ", i);
      a = 0;
      for (k = 0; k < j; k++)
        {
          for (l = 0; l < (first->c_index-a); l++)
            mexPrintf("   ");
          mexPrintf("%-2d ", first->u_index);
          a = first->c_index+1;
          first = first->NZE_R_N;
        }
      for (k = a; k < Size*periods; k++)
        mexPrintf("   ");
      mexPrintf("%-2d ", b[i]);

      first = FNZE_R[i];
      j = NbNZRow[i];
      mexPrintf(" | %-2d ", i);
      a = 0;
      for (k = 0; k < j; k++)
        {
          for (l = 0; l < (first->c_index-a); l++)
            mexPrintf("        ");
          mexPrintf("%8.4f", double (u[first->u_index]));
          a = first->c_index+1;
          first = first->NZE_R_N;
        }
      for (k = a; k < Size*periods; k++)
        mexPrintf("        ");
      mexPrintf("%8.4f", double (u[b[i]]));
      mexPrintf("\n");
    }
}

void
dynSparseMatrix::Insert(const int r, const int c, const int u_index, const int lag_index)
{
  NonZeroElem *firstn, *first, *firsta, *a;
  firstn = mem_mngr.mxMalloc_NZE();
  first = FNZE_R[r];
  firsta = NULL;
  while (first->c_index < c && (a = first->NZE_R_N))
    {
      firsta = first;
      first = a;
    }
  firstn->u_index = u_index;
  firstn->r_index = r;
  firstn->c_index = c;
  firstn->lag_index = lag_index;
  if (first->c_index > c)
    {
      if (first == FNZE_R[r])
        FNZE_R[r] = firstn;
      if (firsta != NULL)
        firsta->NZE_R_N = firstn;
      firstn->NZE_R_N = first;
    }
  else
    {
      first->NZE_R_N = firstn;
      firstn->NZE_R_N = NULL;
    }
  NbNZRow[r]++;
  first = FNZE_C[c];
  firsta = NULL;
  while (first->r_index < r && (a = first->NZE_C_N))
    {
      firsta = first;
      first = a;
    }
  if (first->r_index > r)
    {
      if (first == FNZE_C[c])
        FNZE_C[c] = firstn;
      if (firsta != NULL)
        firsta->NZE_C_N = firstn;
      firstn->NZE_C_N = first;
    }
  else
    {
      first->NZE_C_N = firstn;
      firstn->NZE_C_N = NULL;
    }

  NbNZCol[c]++;
}

void
dynSparseMatrix::Read_SparseMatrix(string file_name, const int Size, int periods, int y_kmin, int y_kmax, bool two_boundaries, int stack_solve_algo, int solve_algo)
{
  unsigned int eq, var;
  int lag;
  filename = file_name;
  mem_mngr.fixe_file_name(file_name);
  /*mexPrintf("steady_state=%d, size=%d, solve_algo=%d, stack_solve_algo=%d, two_boundaries=%d\n",steady_state, Size, solve_algo, stack_solve_algo, two_boundaries);
  mexEvalString("drawnow;");*/
  if (!SaveCode.is_open())
    {
      if (steady_state)
        SaveCode.open((file_name + "_static.bin").c_str(), ios::in | ios::binary);
      else
        SaveCode.open((file_name + "_dynamic.bin").c_str(), ios::in | ios::binary);
      if (!SaveCode.is_open())
        {
          ostringstream tmp;
          if (steady_state)
            tmp << " in Read_SparseMatrix, " << file_name << "_static.bin cannot be opened\n";
          else
            tmp << " in Read_SparseMatrix, " << file_name << "_dynamic.bin cannot be opened\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  IM_i.clear();
  if (two_boundaries)
    {
      if (stack_solve_algo == 5)
        {
          for (int i = 0; i < u_count_init-Size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[make_pair(make_pair(eq, var), lag)] = val;
            }
          for (int j = 0; j < Size; j++)
            IM_i[make_pair(make_pair(j, Size*(periods+y_kmax)), 0)] = j;
        }
      else if (stack_solve_algo >= 0 && stack_solve_algo <= 4)
        {
          for (int i = 0; i < u_count_init-Size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[make_pair(make_pair(var - lag*Size, -lag), eq)] = val;
            }
          for (int j = 0; j < Size; j++)
            IM_i[make_pair(make_pair(Size*(periods+y_kmax), 0), j)] = j;
        }
      else if (stack_solve_algo == 7)
        {
          for (int i = 0; i < u_count_init-Size; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[make_pair(make_pair(eq, lag), var - lag * Size)] = val;
            }
          for (int j = 0; j < Size; j++)
            IM_i[make_pair(make_pair(Size*(periods+y_kmax), 0), j)] = j;
        }

    }
  else
    {
      if ((stack_solve_algo == 5 && !steady_state) || (solve_algo == 5 && steady_state))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[make_pair(make_pair(eq, var), lag)] = val;
            }
        }
      else if (((stack_solve_algo >= 0 || stack_solve_algo <= 4) && !steady_state) || ((solve_algo >= 6 || solve_algo <= 8) && steady_state))
        {
          for (int i = 0; i < u_count_init; i++)
            {
              int val;
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&val), sizeof(val));
              IM_i[make_pair(make_pair(var - lag*Size, -lag), eq)] = val;
            }
        }
    }
  index_vara = (int *) mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int));
  for (int j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
  if (periods+y_kmin+y_kmax > 1)
    for (int i = 1; i < periods+y_kmin+y_kmax; i++)
      {
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int j = 0; j < Size; j++)
          index_vara[j+Size*i] = index_vara[j+Size*(i-1)] + y_size;
      }
  index_equa = (int *) mxMalloc(Size*sizeof(int));
  for (int j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
}

void
dynSparseMatrix::Simple_Init(int Size, map<pair<pair<int, int>, int>, int> &IM, bool &zero_solution)
{
  int i, eq, var, lag;
  map<pair<pair<int, int>, int>, int>::iterator it4;
  NonZeroElem *first;
  pivot = (int *) mxMalloc(Size*sizeof(int));
  pivot_save = (int *) mxMalloc(Size*sizeof(int));
  pivotk = (int *) mxMalloc(Size*sizeof(int));
  pivotv = (double *) mxMalloc(Size*sizeof(double));
  pivotva = (double *) mxMalloc(Size*sizeof(double));
  b = (int *) mxMalloc(Size*sizeof(int));
  line_done = (bool *) mxMalloc(Size*sizeof(bool));

  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op = NULL;
  g_nop_all = 0;
  i = Size*sizeof(NonZeroElem *);
  FNZE_R = (NonZeroElem **) mxMalloc(i);
  FNZE_C = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_R = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_C = (NonZeroElem **) mxMalloc(i);
  i = Size*sizeof(int);
  NbNZRow = (int *) mxMalloc(i);
  NbNZCol = (int *) mxMalloc(i);
  it4 = IM.begin();
  eq = -1;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (i = 0; i < Size; i++)
    {
      line_done[i] = 0;
      FNZE_C[i] = NULL;
      FNZE_R[i] = NULL;
      temp_NZE_C[i] = 0;
      temp_NZE_R[i] = 0;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int u_count1 = Size;
  while (it4 != IM.end())
    {
      var = it4->first.first.second;
      eq = it4->first.first.first;
      lag = it4->first.second;
      if (lag == 0)   /*Build the index for sparse matrix containing the jacobian : u*/
        {
          NbNZRow[eq]++;
          NbNZCol[var]++;
          first = mem_mngr.mxMalloc_NZE();
          first->NZE_C_N = NULL;
          first->NZE_R_N = NULL;
          first->u_index = u_count1;
          first->r_index = eq;
          first->c_index = var;
          first->lag_index = lag;
          if (FNZE_R[eq] == NULL)
            FNZE_R[eq] = first;
          if (FNZE_C[var] == NULL)
            FNZE_C[var] = first;
          if (temp_NZE_R[eq] != NULL)
            temp_NZE_R[eq]->NZE_R_N = first;
          if (temp_NZE_C[var] != NULL)
            temp_NZE_C[var]->NZE_C_N = first;
          temp_NZE_R[eq] = first;
          temp_NZE_C[var] = first;
          u_count1++;
        }
      it4++;
    }
  double cum_abs_sum = 0;
#if USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:cum_abs_sum)
#endif
  for (int i = 0; i < Size; i++)
    {
      b[i] = i;
      cum_abs_sum += fabs(u[i]);
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
  u_count = u_count1;
}

void
dynSparseMatrix::Init_Matlab_Sparse_Simple(int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, bool &zero_solution, mxArray *x0_m)
{
  int eq, var;
  double *b = mxGetPr(b_m);
  if (!b)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't retrieve b vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *Ai = mxGetIr(A_m);
  if (!Ai)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *Aj = mxGetJc(A_m);
  if (!Aj)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't allocate Aj index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *A = mxGetPr(A_m);
  if (!A)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't retrieve A matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }
  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:cum_abs_sum)
#endif
  for (int i = 0; i < Size; i++)
    {
      b[i] = u[i];
      cum_abs_sum += fabs(b[i]);
      x0[i] = y[i];
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  Aj[0] = 0;
  last_var = 0;
  it4 = IM.begin();
  while (it4 != IM.end())
    {
      var = it4->first.first.first;
      if (var != last_var)
        {
          Aj[1+last_var ] = NZE;
          last_var = var;
        }
      eq = it4->first.second;
      int index = it4->second;
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << index << ") out of range for u vector max = " << Size+Size*Size << " allocated = " << u_count_alloc << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (NZE >= max_nze)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, exceeds the capacity of A_m sparse matrix\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      A[NZE] = u[index];
      Ai[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << eq << ") out of range for b vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (var < 0 || var >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << var << ") out of range for index_vara vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << index_vara[var] << ") out of range for y vector max=" << y_size << " (0)\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      it4++;
    }
  Aj[Size] = NZE;
}


void
dynSparseMatrix::Init_UMFPACK_Sparse_Simple(int Size, map<pair<pair<int, int>, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, bool &zero_solution, mxArray *x0_m)
{
  int eq, var;
  //double *b = mxGetPr(b_m);
  *b = (double*)mxMalloc(Size * sizeof(double));
  if (!(*b))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't retrieve b vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse_Simple, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ap = (SuiteSparse_long*)mxMalloc((Size+1) * sizeof(SuiteSparse_long));
  if (!(*Ap))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't allocate Ap index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  size_t prior_nz = IM.size();
  *Ai = (SuiteSparse_long*)mxMalloc(prior_nz * sizeof(SuiteSparse_long));
  if (!(*Ai))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ax = (double*)mxMalloc(prior_nz * sizeof(double));
  if (!(*Ax))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't retrieve Ax matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }

  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < Size; i++)
    {
      int eq = index_vara[i];
      ya[eq+it_*y_size] = y[eq+it_*y_size];
    }
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;

#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:cum_abs_sum)
#endif
  for (int i = 0; i < Size; i++)
    {
      (*b)[i] = u[i];
      cum_abs_sum += fabs((*b)[i]);
      x0[i] = y[i];
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  (*Ap)[0] = 0;
  last_var = 0;
  it4 = IM.begin();
  while (it4 != IM.end())
    {
      var = it4->first.first.first;
      if (var != last_var)
        {
          (*Ap)[1+last_var ] = NZE;
          last_var = var;
        }
      eq = it4->first.second;
      int index = it4->second;
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << index << ") out of range for u vector max = " << Size+Size*Size << " allocated = " << u_count_alloc << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (NZE >= max_nze)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, exceeds the capacity of A_m sparse matrix\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      (*Ax)[NZE] = u[index];
      (*Ai)[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << eq << ") out of range for b vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (var < 0 || var >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << var << ") out of range for index_vara vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        {
          ostringstream tmp;
          tmp << " in Init_Matlab_Sparse_Simple, index (" << index_vara[var] << ") out of range for y vector max=" << y_size << " (0)\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      it4++;
    }
  (*Ap)[Size] = NZE;
}


void
dynSparseMatrix::Init_UMFPACK_Sparse(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, mxArray *x0_m)
{
  int t, eq, var, lag, ti_y_kmin, ti_y_kmax;
  int n = periods * Size;
  *b = (double*)mxMalloc(n * sizeof(double));
  if (!(*b))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't retrieve b vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse_Simple, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ap = (SuiteSparse_long*)mxMalloc((n+1) * sizeof(SuiteSparse_long));
  if (!(*Ap))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't allocate Ap index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  size_t prior_nz = IM.size() * periods;
  *Ai = (SuiteSparse_long*)mxMalloc(prior_nz * sizeof(SuiteSparse_long));
  if (!(*Ai))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ax = (double*)mxMalloc(prior_nz * sizeof(double));
  if (!(*Ax))
    {
      ostringstream tmp;
      tmp << " in Init_UMFPACK_Sparse, can't retrieve Ax matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }
  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;

#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < periods*Size; i++)
    {
      (*b)[i] = 0;
      x0[i] = y[index_vara[Size*y_kmin+i]];
    }
  (*Ap)[0] = 0;
  /*int min_lag = 0;
  int max_lag = 0;*/
  for (t = 0; t < periods; t++)
    {
      last_var = -1;
      it4 = IM.begin();
      while (it4 != IM.end())
        {
          var = it4->first.first.first;
          if (var != last_var)
            {
              (*Ap)[1+last_var + t * Size] = NZE;
              last_var = var;
            }
          eq = it4->first.second+Size*t;
          lag = -it4->first.first.second;
          /*if (t==0)
            {
              if (min_lag > lag)
                min_lag = lag;
              if (max_lag < lag)
                max_lag = lag;
            }*/
          int index = it4->second+ (t-lag) * u_count_init;
          if (var < (periods+y_kmax)*Size)
            {
              ti_y_kmin = -min(t, y_kmin);
              ti_y_kmax = min(periods-(t +1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods-(t+1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin)   /*Build the index for sparse matrix containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
                    {
                      ostringstream tmp;
                      tmp << " in Init_UMFPACK_Sparse, index (" << index << ") out of range for u vector max = " << Size+Size*Size << " allocated = " << u_count_alloc << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (NZE >= max_nze)
                    {
                      ostringstream tmp;
                      tmp << " in Init_UMFPACK_Sparse, exceeds the capacity of A_m sparse matrix\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  (*Ax)[NZE] = u[index];
                  (*Ai)[NZE] = eq - lag * Size;
                  NZE++;
                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= Size * periods)
                    {
                      ostringstream tmp;
                      tmp << " in Init_UMFPACK_Sparse, index (" << eq << ") out of range for b vector\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (var+Size*(y_kmin+t+lag) < 0 || var+Size*(y_kmin+t+lag) >= Size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_UMFPACK_Sparse, index (" << var+Size*(y_kmin+t+lag) << ") out of range for index_vara vector\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (index_vara[var+Size*(y_kmin+t+lag)] < 0 || index_vara[var+Size*(y_kmin+t+lag)] >= y_size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_UMFPACK_Sparse, index (" << index_vara[var+Size*(y_kmin+t+lag)] << ") out of range for y vector max=" << y_size*(periods+y_kmin+y_kmax) << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  (*b)[eq]  += u[index+lag*u_count_init]*y[index_vara[var+Size*(y_kmin+t+lag)]];
                }
            }
          else           /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                {
                  ostringstream tmp;
                  tmp << " in Init_UMFPACK_Sparse, index (" << index << ") out of range for u vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
              if (eq < 0 || eq >= (Size*periods))
                {
                  ostringstream tmp;
                  tmp << " in Init_UMFPACK_Sparse, index (" << eq << ") out of range for b vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
#endif
              (*b)[eq]  += u[index];
            }
          it4++;
        }
    }
  (*Ap)[Size*periods] = NZE;

#ifdef DEBUG
  mexPrintf("*Ax = [");
  for (int i = 0; i < NZE; i++)
    mexPrintf("%f ",(*Ax)[i]);
  mexPrintf("]\n");

  mexPrintf("*Ap = [");
  for (int i = 0; i < n+1; i++)
    mexPrintf("%d ",(*Ap)[i]);
  mexPrintf("]\n");

  mexPrintf("*Ai = [");
  for (int i = 0; i < NZE; i++)
    mexPrintf("%d ",(*Ai)[i]);
  mexPrintf("]\n");
#endif
}

void
dynSparseMatrix::Init_CUDA_Sparse_Simple(int Size, map<pair<pair<int, int>, int>, int> &IM, SuiteSparse_long **Ap, SuiteSparse_long **Ai, double **Ax, double **b, double **x0, bool &zero_solution, mxArray *x0_m)
{
  int eq, var;

  *b = (double*)mxMalloc(Size * sizeof(double));
  if (!(*b))
    {
      ostringstream tmp;
      tmp << " in Init_CUDA_Sparse, can't retrieve b vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *Host_x0 = mxGetPr(x0_m);
  if (!Host_x0)
    {
      ostringstream tmp;
      tmp << " in Init_CUDA_Sparse_Simple, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ap = (SuiteSparse_long*)mxMalloc((Size+1) * sizeof(SuiteSparse_long));
  if (!(*Ap))
    {
      ostringstream tmp;
      tmp << " in Init_CUDA_Sparse, can't allocate Ap index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  size_t prior_nz = IM.size();
  *Ai = (SuiteSparse_long*)mxMalloc(prior_nz * sizeof(SuiteSparse_long));
  if (!(*Ai))
    {
      ostringstream tmp;
      tmp << " in Init_CUDA_Sparse, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *Ax = (double*)mxMalloc(prior_nz * sizeof(double));
  if (!(*Ax))
    {
      ostringstream tmp;
      tmp << " in Init_CUDA_Sparse, can't retrieve Ax matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }

  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < Size; i++)
    {
      int eq = index_vara[i];
      ya[eq+it_*y_size] = y[eq+it_*y_size];
    }

#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;

#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:cum_abs_sum)
#endif
  for (int i = 0; i < Size; i++)
    {
      (*b)[i] = u[i];
      cum_abs_sum += fabs((*b)[i]);
      (*x0)[i] = y[i];
    }
  if (cum_abs_sum < 1e-20)
    zero_solution = true;
  else
    zero_solution = false;

  (*Ap)[0] = 0;
  last_var = -1;
  it4 = IM.begin();
  while (it4 != IM.end())
    {
      var = it4->first.first.first;
      if (var != last_var)
        {
          (*Ap)[1+last_var ] = NZE;
          last_var = var;
        }
      eq = it4->first.second;
      int index = it4->second;
#ifdef DEBUG
      if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
        {
          ostringstream tmp;
          tmp << " in Init_CUDA_Sparse_Simple, index (" << index << ") out of range for u vector max = " << Size+Size*Size << " allocated = " << u_count_alloc << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (NZE >= max_nze)
        {
          ostringstream tmp;
          tmp << " in Init_CUDA_Sparse_Simple, exceeds the capacity of A_m sparse matrix\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      (*Ax)[NZE] = u[index];
      (*Ai)[NZE] = eq;
      NZE++;
#ifdef DEBUG
      if (eq < 0 || eq >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_CUDA_Sparse_Simple, index (" << eq << ") out of range for b vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (var < 0 || var >= Size)
        {
          ostringstream tmp;
          tmp << " in Init_CUDA_Sparse_Simple, index (" << var << ") out of range for index_vara vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (index_vara[var] < 0 || index_vara[var] >= y_size)
        {
          ostringstream tmp;
          tmp << " in Init_CUDA_Sparse_Simple, index (" << index_vara[var] << ") out of range for y vector max=" << y_size << " (0)\n";
          throw FatalExceptionHandling(tmp.str());
        }
#endif
      it4++;
    }
  (*Ap)[Size] = NZE;
}

#ifdef CUDA
void
dynSparseMatrix::Init_CUDA_Sparse(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM, int **Ap, int **Ai, double **Ax, int **Ap_tild, int **Ai_tild, double **A_tild, double **b, double **x0, mxArray *x0_m, int *nnz, int *nnz_tild, int preconditioner)
{
  //cudaError_t cuda_error;
  int t, eq, var, lag, ti_y_kmin, ti_y_kmax;
  int n = periods * Size;
  size_t prior_nz = IM.size() * periods;
  size_t preconditioner_size = 0;
  map<pair<int, int>, int> jacob_struct;

  /* ask cuda how many devices it can find */
  int device_count;
  cudaGetDeviceCount(&device_count);

  cudaSetDevice(CUDA_device);


  double *Host_b = (double*)mxMalloc(n * sizeof(double));
  cudaChk(cudaMalloc((void**)b, n * sizeof(double)), " in Init_Cuda_Sparse, not enought memory to allocate b vector on the graphic card\n");

  double *Host_x0 = mxGetPr(x0_m);
  if (!Host_x0)
    {
      ostringstream tmp;
      tmp << " in Init_Cuda_Sparse, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  cudaChk(cudaMalloc((void**)x0, n * sizeof(double)), " in Init_Cuda_Sparse, not enought memory to allocate x0 vector on the graphic card\n");

  int* Host_Ap = (int*)mxMalloc((n+1) * sizeof(int));


  int* Host_Ai = (int*)mxMalloc(prior_nz * sizeof(int));


  double* Host_Ax = (double*)mxMalloc(prior_nz * sizeof(double));

  int* Host_Ai_tild, * Host_Ap_tild;
  if (preconditioner == 3)
    {
      Host_Ap_tild = (int*) mxMalloc((n+1)*sizeof(int));
      Host_Ai_tild = (int*) mxMalloc(prior_nz*sizeof(int));
      Host_Ap_tild[0] = 0;
    }


  if (preconditioner == 0)
    preconditioner_size = n;
  else if (preconditioner == 1 || preconditioner == 2 || preconditioner == 3)
    preconditioner_size = prior_nz;

  double *Host_A_tild = (double*)mxMalloc(preconditioner_size * sizeof(double));


  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0, NZE_tild = 0;
  int last_eq = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < periods*Size; i++)
    {
      Host_b[i] = 0;
      Host_x0[i] = y[index_vara[Size*y_kmin+i]];
    }

  //Ordered in CSR and not in CSC

  Host_Ap[0] = 0;
  for (t = 0; t < periods; t++)
    {
      last_eq = -1;
      it4 = IM.begin();
      while (it4 != IM.end())
        {
          eq = it4->first.first.first;
          if (eq != last_eq)
            {
#ifdef DEBUG
              if (1+last_eq + t * Size > (n + 1))
                {
                  ostringstream tmp;
                  tmp << " in Init_CUDA_Sparse, 1+last_eq + t * Size (" << 1+last_eq + t * Size << ") out of range for Host_Ap vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
#endif
              Host_Ap[1+last_eq + t * Size] = NZE;
              if (preconditioner == 3 && t == 0)
                 Host_Ap_tild[1+last_eq ] = NZE_tild;
              last_eq = eq;
            }
          var = it4->first.second+Size*t;
          lag = it4->first.first.second;
          int index = it4->second+ (t /*+ lag*/) * u_count_init;
          if (eq < (periods+y_kmax)*Size)
            {
              ti_y_kmin = -min(t, y_kmin);
              ti_y_kmax = min(periods-(t + 1), y_kmax);
              if ((lag <= ti_y_kmax && lag >= ti_y_kmin) || preconditioner == 3)  /*Build the index for sparse matrix containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (index < 0 || index >= u_count_alloc || index > (periods-1)* IM.size() + Size * Size + periods * Size)
                    {
                      ostringstream tmp;
                      tmp << " in Init_CUDA_Sparse, index (" << index << ") out of range for u vector max = " << (periods-1)* IM.size() + Size * Size + periods * Size << " allocated = " << u_count_alloc << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (NZE >= prior_nz)
                    {
                      ostringstream tmp;
                      tmp << " in Init_CUDA_Sparse, exceeds the capacity of A_i or A_x sparse matrix\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  bool to_store = true;
                  if (preconditioner == 0)
                    {
                      if (lag == 0 && it4->first.second == eq)
                        Host_A_tild[var] = u[index];
                    }
                  else if (preconditioner == 1 || preconditioner == 2)
                    Host_A_tild[NZE] = u[index];
                  else if (preconditioner == 3)
                    {
                      if (lag > ti_y_kmax || lag < ti_y_kmin)
                        {
                          Host_b[eq + t * Size]  += u[index]*y[index_vara[var+Size*(y_kmin+lag)]];
                          to_store = false;
                        }
                      if (t == 0)
                        {
                           map<pair<int, int>, int>::const_iterator it = jacob_struct.find(make_pair(eq + t * Size, var));
                           if (it != jacob_struct.end())
                             Host_A_tild[it->second] += u[index];
                           else
                            {
                              jacob_struct[make_pair(eq, var)] = NZE_tild;
                              Host_A_tild[NZE_tild] = u[index];
                              Host_Ai_tild[NZE_tild] = var;
                              NZE_tild++;
                            }
                        }
                    }
                  if (to_store)
                    {
                      Host_Ax[NZE] = u[index];
                      Host_Ai[NZE] = var + lag * Size;
                      NZE++;
                    }
                }
              else
                {
#ifdef DEBUG
                  if (var < 0 || var >= Size * periods)
                    {
                      ostringstream tmp;
                      tmp << " in Init_CUDA_Sparse, index (" << var << ") out of range for b vector\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (var+Size*(y_kmin+t+lag) < 0 || var+Size*(y_kmin+lag) >= Size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_CUDA_Sparse, index (" << var+Size*(y_kmin+lag) << ") out of range for index_vara vector max=" << Size*(periods+y_kmin+y_kmax) << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (index_vara[var+Size*(y_kmin+lag)] < 0 || index_vara[var+Size*(y_kmin+lag)] >= y_size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_CUDA_Sparse, index (" << index_vara[var+Size*(y_kmin+lag)] << ") out of range for y vector max=" << y_size*(periods+y_kmin+y_kmax) << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  Host_b[eq + t * Size]  += u[index]*y[index_vara[var+Size*(y_kmin+lag)]];
                }
            }
          else           // ...and store it in the u vector
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                {
                  ostringstream tmp;
                  tmp << " in Init_CUDA_Sparse, index (" << index << ") out of range for u vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
              if (var < 0 || var >= (Size*periods))
                {
                  ostringstream tmp;
                  tmp << " in Init_CUDA_Sparse, index (" << var << ") out of range for b vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
#endif
              Host_b[var]  += u[index];
            }
          it4++;
        }
    }
  Host_Ap[Size*periods] = NZE;
  if (preconditioner == 3)
    {
      int* tmp_Ap_tild = (int*) mxMalloc((Size + 1) * sizeof(int) );
      int* tmp_Ai_tild = (int*) mxMalloc(NZE_tild * sizeof(int) );
      double* tmp_A_tild = (double*) mxMalloc(NZE_tild * sizeof(double) );
      memcpy(tmp_Ap_tild, Host_Ap_tild, (Size + 1) * sizeof(int));
      memcpy(tmp_Ai_tild, Host_Ai_tild, NZE_tild * sizeof(int));
      memcpy(tmp_A_tild, Host_A_tild, NZE_tild * sizeof(double));
      //int NZE_tild_old = NZE_tild;
      NZE_tild = 0;
      Host_Ap_tild[0] = NZE_tild;

      for (int i = 0; i < Size; i++)
        {
          for(int j = tmp_Ap_tild[i]; j < tmp_Ap_tild[i+1]; j++)
            if (abs(tmp_A_tild[j]) > 1.0e-20 )
              {
                Host_A_tild[NZE_tild] = tmp_A_tild[j];
                Host_Ai_tild[NZE_tild] = tmp_Ai_tild[j];
                NZE_tild++;
              }
          Host_Ap_tild[i+1] = NZE_tild;
        }
      mxFree(tmp_Ap_tild);
      mxFree(tmp_Ai_tild);
      mxFree(tmp_A_tild);
    }

  *nnz = NZE;
  *nnz_tild = NZE_tild;
  if (preconditioner == 1 || preconditioner == 2 || preconditioner == 3)
    preconditioner_size = NZE;


#ifdef DEBUG
  mexPrintf("Host_Ax = [");
  for (int i = 0; i < NZE; i++)
    mexPrintf("%f ",Host_Ax[i]);
  mexPrintf("]\n");

  mexPrintf("Host_Ap = [");
  for (int i = 0; i < n+1; i++)
    mexPrintf("%d ",Host_Ap[i]);
  mexPrintf("]\n");

  mexPrintf("Host_Ai = [");
  for (int i = 0; i < NZE; i++)
    mexPrintf("%d ",Host_Ai[i]);
  mexPrintf("]\n");
#endif
  cudaChk(cudaMalloc((void**)Ai, NZE * sizeof(int)), " in Init_Cuda_Sparse, can't allocate Ai index vector on the graphic card\n");
  cudaChk(cudaMalloc((void**)Ax, NZE * sizeof(double)), "  in Init_Cuda_Sparse, can't allocate Ax on the graphic card\n");
  cudaChk(cudaMalloc((void**)Ap, (n+1) * sizeof(int)), " in Init_Cuda_Sparse, can't allocate Ap index vector on the graphic card\n");
  if (preconditioner == 3)
    {
      cudaChk(cudaMalloc((void**)Ai_tild, NZE_tild * sizeof(int)), " in Init_Cuda_Sparse, can't allocate Ai_tild index vector on the graphic card\n");
      cudaChk(cudaMalloc((void**)Ap_tild, (n+1) * sizeof(int)), " in Init_Cuda_Sparse, can't allocate Ap_tild index vector on the graphic card\n");
    }
  cudaChk(cudaMalloc((void**)A_tild, preconditioner_size * sizeof(double)), "  in Init_Cuda_Sparse, can't allocate A_tild on the graphic card\n");

  cudaChk(cudaMemcpy(*x0,     Host_x0,     n *                   sizeof(double), cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy x0 = Host_x0 failed");
  cudaChk(cudaMemcpy(*b,      Host_b,      n *                   sizeof(double), cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy b = Host_b failed");
  cudaChk(cudaMemcpy(*Ap,     Host_Ap,     (n + 1) *             sizeof(int),    cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy Ap = Host_Ap failed");
  cudaChk(cudaMemcpy(*Ai,     Host_Ai,     NZE *                 sizeof(int),    cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy Ai = Host_Ai failed");
  cudaChk(cudaMemcpy(*Ax,     Host_Ax,     NZE *                 sizeof(double), cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy Ax = Host_Ax failed");
  if (preconditioner == 3)
    {
      cudaChk(cudaMemcpy(*Ap_tild,     Host_Ap_tild,     (n + 1) *             sizeof(int),    cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy Ap_tild = Host_Ap_tild failed");
      cudaChk(cudaMemcpy(*Ai_tild,     Host_Ai_tild,     NZE_tild *                 sizeof(int),    cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy Ai_tild = Host_Ai_til failed");
    }
  cudaChk(cudaMemcpy(*A_tild, Host_A_tild, preconditioner_size * sizeof(double), cudaMemcpyHostToDevice), " in Init_CUDA_Sparse, cudaMemcpy A_tild = Host_A_tild failed");
}
#endif


void
PrintM(int n, double* Ax, mwIndex *Ap, mwIndex *Ai)
{
  int nnz = Ap[n];
  double *A = (double*)mxMalloc(n * n * sizeof(double));
  memset(A,0,n * n  * sizeof(double));
  int k = 0;
  for (int i = 0; i< n; i++)
    {
      for (int j = Ap[i]; j < Ap[i + 1]; j++)
        {
          int row = Ai[j];
          A[row *n + i] = Ax[j];
          k++;
        }
    }
  if (nnz != k)
    mexPrintf("Problem nnz(%d) != number of elements(%d)\n", nnz, k);
  mexPrintf("----------------------\n");
  //mexEvalString("drawnow;");
  for (int i = 0; i < n ; i++)
    {
      for (int j = 0; j < n; j++)
        mexPrintf("%-6.3f ",A[i * n + j]);
      mexPrintf("\n");
    }
  mxFree(A);
}

void
dynSparseMatrix::Init_Matlab_Sparse(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, mxArray *x0_m)
{
  int t, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double *b = mxGetPr(b_m);

  if (!b)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't retrieve b vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *x0 = mxGetPr(x0_m);
  if (!x0)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse_Simple, can't retrieve x0 vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *Aj = mxGetJc(A_m);
  if (!Aj)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't allocate Aj index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *Ai = mxGetIr(A_m);
  if (!Ai)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *A = mxGetPr(A_m);
  if (!A)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't retrieve A matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }

  map<pair<pair<int, int>, int>, int>::iterator it4;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      x0[i] = y[index_vara[Size*y_kmin+i]];
    }
  Aj[0] = 0;
  for (t = 0; t < periods; t++)
    {
      last_var = 0;
      it4 = IM.begin();
      while (it4 != IM.end())
        {
          var = it4->first.first.first;
          if (var != last_var)
            {
              Aj[1+last_var + t * Size] = NZE;
              last_var = var;
            }
          eq = it4->first.second+Size*t;
          lag = -it4->first.first.second;
          int index = it4->second+ (t-lag) * u_count_init;
          if (var < (periods+y_kmax)*Size)
            {
              ti_y_kmin = -min(t, y_kmin);
              ti_y_kmax = min(periods-(t +1), y_kmax);
              int ti_new_y_kmax = min(t, y_kmax);
              int ti_new_y_kmin = -min(periods-(t+1), y_kmin);
              if (lag <= ti_new_y_kmax && lag >= ti_new_y_kmin)   /*Build the index for sparse matrix containing the jacobian : u*/
                {
#ifdef DEBUG
                  if (index < 0 || index >= u_count_alloc || index > Size + Size*Size)
                    {
                      ostringstream tmp;
                      tmp << " in Init_Matlab_Sparse, index (" << index << ") out of range for u vector max = " << Size+Size*Size << " allocated = " << u_count_alloc << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (NZE >= max_nze)
                    {
                      ostringstream tmp;
                      tmp << " in Init_Matlab_Sparse, exceeds the capacity of A_m sparse matrix\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  A[NZE] = u[index];
                  Ai[NZE] = eq - lag * Size;
                  NZE++;
                }
              if (lag > ti_y_kmax || lag < ti_y_kmin)
                {
#ifdef DEBUG
                  if (eq < 0 || eq >= Size * periods)
                    {
                      ostringstream tmp;
                      tmp << " in Init_Matlab_Sparse, index (" << eq << ") out of range for b vector\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (var+Size*(y_kmin+t+lag) < 0 || var+Size*(y_kmin+t+lag) >= Size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_Matlab_Sparse, index (" << var+Size*(y_kmin+t+lag) << ") out of range for index_vara vector\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  if (index_vara[var+Size*(y_kmin+t+lag)] < 0 || index_vara[var+Size*(y_kmin+t+lag)] >= y_size*(periods+y_kmin+y_kmax))
                    {
                      ostringstream tmp;
                      tmp << " in Init_Matlab_Sparse, index (" << index_vara[var+Size*(y_kmin+t+lag)] << ") out of range for y vector max=" << y_size*(periods+y_kmin+y_kmax) << "\n";
                      throw FatalExceptionHandling(tmp.str());
                    }
#endif
                  b[eq]  += u[index+lag*u_count_init]*y[index_vara[var+Size*(y_kmin+t+lag)]];
                }
            }
          else           /* ...and store it in the u vector*/
            {
#ifdef DEBUG
              if (index < 0 || index >= u_count_alloc)
                {
                  ostringstream tmp;
                  tmp << " in Init_Matlab_Sparse, index (" << index << ") out of range for u vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
              if (eq < 0 || eq >= (Size*periods))
                {
                  ostringstream tmp;
                  tmp << " in Init_Matlab_Sparse, index (" << eq << ") out of range for b vector\n";
                  throw FatalExceptionHandling(tmp.str());
                }
#endif
              b[eq]  += u[index];
            }
          it4++;
        }
    }
  Aj[Size*periods] = NZE;
}

void
dynSparseMatrix::Init_GE(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM)
{
  int t, i, eq, var, lag, ti_y_kmin, ti_y_kmax;
  double tmp_b = 0.0;
  map<pair<pair<int, int>, int>, int>::iterator it4;
  NonZeroElem *first;
  pivot = (int *) mxMalloc(Size*periods*sizeof(int));
  pivot_save = (int *) mxMalloc(Size*periods*sizeof(int));
  pivotk = (int *) mxMalloc(Size*periods*sizeof(int));
  pivotv = (double *) mxMalloc(Size*periods*sizeof(double));
  pivotva = (double *) mxMalloc(Size*periods*sizeof(double));
  b = (int *) mxMalloc(Size*periods*sizeof(int));
  line_done = (bool *) mxMalloc(Size*periods*sizeof(bool));
  mem_mngr.init_CHUNK_BLCK_SIZE(u_count);
  g_save_op = NULL;
  g_nop_all = 0;
  i = (periods+y_kmax+1)*Size*sizeof(NonZeroElem *);
  FNZE_R = (NonZeroElem **) mxMalloc(i);
  FNZE_C = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_R = (NonZeroElem **) mxMalloc(i);
  NonZeroElem **temp_NZE_C = (NonZeroElem **) mxMalloc(i);
  i = (periods+y_kmax+1)*Size*sizeof(int);
  NbNZRow = (int *) mxMalloc(i);
  NbNZCol = (int *) mxMalloc(i);

#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      line_done[i] = 0;
    }
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < (periods+y_kmax+1)*Size; i++)
    {
      FNZE_C[i] = NULL;
      FNZE_R[i] = NULL;
      temp_NZE_C[i] = NULL;
      temp_NZE_R[i] = NULL;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }
  int nnz = 0;
  //pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag) schedule(dynamic)
  for (t = 0; t < periods; t++)
    {
      ti_y_kmin = -min(t, y_kmin);
      ti_y_kmax = min(periods-(t+1), y_kmax);
      it4 = IM.begin();
      eq = -1;
      //pragma omp ordered
      while (it4 != IM.end())
        {
          var = it4->first.first.second;
          if (eq != it4->first.first.first+Size*t)
            tmp_b = 0;
          eq = it4->first.first.first+Size*t;
          lag = it4->first.second;
          if (var < (periods+y_kmax)*Size)
            {
              lag = it4->first.second;
              if (lag <= ti_y_kmax && lag >= ti_y_kmin)   /*Build the index for sparse matrix containing the jacobian : u*/
                {
                  nnz++;
                  var += Size*t;
                  NbNZRow[eq]++;
                  NbNZCol[var]++;
                  first = mem_mngr.mxMalloc_NZE();
                  first->NZE_C_N = NULL;
                  first->NZE_R_N = NULL;
                  first->u_index = it4->second+u_count_init*t;
                  first->r_index = eq;
                  first->c_index = var;
                  first->lag_index = lag;
                  if (FNZE_R[eq] == NULL)
                    FNZE_R[eq] = first;
                  if (FNZE_C[var] == NULL)
                    FNZE_C[var] = first;
                  if (temp_NZE_R[eq] != NULL)
                    temp_NZE_R[eq]->NZE_R_N = first;
                  if (temp_NZE_C[var] != NULL)
                    temp_NZE_C[var]->NZE_C_N = first;
                  temp_NZE_R[eq] = first;
                  temp_NZE_C[var] = first;
                }
              else       /*Build the additive terms ooutside the simulation periods related to the first lags and the last leads...*/
                {
                  if (lag < ti_y_kmin)
                    {
                      tmp_b += u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];
                    }
                  else
                    {
                      tmp_b += u[it4->second+u_count_init*t]*y[index_vara[var+Size*(y_kmin+t)]];

                    }
                }
            }
          else           /* ...and store it in the u vector*/
            {
              b[eq] = it4->second+u_count_init*t;
              u[b[eq]] += tmp_b;
              tmp_b = 0;
            }
          it4++;
        }
    }
  //mexPrintf("nnz/n=%f\n", double(nnz)/double(periods*Size));
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
}

int
dynSparseMatrix::Get_u()
{
  if (!u_liste.empty())
    {
      int i = u_liste.back();
      u_liste.pop_back();
      return i;
    }
  else
    {
      if (u_count < u_count_alloc)
        {
          int i = u_count;
          u_count++;
          return i;
        }
      else
        {
          u_count_alloc += 5*u_count_alloc_save;
          u = (double *) mxRealloc(u, u_count_alloc*sizeof(double));
          if (!u)
            {
              ostringstream tmp;
              tmp << " in Get_u, memory exhausted (realloc(" << u_count_alloc*sizeof(double) << "))\n";
              throw FatalExceptionHandling(tmp.str());
            }
          int i = u_count;
          u_count++;
          return i;
        }
    }
}

void
dynSparseMatrix::Delete_u(int pos)
{
  u_liste.push_back(pos);
}

void
dynSparseMatrix::Clear_u()
{
  u_liste.clear();
}

void
dynSparseMatrix::Print_u()
{
  for (unsigned int i = 0; i < u_liste.size(); i++)
    mexPrintf("%d ", u_liste[i]);
}

void
dynSparseMatrix::End_GE(int Size)
{
  mem_mngr.Free_All();
  mxFree(FNZE_R);
  mxFree(FNZE_C);
  mxFree(NbNZRow);
  mxFree(NbNZCol);
  mxFree(b);
  mxFree(line_done);
  mxFree(pivot);
  mxFree(pivot_save);
  mxFree(pivotk);
  mxFree(pivotv);
  mxFree(pivotva);
}

bool
dynSparseMatrix::compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size)
{
  long int i, j, nop = nop4/2;
  double r = 0.0;
  bool OK = true;
  t_save_op_s *save_op_s, *save_opa_s, *save_opaa_s;
  int *diff1, *diff2;
  diff1 = (int *) mxMalloc(nop*sizeof(int));
  diff2 = (int *) mxMalloc(nop*sizeof(int));
  int max_save_ops_first = -1;
  j = i = 0;
  while (i < nop4 && OK)
    {
      save_op_s = (t_save_op_s *) &(save_op[i]);
      save_opa_s = (t_save_op_s *) &(save_opa[i]);
      save_opaa_s = (t_save_op_s *) &(save_opaa[i]);
      diff1[j] = save_op_s->first-save_opa_s->first;
      if (max_save_ops_first < save_op_s->first+diff1[j]*(periods-beg_t))
        {
          max_save_ops_first = save_op_s->first+diff1[j]*(periods-beg_t);
        }
      switch (save_op_s->operat)
        {
        case IFLD:
        case IFDIV:
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first));
          i += 2;
          break;
        case IFLESS:
        case IFSUB:
          diff2[j] = save_op_s->second-save_opa_s->second;
          OK = (save_op_s->operat == save_opa_s->operat && save_opa_s->operat == save_opaa_s->operat
                && diff1[j] == (save_opa_s->first-save_opaa_s->first)
                && diff2[j] == (save_opa_s->second-save_opaa_s->second));
          i += 3;
          break;
        default:
          ostringstream tmp;
          tmp << " in compare, unknown operator = " << save_op_s->operat << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      j++;
    }
  // the same pivot for all remaining periods
  if (OK)
    {
      for (int i = beg_t; i < periods; i++)
        {
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
          for (int j = 0; j < Size; j++)
            pivot[i*Size+j] = pivot[(i-1)*Size+j]+Size;
        }
      if (max_save_ops_first >= u_count_alloc)
        {
          u_count_alloc += max_save_ops_first;
          u = (double *) mxRealloc(u, u_count_alloc*sizeof(double));
          if (!u)
            {
              ostringstream tmp;
              tmp << " in compare, memory exhausted (realloc(" << u_count_alloc*sizeof(double) << "))\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      /*#ifdef USE_OMP
      #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
      #endif*/
      for (int t = 1; t < periods-beg_t-y_kmax; t++)
        {
          int i = j = 0;
          double *up;
          while (i < nop4)
            {
              t_save_op_s *save_op_s = (t_save_op_s *) (&(save_op[i]));
              up = &u[save_op_s->first+t*diff1[j]];
              switch (save_op_s->operat)
                {
                case IFLD:
                  r = *up;
                  i += 2;
                  break;
                case IFDIV:
                  *up /= r;
                  i += 2;
                  break;
                case IFSUB:
                  *up -= u[save_op_s->second+t*diff2[j]]*r;;
                  i += 3;
                  break;
                case IFLESS:
                  *up = -u[save_op_s->second+t*diff2[j]]*r;
                  i += 3;
                  break;
                }
              j++;
            }
        }
      int t1 = max(1, periods-beg_t-y_kmax);
      int periods_beg_t = periods-beg_t;
      /*#ifdef USE_OMP
      #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
      #endif*/
      for (int t = t1; t < periods_beg_t; t++)
        {
          int i = j = 0;
          int gap = periods_beg_t-t;
          /*#ifdef USE_OMP
          #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
          #endif*/
          while (i < nop4)
            {
              t_save_op_s *save_op_s = (t_save_op_s *) (&(save_op[i]));
              if (save_op_s->lag < gap)
                {
                  double *up = &u[save_op_s->first+t*diff1[j]];
                  switch (save_op_s->operat)
                    {
                    case IFLD:
                      r = *up;
                      i += 2;
                      break;
                    case IFDIV:
                      *up /= r;
                      i += 2;
                      break;
                    case IFSUB:
                      *up -= u[save_op_s->second+t*diff2[j]]*r;
                      i += 3;
                      break;
                    case IFLESS:
                      *up = -u[save_op_s->second+t*diff2[j]]*r;
                      i += 3;
                      break;
                    }
                }
              else
                {
                  switch (save_op_s->operat)
                    {
                    case IFLD:
                    case IFDIV:
                      i += 2;
                      break;
                    case IFSUB:
                    case IFLESS:
                      i += 3;
                      break;
                    }
                }
              j++;
            }
        }
    }
  mxFree(diff1);
  mxFree(diff2);
  return OK;
}

int
dynSparseMatrix::complete(int beg_t, int Size, int periods, int *b)
{
  long int i, j, k, nop, nopa, nop1, cal_y, nb_var, pos, max_var, min_var;
  NonZeroElem *first;
  int *save_code;
  int *diff;
  double yy = 0.0, err;

  int size_of_save_code = (1+y_kmax)*Size*(Size+1+4)/2*4;
  save_code = (int *) mxMalloc(size_of_save_code*sizeof(int));
  int size_of_diff = (1+y_kmax)*Size*(Size+1+4);
  diff = (int *) mxMalloc(size_of_diff*sizeof(int));
  cal_y = y_size*y_kmin;

  i = (beg_t+1)*Size-1;
  nop = 0;
  for (j = i; j > i-Size; j--)
    {
      pos = pivot[j];
      nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      save_code[nop] = IFLDZ;
      save_code[nop+1] = 0;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
      for (k = 0; k < nb_var; k++)
        {
          save_code[nop] = IFMUL;
          save_code[nop+1] = index_vara[first->c_index]+cal_y;
          save_code[nop+2] = first->u_index;
          save_code[nop+3] = first->lag_index;
#ifdef DEBUG
          if ((nop+3) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
          nop += 4;
          first = first->NZE_R_N;
        }
      save_code[nop] = IFADD;
      save_code[nop+1] = b[pos];
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
      save_code[nop] = IFSTP;
      save_code[nop+1] = index_vara[j]+y_size*y_kmin;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
#ifdef DEBUG
      if ((nop+2) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
#endif
      nop += 4;
    }
  i = beg_t*Size-1;
  nop1 = nopa = 0;
  for (j = i; j > i-Size; j--)
    {
      pos = pivot[j];
      nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      diff[nopa] = 0;
      diff[nopa+1] = 0;
      nopa += 2;
      nop1 += 4;
      for (k = 0; k < nb_var; k++)
        {
          diff[nopa] = save_code[nop1+1]-(index_vara[first->c_index]+cal_y);
          diff[nopa+1] = save_code[nop1+2]-(first->u_index);
#ifdef DEBUG
          if ((nop1+2) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
          if ((nopa+1) >= size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
          nopa += 2;
          nop1 += 4;
          first = first->NZE_R_N;
        }
      diff[nopa] = save_code[nop1+1]-(b[pos]);
      diff[nopa+1] = 0;
#ifdef DEBUG
      if ((nop1+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
      diff[nopa] = save_code[nop1+1]-(index_vara[j]+y_size*y_kmin);
      diff[nopa+1] = 0;
#ifdef DEBUG
      if ((nop1+4) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
#endif
      nopa += 2;
      nop1 += 4;
    }
  max_var = (periods+y_kmin)*y_size;
  min_var = y_kmin*y_size;
  /*#ifdef USE_OMP
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  #endif*/
  for (int t = periods+y_kmin-1; t >= beg_t+y_kmin; t--)
    {
      int j = 0, k;
      int ti = t-y_kmin-beg_t;
      for (int i = 0; i < nop; i += 4)
        {
          switch (save_code[i])
            {
            case IFLDZ:
              yy = 0;
              break;
            case IFMUL:
              k = save_code[i+1]+ti*diff[j];
              if (k < max_var && k > min_var)
                {
                  yy += y[k]*u[save_code[i+2]+ti*diff[j+1]];
                }
              break;
            case IFADD:
              yy = -(yy+u[save_code[i+1]+ti*diff[j]]);
              break;
            case IFSTP:
              k = save_code[i+1]+ti*diff[j];
              err = yy - y[k];
              y[k] += slowc*(err);
              break;
            }
          j += 2;
        }
    }
  mxFree(save_code);
  mxFree(diff);
  return (beg_t);
}

void
dynSparseMatrix::bksub(int tbreak, int last_period, int Size, double slowc_l)
{
  NonZeroElem *first;
  int i, j, k;
  double yy;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    y[i] = ya[i];
  if (symbolic && tbreak)
    last_period = complete(tbreak, Size, periods, b);
  else
    last_period = periods;
  for (int t = last_period+y_kmin-1; t >= y_kmin; t--)
    {
      int ti = (t-y_kmin)*Size;
      int cal = y_kmin*Size;
      int cal_y = y_size*y_kmin;
      for (i = ti-1; i >= ti-Size; i--)
        {
          j = i+cal;
          int pos = pivot[i+Size];
          int nb_var = At_Row(pos, &first);
          first = first->NZE_R_N;
          nb_var--;
          int eq = index_vara[j]+y_size;
          yy = 0;
          for (k = 0; k < nb_var; k++)
            {
              yy += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
              first = first->NZE_R_N;
            }
          yy = -(yy+y[eq]+u[b[pos]]);
          direction[eq] = yy;
          y[eq] += slowc_l*yy;
        }
    }
}

void
dynSparseMatrix::simple_bksub(int it_, int Size, double slowc_l)
{
  int i, k;
  double yy;
  NonZeroElem *first;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size; i++)
    y[i+it_*y_size] = ya[i+it_*y_size];
  for (i = Size-1; i >= 0; i--)
    {
      int pos = pivot[i];
      int nb_var = At_Row(pos, &first);
      first = first->NZE_R_N;
      nb_var--;
      int eq = index_vara[i];
      yy = 0;
      for (k = 0; k < nb_var; k++)
        {
          yy += y[index_vara[first->c_index]+it_*y_size]*u[first->u_index];
          first = first->NZE_R_N;
        }
      yy = -(yy+y[eq+it_*y_size]+u[b[pos]]);
      direction[eq+it_*y_size] = yy;
      y[eq+it_*y_size] += slowc_l*yy;
    }
}

void
dynSparseMatrix::CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods)
{
  const double epsilon = 1e-7;
  fstream SaveResult;
  ostringstream out;
  out << "Result" << iter;
  SaveResult.open(out.str().c_str(), ios::in);
  if (!SaveResult.is_open())
    {
      ostringstream tmp;
      tmp << " in CheckIt, Result file cannot be opened\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mexPrintf("Reading Result...");
  int row, col;
  SaveResult >> row;
  mexPrintf("row=%d\n", row);
  SaveResult >> col;
  mexPrintf("col=%d\n", col);
  double G1a;
  mexPrintf("Allocated\n");
  NonZeroElem *first;
  for (int j = 0; j < col; j++)
    {
      mexPrintf("j=%d ", j);
      int nb_equ = At_Col(j, &first);
      mexPrintf("nb_equ=%d\n", nb_equ);
      int line;
      if (first)
        line = first->r_index;
      else
        line = -9999999;
      for (int i = 0; i < row; i++)
        {
          SaveResult >> G1a;
          if (line == i)
            {
              if (abs(u[first->u_index]/G1a-1) > epsilon)
                mexPrintf("Problem at r=%d c=%d u[first->u_index]=%5.14f G1a[i][j]=%5.14f %f\n", i, j, u[first->u_index], G1a, u[first->u_index]/G1a-1);
              first = first->NZE_C_N;
              if (first)
                line = first->r_index;
              else
                line = -9999999;
            }
          else
            {
              if (G1a != 0.0)
                mexPrintf("Problem at r=%d c=%d G1a[i][j]=%f\n", i, j, G1a);
            }
        }
    }
  mexPrintf("G1a red done\n");
  SaveResult >> row;
  mexPrintf("row(2)=%d\n", row);
  double *B;
  B = (double *) mxMalloc(row*sizeof(double));
  for (int i = 0; i < row; i++)
    SaveResult >> B[i];
  SaveResult.close();
  mexPrintf("done\n");
  mexPrintf("Comparing...");
  for (int i = 0; i < row; i++)
    {
      if (abs(u[b[i]]+B[i]) > epsilon)
        mexPrintf("Problem at i=%d u[b[i]]=%f B[i]=%f\n", i, u[b[i]], B[i]);
    }
  mxFree(B);
}

void
dynSparseMatrix::Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int *b)
{
  const double epsilon = 1e-10;
  Init_GE(periods, y_kmin, y_kmax, Size, IM_i);
  NonZeroElem *first;
  int cal_y = y_kmin*Size;
  mexPrintf("     ");
  for (int i = 0; i < Size; i++)
    mexPrintf(" %8d", i);
  mexPrintf("\n");
  for (int t = y_kmin; t < periods+y_kmin; t++)
    {
      mexPrintf("t=%5d", t);
      for (int i = 0; i < Size; i++)
        mexPrintf(" %d %1.6f", t*y_size+index_vara[i], y[t*y_size+index_vara[i]]);
      mexPrintf("\n");
    }
  for (int i = 0; i < Size*periods; i++)
    {
      double res = 0;
      int pos = pivot[i];
      mexPrintf("pos[%d]=%d", i, pos);
      int nb_var = At_Row(pos, &first);
      mexPrintf(" nb_var=%d\n", nb_var);
      for (int j = 0; j < nb_var; j++)
        {
          mexPrintf("(y[%d]=%f)*(u[%d]=%f)(r=%d, c=%d)\n", index_vara[first->c_index]+cal_y, y[index_vara[first->c_index]+cal_y], first->u_index, u[first->u_index], first->r_index, first->c_index);
          res += y[index_vara[first->c_index]+cal_y]*u[first->u_index];
          first = first->NZE_R_N;
        }
      double tmp_ = res;
      res += u[b[pos]];
      if (abs(res) > epsilon)
        mexPrintf("Error for equation %d => res=%f y[%d]=%f u[b[%d]]=%f somme(y*u)=%f\n", pos, res, pos, y[index_vara[pos]], pos, u[b[pos]], tmp_);
    }
}

mxArray *
dynSparseMatrix::substract_A_B(mxArray *A_m, mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int j = 0; j < n_A; j++)
    for (unsigned int i = 0; i < m_A; i++)
      {
        size_t index = j*m_A+i;
        C_d[index] = A_d[index] - B_d[index];
      }
  return C_m;
}

mxArray *
dynSparseMatrix::Sparse_substract_A_SB(mxArray *A_m, mxArray *B_m)
{
  size_t n_B = mxGetN(B_m);
  size_t m_B = mxGetM(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  size_t total_nze_B = B_j[n_B];
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxDuplicateArray(A_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_B = 0;
  unsigned int B_col = 0;
  while (nze_B < total_nze_B)
    {
      while (nze_B >= (unsigned int) B_j[B_col+1] && (nze_B < total_nze_B))
        B_col++;
      C_d[B_col*m_B+B_i[nze_B]] -= B_d[nze_B];
      nze_B++;
    }
  return C_m;
}

mxArray *
dynSparseMatrix::Sparse_substract_SA_SB(mxArray *A_m, mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  size_t total_nze_B = B_j[n_B];
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int A_col = 0, B_col = 0, C_col = 0;
  C_j[C_col] = 0;
  while (nze_A < total_nze_A || nze_B < total_nze_B)
    {
      while (nze_A >= (unsigned int) A_j[A_col+1] && (nze_A < total_nze_A))
        A_col++;
      size_t A_row = A_i[nze_A];
      while (nze_B >= (unsigned int) B_j[B_col+1] && (nze_B < total_nze_B))
        B_col++;
      size_t B_row = B_i[nze_B];
      if (A_col == B_col)
        {
          if (A_row == B_row && (nze_B < total_nze_B && nze_A < total_nze_A))
            {
              C_d[nze_C] = A_d[nze_A++] - B_d[nze_B++];
              C_i[nze_C] = A_row;
              while (C_col < A_col)
                C_j[++C_col] = nze_C;
              C_j[A_col+1] = nze_C++;
              C_col = A_col;
            }
          else if (A_row < B_row || (nze_B >= total_nze_B && nze_A < total_nze_A))
            {
              C_d[nze_C] = A_d[nze_A++];
              C_i[nze_C] = A_row;
              while (C_col < A_col)
                C_j[++C_col] = nze_C;
              C_j[A_col+1] = nze_C++;
              C_col = A_col;
            }
          else
            {
              C_d[nze_C] = -B_d[nze_B++];
              C_i[nze_C] = B_row;
              while (C_col < B_col)
                C_j[++C_col] = nze_C;
              C_j[B_col+1] = nze_C++;
              C_col = B_col;
            }
        }
      else if (A_col < B_col || (nze_B >= total_nze_B && nze_A < total_nze_A))
        {
          C_d[nze_C] = A_d[nze_A++];
          C_i[nze_C] = A_row;
          while (C_col < A_col)
            C_j[++C_col] = nze_C;
          C_j[A_col+1] = nze_C++;
          C_col = A_col;
        }
      else
        {
          C_d[nze_C] = -B_d[nze_B++];
          C_i[nze_C] = B_row;
          while (C_col < B_col)
            C_j[++C_col] = nze_C;
          C_j[B_col+1] = nze_C++;
          C_col = B_col;
        }
    }
  while (C_col < n_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray *
dynSparseMatrix::mult_SAT_B(mxArray *A_m, mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
  //unsigned int nze_A = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int j = 0; j < (int)n_B; j++)
    {
      for (unsigned int i = 0; i < n_A; i++)
        {
          double sum = 0;
          size_t nze_A = A_j[i];
          while (nze_A < (unsigned int) A_j[i+1])
            {
              size_t i_A = A_i[nze_A];
              sum += A_d[nze_A++] * B_d[i_A];
            }
          C_d[j*n_A+i] = sum;
        }
    }
  return C_m;
}

mxArray *
dynSparseMatrix::Sparse_mult_SAT_B(mxArray *A_m, mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  size_t m_B = mxGetM(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0;
  //unsigned int nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (unsigned int j = 0; j < n_B; j++)
    {
      for (unsigned int i = 0; i < n_A; i++)
        {
          double sum = 0;
          size_t nze_A = A_j[i];
          while (nze_A < (unsigned int) A_j[i+1])
            {
              size_t i_A = A_i[nze_A];
              sum += A_d[nze_A++] * B_d[i_A];
            }
          if (fabs(sum) > 1e-10)
            {
              C_d[nze_C] = sum;
              C_i[nze_C] = i;
              while (C_col < j)
                C_j[++C_col] = nze_C;
              nze_C++;
            }
        }
    }
  while (C_col < m_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray *
dynSparseMatrix::Sparse_mult_SAT_SB(mxArray *A_m, mxArray *B_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  size_t n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  size_t nze_B = 0, nze_C = 0, nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  for (unsigned int j = 0; j < n_B; j++)
    {
      for (unsigned int i = 0; i < n_A; i++)
        {
          double sum = 0;
          nze_B = B_j[j];
          nze_A = A_j[i];
          while (nze_A < (unsigned int) A_j[i+1] && nze_B < (unsigned int) B_j[j+1])
            {
              size_t i_A = A_i[nze_A];
              size_t i_B = B_i[nze_B];
              if (i_A == i_B)
                sum += A_d[nze_A++] * B_d[nze_B++];
              else if (i_A < i_B)
                nze_A++;
              else
                nze_B++;
            }
          if (fabs(sum) > 1e-10)
            {
              C_d[nze_C] = sum;
              C_i[nze_C] = i;
              while (C_col < j)
                C_j[++C_col] = nze_C;
              nze_C++;
            }
        }
    }
  while (C_col < n_B)
    C_j[++C_col] = nze_C;
  mxSetNzmax(C_m, nze_C);
  return C_m;
}

mxArray *
dynSparseMatrix::Sparse_transpose(mxArray *A_m)
{
  size_t n_A = mxGetN(A_m);
  size_t m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  size_t total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  mxArray *C_m = mxCreateSparse(n_A, m_A, total_nze_A, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0, nze_A = 0;
  memset(C_j, 0, m_A);
  map<pair<mwIndex, unsigned int>, double> B2;
  for (unsigned int i = 0; i < n_A; i++)
    {
      while (nze_A < (unsigned int) A_j[i+1])
        {
          C_j[A_i[nze_A]+1]++;
          B2[make_pair(A_i[nze_A], i)] = A_d[nze_A];
          nze_A++;
        }
    }
  for (unsigned int i = 0; i < m_A; i++)
    C_j[i+1] += C_j[i];
  for (map<pair<mwIndex, unsigned int>, double>::const_iterator it = B2.begin(); it != B2.end(); it++)
    {
      C_d[nze_C] = it->second;
      C_i[nze_C++] = it->first.second;
    }
  return C_m;
}


#define sign(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
bool
dynSparseMatrix::mnbrak(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc)
{
  const double GOLD=1.618034;
  const double GLIMIT=100.0;
  const double TINY=1.0e-20;

  double tmp;
  mexPrintf("bracketing *ax=%f, *bx=%f\n",*ax, *bx);
  //mexEvalString("drawnow;");
  double ulim,u,r,q,fu;
  if (!compute_complete(*ax, fa))
    return false;
  if (!compute_complete(*bx, fb))
    return false;
  if (*fb > *fa)
    {
      tmp = *ax;
      *ax = *bx;
      *bx = tmp;

      tmp = *fa;
      *fa = *fb;
      *fb = tmp;
    }
  *cx=(*bx)+GOLD*(*bx-*ax);
  if (!compute_complete(*cx, fc))
    return false;
  while (*fb > *fc)
    {
      r=(*bx-*ax)*(*fb-*fc);
      q=(*bx-*cx)*(*fb-*fa);
      u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
        (2.0*sign(fmax(fabs(q-r),TINY),q-r));
      ulim=(*bx)+GLIMIT*(*cx-*bx);
      if ((*bx-u)*(u-*cx) > 0.0)
        {
          if (!compute_complete(u, &fu))
            return false;
          if (fu < *fc)
            {
              *ax=(*bx);
              *bx=u;
              *fa=(*fb);
              *fb=fu;
              return true;
            }
          else if (fu > *fb)
            {
              *cx=u;
              *fc=fu;
              return true;
            }
          u=(*cx)+GOLD*(*cx-*bx);
          if (!compute_complete(u, &fu))
            return false;
        }
      else if ((*cx-u)*(u-ulim) > 0.0)
        {
          if (!compute_complete(u, &fu))
            return false;
          if (fu < *fc)
            {
              *bx = *cx;
              *cx = u;
              u = *cx+GOLD*(*cx-*bx);
              *fb = *fc;
              *fc = fu;
              if (!compute_complete(u, &fu))
                return false;
            }
        }
      else if ((u-ulim)*(ulim-*cx) >= 0.0)
        {
          u=ulim;
          if (!compute_complete(u, &fu))
            return false;
        }
      else
        {
          u=(*cx)+GOLD*(*cx-*bx);
          if (!compute_complete(u, &fu))
            return false;
        }
      *ax = *bx;
      *bx = *cx;
      *cx = u;
      *fa = *fb;
      *fb = *fc;
      *fc = fu;
    }
  return true;
}

bool
dynSparseMatrix::golden(double ax, double bx, double cx, double tol, double solve_tolf, double *xmin)
{
  const double R=0.61803399;
  const double C=(1.0-R);
  mexPrintf("golden\n");
  //mexEvalString("drawnow;");
  double f1,f2,x0,x1,x2,x3;
  int iter= 0, max_iter= 100;
  x0=ax;
  x3=cx;
  if (fabs(cx-bx) > fabs(bx-ax))
    {
      x1=bx;
      x2=bx+C*(cx-bx);
    }
  else
    {
      x2=bx;
      x1=bx-C*(bx-ax);
    }
  if (!compute_complete(x1, &f1))
    return false;
  if (!compute_complete(x2, &f2))
    return false;
  while ((fabs(x3-x0) > tol*(fabs(x1)+fabs(x2)) && (f1 > solve_tolf && f2 > solve_tolf)) && (iter < max_iter) && (abs(x1 - x2) > 1e-4))
    {
      if (f2 < f1)
        {
          x0 = x1;
          x1 = x2;
          x2 = R*x1+C*x3;
          f1 = f2;
          if (!compute_complete(x2, &f2))
            return false;
        }
      else
        {
          x3 = x2;
          x2 = x1;
          x1 = R*x2+C*x0;
          f2 = f1;
          if (!compute_complete(x1, &f1))
            return false;
        }
      iter++;
    }
  if (f1 < f2)
    {
      *xmin=x1;
      return true;
    }
  else
    {
      *xmin=x2;
      return true;
    }
}

void
dynSparseMatrix::Solve_Matlab_Relaxation(mxArray *A_m, mxArray *b_m, unsigned int Size, double slowc_l, bool is_two_boundaries, int  it_)
{
  mxArray *B1, *C1, *A2, *B2, *A3, *b1, *b2;
  double *b_m_d = mxGetPr(b_m);
  if (!b_m_d)
    {
      ostringstream tmp;
      tmp << " in Solve_Matlab_Relaxation, can't retrieve b_m vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *A_m_i = mxGetIr(A_m);
  if (!A_m_i)
    {
      ostringstream tmp;
      tmp << " in Solve_Matlab_Relaxation, can't allocate A_m_i index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *A_m_j = mxGetJc(A_m);
  if (!A_m_j)
    {
      ostringstream tmp;
      tmp << " in Solve_Matlab_Relaxation, can't allocate A_m_j index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  double *A_m_d = mxGetPr(A_m);
  if (!A_m_d)
    {
      ostringstream tmp;
      tmp << " in Solve_Matlab_Relaxation, can't retrieve A matrix\n";
      throw FatalExceptionHandling(tmp.str());
    }
  size_t max_nze = A_m_j[Size*periods];
  unsigned int nze = 0;
  size_t var = A_m_j[nze];
  B1 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *B1_i = mxGetIr(B1);
  mwIndex *B1_j = mxGetJc(B1);
  double *B1_d = mxGetPr(B1);
  unsigned int B1_nze = 0;
  unsigned int B1_var = 0;
  B1_i[B1_nze] = 0;
  B1_j[B1_var] = 0;
  C1 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *C1_i = mxGetIr(C1);
  mwIndex *C1_j = mxGetJc(C1);
  double *C1_d = mxGetPr(C1);
  unsigned int C1_nze = 0;
  unsigned int C1_var = 0;
  C1_i[C1_nze] = 0;
  C1_j[C1_var] = 0;
  A2 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *A2_i = mxGetIr(A2);
  mwIndex *A2_j = mxGetJc(A2);
  double *A2_d = mxGetPr(A2);
  unsigned int A2_nze = 0;
  unsigned int A2_var = 0;
  A2_i[A2_nze] = 0;
  A2_j[A2_var] = 0;
  B2 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *B2_i = mxGetIr(B2);
  mwIndex *B2_j = mxGetJc(B2);
  double *B2_d = mxGetPr(B2);
  unsigned int B2_nze = 0;
  unsigned int B2_var = 0;
  B2_i[B2_nze] = 0;
  B2_j[B2_var] = 0;
  A3 = mxCreateSparse(Size, Size, Size*Size, mxREAL);
  mwIndex *A3_i = mxGetIr(A3);
  mwIndex *A3_j = mxGetJc(A3);
  double *A3_d = mxGetPr(A3);
  unsigned int A3_nze = 0;
  unsigned int A3_var = 0;
  A3_i[A3_nze] = 0;
  A3_j[A3_var] = 0;
  b1 = mxCreateDoubleMatrix(Size, 1, mxREAL);
  double *b1_d = mxGetPr(b1);
  b2 = mxCreateDoubleMatrix(Size, 1, mxREAL);
  double *b2_d = mxGetPr(b2);
  size_t eq = 0;
  /*B1 C1
    A2 B2
    A3*/
  while (var < 2*Size && nze < max_nze)
    {
      if ((unsigned int) A_m_j[var+1] <= nze)
        {
          if (var < Size)
            b1_d[var] = b_m_d[var];
          else
            b2_d[var - Size] = b_m_d[var];
          var++;
        }
      eq = A_m_i[nze];
      if (var < Size)
        {
          if (eq < Size)
            {
              while (B1_var < var)
                B1_j[++B1_var] = B1_nze;
              B1_i[B1_nze] = eq;
              B1_d[B1_nze] = A_m_d[nze];
              B1_nze++;
            }
          else
            {
              while (A2_var < var)
                A2_j[++A2_var] = A2_nze;
              A2_i[A2_nze] = eq - Size;
              A2_d[A2_nze] = A_m_d[nze];
              A2_nze++;
            }
        }
      else if (var < 2*Size)
        {
          if (eq < Size)
            {
              while (C1_var < var - Size)
                C1_j[++C1_var] = C1_nze;
              C1_i[C1_nze] = eq;
              C1_d[C1_nze] = A_m_d[nze];
              C1_nze++;
            }
          else if (eq < 2*Size)
            {
              while (B2_var < var - Size)
                B2_j[++B2_var] = B2_nze;
              B2_i[B2_nze] = eq - Size;
              B2_d[B2_nze] = A_m_d[nze];
              B2_nze++;
            }
          else
            {
              while (A3_var < var - Size)
                A3_j[++A3_var] = A3_nze;
              A3_i[A3_nze] = eq - 2*Size;
              A3_d[A3_nze] = A_m_d[nze];
              A3_nze++;
            }
        }
      nze++;
    }
  while (B1_var < Size)
    B1_j[++B1_var] = B1_nze;
  while (C1_var < Size)
    C1_j[++C1_var] = C1_nze;
  while (A2_var < Size)
    A2_j[++A2_var] = A2_nze;
  while (B2_var < Size)
    B2_j[++B2_var] = B2_nze;
  while (A3_var < Size)
    A3_j[++A3_var] = A3_nze;
  mxArray *d1 = NULL;
  vector<pair<mxArray *, mxArray *> > triangular_form;
  double sumc = 0, C_sumc = 1000;
  mxArray *B1_inv = NULL;
  mxArray *B1_inv_t = NULL;
  for (int t = 1; t <= periods; t++)
    {
      if (abs(sumc / C_sumc -1) > 1e-10*res1)
        {
          C_sumc = sumc;
          if (B1_inv)
            mxDestroyArray(B1_inv);
          mexCallMATLAB(1, &B1_inv, 1, &B1, "inv");
          mwIndex *B_inv_j = mxGetJc(B1_inv);
          size_t B_inv_nze = B_inv_j[Size];
          double *B_inv_d = mxGetPr(B1_inv);
          sumc = 0;
          for (unsigned int i = 0; i < B_inv_nze; i++)
            sumc += fabs(B_inv_d[i]);
        }
      B1_inv_t = Sparse_transpose(B1_inv);
      mxArray *S1 = Sparse_mult_SAT_SB(B1_inv_t, C1);

      d1 = mult_SAT_B(B1_inv_t, b1);
      if (t < periods)
        //Computation for the next lines
        {
          mxDestroyArray(B1_inv_t);
          mxArray *A2_t = Sparse_transpose(A2);
          mxDestroyArray(A2);

          mxArray *tmp = Sparse_mult_SAT_SB(A2_t, S1);
          mxDestroyArray(B1);
          B1 = Sparse_substract_SA_SB(B2, tmp);
          mxDestroyArray(tmp);

          tmp = mult_SAT_B(A2_t, d1);
          b1 = substract_A_B(b2, tmp);
          mxDestroyArray(tmp);

          triangular_form.push_back(make_pair(S1, d1));
          mxDestroyArray(A2_t);
        }
      A2 = mxDuplicateArray(A3);

      //I  S1
      //0  B1 C1  =>B1 =
      //   A2 B2  => A2 = A3
      //      A3
      C1_nze = B2_nze = A3_nze = 0;
      C1_var = B2_var = A3_var = 0;

      if (nze < max_nze)
        nze--;
      while (var < (t+2)*Size && nze < max_nze)
        {
          if ((unsigned int) A_m_j[var+1] <= nze)
            {
              b2_d[var - (t+1) * Size] = b_m_d[var];
              var++;
            }
          eq = A_m_i[nze];
          if (eq < (t+1) * Size)
            {
              C1_d[C1_nze] = A_m_d[nze];
              C1_nze++;
            }
          else if (eq < (t+2)*Size)
            {
              B2_d[B2_nze] = A_m_d[nze];
              B2_nze++;
            }
          else
            {
              A3_d[A3_nze] = A_m_d[nze];
              A3_nze++;
            }
          nze++;
        }
    }
  double *d1_d = mxGetPr(d1);
  for (unsigned i = 0; i < Size; i++)
    {
      int eq = index_vara[i+Size*(y_kmin+periods-1)];
      double yy = -(d1_d[i] + y[eq]);
      direction[eq] = yy;
      y[eq] += slowc_l * yy;
    }

  pair<mxArray *, mxArray *> tf;
  for (int t = periods-2; t >= 0; t--)
    {
      mxArray *tmp;
      tf = triangular_form.back();
      triangular_form.pop_back();
      mxArray *tf_first_t = Sparse_transpose(tf.first);
      mxDestroyArray(tf.first);
      tmp = mult_SAT_B(tf_first_t, d1);
      d1 = substract_A_B(tf.second, tmp);
      d1_d = mxGetPr(d1);
      mxDestroyArray(tmp);
      for (unsigned i = 0; i < Size; i++)
        {
          int eq = index_vara[i+Size*(y_kmin+t)];
          double yy = -(d1_d[i] + y[eq]);
          direction[eq] = yy;
          y[eq] += slowc_l * yy;
        }
      mxDestroyArray(tf_first_t);
      mxDestroyArray(tf.second);
    }
  mxDestroyArray(B1);
  mxDestroyArray(C1);
  mxDestroyArray(A2);
  mxDestroyArray(B2);
  mxDestroyArray(A3);
  mxDestroyArray(b1);
  mxDestroyArray(b2);
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
}

void
dynSparseMatrix::Solve_Matlab_LU_UMFPack(mxArray *A_m, mxArray *b_m, int Size, double slowc_l, bool is_two_boundaries, int  it_)
{
  size_t n = mxGetM(A_m);
  mxArray *z;
  mxArray *rhs[2];
  rhs[0] = A_m;
  rhs[1] = b_m;
  mexCallMATLAB(1, &z, 2, rhs, "mldivide");
  double *res = mxGetPr(z);
  if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(res[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc_l * yy;
      }
  else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i];
        double yy = -(res[i] + y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc_l * yy;
      }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(z);
}

void
dynSparseMatrix::End_Matlab_LU_UMFPack()
{
  if (Symbolic)
    umfpack_dl_free_symbolic (&Symbolic) ;
  if (Numeric)
    umfpack_dl_free_numeric (&Numeric) ;
}


void
dynSparseMatrix::End_Solver()
{
  if (((stack_solve_algo == 0 || stack_solve_algo == 4) && !steady_state) || (solve_algo == 6 && steady_state))
    End_Matlab_LU_UMFPack();
}

void
dynSparseMatrix::Print_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, int n)
{
  int k = 0;
  for (int i = 0; i < n; i++)
      for (int j = Ap[i]; j < Ap[i+1]; j++)
        mexPrintf("(%d, %d)    %f\n", Ai[j]+1, i+1, Ax[k++]);
}

void
dynSparseMatrix::Solve_LU_UMFPack(SuiteSparse_long *Ap, SuiteSparse_long *Ai, double *Ax, double *b, int n, int Size, double slowc_l, bool is_two_boundaries, int  it_)
{
  SuiteSparse_long status, sys = 0;
#ifndef _MSC_VER
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], res [n];
#else
  double *Control, *Info, *res;
  Control = (double*)mxMalloc(UMFPACK_CONTROL * sizeof(double));
  Info = (double*)mxMalloc(UMFPACK_INFO * sizeof(double));
  res = (double*)mxMalloc(n * sizeof(double));
#endif

  umfpack_dl_defaults(Control);
  Control [UMFPACK_PRL] = 5;
  status = 0;
  if (iter == 0)
    {
      status = umfpack_dl_symbolic(n, n, Ap, Ai, Ax, &Symbolic, Control, Info);
      if (status < 0)
        {
          umfpack_dl_report_info(Control, Info);
          umfpack_dl_report_status(Control, status);
          ostringstream  Error;
          Error << " umfpack_dl_symbolic failed\n";
          throw FatalExceptionHandling(Error.str());
        }
    }
  if (iter > 0)
    umfpack_dl_free_numeric(&Numeric) ;
  status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info);
  if (status < 0)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      ostringstream  Error;
      Error << " umfpack_dl_numeric failed\n";
      throw FatalExceptionHandling(Error.str());
    }
  status = umfpack_dl_solve(sys, Ap, Ai, Ax, res, b, Numeric, Control, Info);
  if (status != UMFPACK_OK)
    {
      umfpack_dl_report_info(Control, Info);
      umfpack_dl_report_status(Control, status);
      ostringstream  Error;
      Error << " umfpack_dl_solve failed\n";
      throw FatalExceptionHandling(Error.str());
    }

  if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(res[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc_l * yy;
      }
  else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i];
        double yy = -(res[i] + y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc_l * yy;
      }

  mxFree(Ap);
  mxFree(Ai);
  mxFree(Ax);
  mxFree(b);
#ifdef _MSC_VER
  mxFree(Control);
  mxFree(Info);
  mxFree(res);
#endif
}


void
dynSparseMatrix::Solve_LU_UMFPack(mxArray *A_m, mxArray *b_m, int Size, double slowc_l, bool is_two_boundaries, int  it_)
{
  SuiteSparse_long n = mxGetM(A_m);

  SuiteSparse_long *Ap = (SuiteSparse_long*)mxGetJc (A_m);

  SuiteSparse_long *Ai = (SuiteSparse_long*)mxGetIr(A_m);
  double*  Ax = mxGetPr(A_m);
  double*  B  = mxGetPr(b_m);
  SuiteSparse_long status, sys = 0;
#ifndef _MSC_VER
  double Control [UMFPACK_CONTROL], Info [UMFPACK_INFO], res [n];
#else
  double *Control, *Info, *res;
  Control = (double*)mxMalloc(UMFPACK_CONTROL * sizeof(double));
  Info = (double*)mxMalloc(UMFPACK_INFO * sizeof(double));
  res = (double*)mxMalloc(n * sizeof(double));
#endif
  void *Symbolic, *Numeric ;
  umfpack_dl_defaults (Control) ;

  status = umfpack_dl_symbolic (n, n, Ap, Ai, Ax, &Symbolic, Control, Info) ;
  if (status != UMFPACK_OK)
    umfpack_dl_report_info ((double*) NULL, Info) ;

  status = umfpack_dl_numeric (Ap, Ai, Ax, Symbolic, &Numeric, Control, Info) ;
  if (status != UMFPACK_OK)
    umfpack_dl_report_info ((double*) NULL, Info) ;

  status = umfpack_dl_solve (sys, Ap, Ai, Ax, res, B, Numeric, Control, Info) ;
  if (status != UMFPACK_OK)
    umfpack_dl_report_info ((double*) NULL, Info) ;
  //double *res = mxGetPr(z);
  if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(res[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc_l * yy;
      }
  else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i];
        double yy = -(res[i] + y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc_l * yy;
      }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
#ifdef _MSC_VER
  mxFree(Control);
  mxFree(Info);
  mxFree(res);
#endif

}


#ifdef CUDA
void
printM(int n,double *Ax, int* Ap, int* Ai,  cusparseMatDescr_t descrA, cusparseHandle_t cusparse_handle)
{
  //cudaError_t cuda_error;
  //cusparseStatus_t cusparse_status;
  double * A_dense;
  cudaChk(cudaMalloc((void**) &A_dense, n * n *sizeof(double)), "A_dense cudaMalloc has failed\n");


  cusparseChk(cusparseDcsr2dense(cusparse_handle, n, n, descrA,
                                 Ax, Ap,Ai, A_dense, n), "cusparseDcsr2dense has failed\n");
  double *A_dense_hoste = (double*)mxMalloc(n * n * sizeof(double));
  cudaChk(cudaMemcpy(A_dense_hoste, A_dense, n * n * sizeof(double),cudaMemcpyDeviceToHost), " cudaMemcpy(A_dense_hoste, A_dense) has failed\n");
  mexPrintf("----------------------\n");
  mexPrintf("FillMode=%d, IndexBase=%d, MatType=%d, DiagType=%d\n",cusparseGetMatFillMode(descrA), cusparseGetMatIndexBase(descrA), cusparseGetMatType(descrA), cusparseGetMatDiagType(descrA));
  //mexEvalString("drawnow;");
  for (int i = 0; i < n ; i++)
    {
      for (int j = 0; j < n; j++)
        mexPrintf("%-6.3f ",A_dense_hoste[i + j * n]);
      mexPrintf("\n");
    }
  mxFree(A_dense_hoste);
  cudaChk(cudaFree(A_dense), "cudaFree(A_dense) has failed\n");
}



void
dynSparseMatrix::Solve_CUDA_BiCGStab_Free(double* tmp_vect_host, double* p, double* r, double* v, double* s, double* t, double* y_, double* z, double* tmp_,
                                       int* Ai, double* Ax, int* Ap, double* x0, double* b, double* A_tild, int* A_tild_i, int* A_tild_p/*, double* Lx, int* Li, int* Lp,
                                       double* Ux, int* Ui, int* Up, int* device_n*/, cusparseSolveAnalysisInfo_t infoL, cusparseSolveAnalysisInfo_t infoU,
                                       cusparseMatDescr_t descrL, cusparseMatDescr_t descrU, int preconditioner)
{
  //cudaError_t cuda_error;
  //cusparseStatus_t cusparse_status;
  mxFree(tmp_vect_host);
  cudaChk(cudaFree(p), "  in Solve_Cuda_BiCGStab, can't free p\n");
  cudaChk(cudaFree(r), "  in Solve_Cuda_BiCGStab, can't free r\n");
  cudaChk(cudaFree(v), "  in Solve_Cuda_BiCGStab, can't free v\n");
  cudaChk(cudaFree(s), "  in Solve_Cuda_BiCGStab, can't free s\n");
  cudaChk(cudaFree(t), "  in Solve_Cuda_BiCGStab, can't free t\n");
  cudaChk(cudaFree(y_), "  in Solve_Cuda_BiCGStab, can't free y_\n");
  cudaChk(cudaFree(z), "  in Solve_Cuda_BiCGStab, can't free z\n");
  cudaChk(cudaFree(tmp_), "  in Solve_Cuda_BiCGStab, can't free tmp_\n");
  cudaChk(cudaFree(Ai), "  in Solve_Cuda_BiCGStab, can't free Ai\n");
  cudaChk(cudaFree(Ax), "  in Solve_Cuda_BiCGStab, can't free Ax\n");
  cudaChk(cudaFree(Ap), "  in Solve_Cuda_BiCGStab, can't free Ap\n");
  cudaChk(cudaFree(x0), "  in Solve_Cuda_BiCGStab, can't free x0\n");
  cudaChk(cudaFree(b), "  in Solve_Cuda_BiCGStab, can't free b\n");
  /*if (preconditioner == 0)
    {*/
      cudaChk(cudaFree(A_tild), "  in Solve_Cuda_BiCGStab, can't free A_tild (1)\n");
      cudaChk(cudaFree(A_tild_i), "  in Solve_Cuda_BiCGStab, can't free A_tild_i (1)\n");
      cudaChk(cudaFree(A_tild_p), "  in Solve_Cuda_BiCGStab, can't free A_tild_p (1)\n");
    /*}
  else
    {
      cudaChk(cudaFree(Lx), "  in Solve_Cuda_BiCGStab, can't free Lx\n");
      cudaChk(cudaFree(Li), "  in Solve_Cuda_BiCGStab, can't free Li\n");
      cudaChk(cudaFree(Lp), "  in Solve_Cuda_BiCGStab, can't free Lp\n");
      cudaChk(cudaFree(Ux), "  in Solve_Cuda_BiCGStab, can't free Ux\n");
      cudaChk(cudaFree(Ui), "  in Solve_Cuda_BiCGStab, can't free Ui\n");
      cudaChk(cudaFree(Up), "  in Solve_Cuda_BiCGStab, can't free Up\n");
    }*/
  //cudaChk(cudaFree(device_n), "  in Solve_Cuda_BiCGStab, can't free device_n\n");
  if (preconditioner == 1 || preconditioner == 2 || preconditioner == 3)
    {
      cusparseChk(cusparseDestroySolveAnalysisInfo(infoL),
                  "  in Solve_Cuda_BiCGStab, cusparseDestroySolveAnalysisInfo has failed for infoL\n");
      cusparseChk(cusparseDestroySolveAnalysisInfo(infoU),
                  "  in Solve_Cuda_BiCGStab, cusparseDestroySolveAnalysisInfo has failed for infoU\n");
    }
  cusparseChk(cusparseDestroyMatDescr(descrL),
              " in Solve_Cuda_BiCGStab, matrix descriptor destruction failed for descrL\n");
  cusparseChk(cusparseDestroyMatDescr(descrU),
              " in Solve_Cuda_BiCGStab, matrix descriptor destruction failed for descrU\n");
}
#endif

void
Solve(double* Ax, int* Ap, int* Ai, double *b, int n, bool Lower, double *x)
{
  if (Lower)
    {
      for (int i = 0; i < n; i++)
        {
          double sum = 0;
          for(int j = Ap[i]; j < Ap[i+1]; j++)
            {
              int k = Ai[j];
              if (k < i)
                sum += x[k] * Ax[j];
            }
          x[i] = b[i] - sum;
        }
    }
  else
    {
      for (int i = n-1 ; i >= 0; i--)
        {
          double sum = 0, mul = 1;
          for(int j = Ap[i]; j < Ap[i+1]; j++)
            {
              int k = Ai[j];
              if (k > i)
                sum += x[k] * Ax[j];
              else if (k == i)
                mul = Ax[j];
            }
          x[i] = (b[i] - sum) / mul;
        }
    }
}

void
Check(int n, double* Ax, int* Ap, int* Ai, double* b, double *x, bool Lower)
{
  if (Lower)
    {
      for (int i = 0; i < n; i++)
        {
          double sum = 0;
          for(int j = Ap[i]; j < Ap[i+1]; j++)
            {
              int k = Ai[j];
              if (k < i)
                sum += x[k] * Ax[j];
            }
          double err =  b[i] - sum - x[i];
          if (abs(err) > 1e-10)
            mexPrintf("error at i=%d\n",i);
        }
    }
  else
    {
      for (int i = n-1 ; i >= 0; i--)
        {
          double sum = 0;
          for(int j = Ap[i]; j < Ap[i+1]; j++)
            {
              int k = Ai[j];
              if (k >= i)
                sum += x[k] * Ax[j];
            }
          double err =  b[i] - sum;
          if (abs(err) > 1e-10)
            mexPrintf("error at i=%d\n",i);
        }
    }
}

#ifdef CUDA
int
dynSparseMatrix::Solve_CUDA_BiCGStab(int *Ap, int *Ai, double *Ax, int *Ap_tild, int *Ai_tild, double *A_tild, double *b, double *x0, int n, int Size, double slowc_l, bool is_two_boundaries,
                                  int  it_, int nnz, int nnz_tild, int preconditioner, int max_iterations, int block)
{
  cusparseSolveAnalysisInfo_t info, infoL, infoU;
  cusparseMatDescr_t descrL, descrU;
  const double tol = 1.0e-6;//1.0e-6;
  const double eps = 1.0e-16;
  double *p, *r, *r0, *v, *s, *t, *y_, *z, *tmp_;
  int *A_tild_i, *A_tild_p;
  double *Qx;
  int *Qi, *Qj;
  double *Px;
  int *Pi, *Pj;
  int Q_nnz, P_nnz;
  int W_nnz;
  double bnorm;
  double tmp1, tmp2;
  int refinement_needed = 0, stagnation = 0;
  int max_refinement = min(min(int(floor(double(n)/50)),10),n-max_iterations), max_stagnation = 3;
  int nblocks = ceil(double(n) / double(1024));
  int n_threads;
  if (nblocks == 0)
    n_threads = n;
  else
    n_threads = 1024;
  int periods = n / Size;

  double * tmp_vect_host = (double*)mxMalloc(n * sizeof(double));

  cublasChk(cublasDnrm2(cublas_handle, n,b, 1, &bnorm),
            "  in Solve_Cuda_BiCGStab, cublasDnrm2(b) has failed\n");

  double tolb = tol * bnorm;

  if (bnorm == 0.0)
    {
      // if b = 0 the A.x = 0 => x = 0
      cudaChk(cudaFree(Ai), "  in Solve_Cuda_BiCGStab, can't free Ai\n");
      cudaChk(cudaFree(Ax), "  in Solve_Cuda_BiCGStab, can't free Ax\n");
      cudaChk(cudaFree(Ap), "  in Solve_Cuda_BiCGStab, can't free Ap\n");
      if (preconditioner == 3)
        {
          cudaChk(cudaFree(Ai_tild), "  in Solve_Cuda_BiCGStab, can't free Ai_tild\n");
          cudaChk(cudaFree(Ap_tild), "  in Solve_Cuda_BiCGStab, can't free Ap_tild\n");
        }
      cudaChk(cudaFree(A_tild), "  in Solve_Cuda_BiCGStab, can't free A_tild\n");
      cudaChk(cudaFree(x0), "  in Solve_Cuda_BiCGStab, can't free x0\n");
      cudaChk(cudaFree(b), "  in Solve_Cuda_BiCGStab, can't free b\n");
      if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -y[eq];
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i];
            double yy = -y[eq+it_*y_size];
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc * yy;
          }
      return 0;
    }

  int iteration = 0;
  bool convergence = false;
  double zeros = 0.0, one = 1.0, m_one = -1.0;

  cudaChk(cudaMalloc((void**)&tmp_, n * sizeof(double)), "  in Solve_Cuda_Sparse, can't allocate tmp_ on the graphic card\n");

  cudaChk(cudaMalloc((void**)&r, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate r on the graphic card\n");

  cudaChk(cudaMemcpy(r, b, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy r = b has failed\n");

  //r = b - A * x0
  cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE, n,
                                   n, nnz, &m_one,
                                   CUDA_descr, Ax,
                                   Ap, Ai,
                                   x0, &one,
                                   r), "in Solve_Cuda_BiCGStab, cusparseDcsrmv A * x0 has failed");

  cudaChk(cudaMemcpy(tmp_vect_host, r, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = p_tild has failed\n");
  /*mexPrintf("r\n");
  for (int i = 0; i < n; i++)
    mexPrintf("%f\n",tmp_vect_host[i]);*/

  cudaChk(cudaMalloc((void**)&r0, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate r0 on the graphic card\n");
  cudaChk(cudaMemcpy(r0, r, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy r0 = r has failed\n");

  cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                        r, 1,
                        &tmp1),
            "  in Solve_Cuda_BiCGStab, cublasDnrm2(r) has failed\n");
  double conv_criteria = tmp1;

  convergence = conv_criteria < tolb;
  if (convergence)
    {
      /* the initial value (x0) is solution of A x = b*/
      cudaChk(cudaFree(Ai), "  in Solve_Cuda_BiCGStab, can't free Ai\n");
      cudaChk(cudaFree(Ax), "  in Solve_Cuda_BiCGStab, can't free Ax\n");
      cudaChk(cudaFree(Ap), "  in Solve_Cuda_BiCGStab, can't free Ap\n");
      if (preconditioner == 3)
        {
          cudaChk(cudaFree(Ai_tild), "  in Solve_Cuda_BiCGStab, can't free Ai_tild\n");
          cudaChk(cudaFree(Ap_tild), "  in Solve_Cuda_BiCGStab, can't free Ap_tild\n");
        }
      cudaChk(cudaFree(A_tild), "  in Solve_Cuda_BiCGStab, can't free A_tild\n");
      cudaChk(cudaFree(x0), "  in Solve_Cuda_BiCGStab, can't free x0\n");
      cudaChk(cudaFree(b), "  in Solve_Cuda_BiCGStab, can't free b\n");
      return 0;
    }


  if (preconditioner == 0)
    {
      //Apply the Jacobi preconditioner
      /*VecDiv<<<nblocks, n_threads>>>(r_, A_tild, z_, n);
      cuda_error = cudaMemcpy(zz_, z_, n * sizeof(double), cudaMemcpyDeviceToDevice);*/
    }
  else if (preconditioner == 1)
    {
      //Apply an incomplete LU decomposition of A as preconditioner
      cusparseChk(cusparseCreateSolveAnalysisInfo(&info), "  in Solve_Cuda_BiCGStab, cusparseCreateSolveAnalysisInfo for info has failed\n");

      cusparseChk(cusparseDcsrsv_analysis(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n, nnz, CUDA_descr,
                                          A_tild, Ap, Ai,
                                          info),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrsm_analysis(info) has failed\n");

      cusparseChk(cusparseDcsrilu0(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                   n, CUDA_descr,
                                   A_tild, Ap, Ai,
                                   info),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrilu0 has failed\n");

      //Make a copy of the indexes in A_tild_i and A_tild_p to use it the Bicgstab algorithm
      cudaChk(cudaMalloc((void**)&A_tild_i, nnz * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild_i on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_i, Ai, nnz * sizeof(int), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_i = Ai has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_p, (n + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild_p on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_p, Ap, (n + 1) * sizeof(int), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_p = Ap has failed\n");
    }
  else if (preconditioner == 2)
    {
      //Because the Jacobian matrix A is store in CSC format in matlab
      // we have to transpose it to get a CSR format used by CUDA
      mwIndex* Awi, *Awp;
      double* A_tild_host = (double*)mxMalloc(nnz*sizeof(double));
      Awi = (mwIndex*)mxMalloc(nnz * sizeof(mwIndex));
      Awp = (mwIndex*)mxMalloc((n + 1) * sizeof(mwIndex));
      int* Aii = (int*)mxMalloc(nnz * sizeof(int));
      int* Aip = (int*)mxMalloc((n + 1) * sizeof(int));
      cudaChk(cudaMemcpy(A_tild_host, A_tild, nnz*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_host = A_tild has failed\n");
      cudaChk(cudaMemcpy(Aii, Ai, nnz*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aii = Ai has failed\n");
      cudaChk(cudaMemcpy(Aip, Ap, (n+1)*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aip = Ai has failed\n");
      for (int i = 0; i < nnz; i++)
        Awi[i] = Aii[i];
      for (int i = 0; i < n + 1; i++)
        Awp[i] = Aip[i];
      mxFree(Aii);
      mxFree(Aip);
      mxArray * At_m = mxCreateSparse(n,n,nnz,mxREAL);
      mxSetIr(At_m, Awi);
      mxSetJc(At_m, Awp);
      mxSetPr(At_m, A_tild_host);
      mxArray *A_m;
      mexCallMATLAB(1, &A_m, 1, &At_m, "transpose");
      mxDestroyArray(At_m);

      /*mexPrintf("A_m\n");
      mexCallMATLAB(0, NULL, 1, &A_m, "disp_dense");*/
      /*mxFree(Awi);
      mxFree(Awp);*/

      /*[L1, U1] = ilu(g1a=;*/
      const char *field_names[] = {"type", "droptol", "milu", "udiag", "thresh"};
      const int type = 0;
      const int droptol = 1;
      const int milu = 2;
      const int udiag = 3;
      const int thresh = 4;
      mwSize dims[1] = {(mwSize)1 };
      mxArray *Setup = mxCreateStructArray(1, dims, 5, field_names);
      mxSetFieldByNumber(Setup, 0, type, mxCreateString("ilutp"));
      //mxSetFieldByNumber(Setup, 0, type, mxCreateString("nofill"));
      mxSetFieldByNumber(Setup, 0, droptol, mxCreateDoubleScalar(lu_inc_tol));
      mxSetFieldByNumber(Setup, 0, milu, mxCreateString("off"));
      mxSetFieldByNumber(Setup, 0, udiag, mxCreateDoubleScalar(0));
      mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(1));
      //mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(1));
      mxArray *lhs0[2];
      mxArray *rhs0[2];
      rhs0[0] = A_m;
      rhs0[1] = Setup;
      ostringstream tmp;
      if (mexCallMATLAB(2, lhs0, 2, rhs0, "ilu"))
        {
          tmp << " In BiCGStab, the incomplet LU decomposition (ilu) ahs failed.\n";
          throw FatalExceptionHandling(tmp.str());
        }
      mxDestroyArray(Setup);


 /*     //ILUT preconditionner computed by Matlab (todo: in futur version of cuda replace it by a new equivalent cuda function)
      const char *field_names[] = {"type", "droptol", "milu", "udiag", "thresh"};
      const int type = 0;
      const int droptol = 1;
      const int milu = 2;
      const int udiag = 3;
      const int thresh = 4;
      mwSize dims[1] = {(mwSize)1 };
      mxArray *Setup = mxCreateStructArray(1, dims, 5, field_names);
      mxSetFieldByNumber(Setup, 0, type, mxCreateString("ilutp"));
      mxSetFieldByNumber(Setup, 0, droptol, mxCreateDoubleScalar(lu_inc_tol));
      mxSetFieldByNumber(Setup, 0, milu, mxCreateString("off"));
      mxSetFieldByNumber(Setup, 0, udiag, mxCreateDoubleScalar(0));
      mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(0));
      mxArray *lhs0[2], *rhs0[2];
      rhs0[0] = A_m;
      rhs0[1] = Setup;
      mexCallMATLAB(1, lhs0, 2, rhs0, "ilu");
*/
      // To store the resultng matrix in a CSR format we have to transpose it
      mxArray *Wt = lhs0[0];
      mwIndex* Wtj = mxGetJc(Wt);
      nnz = Wtj[n];
      mxArray* W;
      mexCallMATLAB(1, &W, 1, &Wt, "transpose");
      mxDestroyArray(Wt);
      double* pW = mxGetPr(W);
      mwIndex* Wi = mxGetIr(W);
      mwIndex* Wp = mxGetJc(W);
      int *Wii = (int*)mxMalloc(nnz * sizeof(int));
      int *Wip = (int*)mxMalloc((n + 1) * sizeof(int));
      for (int i = 0; i < nnz; i++)
        Wii[i] = Wi[i];
      for (int i = 0; i < n + 1; i++)
        Wip[i] = Wp[i];

      //mxFree(A_tild_host);

      cudaChk(cudaFree(A_tild), "cudaFree(A_tild) has failed\n");

      cudaChk(cudaMalloc((void**)&A_tild, nnz * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild, pW, nnz * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild = pW has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_i, nnz * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Ai on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_i, Wii, nnz * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_i = A_tild_i_host has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_p, (n + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild_p on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_p, Wip, (n + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_p = A_tild_j_host has failed\n");
      /*mxFree(pW);
      mxFree(Wi);
      mxFree(Wj);*/
      mxDestroyArray(W);
      mxFree(Wii);
      mxFree(Wip);
    }
  else if (preconditioner == 3)
    {
      mwIndex* Aowi, *Aowp;
      double* A_host = (double*)mxMalloc(nnz*sizeof(double));
      Aowi = (mwIndex*)mxMalloc(nnz * sizeof(mwIndex));
      Aowp = (mwIndex*)mxMalloc((n + 1) * sizeof(mwIndex));
      int* Aoii = (int*)mxMalloc(nnz * sizeof(int));
      int* Aoip = (int*)mxMalloc((n + 1) * sizeof(int));
      cudaChk(cudaMemcpy(A_host, Ax, nnz*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_host = A_tild has failed\n");
      cudaChk(cudaMemcpy(Aoii, Ai, nnz*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aii = Ai_tild has failed\n");
      cudaChk(cudaMemcpy(Aoip, Ap, (n+1)*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aip = Ap_tild has failed\n");
      for (int i = 0; i < nnz; i++)
        Aowi[i] = Aoii[i];
      for (int i = 0; i < n + 1; i++)
        Aowp[i] = Aoip[i];
      mxFree(Aoii);
      mxFree(Aoip);
      mxArray * Ao_m = mxCreateSparse(n,n,nnz,mxREAL);
      mxSetIr(Ao_m, Aowi);
      mxSetJc(Ao_m, Aowp);
      mxSetPr(Ao_m, A_host);
      /*mexPrintf("A_m\n");
      mxArray *Aoo;
      mexCallMATLAB(1, &Aoo, 1, &Ao_m, "transpose");
      mexCallMATLAB(0, NULL, 1, &Aoo, "disp_dense");
      mxDestroyArray(Ao_m);
      mxDestroyArray(Aoo);*/

      //Because the Jacobian matrix A is store in CSC format in matlab
      // we have to transpose it to get a CSR format used by CUDA
      mwIndex* Awi, *Awp;
      double* A_tild_host = (double*)mxMalloc(nnz_tild*sizeof(double));
      Awi = (mwIndex*)mxMalloc(nnz_tild * sizeof(mwIndex));
      Awp = (mwIndex*)mxMalloc((Size + 1) * sizeof(mwIndex));
      int* Aii = (int*)mxMalloc(nnz_tild * sizeof(int));
      int* Aip = (int*)mxMalloc((Size + 1) * sizeof(int));
      cudaChk(cudaMemcpy(A_tild_host, A_tild, nnz_tild*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_host = A_tild has failed\n");
      cudaChk(cudaMemcpy(Aii, Ai_tild, nnz_tild*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aii = Ai_tild has failed\n");
      cudaChk(cudaMemcpy(Aip, Ap_tild, (Size+1)*sizeof(int), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy Aip = Ap_tild has failed\n");
      for (int i = 0; i < nnz_tild; i++)
        Awi[i] = Aii[i];
      for (int i = 0; i < Size + 1; i++)
        Awp[i] = Aip[i];
      /*for (int i = 0; i < nnz_tild; i++)
        mexPrintf("%20.17f\n",A_tild_host[i]);*/
      mxFree(Aii);
      mxFree(Aip);
      mxArray * At_m = mxCreateSparse(Size,Size,nnz_tild,mxREAL);
      mxSetIr(At_m, Awi);
      mxSetJc(At_m, Awp);
      mxSetPr(At_m, A_tild_host);
      mxArray *A_m;
      mexCallMATLAB(1, &A_m, 1, &At_m, "transpose");
      /*mexPrintf("A_tild_m\n");
      mexCallMATLAB(0, NULL, 1, &A_m, "disp_dense");*/
      mxDestroyArray(At_m);
      mxArray *P, *Q, *L, *U;
      mxArray *lhs0[4];
      mexCallMATLAB(4, lhs0, 1, &A_m, "lu");

      mxArray *P0, *Q0, *L0, *U0;
      L0 = lhs0[0];
      U0 = lhs0[1];
      P0 = lhs0[2];
      Q0 = lhs0[3];
      mexCallMATLAB(1, &P, 1, &P0, "transpose");
      mexCallMATLAB(1, &Q, 1, &Q0, "transpose");
      mexCallMATLAB(1, &L, 1, &L0, "transpose");
      mexCallMATLAB(1, &U, 1, &U0, "transpose");
      mxDestroyArray(P0);
      mxDestroyArray(Q0);
      mxDestroyArray(L0);
      mxDestroyArray(U0);
      /*L = lhs0[0];
      U = lhs0[1];
      P = lhs0[2];
      Q = lhs0[3];*/

      /*mexPrintf("L\n");
      mexCallMATLAB(0, NULL, 1, &L, "disp_dense");

      mexPrintf("U\n");
      mexCallMATLAB(0, NULL, 1, &U, "disp_dense");

      mexPrintf("P\n");
      mexCallMATLAB(0, NULL, 1, &P, "disp_dense");

      mexPrintf("Q\n");
      mexCallMATLAB(0, NULL, 1, &Q, "disp_dense");*/

      mwIndex* Qiw_host = mxGetIr(Q);
      mwIndex* Qjw_host = mxGetJc(Q);
      double*  Qx_host = mxGetPr(Q);
      Q_nnz = Qjw_host[Size];
      mexPrintf("Q_nnz=%d\n",Q_nnz);
      int *Qi_host = (int*)mxMalloc(Q_nnz * periods * sizeof(int));
      double *Q_x_host = (double*)mxMalloc(Q_nnz * periods * sizeof(double));
      int *Qj_host = (int*)mxMalloc((n + 1) * sizeof(int));
      for (int t = 0; t < periods; t++)
        {
          for (int i = 0; i < Q_nnz; i++)
            {
              Qi_host[i + t * Q_nnz] = Qiw_host[i] + t * Size;
              Q_x_host[i + t * Q_nnz] = Qx_host[i];
            }
          for (int i = 0; i < Size; i++)
            {
              Qj_host[i + t * Size] = Qjw_host[i] + t * Q_nnz;
            }
        }
      Qj_host[periods * Size] = periods * Q_nnz;


      /*mwIndex *Qtiw_host  = (mwIndex*) mxMalloc(Q_nnz * periods * sizeof(mwIndex));
      double *Qt_x_host = (double*)mxMalloc(Q_nnz * periods * sizeof(double));
      mwIndex *Qtjw_host = (mwIndex*)mxMalloc((n + 1) * sizeof(mwIndex));
      mexPrintf("n = %d\n",n);
      for (int i = 0; i < n + 1; i++)
        Qtjw_host[i] = Qj_host[i];
      for (int i = 0; i < Q_nnz * periods; i++)
        {
          Qtiw_host[i] = Qi_host[i];
          Qt_x_host[i] = Q_x_host[i];
        }
      mxArray* Qt_m = mxCreateSparse(n,n,Q_nnz * periods,mxREAL);
      mxSetIr(Qt_m, Qtiw_host);
      mxSetJc(Qt_m, Qtjw_host);
      mxSetPr(Qt_m, Qt_x_host);
      mexPrintf("Qt_m\n");
      mexCallMATLAB(0, NULL, 1, &Qt_m, "disp_dense");*/


      /*mexPrintf("Qtjw_host[periods * Size=%d]=%d\n", periods * Size, Qtjw_host[periods * Size]);
      for (int i = 0; i < n; i++)
        for (int j = Qtjw_host[i]; j < Qtjw_host[i+1]; j++)
           mexPrintf("(i=%d, j=%d) = %f\n", i, Qtiw_host[j], Qt_x_host[j]);*/
      //mxDestroyArray(Qt_m);


      cudaChk(cudaMalloc((void**)&Qx, Q_nnz * periods * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate Qx on the graphic card\n");
      cudaChk(cudaMemcpy(Qx, Q_x_host, Q_nnz * periods * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Qx = Qx_host has failed\n");
      cudaChk(cudaMalloc((void**)&Qi, Q_nnz * periods * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Qi on the graphic card\n");
      cudaChk(cudaMemcpy(Qi, Qi_host, Q_nnz * periods * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Qi = Qi_host has failed\n");
      cudaChk(cudaMalloc((void**)&Qj, (Size * periods + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Qj on the graphic card\n");
      cudaChk(cudaMemcpy(Qj, Qj_host, (Size * periods + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Qj = Qj_host has failed\n");
      mxFree(Qi_host);
      mxFree(Qj_host);
      mxFree(Q_x_host);
      mxDestroyArray(Q);


      mwIndex* Piw_host = mxGetIr(P);
      mwIndex* Pjw_host = mxGetJc(P);
      double*  Px_host = mxGetPr(P);
      P_nnz = Pjw_host[Size];
      int *Pi_host = (int*)mxMalloc(P_nnz * periods * sizeof(int));
      double *P_x_host = (double*)mxMalloc(P_nnz * periods * sizeof(double));
      int *Pj_host = (int*)mxMalloc((n + 1) * sizeof(int));
      for (int t = 0; t < periods; t++)
        {
          for (int i = 0; i < P_nnz; i++)
            {
              Pi_host[i + t * P_nnz] = Piw_host[i] + t * Size;
              P_x_host[i + t * P_nnz] = Px_host[i];
            }
          for (int i = 0; i < Size; i++)
            Pj_host[i + t * Size] = Pjw_host[i] + t * P_nnz;
        }
      Pj_host[periods * Size] = periods * P_nnz;

      /*mwIndex *Ptiw_host  = (mwIndex*) mxMalloc(P_nnz * periods * sizeof(mwIndex));
      double *Pt_x_host = (double*)mxMalloc(P_nnz * periods * sizeof(double));
      mwIndex *Ptjw_host = (mwIndex*)mxMalloc((n + 1) * sizeof(mwIndex));
      for (int i = 0; i < n + 1; i++)
        Ptjw_host[i] = Pj_host[i];
      for (int i = 0; i < P_nnz * periods; i++)
        {
          Ptiw_host[i] = Pi_host[i];
          Pt_x_host[i] = P_x_host[i];
        }
      mxArray* Pt_m = mxCreateSparse(n,n,P_nnz * periods,mxREAL);
      mxSetIr(Pt_m, Ptiw_host);
      mxSetJc(Pt_m, Ptjw_host);
      mxSetPr(Pt_m, Pt_x_host);
      mexPrintf("Pt_m\n");
      mexCallMATLAB(0, NULL, 1, &Pt_m, "disp_dense");
      mxDestroyArray(Pt_m);*/


      cudaChk(cudaMalloc((void**)&Px, P_nnz * periods * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate Px on the graphic card\n");
      cudaChk(cudaMemcpy(Px, P_x_host, P_nnz * periods * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Px = Px_host has failed\n");
      cudaChk(cudaMalloc((void**)&Pi, P_nnz * periods * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pi on the graphic card\n");
      cudaChk(cudaMemcpy(Pi, Pi_host, P_nnz * periods * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Pi = Pi_host has failed\n");
      cudaChk(cudaMalloc((void**)&Pj, (Size * periods + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pj on the graphic card\n");
      cudaChk(cudaMemcpy(Pj, Pj_host, (Size * periods + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Pj = Pj_host has failed\n");
      mxFree(Pi_host);
      mxFree(Pj_host);
      mxFree(P_x_host);
      mxDestroyArray(P);

      /*mwIndex* Piw_host = mxGetIr(P);
      mwIndex* Pjw_host = mxGetJc(P);
      double*  Px_host = mxGetPr(P);
      P_nnz = Pjw_host[Size];
      int *Pi_host = (int*)mxMalloc(P_nnz * sizeof(int));
      int *Pj_host = (int*)mxMalloc((Size + 1) * sizeof(int));
      for (int i = 0; i < P_nnz; i++)
        Pi_host[i] = Piw_host[i];
      for (int i = 0; i < Size + 1; i++)
        Pj_host[i] = Pjw_host[i];

      cudaChk(cudaMalloc((void**)&Px, P_nnz * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate Px on the graphic card\n");
      cudaChk(cudaMemcpy(Px, Px_host, P_nnz * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Px = Px_host has failed\n");
      cudaChk(cudaMalloc((void**)&Pi, P_nnz * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pi on the graphic card\n");
      cudaChk(cudaMemcpy(Pi, Pi_host, P_nnz * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Pi = Pi_host has failed\n");
      cudaChk(cudaMalloc((void**)&Pj, (Size + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pj on the graphic card\n");
      cudaChk(cudaMemcpy(Pj, Pj_host, (Size + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy Pj = Pj_host has failed\n");
      mxFree(Pi_host);
      mxFree(Pj_host);
      mxDestroyArray(P);*/

      /*mexPrintf("L\n");
      mexCallMATLAB(0, NULL, 1, &L, "disp_dense");

      mexPrintf("U\n");
      mexCallMATLAB(0, NULL, 1, &U, "disp_dense");*/

      mwIndex* Liw_host = mxGetIr(L);
      mwIndex* Ljw_host = mxGetJc(L);
      double*  Lx_host = mxGetPr(L);
      int L_nnz = Ljw_host[Size];

      mwIndex* Uiw_host = mxGetIr(U);
      mwIndex* Ujw_host = mxGetJc(U);
      double*  Ux_host = mxGetPr(U);
      int U_nnz = Ujw_host[Size];

      double *pW = (double*)mxMalloc((L_nnz + U_nnz - Size) * periods * sizeof(double));
      int *Wi = (int*)mxMalloc((L_nnz + U_nnz - Size) * periods * sizeof(int));
      int *Wj = (int*)mxMalloc((n + 1) * sizeof(int));
      Wj[0] = 0;
      W_nnz = 0;
      for (int t = 0; t < periods; t++)
        for (int i = 0; i < Size ; i++)
          {
            for (mwIndex l  = Ujw_host[i]; l < Ujw_host[i+1]; l++)
              {
                Wi[W_nnz] = Uiw_host[l] + t * Size;
                pW[W_nnz] = Ux_host[l];
                //mexPrintf("Wj[%d] = %d, Wi[%d] = Uiw_host[%d] + t * Size = %d, pW[%d]=%f\n", i + t * Size, Wj[i + t * Size], W_nnz, l, Uiw_host[l] + t * Size, W_nnz, Ux_host[l]);
                W_nnz++;
              }
            for (mwIndex l  = Ljw_host[i]; l < Ljw_host[i+1]; l++)
              {
                if (Liw_host[l] > i)
                  {
                    Wi[W_nnz] = Liw_host[l] + t * Size;
                    pW[W_nnz] = Lx_host[l];
                    //mexPrintf("Wj[%d] = %d, Wi[%d] = Liw_host[%d] + t * Size = %d, pW[%d]=%f\n", i  + t * Size, Wj[i + t * Size], W_nnz, l, Liw_host[l] + t * Size, W_nnz, Lx_host[l]);
                    W_nnz++;
                  }
              }
            Wj[i + 1 + t * Size] = W_nnz;
          }
      //mexPrintf("Wj[%d] = %d, n=%d\n", periods * Size, Wj[periods * Size], n);
      cudaChk(cudaMalloc((void**)&A_tild, W_nnz * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate Px on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild, pW, W_nnz * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild = pW has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_i, W_nnz * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pi on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_i, Wi, W_nnz * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_i = Wi has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_p, (n + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Pj on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_p, Wj, (n + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_p = Wj has failed\n");

      /*mwIndex *Wwi = (mwIndex*)mxMalloc(W_nnz * sizeof(mwIndex));
      mwIndex *Wwj = (mwIndex*)mxMalloc((n + 1) * sizeof(mwIndex));
      for (int i = 0; i < W_nnz; i++)
        Wwi[i] = Wi[i];
      for (int i = 0; i < n + 1; i++)
        Wwj[i] = Wj[i];
      mxFree(Wi);
      mxFree(Wj);
      mxArray* Ao_tild = mxCreateSparse(n,n,W_nnz,mxREAL);
      mxSetIr(Ao_tild, Wwi);
      mxSetJc(Ao_tild, Wwj);
      mxSetPr(Ao_tild, pW);
      mexPrintf("Ao_tild\n");
      mexCallMATLAB(0, NULL, 1, &Ao_tild, "disp_dense");
      mxDestroyArray(Ao_tild);*/


      /*ostringstream tmp;
      tmp << "debugging";
      mexWarnMsgTxt(tmp.str().c_str());
      return 4;*/

      /* /**Apply the permutation matrices (P and Q) to the b vector of system to solve :
       b_tild = P-1 . b  = P' . b */
      /*cudaChk(cudaMalloc((void**)&b_tild, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate b_tild on the graphic card\n");
      cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
                                 n, n, nnz, &one, CUDA_descr,
                                 Px, Pj, Pi,
                                 b, &zeros,
                                 b_tild),
                  "  in Solve_Cuda_BiCGStab, b_tild = cusparseDcsrmv(P', b) has failed\n");

      cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_TRANSPOSE,
                                 n, n, nnz, &one, CUDA_descr,
                                 Px, Pj, Pi,
                                 b, &zeros,
                                 b),
                  "  in Solve_Cuda_BiCGStab, b = cusparseDcsrmv(P', b) has failed\n");
      */
      /*mexPrintf("Wt = lu(A_m)\n");
      mexCallMATLAB(0, NULL, 1, &Wt, "disp_dense");*/
      /*ostringstream tmp;
      tmp << "debugging";
      mexWarnMsgTxt(tmp.str().c_str());
      return 4;*/
      // To store the resultng matrix in a CSR format we have to transpose it
      /*mwIndex* Wtj = mxGetJc(Wt);
      nnz = Wtj[n];
      mxArray* W;
      mexCallMATLAB(1, &W, 1, &Wt, "transpose");
      mxDestroyArray(Wt);
      pW = mxGetPr(W);
      Wwi = mxGetIr(W);
      mwIndex* Wp = mxGetJc(W);
      int *Wii = (int*)mxMalloc(nnz * sizeof(int));
      int *Wip = (int*)mxMalloc((n + 1) * sizeof(int));
      for (int i = 0; i < nnz; i++)
        Wii[i] = Wi[i];
      for (int i = 0; i < n + 1; i++)
        Wip[i] = Wp[i];

      //mxFree(A_tild_host);

      cudaChk(cudaFree(Ai_tild), "  in Solve_Cuda_BiCGStab, cudaFree(Ai_tild) has failed\n");
      cudaChk(cudaFree(Ap_tild), "  in Solve_Cuda_BiCGStab, cudaFree(Ap_tild) has failed\n");
      cudaChk(cudaFree(A_tild), "  in Solve_Cuda_BiCGStab, cudaFree(A_tild) has failed\n");

      cudaChk(cudaMalloc((void**)&A_tild, nnz * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild, pW, nnz * sizeof(double), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild = pW has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_i, nnz * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate Ai on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_i, Wii, nnz * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_i = A_tild_i_host has failed\n");
      cudaChk(cudaMalloc((void**)&A_tild_p, (n + 1) * sizeof(int)), "  in Solve_Cuda_BiCGStab, can't allocate A_tild_p on the graphic card\n");
      cudaChk(cudaMemcpy(A_tild_p, Wip, (n + 1) * sizeof(int), cudaMemcpyHostToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy A_tild_p = A_tild_j_host has failed\n");
      mxDestroyArray(W);
      mxFree(Wii);
      mxFree(Wip);*/
    }
  if (preconditioner == 1 || preconditioner == 2 || preconditioner == 3)
    {
      cusparseChk(cusparseCreateMatDescr(&descrL),
                  "  in Solve_Cuda_BiCGStab, cusparseCreateMatDescr has failed for descrL\n");
      cusparseChk(cusparseSetMatIndexBase(descrL, CUSPARSE_INDEX_BASE_ZERO),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatIndexBase has failed for descrL\n");
      cusparseChk(cusparseSetMatType(descrL, CUSPARSE_MATRIX_TYPE_GENERAL),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatType has failed for descrL\n");
      cusparseChk(cusparseSetMatFillMode(descrL, CUSPARSE_FILL_MODE_LOWER),
                  "  in Solve_Cuda_BiCGStab, cusparseSetFillMod has failed for descrL\n");
      cusparseChk(cusparseSetMatDiagType(descrL, CUSPARSE_DIAG_TYPE_UNIT),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatDiagType has failed for descrL\n");

      cusparseChk(cusparseCreateMatDescr(&descrU),
                  "  in Solve_Cuda_BiCGStab, cusparseCreateMatDescr has failed for descrU\n");
      cusparseChk(cusparseSetMatIndexBase(descrU, CUSPARSE_INDEX_BASE_ZERO),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatIndexBase has failed for descrU\n");
      cusparseChk(cusparseSetMatType(descrU, CUSPARSE_MATRIX_TYPE_GENERAL),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatType has failed for descrU\n");
      cusparseChk(cusparseSetMatFillMode(descrU, CUSPARSE_FILL_MODE_UPPER),
                  "  in Solve_Cuda_BiCGStab, cusparseSetFillMod has failed for descrU\n");
      cusparseChk(cusparseSetMatDiagType(descrU, CUSPARSE_DIAG_TYPE_NON_UNIT),
                  "  in Solve_Cuda_BiCGStab, cusparseSetMatDiagType has failed for descrU\n");

      int host_nnz_tild;
      if  (preconditioner == 3)
        host_nnz_tild = W_nnz;
      else
        host_nnz_tild = nnz;

      if (preconditioner == 1)
        cusparseChk(cusparseDestroySolveAnalysisInfo(info),
                    "  in Solve_Cuda_BiCGStab, cusparseDestroySolveAnalysisInfo has failed for info\n");

      cusparseChk(cusparseCreateSolveAnalysisInfo(&infoL),
                  "  in Solve_Cuda_BiCGStab, cusparseCreateSolveAnalysisInfo has failed for infoL\n");
      cusparseChk(cusparseDcsrsv_analysis(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n, host_nnz_tild, descrL,
                                          A_tild, A_tild_p, A_tild_i,
                                          infoL),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrsm_analysis for infoL has failed\n");

      cusparseChk(cusparseCreateSolveAnalysisInfo(&infoU),
                  "  in Solve_Cuda_BiCGStab, cusparseCreateSolveAnalysisInfo has failed for infoU\n");
      cusparseChk(cusparseDcsrsv_analysis(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                          n, host_nnz_tild, descrU,
                                          A_tild, A_tild_p, A_tild_i,
                                          infoU),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrsm_analysis for infoU has failed\n");
    }

  cudaChk(cudaMalloc((void**)&v, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate v on the graphic card\n");
  cudaChk(cudaMalloc((void**)&p, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate p on the graphic card\n");
  //cudaChk(cudaMemset(p, 0, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, cudaMemset p = 0 has failed\n");
  cudaChk(cudaMalloc((void**)&s, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate s on the graphic card\n");
  cudaChk(cudaMalloc((void**)&t, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate t on the graphic card\n");
  cudaChk(cudaMalloc((void**)&y_, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate y_ on the graphic card\n");
  cudaChk(cudaMalloc((void**)&z, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate z on the graphic card\n");

  double rho = 1.0, alpha = 1.0, omega = 1.0;


  //residual = P*B*Q - L*U;
  //norm(Z,1) should be close to 0


  while (iteration < 50/*max_iterations*/ && !convergence)
    {
      double rho_prev = rho;
      /**store in s previous value of r*/
      cudaChk(cudaMemcpy(s, r, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy s = r has failed\n");

      /**rho = r0 . r*/
      cublasChk(cublasDdot(cublas_handle, n, // numerator
                           r0, 1,
                           r, 1,
                           &rho),
                "  in Solve_Cuda_BiCGStab, rho = cublasDdot(r0, r) has failed\n");

      mexPrintf("rho=%f\n",rho);

      double beta;

      if (iteration == 0)
        {
          cudaChk(cudaMemcpy(p, r, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy p = r has failed\n");
        }
      else
        {
          /**beta = (rho / rho_prev) . (alpha / omega);*/
          beta = rho / rho_prev * alpha / omega;

          /**p = r + beta * (p - omega * v)*/
          // tmp_ = p - omega * v
          VecAdd<<<nblocks, n_threads>>>(tmp_, p, -omega, v, n);
          //p = r + beta * tmp_
          VecAdd<<<nblocks, n_threads>>>(p, r, beta, tmp_, n);
        }

      /**y_ solution of A_tild * y_ = p <=> L . U . y_ = p*/
      //  L tmp_ = p => tmp_ = L^-1 p, with tmp_ = U . y_

      if (preconditioner == 3)
        {
          double *p_tild;

          cudaChk(cudaMemcpy(tmp_vect_host, p, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = p has failed\n");
          /*mexPrintf("p\n");
          for (int i = 0; i < n; i++)
             mexPrintf("%f\n",tmp_vect_host[i]);*/

          cudaChk(cudaMalloc((void**)&p_tild, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate b_tild on the graphic card\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, P_nnz * periods, &one, CUDA_descr,
                                     Px, Pj, Pi,
                                     p, &zeros,
                                     p_tild),
                      "  in Solve_Cuda_BiCGStab, p_tild = cusparseDcsrmv(P', p) has failed\n");

          /*mexPrintf("P\n");
          printM(n, Px, Pj, Pi, CUDA_descr, cusparse_handle);*/

          cudaChk(cudaMemcpy(tmp_vect_host, p_tild, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = p_tild has failed\n");
          /*mexPrintf("p_tild\n");
          for (int i = 0; i < n; i++)
             mexPrintf("%f\n",tmp_vect_host[i]);*/

          cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                           n, &one,
                                           descrL,
                                           A_tild, A_tild_p, A_tild_i,
                                           infoL, p_tild,
                                           tmp_),
                      "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for L . tmp_ = p_tild has failed\n");
          cudaChk(cudaFree(p_tild), "  in Solve_Cuda_BiCGStab, can't free p_tild\n");

          cudaChk(cudaMemcpy(tmp_vect_host, tmp_, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = v has failed\n");
          /*mexPrintf("tmp_\n");
          for (int i = 0; i < n; i++)
             mexPrintf("%f\n",tmp_vect_host[i]);*/
        }
      else
        cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                         n, &one,
                                         descrL,
                                         A_tild, A_tild_p, A_tild_i,
                                         infoL, p,
                                         tmp_),
                    "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for L . tmp_ = p has failed\n");

      //  U . y_ = L^-1 p <=> U . y_ = tmp_ => y_ = U^-1 L^-1 p
      cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       n, &one,
                                       descrU,
                                       A_tild, A_tild_p, A_tild_i,
                                       infoU, tmp_,
                                       y_),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for U . y_ = tmp_ has failed\n");

      /*cudaChk(cudaMemcpy(tmp_vect_host, y_, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = v has failed\n");
      mexPrintf("y_\n");
      for (int i = 0; i < n; i++)
        mexPrintf("%f\n",tmp_vect_host[i]);*/

      if (preconditioner == 3)
        {
          double *y_tild;
          cudaChk(cudaMalloc((void**)&y_tild, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate b_tild on the graphic card\n");
          cudaChk(cudaMemcpy(y_tild, y_, n  * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy y_tild = y_ has failed\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, Q_nnz * periods, &one, CUDA_descr,
                                     Qx, Qj, Qi,
                                     y_tild, &zeros,
                                     y_),
                      "  in Solve_Cuda_BiCGStab, y_ = cusparseDcsrmv(Q', y_tild) has failed\n");
          cudaChk(cudaFree(y_tild), "  in Solve_Cuda_BiCGStab, can't free y_tild\n");
        }
      /*cudaChk(cudaMemcpy(tmp_vect_host, y_, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = v has failed\n");
      mexPrintf("y_\n");
      for (int i = 0; i < n; i++)
        mexPrintf("%f\n",tmp_vect_host[i]);*/
      /**v = A*y_*/
      cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 n, n, nnz, &one, CUDA_descr,
                                 Ax, Ap, Ai,
                                 y_, &zeros,
                                 v),
                  "  in Solve_Cuda_BiCGStab, v = cusparseDcsrmv(A, y_) has failed\n");
      cudaChk(cudaMemcpy(tmp_vect_host, v, n*sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = v has failed\n");
      /*mexPrintf("v\n");
      for (int i = 0; i < n; i++)
        mexPrintf("%f\n",tmp_vect_host[i]);*/



      /**alpha = rho / (rr0 . v) with rr0 = r0*/
      cublasChk(cublasDdot(cublas_handle, n, // numerator
                           r0, 1,
                           v, 1,
                           &tmp1),
                "  in Solve_Cuda_BiCGStab, cublasDdot(r0, v) has failed\n");

      alpha = rho / tmp1;
      mexPrintf("rho = %f, tmp1 = %f\n", rho, tmp1);
      mexPrintf("alpha = %f\n", alpha);

      if (alpha == 0 || isinf(alpha) || isnan(alpha))
        {
          Solve_CUDA_BiCGStab_Free(tmp_vect_host, p, r, v, s, t, y_, z, tmp_, Ai, Ax, Ap, x0, b, A_tild, A_tild_i, A_tild_p, infoL, infoU, descrL, descrU, preconditioner);
          ostringstream tmp;
          tmp << "one of the scalar quantities (alpha=" << alpha << ") calculated during BICGSTAB became too small or too large to continue computing, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
          return 4;
        }

      /** Check for potential stagnation*/
      cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                            y_, 1,
                            &tmp1),
                "  in Solve_Cuda_BiCGStab, cublasDnrm2(y_) has failed\n");
      cublasChk(cublasDnrm2(cublas_handle, n, // denominator
                            x0, 1,
                            &tmp2),
                "  in Solve_Cuda_BiCGStab, cublasDnrm2(y_) has failed\n");
      mexPrintf("abs(alpha)*tmp1  = %f, alpha = %f, tmp1 = %f, tmp2 = %f, eps = %f\n",abs(alpha)*tmp1 , alpha, tmp1, tmp2, eps);
      if (abs(alpha)*tmp1  < eps * tmp2)
        stagnation++;
      else
        stagnation = 0;

      /**x = x + alpha * y_*/
      VecInc<<<nblocks, n_threads>>>(x0, alpha, y_, n);

      /**s = r_prev - alpha *v with r_prev = s*/
      VecInc<<<nblocks, n_threads>>>(s, -alpha, v, n);

      /**Has BiCGStab converged?*/
      cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                            s, 1,
                            &tmp1),
                "  in Solve_Cuda_BiCGStab, cublasDnrm2(s) has failed\n");
      conv_criteria = tmp1;
      mexPrintf("conv_criteria = %f, tolb = %f\n", conv_criteria, tolb);
      convergence = conv_criteria < tolb;

      if (convergence || stagnation >= max_stagnation || refinement_needed)
        {
          /**s = b - A * x0*/
          cudaChk(cudaMemcpy(s, b, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy s = b has failed\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, nnz, &m_one, CUDA_descr,
                                     Ax, Ap, Ai,
                                     x0, &one,
                                     s),
                      "  in Solve_Cuda_BiCGStab, s = b - cusparseDcsrmv(A, x0) has failed\n");
          cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                                s, 1,
                                &tmp1),
                    "  in Solve_Cuda_BiCGStab, cublasDnrm2(s) has failed\n");
          conv_criteria = tmp1;
          convergence = conv_criteria < tolb;
          if (convergence)
            {
              break;
            }
          else
            {
              if (stagnation >= max_stagnation && refinement_needed == 0)
                stagnation = 0;
              refinement_needed++;
              if (refinement_needed > max_refinement)
                {
                  Solve_CUDA_BiCGStab_Free(tmp_vect_host, p, r, v, s, t, y_, z, tmp_, Ai, Ax, Ap, x0, b, A_tild, A_tild_i, A_tild_p, infoL, infoU, descrL, descrU, preconditioner);
                  ostringstream tmp;
                  tmp << "Error in bytecode: BiCGStab stagnated (Two consecutive iterates were the same.), in block " << block+1;
                  mexWarnMsgTxt(tmp.str().c_str());
                  return 3;
                }
            }
        }

      /**z solution of A_tild * z = s*/
      //  L tmp_ = s => tmp_ = L^-1 s, with tmp_ = U . z
      if (preconditioner == 3)
        {
          double *s_tild;
          cudaChk(cudaMalloc((void**)&s_tild, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate b_tild on the graphic card\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, P_nnz * periods, &one, CUDA_descr,
                                     Px, Pj, Pi,
                                     s, &zeros,
                                     s_tild),
                      "  in Solve_Cuda_BiCGStab, s_tild = cusparseDcsrmv(P', s) has failed\n");
          cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                           n, &one,
                                           descrL,
                                           A_tild, A_tild_p, A_tild_i,
                                           infoL, s_tild,
                                           tmp_),
                      "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for L . tmp_ = s_tild has failed\n");
          cudaChk(cudaFree(s_tild), "  in Solve_Cuda_BiCGStab, can't free s_tild\n");
        }
      else
        cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                         n, &one,
                                         descrL,
                                         //Lx, Lp, Li,
                                         A_tild, A_tild_p, A_tild_i,
                                         infoL, s,
                                         tmp_),
                    "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for L . tmp_ = s has failed\n");
      //  U . z = L^-1 s <=> U . z = tmp_ => z = U^-1 L^-1 s
      cusparseChk(cusparseDcsrsv_solve(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                       n, &one,
                                       descrU,
                                       //Ux, Up, Ui,
                                       A_tild, A_tild_p, A_tild_i,
                                       infoU, tmp_,
                                       z),
                  "  in Solve_Cuda_BiCGStab, cusparseDcsrsv_solve for U . z = tmp_ has failed\n");
      if (preconditioner == 3)
        {
          double *z_tild;
          cudaChk(cudaMalloc((void**)&z_tild, n * sizeof(double)), "  in Solve_Cuda_BiCGStab, can't allocate z_tild on the graphic card\n");
          cudaChk(cudaMemcpy(z_tild, z, n  * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy z_tild = z has failed\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, Q_nnz * periods, &one, CUDA_descr,
                                     Qx, Qj, Qi,
                                     z_tild, &zeros,
                                     z),
                      "  in Solve_Cuda_BiCGStab, z = cusparseDcsrmv(Q, z_tild) has failed\n");
          cudaChk(cudaFree(z_tild), "  in Solve_Cuda_BiCGStab, can't free x_tild\n");
        }
      /**t = A * z*/
      cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                 n, n, nnz, &one, CUDA_descr,
                                 Ax, Ap, Ai,
                                 z, &zeros,
                                 t),
                  "  in Solve_Cuda_BiCGStab, t = cusparseDcsrmv(A, z) has failed\n");

      /** omega = (t' s) / (t' t)*/
      cublasChk(cublasDdot(cublas_handle, n, // numerator
                           t, 1,
                           s, 1,
                           &tmp1),
                "  in Solve_Cuda_BiCGStab, cublasDdot(t, s) has failed\n");

      cublasChk(cublasDdot(cublas_handle, n, // numerator
                           t, 1,
                           t, 1,
                           &tmp2),
                "  in Solve_Cuda_BiCGStab, cublasDdot(t, t) has failed\n");

      omega = tmp1 / tmp2;

      if (omega == 0 || isinf(omega) || isnan(omega))
        {
          Solve_CUDA_BiCGStab_Free(tmp_vect_host, p, r, v, s, t, y_, z, tmp_, Ai, Ax, Ap, x0, b, A_tild, A_tild_i, A_tild_p, infoL, infoU, descrL, descrU, preconditioner);
          ostringstream tmp;
          mexEvalString("diary off;");
          tmp << "one of the scalar quantities (omega=" << omega << ") calculated during BICGSTAB became too small or too large to continue computing, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
          return 4;
        }

      /**x = x +  omega * z*/
      VecInc<<<nblocks, n_threads>>>(x0, omega, z, n);

      /**r = s - omega * t*/
      VecAdd<<<nblocks, n_threads>>>(r, s, -omega, t, n);

      /**Has BiCGStab converged?*/
      cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                            r, 1,
                            &tmp1),
                "  in Solve_Cuda_BiCGStab, cublasDnrm2(r) has failed\n");
      conv_criteria = tmp1;

      convergence = conv_criteria < tolb;

      if (convergence || stagnation >= max_stagnation || refinement_needed)
        {
          /**r = b - A * x0*/
          cudaChk(cudaMemcpy(r, b, n * sizeof(double), cudaMemcpyDeviceToDevice), "  in Solve_Cuda_BiCGStab, cudaMemcpy r = b has failed\n");
          cusparseChk(cusparseDcsrmv(cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                                     n, n, nnz, &m_one, CUDA_descr,
                                     Ax, Ap, Ai,
                                     x0, &one,
                                     r),
                      "  in Solve_Cuda_BiCGStab, r = b - cusparseDcsrmv(A, x0) has failed\n");
          cublasChk(cublasDnrm2(cublas_handle, n, // numerator
                                r, 1,
                                &tmp1),
                    "  in Solve_Cuda_BiCGStab, cublasDnrm2(r) has failed\n");
          conv_criteria = tmp1;
          convergence = conv_criteria < tolb;
          if (convergence)
            {
              mexPrintf("convergence achieved\n");
              break;
            }
          else
            {
              if (stagnation >= max_stagnation && refinement_needed == 0)
                stagnation = 0;
              refinement_needed++;
              if (refinement_needed > max_refinement)
                {
                  Solve_CUDA_BiCGStab_Free(tmp_vect_host, p, r, v, s, t, y_, z, tmp_, Ai, Ax, Ap, x0, b, A_tild, A_tild_i, A_tild_p, /*Lx, Li, Lp, Ux, Ui, Up, device_n, */infoL, infoU, descrL, descrU, preconditioner);
                  ostringstream tmp;
                  mexEvalString("diary off;");
                  tmp << "Error in bytecode: BiCGStab stagnated (Two consecutive iterates were the same.), in block " << block+1;
                  mexWarnMsgTxt(tmp.str().c_str());
                  return 3;
                }
            }
        }

      iteration++;
    }
  cudaChk(cudaMemcpy(tmp_vect_host, x0, n * sizeof(double), cudaMemcpyDeviceToHost), "  in Solve_Cuda_BiCGStab, cudaMemcpy tmp_vect_host = x0 has failed\n");

  if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(tmp_vect_host[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc * yy;
      }
  else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i];
        double yy = -(tmp_vect_host[i] + y[eq+it_*y_size]);
        direction[eq] = yy;
        y[eq+it_*y_size] += slowc * yy;
      }
  Solve_CUDA_BiCGStab_Free(tmp_vect_host, p, r, v, s, t, y_, z, tmp_, Ai, Ax, Ap, x0, b, A_tild, A_tild_i, A_tild_p, infoL, infoU, descrL, descrU, preconditioner);

  if (iteration >= max_iterations)
    {
      ostringstream tmp;
      mexEvalString("diary off;");
      tmp << "Error in bytecode: No convergence inside BiCGStab, in block " << block+1;
      mexWarnMsgTxt(tmp.str().c_str());
      return 1;
    }
  else
    return 0;
}
#endif

void
dynSparseMatrix::Solve_Matlab_GMRES(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m)
{
#ifdef OCTAVE_MEX_FILE
  ostringstream tmp;
  if (steady_state)
    tmp << " GMRES method is not implemented in Octave. You cannot use solve_algo=7, change solve_algo.\n";
  else
    tmp << " GMRES method is not implemented in Octave. You cannot use stack_solve_algo=2, change stack_solve_algo.\n";
  throw FatalExceptionHandling(tmp.str());
#endif
  size_t n = mxGetM(A_m);
  const char *field_names[] = {"droptol", "type"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, 2, field_names);
  mxSetFieldByNumber(Setup, 0, 0, mxCreateDoubleScalar(lu_inc_tol));
  mxSetFieldByNumber(Setup, 0, 1, mxCreateString("ilutp"));
  mxArray *lhs0[2];
  mxArray *rhs0[2];
  rhs0[0] = A_m;
  rhs0[1] = Setup;
  if (mexCallMATLAB(2, lhs0, 2, rhs0, "ilu"))
    throw FatalExceptionHandling("In GMRES, the incomplet LU decomposition (ilu) ahs failed.");
  mxArray *L1 = lhs0[0];
  mxArray *U1 = lhs0[1];
  /*[za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);*/
  mxArray *rhs[8];
  rhs[0] = A_m;
  rhs[1] = b_m;
  rhs[2] = mxCreateDoubleScalar(Size);
  rhs[3] = mxCreateDoubleScalar(1e-6);
  rhs[4] = mxCreateDoubleScalar((double)n);
  rhs[5] = L1;
  rhs[6] = U1;
  rhs[7] = x0_m;
  mxArray *lhs[2];
  mexCallMATLAB(2, lhs, 8, rhs, "gmres");
  mxArray *z = lhs[0];
  mxArray *flag = lhs[1];
  double *flag1 = mxGetPr(flag);
  mxDestroyArray(rhs0[1]);
  mxDestroyArray(rhs[2]);
  mxDestroyArray(rhs[3]);
  mxDestroyArray(rhs[4]);
  mxDestroyArray(rhs[5]);
  mxDestroyArray(rhs[6]);
  if (*flag1 > 0)
    {
      ostringstream tmp;
      if (*flag1 == 1)
        {
          tmp << "Error in bytecode: No convergence inside GMRES, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      else if (*flag1 == 2)
        {
          tmp << "Error in bytecode: Preconditioner is ill-conditioned, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      else if (*flag1 == 3)
        {
          tmp << "Error in bytecode: GMRES stagnated (Two consecutive iterates were the same.), in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      lu_inc_tol /= 10;
    }
  else
    {
      double *res = mxGetPr(z);
      if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(z);
  mxDestroyArray(flag);
}

void
dynSparseMatrix::Solve_Matlab_BiCGStab(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m, int preconditioner)
{
  /* precond = 0  => Jacobi
     precond = 1  => Incomplet LU decomposition*/
  size_t n = mxGetM(A_m);
  mxArray *L1, *U1, *Diag;

  mxArray *rhs0[4];
  if (preconditioner == 0)
    {
      mxArray *lhs0[1];
      rhs0[0] = A_m;
      rhs0[1] = mxCreateDoubleScalar(0);
      mexCallMATLAB(1, lhs0, 2, rhs0, "spdiags");
      mxArray* tmp = lhs0[0];
      double* tmp_val = mxGetPr(tmp);
      Diag = mxCreateSparse(n, n, n, mxREAL);
      mwIndex *Diag_i = mxGetIr(Diag);
      mwIndex *Diag_j = mxGetJc(Diag);
      double *Diag_val = mxGetPr(Diag);
      for (size_t i = 0; i < n; i++)
        {
          Diag_val[i] = tmp_val[i];
          Diag_j[i] = i;
          Diag_i[i] = i;
        }
      Diag_j[n] = n;
    }
  else if (preconditioner == 1)
    {
      /*[L1, U1] = ilu(g1a=;*/
      const char *field_names[] = {"type", "droptol", "milu", "udiag", "thresh"};
      const int type = 0;
      const int droptol = 1;
      const int milu = 2;
      const int udiag = 3;
      const int thresh = 4;
      mwSize dims[1] = {(mwSize)1 };
      mxArray *Setup = mxCreateStructArray(1, dims, 5, field_names);
      mxSetFieldByNumber(Setup, 0, type, mxCreateString("ilutp"));
      mxSetFieldByNumber(Setup, 0, droptol, mxCreateDoubleScalar(lu_inc_tol));
      mxSetFieldByNumber(Setup, 0, milu, mxCreateString("off"));
      mxSetFieldByNumber(Setup, 0, udiag, mxCreateDoubleScalar(0));
      mxSetFieldByNumber(Setup, 0, thresh, mxCreateDoubleScalar(1));
      mxArray *lhs0[2];
      mxArray *rhs0[2];
      rhs0[0] = A_m;
      rhs0[1] = Setup;
      if (mexCallMATLAB(2, lhs0, 2, rhs0, "ilu"))
        {
          ostringstream tmp;
          tmp  << " In BiCGStab, the incomplet LU decomposition (ilu) ahs failed.\n";
          throw FatalExceptionHandling(tmp.str());
        }
      L1 = lhs0[0];
      U1 = lhs0[1];
      mxDestroyArray(Setup);
    }
  double flags = 2;
  mxArray *z;
  if (steady_state)  /*Octave BicStab algorihtm involves a 0 division in case of a preconditionner equal to the LU decomposition of A matrix*/
    {
      mxArray *res = mult_SAT_B(Sparse_transpose(A_m), x0_m);
      double *resid = mxGetPr(res);
      double *b = mxGetPr(b_m);
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < (int)n; i++)
        resid[i] = b[i] - resid[i];
      mxArray *rhs[2];
      mxArray *lhs[1];
      rhs[0] = L1;
      rhs[1] = res;
      mexCallMATLAB(1, lhs, 2, rhs, "mldivide");
      rhs[0] = U1;
      rhs[1] = lhs[0];
      mexCallMATLAB(1, lhs, 2, rhs, "mldivide");
      z = lhs[0];
      double *phat = mxGetPr(z);
      double *x0 = mxGetPr(x0_m);
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < (int)n; i++)
        phat[i] = x0[i] + phat[i];

      /*Check the solution*/
      res = mult_SAT_B(Sparse_transpose(A_m), z);
      resid = mxGetPr(res);
      double cum_abs = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:cum_abs)
#endif
      for (int i = 0; i < (int)n; i++)
        {
          resid[i] = b[i] - resid[i];
          cum_abs += fabs(resid[i]);
        }
      if (cum_abs > 1e-7)
        flags = 2;
      else
        flags = 0;
      mxDestroyArray(res);
    }
  //else

  if (flags == 2)
    {
      if (preconditioner == 0)
        {
          /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
          mxArray *rhs[5];
          rhs[0] = A_m;
          rhs[1] = b_m;
          rhs[2] = mxCreateDoubleScalar(1e-6);
          rhs[3] = mxCreateDoubleScalar((double)n);
          rhs[4] = Diag;
          //rhs[5] = x0_m;
          mxArray *lhs[2];
          mexCallMATLAB(2, lhs, 5, rhs, "bicgstab");
          z = lhs[0];
          mxArray *flag = lhs[1];
          double *flag1 = mxGetPr(flag);
          flags = flag1[0];
          mxDestroyArray(flag);
          mxDestroyArray(rhs[2]);
          mxDestroyArray(rhs[3]);
          mxDestroyArray(rhs[4]);
        }
      else if (preconditioner == 1)
        {
          /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
          mxArray *rhs[7];
          rhs[0] = A_m;
          rhs[1] = b_m;
          rhs[2] = mxCreateDoubleScalar(1e-6);
          rhs[3] = mxCreateDoubleScalar((double)n);
          rhs[4] = L1;
          rhs[5] = U1;
          rhs[6] = x0_m;
          mxArray *lhs[2];
          mexCallMATLAB(2, lhs, 7, rhs, "bicgstab");
          z = lhs[0];
          mxArray *flag = lhs[1];
          double *flag1 = mxGetPr(flag);
          flags = flag1[0];
          mxDestroyArray(flag);
          mxDestroyArray(rhs[2]);
          mxDestroyArray(rhs[3]);
          mxDestroyArray(rhs[4]);
          mxDestroyArray(rhs[5]);
        }
    }


  if (flags > 0)
    {
      ostringstream tmp;
      if (flags == 1)
        {
          tmp << "Error in bytecode: No convergence inside BiCGStab, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      else if (flags == 2)
        {
          tmp << "Error in bytecode: Preconditioner is ill-conditioned, in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      else if (flags == 3)
        {
          tmp << "Error in bytecode: BiCGStab stagnated (Two consecutive iterates were the same.), in block " << block+1;
          mexWarnMsgTxt(tmp.str().c_str());
        }
      lu_inc_tol /= 10;
    }
  else
    {
      double *res = mxGetPr(z);
      if (is_two_boundaries)
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i];
            double yy = -(res[i] + y[eq+it_*y_size]);
            direction[eq] = yy;
            y[eq+it_*y_size] += slowc * yy;
          }
    }
  mxDestroyArray(A_m);
  mxDestroyArray(b_m);
  mxDestroyArray(z);
}

void
dynSparseMatrix::Singular_display(int block, int Size)
{
  bool zero_solution;
  Simple_Init(Size, IM_i, zero_solution);
  NonZeroElem *first;
  mxArray *rhs[1];
  rhs[0] = mxCreateDoubleMatrix(Size, Size, mxREAL);
  double *pind;
  pind = mxGetPr(rhs[0]);
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int j = 0; j < Size * Size; j++)
    pind[j] = 0.0;
  for (int ii = 0; ii < Size; ii++)
    {
      int nb_eq = At_Col(ii, &first);
      for (int j = 0; j < nb_eq; j++)
        {
          int k = first->u_index;
          int jj = first->r_index;
          pind[ii * Size + jj ] = u[k];
          first = first->NZE_C_N;
        }
    }
  mxArray *lhs[3];
  mexCallMATLAB(3, lhs, 1, rhs, "svd");
  mxArray* SVD_u = lhs[0];
  mxArray* SVD_s = lhs[1];
  //mxArray* SVD_v = lhs[2];
  double *SVD_ps = mxGetPr(SVD_s);
  double *SVD_pu = mxGetPr(SVD_u);
  for (int i = 0; i < Size; i++)
    {
      if (abs(SVD_ps[i * (1 + Size)]) < 1e-12)
        {
          mexPrintf(" The following equations form a linear combination:\n    ");
          double max_u = 0;
          for (int j = 0; j < Size; j++)
            if (abs(SVD_pu[j + i * Size]) > abs(max_u))
              max_u = SVD_pu[j + i * Size];
          vector<int> equ_list;
          for (int j = 0; j < Size; j++)
            {
              double rr = SVD_pu[j + i * Size] / max_u;
              if ( rr < -1e-10)
                {
                  equ_list.push_back(j);
                  if (rr != -1)
                    mexPrintf(" - %3.2f*Dequ_%d_dy",abs(rr),j+1);
                  else
                    mexPrintf(" - Dequ_%d_dy",j+1);
                }
              else if (rr > 1e-10)
                {
                  equ_list.push_back(j);
                  if (j > 0)
                    if (rr != 1)
                      mexPrintf(" + %3.2f*Dequ_%d_dy",rr,j+1);
                    else
                      mexPrintf(" + Dequ_%d_dy",j+1);
                  else if (rr != 1)
                    mexPrintf(" %3.2f*Dequ_%d_dy",rr,j+1);
                  else
                    mexPrintf(" Dequ_%d_dy",j+1);
                }
            }
          mexPrintf(" = 0\n");
          /*mexPrintf(" with:\n");
          it_code = get_begin_block(block);
          for (int j=0; j < Size; j++)
            {
              if (find(equ_list.begin(), equ_list.end(), j) != equ_list.end())
                mexPrintf("  equ_%d: %s\n",j, print_expression(it_code_expr, false, Size, block, steady_state, 0, 0, it_code, true).c_str());
            }*/
        }
    }
  mxDestroyArray(lhs[0]);
  mxDestroyArray(lhs[1]);
  mxDestroyArray(lhs[2]);
  ostringstream tmp;
  if (block > 1)
    tmp << " in Solve_ByteCode_Sparse_GaussianElimination, singular system in block " << block+1 << "\n";
  else
    tmp << " in Solve_ByteCode_Sparse_GaussianElimination, singular system\n";
  throw FatalExceptionHandling(tmp.str());
}


bool
dynSparseMatrix::Solve_ByteCode_Sparse_GaussianElimination(int Size, int blck, int it_)
{
  bool one;
  int pivj = 0, pivk = 0;
  double *piv_v;
  int *pivj_v, *pivk_v, *NR;
  int l, N_max;
  NonZeroElem *first, *firsta, *first_suba;
  double piv_abs;
  NonZeroElem **bc;
  bc = (NonZeroElem **) mxMalloc(Size*sizeof(*bc));
  piv_v = (double *) mxMalloc(Size*sizeof(double));
  pivj_v = (int *) mxMalloc(Size*sizeof(int));
  pivk_v = (int *) mxMalloc(Size*sizeof(int));
  NR = (int *) mxMalloc(Size*sizeof(int));

  for (int i = 0; i < Size; i++)
    {
      /*finding the max-pivot*/
      double piv = piv_abs = 0;
      int nb_eq = At_Col(i, &first);
      l = 0;
      N_max = 0;
      one = false;
      piv_abs = 0;
      for (int j = 0; j < nb_eq; j++)
        {
          if (!line_done[first->r_index])
            {
              int k = first->u_index;
              int jj = first->r_index;
              int NRow_jj = NRow(jj);

              piv_v[l] = u[k];
              double piv_fabs = fabs(u[k]);
              pivj_v[l] = jj;
              pivk_v[l] = k;
              NR[l] = NRow_jj;
              if (NRow_jj == 1 && !one)
                {
                  one = true;
                  piv_abs = piv_fabs;
                  N_max = NRow_jj;
                }
              if (!one)
                {
                  if (piv_fabs > piv_abs)
                    piv_abs = piv_fabs;
                  if (NRow_jj > N_max)
                    N_max = NRow_jj;
                }
              else
                {
                  if (NRow_jj == 1)
                    {
                      if (piv_fabs > piv_abs)
                        piv_abs = piv_fabs;
                      if (NRow_jj > N_max)
                        N_max = NRow_jj;
                    }
                }
              l++;
            }
          first = first->NZE_C_N;
        }
      if (piv_abs < eps)
        {
          mxFree(piv_v);
          mxFree(pivj_v);
          mxFree(pivk_v);
          mxFree(NR);
          mxFree(bc);
          if (steady_state)
            {
              if (blck > 1)
                mexPrintf("Error: singular system in Simulate_NG in block %d\n", blck+1);
              else
                mexPrintf("Error: singular system in Simulate_NG\n");
              return true;
            }
          else
            {
              ostringstream tmp;
              if (blck > 1)
                tmp << " in Solve_ByteCode_Sparse_GaussianElimination, singular system in block " << blck+1 << "\n";
              else
                tmp << " in Solve_ByteCode_Sparse_GaussianElimination, singular system\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      double markovitz = 0, markovitz_max = -9e70;
      if (!one)
        {
          for (int j = 0; j < l; j++)
            {
              if (N_max > 0 && NR[j] > 0)
                {
                  if (fabs(piv_v[j]) > 0)
                    {
                      if (markowitz_c > 0)
                        markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j])/piv_abs;
              if (markovitz > markovitz_max)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j];   //Line number
                  pivk = pivk_v[j];   //positi
                  markovitz_max = markovitz;
                }
            }
        }
      else
        {
          for (int j = 0; j < l; j++)
            {
              if (N_max > 0 && NR[j] > 0)
                {
                  if (fabs(piv_v[j]) > 0)
                    {
                      if (markowitz_c > 0)
                        markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                    }
                  else
                    markovitz = 0;
                }
              else
                markovitz = fabs(piv_v[j])/piv_abs;
              if (NR[j] == 1)
                {
                  piv = piv_v[j];
                  pivj = pivj_v[j];   //Line number
                  pivk = pivk_v[j];   //positi
                  markovitz_max = markovitz;
                }
            }
        }
      pivot[i] = pivj;
      pivotk[i] = pivk;
      pivotv[i] = piv;
      line_done[pivj] = true;

      /*divide all the non zeros elements of the line pivj by the max_pivot*/
      int nb_var = At_Row(pivj, &first);
      for (int j = 0; j < nb_var; j++)
        {
          u[first->u_index] /= piv;
          first = first->NZE_R_N;
        }
      u[b[pivj]] /= piv;
      /*substract the elements on the non treated lines*/
      nb_eq = At_Col(i, &first);
      NonZeroElem *first_piva;
      int nb_var_piva = At_Row(pivj, &first_piva);
      int nb_eq_todo = 0;
      for (int j = 0; j < nb_eq && first; j++)
        {
          if (!line_done[first->r_index])
            bc[nb_eq_todo++] = first;
          first = first->NZE_C_N;
        }
      //pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
      for (int j = 0; j < nb_eq_todo; j++)
        {
          first = bc[j];
          int row = first->r_index;
          double first_elem = u[first->u_index];

          int nb_var_piv = nb_var_piva;
          NonZeroElem *first_piv = first_piva;
          NonZeroElem *first_sub;
          int nb_var_sub = At_Row(row, &first_sub);
          int l_sub = 0, l_piv = 0;
          int sub_c_index = first_sub->c_index, piv_c_index = first_piv->c_index;
          while (l_sub < nb_var_sub || l_piv < nb_var_piv)
            {
              if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                {
                  first_sub = first_sub->NZE_R_N;
                  if (first_sub)
                    sub_c_index = first_sub->c_index;
                  else
                    sub_c_index = Size;
                  l_sub++;
                }
              else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                {
                  int tmp_u_count = Get_u();
                  Insert(row, first_piv->c_index, tmp_u_count, 0);
                  u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                  first_piv = first_piv->NZE_R_N;
                  if (first_piv)
                    piv_c_index = first_piv->c_index;
                  else
                    piv_c_index = Size;
                  l_piv++;
                }
              else
                {
                  if (i == sub_c_index)
                    {
                      firsta = first;
                      first_suba = first_sub->NZE_R_N;
                      Delete(first_sub->r_index, first_sub->c_index);
                      first = firsta->NZE_C_N;
                      first_sub = first_suba;
                      if (first_sub)
                        sub_c_index = first_sub->c_index;
                      else
                        sub_c_index = Size;
                      l_sub++;
                      first_piv = first_piv->NZE_R_N;
                      if (first_piv)
                        piv_c_index = first_piv->c_index;
                      else
                        piv_c_index = Size;
                      l_piv++;
                    }
                  else
                    {
                      u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                      first_sub = first_sub->NZE_R_N;
                      if (first_sub)
                        sub_c_index = first_sub->c_index;
                      else
                        sub_c_index = Size;
                      l_sub++;
                      first_piv = first_piv->NZE_R_N;
                      if (first_piv)
                        piv_c_index = first_piv->c_index;
                      else
                        piv_c_index = Size;
                      l_piv++;
                    }
                }
            }
          u[b[row]] -= u[b[pivj]]*first_elem;
        }
    }
  double slowc_lbx = slowc;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size; i++)
    ya[i+it_*y_size] = y[i+it_*y_size];

  slowc_save = slowc;
  simple_bksub(it_, Size, slowc_lbx);
  End_GE(Size);
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  mxFree(bc);
  return false;
}

void
dynSparseMatrix::Solve_ByteCode_Symbolic_Sparse_GaussianElimination(int Size, bool symbolic, int Block_number)
{
  /*Triangularisation at each period of a block using a simple gaussian Elimination*/
  t_save_op_s *save_op_s;
  int *save_op = NULL, *save_opa = NULL, *save_opaa = NULL;
  long int nop = 0, nopa = 0;
  bool record = false;
  double *piv_v;
  double piv_abs;
  int *pivj_v, *pivk_v, *NR;
  int pivj = 0, pivk = 0;
  NonZeroElem *first;
  int tmp_u_count, lag;
  int tbreak = 0, last_period = periods;

  piv_v = (double *) mxMalloc(Size*sizeof(double));
  pivj_v = (int *) mxMalloc(Size*sizeof(int));
  pivk_v = (int *) mxMalloc(Size*sizeof(int));
  NR = (int *) mxMalloc(Size*sizeof(int));
  //clock_t time00 = clock();
  NonZeroElem **bc;
  bc = (NonZeroElem **) mxMalloc(Size*sizeof(first));

  for (int t = 0; t < periods; t++)
    {
      /*clock_t time11 = clock();
      mexPrintf("t=%d, record = %d\n",t, record);*/
#ifdef OCTAVE_MEX_FILE
      OCTAVE_QUIT;
#else
    	if ( utIsInterruptPending() )
		    throw UserExceptionHandling();
#endif

      if (record && symbolic)
        {
          /*if (save_op)
            {
              mxFree(save_op);
              save_op = NULL;
            }*/
          save_op = (int *) mxMalloc(nop*sizeof(int));
          nopa = nop;
        }
      nop = 0;
      Clear_u();
      int ti = t*Size;
      for (int i = ti; i < Size+ti; i++)
        {
          /*finding the max-pivot*/
          double piv = piv_abs = 0;
          int nb_eq = At_Col(i, 0, &first);
          if ((symbolic && t <= start_compare) || !symbolic)
            {
              int l = 0, N_max = 0;
              bool one = false;
              piv_abs = 0;
              for (int j = 0; j < nb_eq; j++)
                {
                  if (!line_done[first->r_index])
                    {
                      int k = first->u_index;
                      int jj = first->r_index;
                      int NRow_jj = NRow(jj);
                      piv_v[l] = u[k];
                      double piv_fabs = fabs(u[k]);
                      pivj_v[l] = jj;
                      pivk_v[l] = k;
                      NR[l] = NRow_jj;
                      if (NRow_jj == 1 && !one)
                        {
                          one = true;
                          piv_abs = piv_fabs;
                          N_max = NRow_jj;
                        }
                      if (!one)
                        {
                          if (piv_fabs > piv_abs)
                            piv_abs = piv_fabs;
                          if (NRow_jj > N_max)
                            N_max = NRow_jj;
                        }
                      else
                        {
                          if (NRow_jj == 1)
                            {
                              if (piv_fabs > piv_abs)
                                piv_abs = piv_fabs;
                              if (NRow_jj > N_max)
                                N_max = NRow_jj;
                            }
                        }
                      l++;
                    }
                  first = first->NZE_C_N;
                }
              double markovitz = 0, markovitz_max = -9e70;
              int NR_max = 0;
              if (!one)
                {
                  for (int j = 0; j < l; j++)
                    {
                      if (N_max > 0 && NR[j] > 0)
                        {
                          if (fabs(piv_v[j]) > 0)
                            {
                              if (markowitz_c > 0)
                                markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                              else
                                markovitz = fabs(piv_v[j])/piv_abs;
                            }
                          else
                            markovitz = 0;
                        }
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                      if (markovitz > markovitz_max)
                        {
                          piv = piv_v[j];
                          pivj = pivj_v[j];   //Line number
                          pivk = pivk_v[j];   //positi
                          markovitz_max = markovitz;
                          NR_max = NR[j];
                        }
                    }
                }
              else
                {
                  for (int j = 0; j < l; j++)
                    {
                      if (N_max > 0 && NR[j] > 0)
                        {
                          if (fabs(piv_v[j]) > 0)
                            {
                              if (markowitz_c > 0)
                                markovitz = exp(log(fabs(piv_v[j])/piv_abs)-markowitz_c*log(double (NR[j])/double (N_max)));
                              else
                                markovitz = fabs(piv_v[j])/piv_abs;
                            }
                          else
                            markovitz = 0;
                        }
                      else
                        markovitz = fabs(piv_v[j])/piv_abs;
                      if (NR[j] == 1)
                        {
                          piv = piv_v[j];
                          pivj = pivj_v[j];   //Line number
                          pivk = pivk_v[j];   //positi
                          markovitz_max = markovitz;
                          NR_max = NR[j];
                        }
                    }
                }
              if (fabs(piv) < eps)
                mexPrintf("==> Error NR_max=%d, N_max=%d and piv=%f, piv_abs=%f, markovitz_max=%f\n", NR_max, N_max, piv, piv_abs, markovitz_max);
              if (NR_max == 0)
                mexPrintf("==> Error NR_max=0 and piv=%f, markovitz_max=%f\n", piv, markovitz_max);
              pivot[i] = pivj;
              pivot_save[i] = pivj;
              pivotk[i] = pivk;
              pivotv[i] = piv;
            }
          else
            {
              pivj = pivot[i-Size]+Size;
              pivot[i] = pivj;
              At_Pos(pivj, i, &first);
              pivk = first->u_index;
              piv = u[pivk];
              piv_abs = fabs(piv);
            }
          line_done[pivj] = true;

          if (record && symbolic)
            {
              if (nop+1 >= nopa)
                {
                  nopa = long (mem_increasing_factor*(double) nopa);
                  save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                }
              save_op_s = (t_save_op_s *) (&(save_op[nop]));
              save_op_s->operat = IFLD;
              save_op_s->first = pivk;
              save_op_s->lag = 0;
              nop += 2;
              if (piv_abs < eps)
                {
                  ostringstream tmp;
                  if (Block_number > 1)
                    tmp << " in Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system in block " << Block_number+1 << "\n";
                  else
                    tmp << " in Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system\n";
                  throw FatalExceptionHandling(tmp.str());
                }
              /*divide all the non zeros elements of the line pivj by the max_pivot*/
              int nb_var = At_Row(pivj, &first);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  if (nop+j*2+1 >= nopa)
                    {
                      nopa = long (mem_increasing_factor*(double) nopa);
                      save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                    }
                  save_op_s = (t_save_op_s *) (&(save_op[nop+j*2]));
                  save_op_s->operat = IFDIV;
                  save_op_s->first = first->u_index;
                  save_op_s->lag = first->lag_index;
                  first = first->NZE_R_N;
                }
              nop += nb_var*2;
              u[b[pivj]] /= piv;
              if (nop+1 >= nopa)
                {
                  nopa = long (mem_increasing_factor*(double) nopa);
                  save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                }
              save_op_s = (t_save_op_s *) (&(save_op[nop]));
              save_op_s->operat = IFDIV;
              save_op_s->first = b[pivj];
              save_op_s->lag = 0;
              nop += 2;
              /*substract the elements on the non treated lines*/
              nb_eq = At_Col(i, &first);
              NonZeroElem *first_piva;
              int nb_var_piva = At_Row(pivj, &first_piva);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
//#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) shared(nb_var_piva, first_piva, nopa, save_op) reduction(+:nop)
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  t_save_op_s *save_op_s_l;
                  NonZeroElem *first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  if (nop+1 >= nopa)
                    {
                      nopa = long (mem_increasing_factor*(double) nopa);
                      save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                    }
                  save_op_s_l = (t_save_op_s *) (&(save_op[nop]));
                  save_op_s_l->operat = IFLD;
                  save_op_s_l->first = first->u_index;
                  save_op_s_l->lag = abs(first->lag_index);
                  nop += 2;

                  int nb_var_piv = nb_var_piva;
                  NonZeroElem *first_piv = first_piva;
                  NonZeroElem *first_sub;
                  int nb_var_sub = At_Row(row, &first_sub);
                  int l_sub = 0;
                  int l_piv = 0;
                  int sub_c_index = first_sub->c_index;
                  int piv_c_index = first_piv->c_index;
                  int tmp_lag = first_sub->lag_index;
                  while (l_sub < (nb_var_sub/*=NRow(row)*/) || l_piv < nb_var_piv)
                    {
                      if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                        {
                          //There is no nonzero element at row pivot for this column=> Nothing to do for the current element got to next column
                          first_sub = first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index = first_sub->c_index;
                          else
                            sub_c_index = Size*periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          // There is an nonzero element at row pivot but not at the current row=> insert a negative element in the current row
                          tmp_u_count = Get_u();
                          lag = first_piv->c_index/Size-row/Size;
                          //#pragma omp critical
                            {
                              Insert(row, first_piv->c_index, tmp_u_count, lag);
                            }
                          u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                          if (nop+2 >= nopa)
                            {
                              nopa = long (mem_increasing_factor*(double) nopa);
                              save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                            }
                          save_op_s_l = (t_save_op_s *) (&(save_op[nop]));
                          save_op_s_l->operat = IFLESS;
                          save_op_s_l->first = tmp_u_count;
                          save_op_s_l->second = first_piv->u_index;
                          save_op_s_l->lag = max(first_piv->lag_index, abs(tmp_lag));
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = Size*periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem *firsta = first;
                              NonZeroElem *first_suba = first_sub->NZE_R_N;
                              //#pragma omp critical
                                {
                                  Delete(first_sub->r_index, first_sub->c_index);
                                }
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                              if (nop+3 >= nopa)
                                {
                                  nopa = long (mem_increasing_factor*(double) nopa);
                                  save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                                }
                              save_op_s_l = (t_save_op_s *) (&(save_op[nop]));
                              save_op_s_l->operat = IFSUB;
                              save_op_s_l->first = first_sub->u_index;
                              save_op_s_l->second = first_piv->u_index;
                              save_op_s_l->lag = max(abs(tmp_lag), first_piv->lag_index);
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]]*first_elem;

                  if (nop+3 >= nopa)
                    {
                      nopa = long (mem_increasing_factor*(double) nopa);
                      save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                    }
                  save_op_s_l = (t_save_op_s *) (&(save_op[nop]));
                  save_op_s_l->operat = IFSUB;
                  save_op_s_l->first = b[row];
                  save_op_s_l->second = b[pivj];
                  save_op_s_l->lag = abs(tmp_lag);
                  nop += 3;
                }
            }
          else if(symbolic)
            {
              nop += 2;
              if (piv_abs < eps)
                {
                  ostringstream tmp;
                  if (Block_number > 1)
                    tmp << " in Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system in block " << Block_number+1 << "\n";
                  else
                    tmp << " in Solve_ByteCode_Symbolic_Sparse_GaussianElimination, singular system\n";
                  throw FatalExceptionHandling(tmp.str());
                }
              /*divide all the non zeros elements of the line pivj by the max_pivot*/
              int nb_var = At_Row(pivj, &first);
              for (int j = 0; j < nb_var; j++)
                {
                  u[first->u_index] /= piv;
                  first = first->NZE_R_N;
                }
              nop += nb_var*2;
              u[b[pivj]] /= piv;
              nop += 2;
              /*substract the elements on the non treated lines*/
              nb_eq = At_Col(i, &first);
              NonZeroElem *first_piva;
              int nb_var_piva = At_Row(pivj, &first_piva);

              int nb_eq_todo = 0;
              for (int j = 0; j < nb_eq && first; j++)
                {
                  if (!line_done[first->r_index])
                    bc[nb_eq_todo++] = first;
                  first = first->NZE_C_N;
                }
//#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) shared(nb_var_piva, first_piva, nopa, save_op) reduction(+:nop)
              for (int j = 0; j < nb_eq_todo; j++)
                {
                  NonZeroElem *first = bc[j];
                  int row = first->r_index;
                  double first_elem = u[first->u_index];
                  nop += 2;
                  int nb_var_piv = nb_var_piva;
                  NonZeroElem *first_piv = first_piva;
                  NonZeroElem *first_sub;
                  int nb_var_sub = At_Row(row, &first_sub);
                  int l_sub = 0;
                  int l_piv = 0;
                  int sub_c_index = first_sub->c_index;
                  int piv_c_index = first_piv->c_index;
                  while (l_sub < (nb_var_sub /*= NRow(row)*/) || l_piv < nb_var_piv)
                    {
                      if (l_sub < nb_var_sub && (sub_c_index < piv_c_index || l_piv >= nb_var_piv))
                        {
                          //There is no nonzero element at row pivot for this column=> Nothing to do for the current element got to next column
                          first_sub = first_sub->NZE_R_N;
                          if (first_sub)
                            sub_c_index = first_sub->c_index;
                          else
                            sub_c_index = Size*periods;
                          l_sub++;
                        }
                      else if (sub_c_index > piv_c_index || l_sub >= nb_var_sub)
                        {
                          // There is an nonzero element at row pivot but not at the current row=> insert a negative element in the current row
                          tmp_u_count = Get_u();
                          lag = first_piv->c_index/Size-row/Size;
                          //#pragma omp critical
                           {
                             Insert(row, first_piv->c_index, tmp_u_count, lag);
                           }
                          u[tmp_u_count] = -u[first_piv->u_index]*first_elem;
                          nop += 3;
                          first_piv = first_piv->NZE_R_N;
                          if (first_piv)
                            piv_c_index = first_piv->c_index;
                          else
                            piv_c_index = Size*periods;
                          l_piv++;
                        }
                      else /*first_sub->c_index==first_piv->c_index*/
                        {
                          if (i == sub_c_index)
                            {
                              NonZeroElem *firsta = first;
                              NonZeroElem *first_suba = first_sub->NZE_R_N;
                              //#pragma omp critical
                                {
                                  Delete(first_sub->r_index, first_sub->c_index);
                                }
                              first = firsta->NZE_C_N;
                              first_sub = first_suba;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                          else
                            {
                              u[first_sub->u_index] -= u[first_piv->u_index]*first_elem;
                              nop += 3;
                              first_sub = first_sub->NZE_R_N;
                              if (first_sub)
                                sub_c_index = first_sub->c_index;
                              else
                                sub_c_index = Size*periods;
                              l_sub++;
                              first_piv = first_piv->NZE_R_N;
                              if (first_piv)
                                piv_c_index = first_piv->c_index;
                              else
                                piv_c_index = Size*periods;
                              l_piv++;
                            }
                        }
                    }
                  u[b[row]] -= u[b[pivj]]*first_elem;
                  nop += 3;
                }
            }
        }
      if (symbolic)
        {
		  if (t > int(periods*0.35))
            {
              symbolic = false;
              mxFree(save_opaa);
              mxFree(save_opa);
              mxFree(save_op);
            }
          else if (record && (nop == nop1))
            {
              if (t > int(periods*0.35))
                {
                  symbolic = false;
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = NULL;
                    }
                  if (save_opa)
                    {
                      mxFree(save_opa);
                      save_opa = NULL;
                    }
                  if (save_op)
                    {
                      mxFree(save_op);
                      save_op = NULL;
                    }
                }
              else if (save_opa && save_opaa)
                {
                  if (compare(save_op, save_opa, save_opaa, t, periods, nop, Size))
                    {
                      tbreak = t;
                      tbreak_g = tbreak;
                      //mexPrintf("time=%f\n",(1000.0*(double (clock())-double (time11)))/double (CLOCKS_PER_SEC));
                      break;
                    }
                }
              if (save_opa)
                {
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = NULL;
                    }
                  save_opaa = save_opa;
                }
              save_opa = save_op;
            }
          else
            {
              if (nop == nop1)
                record = true;
              else
                {
                  record = false;
                  if (save_opa)
                    {
                      mxFree(save_opa);
                      save_opa = NULL;
                    }
                  if (save_opaa)
                    {
                      mxFree(save_opaa);
                      save_opaa = NULL;
                    }
                }
            }
          nop2 = nop1;
          nop1 = nop;
        }
      //mexPrintf("time=%f\n",(1000.0*(double (clock())-double (time11)))/double (CLOCKS_PER_SEC));
    }
  mxFree(bc);
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
  /*mexPrintf("tbreak=%d, periods=%d time required=%f\n",tbreak,periods, (1000.0*(double (clock())-double (time00)))/double (CLOCKS_PER_SEC));
  mexEvalString("drawnow;");
  time00 = clock();*/
  nop_all += nop;
  if (symbolic)
    {
      if (save_op)
        mxFree(save_op);
      if (save_opa)
        mxFree(save_opa);
      if (save_opaa)
        mxFree(save_opaa);
    }

  /*The backward substitution*/
  double slowc_lbx = slowc;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  slowc_save = slowc;
  bksub(tbreak, last_period, Size, slowc_lbx);
  /*mexPrintf("remaining operations and bksub time required=%f\n",tbreak,periods, (1000.0*(double (clock())-double (time00)))/double (CLOCKS_PER_SEC));
  mexEvalString("drawnow;");*/
  End_GE(Size);
}


void
dynSparseMatrix::Grad_f_product(int n, mxArray *b_m, double* vectr, mxArray *A_m, SuiteSparse_long *Ap, SuiteSparse_long *Ai, double* Ax, double* b_)
{
  if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
    {
      NonZeroElem *first;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) private(first)
#endif
      for (int i = 0; i < n; i++)
        {
          double sum = 0;
          first = FNZE_R[i];
          if (first)
            for (int k = 0; k < NbNZRow[i]; k++)
              {
                sum += u[first->u_index] * u[b[first->c_index]];
                first = first->NZE_R_N;
              }
          vectr[i] = sum;
        }
    }
  else
    {
      if (!((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4) && !steady_state)))
        {
          mwIndex *Ai = mxGetIr(A_m);
          if (!Ai)
            {
              ostringstream tmp;
              tmp << " in Init_Matlab_Sparse_Simple, can't allocate Ai index vector\n";
              throw FatalExceptionHandling(tmp.str());
            }
          mwIndex *Aj = mxGetJc(A_m);
          if (!Aj)
            {
              ostringstream tmp;
              tmp << " in Init_Matlab_Sparse_Simple, can't allocate Aj index vector\n";
              throw FatalExceptionHandling(tmp.str());
            }
          double *A = mxGetPr(A_m);
          if (!A)
            {
              ostringstream tmp;
              tmp << " in Init_Matlab_Sparse_Simple, can't retrieve A matrix\n";
              throw FatalExceptionHandling(tmp.str());
            }
          b_ = mxGetPr(b_m);
          if (!b_)
            {
              ostringstream tmp;
              tmp << " in Init_Matlab_Sparse_Simple, can't retrieve b matrix\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      memset(vectr, 0, n * sizeof(double));
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) /*shared(vectr)*/
#endif
      for (int i = 0; i < n; i++)
        for (SuiteSparse_long j = Ap[i]; j < Ap[i+1]; j++)
          vectr[Ai[j]] += Ax[j] * b_[i];
    }
}

void
dynSparseMatrix::Check_and_Correct_Previous_Iteration(int block_num, int y_size, int size, double crit_opt_old)
{
  double top = 1.0;
  double bottom = 0.1;
  //mexPrintf("res2=%f > g0=%f, res1=%f, iter=%d it_=%d\n", res2, g0, res1, iter, it_);
  if (isnan(res1) || isinf(res1) || (res2 > g0 && iter > 0))
    {
      while ((isnan(res1) || isinf(res1)))
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
            }
          /*mexPrintf("reducing solwc_save = %e, it_=%d, y_size=%d, size=%d, y[%d]=%e, ya[%d]=%e,\n y[%d]=%e, ya[%d]=%e\n",slowc_save, it_, y_size, size-1, index_vara[0]+it_*y_size, y[index_vara[0]+it_*y_size], index_vara[0]+it_*y_size, ya[index_vara[0]+it_*y_size]
                                                                                                       , index_vara[size-1]+it_*y_size, y[index_vara[size-1]+it_*y_size], index_vara[size-1]+it_*y_size, ya[index_vara[size-1]+it_*y_size]);*/
           //mexPrintf("->slowc_save=%f\n",slowc_save);
           compute_complete(true, res1, res2, max_res, max_res_idx);
        }

      while (res2 > g0 && slowc_save > 1e-1)
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.5;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
            }
          /*mexPrintf("reducing solwc_save = %e, it_=%d, y_size=%d, size=%d, y[%d]=%e, ya[%d]=%e,\n y[%d]=%e, ya[%d]=%e\n",slowc_save, it_, y_size, size-1, index_vara[0]+it_*y_size, y[index_vara[0]+it_*y_size], index_vara[0]+it_*y_size, ya[index_vara[0]+it_*y_size]                                                                                            , index_vara[size-1]+it_*y_size, y[index_vara[size-1]+it_*y_size], index_vara[size-1]+it_*y_size, ya[index_vara[size-1]+it_*y_size]);*/
          //mexPrintf("->slowc_save=%f\n",slowc_save);
          compute_complete(true, res1, res2, max_res, max_res_idx);
        }
      double ax = slowc_save-0.001, bx = slowc_save+0.001, cx = slowc_save, fa, fb, fc, xmin;
      if (false/*slowc_save > 2e-1*/)
        if (mnbrak(&ax, &bx, &cx, &fa, &fb, &fc))
          if (golden(ax, bx, cx, 1e-1, solve_tolf, &xmin))
            slowc_save = xmin;
      //mexPrintf("cx=%f\n", cx);
      //mexPrintf("ax= %f, bx=%f, cx=%f, fa=%f, fb=%f, fc=%d\n", ax, bx, cx, fa, fb, fc);

      //if (!(isnan(res1) || isinf(res1))/* && !(isnan(g0) || isinf(g0))*//*|| (res2 > g0 && iter > 1)*/)
      if (false)
        {

          double *p = (double*)mxMalloc(size * sizeof(double));
          Grad_f_product(size, b_m_save, p, A_m_save, Ap_save, Ai_save, Ax_save, b_save);
          double slope=0.0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:slope)
#endif
          for (int i = 1; i < size; i++)
            slope += - direction[i] * p[i];
          /*if (slope > 0)
            mexPrintf("Roundoff in lnsearch\n");
          else*/
            {
              prev_slowc_save = 1;
              double crit_opt = res2/2;
              double max_try_iteration = 100;
              double small_ = 1.0e-4;
              bool try_at_cvg = false;
              while ((try_at_iteration < max_try_iteration) && (!try_at_cvg) && (abs(prev_slowc_save - slowc_save) > 1e-10))
                {
                  crit_opt = res2 / 2;
                  if (slowc_save < 1e-7)
                    {
                      try_at_cvg = true;
                      continue;
                    }
                  else if ((crit_opt <= crit_opt_old + small_ * slowc_save * slope) && !(isnan(res1) || isinf(res1)))
                    {
                      try_at_cvg = true;
                      continue;
                    }
                  else if (try_at_iteration == 0)
                    {
                      prev_slowc_save = slowc_save;
                      //slowc_save = max(- top * slope / ( (crit_opt - crit_opt_old - slope)), bottom);
                      slowc_save /= 1.2;
                    }
                  else
                    {
                      double t1 = crit_opt - slope * slowc_save - crit_opt_old;
                      double t2 = glambda2 - slope * prev_slowc_save - crit_opt_old;
                      double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
                      double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
                      if (a == 0)
                        slowc_save = max(min( - slope/(2 * b) , top * slowc_save), bottom * slowc_save);
                      else
                        {
                          double delta = b*b - 3 * a * slope;
                          if (delta <= 0)
                            slowc_save = top * slowc_save;
                          else if (b <= 0)
                            slowc_save = max(min(-b + sqrt(delta) / (3 * a), top * slowc_save), bottom * slowc_save);
                          else
                            slowc_save = max(min(- slope / (b + sqrt(delta)), top * slowc_save), bottom * slowc_save);
                        }
                    }
                  if (abs(prev_slowc_save - slowc_save) < 1e-10)
                    slowc_save /= 1.1;
                  //mexPrintf("=>slowc_save=%f, prev_slowc_save=%f\n",slowc_save, prev_slowc_save);
                  prev_slowc_save = slowc_save;
                  glambda2 = crit_opt;
                  try_at_iteration++;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
                  for (int i = 0; i < size; i++)
                    {
                      int eq = index_vara[i];
                      y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
                    }
                  compute_complete(true, res1, res2, max_res, max_res_idx);
                }
            }
          mxFree(p);
        }
      //if (print_it)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
        }
      compute_complete(false, res1, res2, max_res, max_res_idx);
    }
  else
    {
      //mexPrintf("slowc_save=%f res1=%f\n",slowc_save, res1);
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq+it_*y_size] = ya[eq+it_*y_size] + slowc_save * direction[eq+it_*y_size];
        }
    }
  slowc_save = slowc;
}

bool
dynSparseMatrix::Simulate_One_Boundary(int block_num, int y_size, int y_kmin, int y_kmax, int size, bool cvg)
{
  //int i;
  mxArray *b_m = NULL, *A_m = NULL, *x0_m = NULL;
  SuiteSparse_long *Ap = NULL, *Ai = NULL;
  double *Ax = NULL, *b = NULL;
  int preconditioner = 1;


  try_at_iteration = 0;
  Clear_u();
  bool singular_system = false;
  u_count_alloc_save = u_count_alloc;

  if (isnan(res1) || isinf(res1))
    {
#ifdef DEBUG
      for (int j = 0; j < y_size; j++)
        {
          bool select = false;
          for (int i = 0; i < size; i++)
            if (j == index_vara[i])
              {
                select = true;
                break;
              }
          if (select)
            mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n", get_variable(eEndogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
          else
            mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n", get_variable(eEndogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
        }
#endif
      if (steady_state)
        {
          if (iter == 0)
            mexPrintf(" the initial values of endogenous variables are too far from the solution.\nChange them!\n");
          else
            mexPrintf(" dynare cannot improve the simulation in block %d at time %d (variable %d)\n", block_num+1, it_+1, index_vara[max_res_idx]+1);
          mexEvalString("drawnow;");
          //return singular_system;
        }
      else
        {
          ostringstream tmp;
          if (iter == 0)
            tmp << " in Simulate_One_Boundary, The initial values of endogenous variables are too far from the solution.\nChange them!\n";
          else
            tmp << " in Simulate_One_Boundary, Dynare cannot improve the simulation in block " << block_num+1 << " at time " << it_+1 << " (variable " << index_vara[max_res_idx]+1 << "%d)\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  if (print_it)
    {
      if (steady_state)
        {
          switch (solve_algo)
            {
            case 0:
              mexPrintf("MODEL STEADY STATE: MATLAB fsolve\n");
              break;
            case 1:
              mexPrintf("MODEL STEADY STATE: MATLAB solve1\n");
              break;
            case 2:
            case 4:
              mexPrintf("MODEL STEADY STATE: block decomposition + MATLAB solve1\n");
              break;
            case 3:
              mexPrintf("MODEL STEADY STATE: MATLAB csolve\n");
              break;
            case 5:
              mexPrintf("MODEL STEADY STATE: (method=ByteCode own solver)\n");
              break;
            case 6:
              mexPrintf("MODEL STEADY STATE: Sparse LU\n");
              break;
            case 7:
                mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=GMRES)\n", preconditioner, true).c_str());
                //mexPrintf("MODEL STEADY STATE: (method=GMRES)\n");
              break;
            case 8:
              mexPrintf(preconditioner_print_out("MODEL STEADY STATE: (method=BiCGStab)\n", preconditioner, true).c_str());
              //mexPrintf("MODEL STEADY STATE: (method=BiCGStab)\n");
              break;
            default:
              mexPrintf("MODEL STEADY STATE: (method=Unknown - %d - )\n", stack_solve_algo);
            }
        }

      mexPrintf("-----------------------------------\n");
      mexPrintf("      Simulate iteration no %d     \n", iter+1);
      mexPrintf("      max. error=%.10e       \n", double (max_res));
      mexPrintf("      sqr. error=%.10e       \n", double (res2));
      mexPrintf("      abs. error=%.10e       \n", double (res1));
      mexPrintf("-----------------------------------\n");
    }
  bool zero_solution;

  if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
    Simple_Init(size, IM_i, zero_solution);
  else
    {
      b_m = mxCreateDoubleMatrix(size, 1, mxREAL);
      if (!b_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_One_Boundary, can't allocate b_m vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      A_m = mxCreateSparse(size, size, min(int (IM_i.size()*2), size * size), mxREAL);
      if (!A_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_One_Boundary, can't allocate A_m matrix\n";
          throw FatalExceptionHandling(tmp.str());
        }
      x0_m = mxCreateDoubleMatrix(size, 1, mxREAL);
      if (!x0_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_One_Boundary, can't allocate x0_m vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      if (!((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 4) && !steady_state)))
        {
          Init_Matlab_Sparse_Simple(size, IM_i, A_m, b_m, zero_solution, x0_m);
          A_m_save = mxDuplicateArray(A_m);
          b_m_save = mxDuplicateArray(b_m);
        }
      else
        {
          Init_UMFPACK_Sparse_Simple(size, IM_i, &Ap, &Ai, &Ax, &b, zero_solution, x0_m);
          if (Ap_save[size] != Ap[size])
            {
              mxFree(Ai_save);
              mxFree(Ax_save);
              Ai_save = (SuiteSparse_long*)mxMalloc(Ap[size] * sizeof(SuiteSparse_long));
              Ax_save = (double*)mxMalloc(Ap[size] * sizeof(double));
            }
          memcpy(Ap_save, Ap, (size + 1) * sizeof(SuiteSparse_long));
          memcpy(Ai_save, Ai, Ap[size] * sizeof(SuiteSparse_long));
          memcpy(Ax_save, Ax, Ap[size] * sizeof(double));
          memcpy(b_save, b, size * sizeof(double));
        }
    }
  if (zero_solution)
    {
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          double yy = -(y[eq+it_*y_size]);
          direction[eq] = yy;
          y[eq+it_*y_size] += slowc * yy;
        }
    }
  else
    {
      if ((solve_algo == 5 && steady_state) || (stack_solve_algo == 5 && !steady_state))
        singular_system = Solve_ByteCode_Sparse_GaussianElimination(size, block_num, it_);
      else if ((solve_algo == 7 && steady_state) || (stack_solve_algo == 2 && !steady_state))
        Solve_Matlab_GMRES(A_m, b_m, size, slowc, block_num, false, it_, x0_m);
      else if ((solve_algo == 8 && steady_state) || (stack_solve_algo == 3 && !steady_state))
        Solve_Matlab_BiCGStab(A_m, b_m, size, slowc, block_num, false, it_, x0_m, preconditioner);
      else if ((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4) && !steady_state))
        Solve_LU_UMFPack(Ap, Ai, Ax, b, size, size, slowc, true, 0);
    }
  return singular_system;
}




bool
dynSparseMatrix::solve_linear(const int block_num, const int y_size, const int y_kmin, const int y_kmax, const int size, const int iter)
{
  bool cvg = false;
  double crit_opt_old = res2/2;
  compute_complete(false, res1, res2, max_res, max_res_idx);
  cvg = (max_res < solve_tolf);
  if (!cvg || isnan(res1) || isinf(res1))
    {
      if (iter)
        Check_and_Correct_Previous_Iteration(block_num, y_size, size, crit_opt_old);
      bool singular_system = Simulate_One_Boundary(block_num, y_size, y_kmin, y_kmax, size, cvg);
      if (singular_system)
        Singular_display(block_num, size);
    }
  return cvg;
}

void
dynSparseMatrix::solve_non_linear(const int block_num, const int y_size, const int y_kmin, const int y_kmax, const int size)

{
  max_res_idx = 0;
  bool cvg = false;
  iter = 0;
  glambda2 = g0 = very_big;
  //try_at_iteration = 0;
  while ((!cvg) && (iter < maxit_))
    {
      cvg = solve_linear(block_num, y_size, y_kmin, y_kmax, size, iter);
      g0 = res2;
      iter++;
    }
  if (!cvg)
    {
      ostringstream tmp;
      if (steady_state)
        tmp << " in Solve Forward complete, convergence not achieved in block " << block_num+1 << ", after " << iter << " iterations\n";
      else
        tmp << " in Solve Forward complete, convergence not achieved in block " << block_num+1 << ", at time " << it_ << ", after " << iter << " iterations\n";
      throw FatalExceptionHandling(tmp.str());
    }
}

void
dynSparseMatrix::Simulate_Newton_One_Boundary(const bool forward)
{
  g1 = (double *) mxMalloc(size*size*sizeof(double));
  r = (double *) mxMalloc(size*sizeof(double));
  //mexPrintf("Simulate_Newton_One_Boundary, block_num=%d, size=%d, steady=%d, forward=%d, iter=%d, is_linear=%d\n", block_num, size, steady_state, forward, iter, is_linear);
  iter = 0;
  if ((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4) && !steady_state))
    {
      Ap_save = (SuiteSparse_long*)mxMalloc((size + 1) * sizeof(SuiteSparse_long));
      Ap_save[size] = 0;
      Ai_save = (SuiteSparse_long*)mxMalloc(1 * sizeof(SuiteSparse_long));
      Ax_save = (double*)mxMalloc(1 * sizeof(double));
      b_save = (double*)mxMalloc((size) * sizeof(SuiteSparse_long));
    }
  if (steady_state)
    {
      it_ = 0;
      if (!is_linear)
        solve_non_linear(block_num, y_size, 0, 0, size);
      else
        solve_linear(block_num, y_size, 0, 0, size, 0);
    }
  else if (forward)
    {
      if (!is_linear)
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            solve_non_linear(block_num, y_size, y_kmin, y_kmax, size);
        }
      else
        {
          for (int it_ = y_kmin; it_ < periods+y_kmin; it_++)
            solve_linear(block_num, y_size, y_kmin, y_kmax, size, 0);
        }
    }
  else
    {
      if (!is_linear)
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            solve_non_linear(block_num, y_size, y_kmin, y_kmax, size);
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            solve_linear(block_num, y_size, y_kmin, y_kmax, size, 0);
        }
    }
  if ((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1 || stack_solve_algo == 4) && !steady_state))
    {
      mxFree(Ap_save);
      mxFree(Ai_save);
      mxFree(Ax_save);
      mxFree(b_save);
    }
  mxFree(g1);
  mxFree(r);
}

string
dynSparseMatrix::preconditioner_print_out(string s, int preconditioner, bool ss)
{
  int n = s.length();
  string tmp = ", preconditioner=";
  switch(preconditioner)
    {
    case 0:
      if (ss)
        tmp.append("Jacobi on static jacobian");
      else
        tmp.append("Jacobi on dynamic jacobian");
      break;
    case 1:
      if (ss)
        tmp.append("incomplet lutp on static jacobian");
      else 
        tmp.append("incomplet lu0 on dynamic jacobian");
      break;
    case 2:
      tmp.append("incomplet lutp on dynamic jacobian");
      break;
    case 3:
      tmp.append("lu on static jacobian");
      break;
    }
  s.insert(n - 2, tmp);
  return s;
}

void
dynSparseMatrix::Simulate_Newton_Two_Boundaries(int blck, int y_size, int y_kmin, int y_kmax, int Size, int periods, bool cvg, int minimal_solving_periods, int stack_solve_algo, unsigned int endo_name_length, char *P_endo_names)
{
  double top = 0.5;
  double bottom = 0.1;
#ifdef CUDA
  int nnz, nnz_tild;
  int *Ap_i, *Ai_i;
  int *Ap_i_tild, *Ai_i_tild;
  double *x0, *A_tild;

#endif
  int preconditioner = 2;
  if (start_compare == 0)
    start_compare = y_kmin;
  u_count_alloc_save = u_count_alloc;
  clock_t t1 = clock();
  nop1 = 0;
  mxArray *b_m = NULL, *A_m = NULL, *x0_m = NULL;
  double *Ax = NULL, *b;
  SuiteSparse_long *Ap = NULL, *Ai = NULL;




  if (iter > 0)
    {
      if (print_it)
        {
          mexPrintf("Sim : %f ms\n", (1000.0*(double (clock())-double (time00)))/double (CLOCKS_PER_SEC));
          mexEvalString("drawnow;");
        }
      time00 = clock();
    }
  if (isnan(res1) || isinf(res1) || (res2 > 12*g0 && iter > 0))
    {
      if (iter == 0 || fabs(slowc_save) < 1e-8)
        {
          for (int j = 0; j < y_size; j++)
            {
              ostringstream res;
              for (unsigned int i = 0; i < endo_name_length; i++)
                if (P_endo_names[CHAR_LENGTH*(j+i*y_size)] != ' ')
                  res << P_endo_names[CHAR_LENGTH*(j+i*y_size)];
              bool select = false;
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
                    select = true;
                    break;
                  }
              if (select)
                mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n", res.str().c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n", res.str().c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
            }
          ostringstream Error;
          if (iter == 0)
            Error << " in Simulate_Newton_Two_Boundaries, the initial values of endogenous variables are too far from the solution.\nChange them!\n";
          else
            Error << " in Simulate_Newton_Two_Boundaries, dynare cannot improve the simulation in block " << blck+1 << " at time " << it_+1 << " (variable " << index_vara[max_res_idx]+1 << ")\n";
          throw FatalExceptionHandling(Error.str());
        }
      if (!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0)) && (stack_solve_algo == 4 || stack_solve_algo == 5))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max(-gp0 / (2 * (res2 - g0 - gp0)), bottom);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min(-b + sqrt(b*b - 3 * a * gp0) / (3 * a), top * slowc_save), bottom * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration++;
          if (slowc_save <= bottom)
            {
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
              for (int i = 0; i < y_size*(periods+y_kmin); i++)
                y[i] = ya[i]+direction[i];
              g0 = res2;
              gp0 = -res2;
              try_at_iteration = 0;
              iter--;
              return;
            }
        }
      else
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
        }
      if (print_it)
        {
          if (isnan(res1) || isinf(res1))
            mexPrintf("The model cannot be evaluated, trying to correct it using slowc=%f\n", slowc_save);
          else
            mexPrintf("Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
        }
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < y_size*(periods+y_kmin); i++)
        y[i] = ya[i]+slowc_save*direction[i];
      iter--;
      return;
    }
  u_count += u_count_init;
  if (stack_solve_algo == 5)
    {
      if (alt_symbolic && alt_symbolic_count < alt_symbolic_count_max)
        {
          mexPrintf("Pivoting method will be applied only to the first periods.\n");
          alt_symbolic = false;
          symbolic = true;
          markowitz_c = markowitz_c_s;
          alt_symbolic_count++;
        }
      if (((res1/res1a-1) > -0.3) && symbolic && iter > 0)
        {
          if (restart > 2)
            {
              mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied to all periods.\n");
              symbolic = false;
              alt_symbolic = true;
              markowitz_c_s = markowitz_c;
              markowitz_c = 0;
            }
          else
            {
              mexPrintf("Divergence or slowdown occured during simulation.\nIn the next iteration, pivoting method will be applied for a longer period.\n");
              start_compare = min(tbreak_g, periods);
              restart++;
            }
        }
      else
        {
          start_compare = max(y_kmin, minimal_solving_periods);
          restart = 0;
        }
    }
  res1a = res1;

  if (print_it)
    {
      if (iter == 0)
        {
          switch (stack_solve_algo)
            {
            case 0:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU)\n");
              break;
            case 1:
              mexPrintf("MODEL SIMULATION: (method=Relaxation)\n");
              break;
            case 2:
              mexPrintf(preconditioner_print_out("MODEL SIMULATION: (method=GMRES)\n", preconditioner, false).c_str());
              break;
            case 3:
              mexPrintf(preconditioner_print_out("MODEL SIMULATION: (method=BiCGStab)\n", preconditioner, false).c_str());
              break;
            case 4:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU & optimal path length)\n");
              break;
            case 5:
              mexPrintf("MODEL SIMULATION: (method=ByteCode own solver)\n");
              break;
            case 7:
              mexPrintf(preconditioner_print_out("MODEL SIMULATION: (method=GPU BiCGStab)\n", preconditioner, false).c_str());
              break;
            default:
              mexPrintf("MODEL SIMULATION: (method=Unknown - %d - )\n", stack_solve_algo);
            }
        }
      mexPrintf("-----------------------------------\n");
      mexPrintf("      Simulate iteration no %d     \n", iter+1);
      mexPrintf("      max. error=%.10e       \n", double (max_res));
      mexPrintf("      sqr. error=%.10e       \n", double (res2));
      mexPrintf("      abs. error=%.10e       \n", double (res1));
      mexPrintf("-----------------------------------\n");
      mexEvalString("drawnow;");
    }
  if (cvg)
    {
      return;
    }
  else
    {
      if (stack_solve_algo == 5)
        Init_GE(periods, y_kmin, y_kmax, Size, IM_i);
      else
        {
          b_m = mxCreateDoubleMatrix(periods*Size, 1, mxREAL);
          if (!b_m)
            {
              ostringstream tmp;
              tmp << " in Simulate_Newton_Two_Boundaries, can't allocate b_m vector\n";
              throw FatalExceptionHandling(tmp.str());
            }
          x0_m = mxCreateDoubleMatrix(periods*Size, 1, mxREAL);
          if (!x0_m)
            {
              ostringstream tmp;
              tmp << " in Simulate_Newton_Two_Boundaries, can't allocate x0_m vector\n";
              throw FatalExceptionHandling(tmp.str());
            }
          if (stack_solve_algo != 0 && stack_solve_algo != 4 && stack_solve_algo != 7)
            {
              A_m = mxCreateSparse(periods*Size, periods*Size, IM_i.size()* periods*2, mxREAL);
              if (!A_m)
                {
                  ostringstream tmp;
                  tmp << " in Simulate_Newton_Two_Boundaries, can't allocate A_m matrix\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
          if (stack_solve_algo == 0 || stack_solve_algo == 4)
            Init_UMFPACK_Sparse(periods, y_kmin, y_kmax, Size, IM_i, &Ap, &Ai, &Ax, &b, x0_m);
#ifdef CUDA
          else if (stack_solve_algo == 7)
            Init_CUDA_Sparse(periods, y_kmin, y_kmax, Size, IM_i, &Ap_i, &Ai_i, &Ax, &Ap_i_tild, &Ai_i_tild, &A_tild, &b, &x0, x0_m, &nnz, &nnz_tild, preconditioner);
#endif
          else
            Init_Matlab_Sparse(periods, y_kmin, y_kmax, Size, IM_i, A_m, b_m, x0_m);

        }
      //if (iter > 0)
      /*mexPrintf("--> stack_solve_algo=%d\n", stack_solve_algo);
      mexEvalString("drawnow;");*/

      if (stack_solve_algo == 0 || stack_solve_algo == 4)
        Solve_LU_UMFPack(Ap, Ai, Ax, b, Size * periods, Size, slowc, true, 0);
      else if (stack_solve_algo == 1)
        Solve_Matlab_Relaxation(A_m, b_m, Size, slowc, true, 0);
      else if (stack_solve_algo == 2)
        Solve_Matlab_GMRES(A_m, b_m, Size, slowc, blck, true, 0, x0_m);
      else if (stack_solve_algo == 3)
        Solve_Matlab_BiCGStab(A_m, b_m, Size, slowc, blck, true, 0, x0_m, 1);
      else if (stack_solve_algo == 5)
        Solve_ByteCode_Symbolic_Sparse_GaussianElimination(Size, symbolic, blck);
#ifdef CUDA
      else if (stack_solve_algo == 7)
        Solve_CUDA_BiCGStab(Ap_i, Ai_i, Ax, Ap_i_tild, Ai_i_tild, A_tild, b, x0, Size * periods, Size, slowc, true, 0, nnz, nnz_tild, preconditioner, Size * periods, blck);
#endif
    }
  if (print_it)
    {
      clock_t t2 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double (t2) - double (t1))/double (CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }
  if ((!steady_state && (stack_solve_algo == 4 /*|| stack_solve_algo == 0*/))/* || steady_state*/)
    {
      clock_t t2 = clock();
      double ax = -0.1, bx = 1.1, cx = 0.5, fa, fb, fc, xmin;

      if (!mnbrak(&ax, &bx, &cx, &fa, &fb, &fc))
        return;
      //mexPrintf("ax= %f, bx=%f, cx=%f, fa=%f, fb=%f, fc=%d\n", ax, bx, cx, fa, fb, fc);
      if (!golden(ax, bx, cx, 1e-1, solve_tolf, &xmin))
        return;
      slowc = xmin;
      clock_t t3 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double (t3) - double (t2))/double (CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }
  time00 = clock();
  if (tbreak_g == 0)
    tbreak_g = periods;
  return;
}

void
dynSparseMatrix::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
{
  u_count = u_count_int * periods;
  u_count_alloc = 2*u_count;
#ifdef DEBUG
  mexPrintf("fixe_u : alloc(%d double)\n", u_count_alloc);
#endif
  (*u) = (double *) mxMalloc(u_count_alloc*sizeof(double));
#ifdef DEBUG
  mexPrintf("*u=%d\n", *u);
#endif
  memset((*u), 0, u_count_alloc*sizeof(double));
  u_count_init = max_lag_plus_max_lead_plus_1;
}

