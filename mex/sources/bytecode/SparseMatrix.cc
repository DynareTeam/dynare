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

//#define _GLIBCXX_USE_C99_FENV_TR1 1
//#include <cfenv>

#include <cstring>
#include <ctime>
#include <sstream>
#include "SparseMatrix.hh"

SparseMatrix::SparseMatrix()
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
}

int
SparseMatrix::NRow(int r)
{
  return NbNZRow[r];
}

int
SparseMatrix::NCol(int c)
{
  return NbNZCol[c];
}

int
SparseMatrix::At_Row(int r, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  return NbNZRow[r];
}

int
SparseMatrix::Union_Row(int row1, int row2)
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
SparseMatrix::At_Pos(int r, int c, NonZeroElem **first)
{
  (*first) = FNZE_R[r];
  while ((*first)->c_index != c)
    (*first) = (*first)->NZE_R_N;
  return NbNZRow[r];
}

int
SparseMatrix::At_Col(int c, NonZeroElem **first)
{
  (*first) = FNZE_C[c];
  return NbNZCol[c];
}

int
SparseMatrix::At_Col(int c, int lag, NonZeroElem **first)
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
SparseMatrix::Delete(const int r, const int c)
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
SparseMatrix::Print(int Size, int *b)
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
SparseMatrix::Insert(const int r, const int c, const int u_index, const int lag_index)
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
SparseMatrix::Read_SparseMatrix(string file_name, const int Size, int periods, int y_kmin, int y_kmax, bool steady_state, bool two_boundaries, int stack_solve_algo, int solve_algo)
{
  unsigned int eq, var;
  int i, j, lag;
  filename = file_name;
  mem_mngr.fixe_file_name(file_name);
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
          for (i = 0; i < u_count_init-Size; i++)
            {
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
              IM_i[make_pair(make_pair(eq, var), lag)] = j;
            }
          for (j = 0; j < Size; j++)
            IM_i[make_pair(make_pair(j, Size*(periods+y_kmax)), 0)] = j;
        }
      else if (stack_solve_algo >= 0 || stack_solve_algo <= 4)
        {
          for (i = 0; i < u_count_init-Size; i++)
            {
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
              IM_i[make_pair(make_pair(var - lag*Size, -lag), eq)] = j;
            }
          for (j = 0; j < Size; j++)
            IM_i[make_pair(make_pair(Size*(periods+y_kmax), 0), j)] = j;
        }

    }
  else
    {
      if ((stack_solve_algo == 5 && !steady_state) || (solve_algo == 5 && steady_state))
        {
          for (i = 0; i < u_count_init; i++)
            {
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
              IM_i[make_pair(make_pair(eq, var), lag)] = j;
            }
        }
      else if (((stack_solve_algo >= 0 || stack_solve_algo <= 4) && !steady_state) || ((solve_algo >= 6 || solve_algo <= 8) && steady_state))
        {
          for (i = 0; i < u_count_init; i++)
            {
              SaveCode.read(reinterpret_cast<char *>(&eq), sizeof(eq));
              SaveCode.read(reinterpret_cast<char *>(&var), sizeof(var));
              SaveCode.read(reinterpret_cast<char *>(&lag), sizeof(lag));
              SaveCode.read(reinterpret_cast<char *>(&j), sizeof(j));
              IM_i[make_pair(make_pair(var - lag*Size, -lag), eq)] = j;
            }
        }
    }
  index_vara = (int *) mxMalloc(Size*(periods+y_kmin+y_kmax)*sizeof(int));
  for (j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_vara[j]), sizeof(*index_vara));
  if (periods+y_kmin+y_kmax > 1)
    for (i = 1; i < periods+y_kmin+y_kmax; i++)
      for (j = 0; j < Size; j++)
        index_vara[j+Size*i] = index_vara[j+Size*(i-1)]+y_size;
  index_equa = (int *) mxMalloc(Size*sizeof(int));
  for (j = 0; j < Size; j++)
    SaveCode.read(reinterpret_cast<char *>(&index_equa[j]), sizeof(*index_equa));
}

void
SparseMatrix::Simple_Init(int Size, map<pair<pair<int, int>, int>, int> &IM, bool &zero_solution)
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
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < Size; i++)
    {
      line_done[i] = 0;
      FNZE_C[i] = 0;
      FNZE_R[i] = 0;
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
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  double cum_abs_sum = 0;
  for (i = 0; i < Size; i++)
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
SparseMatrix::Init_Matlab_Sparse_Simple(int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, bool &zero_solution, mxArray *x0_m)
{
  int i, eq, var;
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
  for (i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  double cum_abs_sum = 0;
  for (i = 0; i < Size; i++)
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
  last_var = -1;
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
SparseMatrix::Init_Matlab_Sparse(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM, mxArray *A_m, mxArray *b_m, mxArray *x0_m)
{
  int t, i, eq, var, lag, ti_y_kmin, ti_y_kmax;
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
  mwIndex *Ai = mxGetIr(A_m);
  if (!Ai)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't allocate Ai index vector\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mwIndex *Aj = mxGetJc(A_m);
  if (!Aj)
    {
      ostringstream tmp;
      tmp << " in Init_Matlab_Sparse, can't allocate Aj index vector\n";
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
  for (i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
#ifdef DEBUG
  unsigned int max_nze = mxGetNzmax(A_m);
#endif
  unsigned int NZE = 0;
  int last_var = 0;
  for (i = 0; i < periods*Size; i++)
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
SparseMatrix::Init_GE(int periods, int y_kmin, int y_kmax, int Size, map<pair<pair<int, int>, int>, int> &IM)
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

  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < periods*Size; i++)
    {
      b[i] = 0;
      line_done[i] = 0;
    }
  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  for (i = 0; i < (periods+y_kmax+1)*Size; i++)
    {
      FNZE_C[i] = 0;
      FNZE_R[i] = 0;
      temp_NZE_C[i] = NULL;
      temp_NZE_R[i] = NULL;
      NbNZRow[i] = 0;
      NbNZCol[i] = 0;
    }

  //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(it4, ti_y_kmin, ti_y_kmax, eq, var, lag) schedule(dynamic)
  for (t = 0; t < periods; t++)
    {
      ti_y_kmin = -min(t, y_kmin);
      ti_y_kmax = min(periods-(t+1), y_kmax);
      it4 = IM.begin();
      eq = -1;
      //#pragma omp ordered
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
  mxFree(temp_NZE_R);
  mxFree(temp_NZE_C);
}

int
SparseMatrix::Get_u()
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
SparseMatrix::Delete_u(int pos)
{
  u_liste.push_back(pos);
}

void
SparseMatrix::Clear_u()
{
  u_liste.clear();
}

void
SparseMatrix::Print_u()
{
  for (unsigned int i = 0; i < u_liste.size(); i++)
    mexPrintf("%d ", u_liste[i]);
}

void
SparseMatrix::End_GE(int Size)
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
SparseMatrix::compare(int *save_op, int *save_opa, int *save_opaa, int beg_t, int periods, long int nop4,  int Size)
{
  long int i, j, nop = nop4/2, t, k;
  double r = 0.0;
  bool OK = true;
  t_save_op_s *save_op_s, *save_opa_s, *save_opaa_s;
  int *diff1, *diff2;
  diff1 = (int *) mxMalloc(nop*sizeof(int));
  diff2 = (int *) mxMalloc(nop*sizeof(int));
  int max_save_ops_first = -1;
  j = k = i = 0;
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
      //#pragma omp parallel for  num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) ordered private(j) schedule(dynamic)
      for (i = beg_t; i < periods; i++)
        {
          for (j = 0; j < Size; j++)
            {
              ///#pragma omp ordered
              pivot[i*Size+j] = pivot[(i-1)*Size+j]+Size;
            }
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
      double *up;
      for (t = 1; t < periods-beg_t-y_kmax; t++)
        {
          i = j = 0;
          while (i < nop4)
            {
              save_op_s = (t_save_op_s *) (&(save_op[i]));
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
      for (t = t1; t < periods_beg_t; t++)
        {
          i = j = 0;
          while (i < nop4)
            {
              save_op_s = (t_save_op_s *) (&(save_op[i]));
              if (save_op_s->lag < (periods_beg_t-t))
                {
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
SparseMatrix::complete(int beg_t, int Size, int periods, int *b)
{
  long int i, j, k, nop, nopa, nop1, cal_y, nb_var, pos, t, ti, max_var, min_var;
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
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
      nop += 4;
      for (k = 0; k < nb_var; k++)
        {
          save_code[nop] = IFMUL;
          save_code[nop+1] = index_vara[first->c_index]+cal_y;
          save_code[nop+2] = first->u_index;
          save_code[nop+3] = first->lag_index;
          if ((nop+3) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
          nop += 4;
          first = first->NZE_R_N;
        }
      save_code[nop] = IFADD;
      save_code[nop+1] = b[pos];
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
      if ((nop+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
      nop += 4;
      save_code[nop] = IFSTP;
      save_code[nop+1] = index_vara[j]+y_size*y_kmin;
      save_code[nop+2] = 0;
      save_code[nop+3] = 0;
      if ((nop+2) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop+2, size_of_save_code);
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
          if ((nop1+2) >= size_of_save_code)
            mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
          if ((nopa+1) >= size_of_diff)
            mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
          nopa += 2;
          nop1 += 4;
          first = first->NZE_R_N;
        }
      diff[nopa] = save_code[nop1+1]-(b[pos]);
      diff[nopa+1] = 0;
      if ((nop1+3) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
      nopa += 2;
      nop1 += 4;
      diff[nopa] = save_code[nop1+1]-(index_vara[j]+y_size*y_kmin);
      diff[nopa+1] = 0;
      if ((nop1+4) >= size_of_save_code)
        mexPrintf("out of save_code[%d] (bound=%d)\n", nop1+2, size_of_save_code);
      if ((nopa+1) >= size_of_diff)
        mexPrintf("out of diff[%d] (bound=%d)\n", nopa+2, size_of_diff);
      nopa += 2;
      nop1 += 4;
    }
  max_var = (periods+y_kmin)*y_size;
  min_var = y_kmin*y_size;
  for (t = periods+y_kmin-1; t >= beg_t+y_kmin; t--)
    {
      j = 0;
      ti = t-y_kmin-beg_t;
      for (i = 0; i < nop; i += 4)
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
SparseMatrix::bksub(int tbreak, int last_period, int Size, double slowc_l)
{
  NonZeroElem *first;
  int i, j, k;
  double yy;
  for (i = 0; i < y_size*(periods+y_kmin); i++)
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
SparseMatrix::simple_bksub(int it_, int Size, double slowc_l)
{
  int i, k;
  double yy;
  NonZeroElem *first;
  for (i = 0; i < y_size; i++)
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
SparseMatrix::CheckIt(int y_size, int y_kmin, int y_kmax, int Size, int periods, int iter)
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
SparseMatrix::Check_the_Solution(int periods, int y_kmin, int y_kmax, int Size, double *u, int *pivot, int *b)
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
SparseMatrix::substract_A_B(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  double *A_d = mxGetPr(A_m);
  unsigned int n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
  for (unsigned int j = 0; j < n_A; j++)
    for (unsigned int i = 0; i < m_A; i++)
      {
        unsigned int index = j*m_A+i;
        C_d[index] = A_d[index] - B_d[index];
      }
  return C_m;
}

mxArray *
SparseMatrix::Sparse_substract_A_SB(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_B = mxGetN(B_m);
  unsigned int m_B = mxGetM(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  unsigned int total_nze_B = B_j[n_B];
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
SparseMatrix::Sparse_substract_SA_SB(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  unsigned int total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  unsigned int n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  unsigned int total_nze_B = B_j[n_B];
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
      int A_row = A_i[nze_A];
      while (nze_B >= (unsigned int) B_j[B_col+1] && (nze_B < total_nze_B))
        B_col++;
      int B_row = B_i[nze_B];
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
SparseMatrix::mult_SAT_B(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  unsigned int n_B = mxGetN(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateDoubleMatrix(m_A, n_B, mxREAL);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_A = 0;
  for (unsigned int j = 0; j < n_B; j++)
    {
      for (unsigned int i = 0; i < n_A; i++)
        {
          double sum = 0;
          nze_A = A_j[i];
          while (nze_A < (unsigned int) A_j[i+1])
            {
              unsigned int i_A = A_i[nze_A];
              sum += A_d[nze_A++] * B_d[i_A];
            }
          C_d[j*n_A+i] = sum;
        }
    }
  return C_m;
}

mxArray *
SparseMatrix::Sparse_mult_SAT_B(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  unsigned int n_B = mxGetN(B_m);
  unsigned int m_B = mxGetM(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0, nze_A = 0;
  unsigned int C_col = 0;
  C_j[C_col] = 0;
  for (unsigned int j = 0; j < n_B; j++)
    {
      for (unsigned int i = 0; i < n_A; i++)
        {
          double sum = 0;
          nze_A = A_j[i];
          while (nze_A < (unsigned int) A_j[i+1])
            {
              unsigned int i_A = A_i[nze_A];
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
SparseMatrix::Sparse_mult_SAT_SB(mxArray *A_m, mxArray *B_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  double *A_d = mxGetPr(A_m);
  unsigned int n_B = mxGetN(B_m);
  mwIndex *B_i = mxGetIr(B_m);
  mwIndex *B_j = mxGetJc(B_m);
  double *B_d = mxGetPr(B_m);
  mxArray *C_m = mxCreateSparse(m_A, n_B, m_A*n_B, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_B = 0, nze_C = 0, nze_A = 0;
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
              unsigned int i_A = A_i[nze_A];
              unsigned int i_B = B_i[nze_B];
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
SparseMatrix::Sparse_transpose(mxArray *A_m)
{
  unsigned int n_A = mxGetN(A_m);
  unsigned int m_A = mxGetM(A_m);
  mwIndex *A_i = mxGetIr(A_m);
  mwIndex *A_j = mxGetJc(A_m);
  unsigned int total_nze_A = A_j[n_A];
  double *A_d = mxGetPr(A_m);
  mxArray *C_m = mxCreateSparse(n_A, m_A, total_nze_A, mxREAL);
  mwIndex *C_i = mxGetIr(C_m);
  mwIndex *C_j = mxGetJc(C_m);
  double *C_d = mxGetPr(C_m);
  unsigned int nze_C = 0, nze_A = 0;
  memset(C_j, 0, m_A);
  map<pair<unsigned int, unsigned int>, double> B2;
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
  for (map<pair<unsigned int, unsigned int>, double>::const_iterator it = B2.begin(); it != B2.end(); it++)
    {
      C_d[nze_C] = it->second;
      C_i[nze_C++] = it->first.second;
    }
  return C_m;
}

void
SparseMatrix::Solve_Matlab_Relaxation(mxArray *A_m, mxArray *b_m, unsigned int Size, double slowc_l, bool is_two_boundaries, int  it_)
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
  unsigned int max_nze = A_m_j[Size*periods];
  unsigned int nze = 0;
  unsigned int var = A_m_j[nze];
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
  unsigned int eq = 0;
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
          unsigned int B_inv_nze = B_inv_j[Size];
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
SparseMatrix::Solve_Matlab_LU_UMFPack(mxArray *A_m, mxArray *b_m, int Size, double slowc_l, bool is_two_boundaries, int  it_)
{
  int n = mxGetM(A_m);
  mxArray *z;
  mxArray *rhs[2];
  rhs[0] = A_m;
  rhs[1] = b_m;
  mexCallMATLAB(1, &z, 2, rhs, "mldivide");
  double *res = mxGetPr(z);
  if (is_two_boundaries)
    for (int i = 0; i < n; i++)
      {
        int eq = index_vara[i+Size*y_kmin];
        double yy = -(res[i] + y[eq]);
        direction[eq] = yy;
        y[eq] += slowc_l * yy;
      }
  else
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
SparseMatrix::Solve_Matlab_GMRES(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, bool steady_state, mxArray *x0_m)
{
  ostringstream tmp;
#ifdef OCTAVE_MEX_FILE
  if (steady_state)
    tmp << " GMRES method is not implemented in Octave. You cannot use solve_algo=7, change solve_algo.\n";
  else
    tmp << " GMRES method is not implemented in Octave. You cannot use stack_solve_algo=2, change stack_solve_algo.\n";
  throw FatalExceptionHandling(tmp.str());
#endif
  int n = mxGetM(A_m);
  /*const char *field_names[] = {"droptol"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, 1, field_names);
  mxSetFieldByNumber(Setup, 0, 0, mxCreateDoubleScalar(lu_inc_tol));*/
  
  const char *field_names[] = {"type","droptol"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, 2, field_names);
  mxSetFieldByNumber(Setup, 0, 1, mxCreateDoubleScalar(lu_inc_tol));
  mxSetFieldByNumber(Setup, 0, 0, mxCreateString("ilutp"));
  
  mxArray *lhs0[2];
  mxArray *rhs0[2];
  rhs0[0] = A_m;
  rhs0[1] = Setup;

  if (mexCallMATLAB(2, lhs0, 2, rhs0, "ilu"))
    {
      tmp << " In GMRES, the incomplet LU decomposition (ilu) ahs failed.\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mxArray *L1 = lhs0[0];
  mxArray *U1 = lhs0[1];
  /*[za,flag1] = gmres(g1a,b,Blck_size,1e-6,Blck_size*periods,L1,U1);*/
  mxArray *rhs[8];
  rhs[0] = A_m;
  rhs[1] = b_m;
  rhs[2] = mxCreateDoubleScalar(Size);
  rhs[3] = mxCreateDoubleScalar(1e-6);
  rhs[4] = mxCreateDoubleScalar(n);
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
        for (int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
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
SparseMatrix::Solve_Matlab_BiCGStab(mxArray *A_m, mxArray *b_m, int Size, double slowc, int block, bool is_two_boundaries, int it_, mxArray *x0_m, bool steady_state)
{
  ostringstream tmp;
  unsigned int n = mxGetM(A_m);
  /*[L1, U1]=luinc(g1a,luinc_tol);*/
  /*const char *field_names[] = {"droptol"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, 1, field_names);
  mxSetFieldByNumber(Setup, 0, 0, mxCreateDoubleScalar(lu_inc_tol));*/
  
  const char *field_names[] = {"type","droptol"};
  mwSize dims[1] = { 1 };
  mxArray *Setup = mxCreateStructArray(1, dims, 2, field_names);
  mxSetFieldByNumber(Setup, 0, 1, mxCreateDoubleScalar(lu_inc_tol));
  mxSetFieldByNumber(Setup, 0, 0, mxCreateString("ilutp"));
  
  mxArray *lhs0[2];
  mxArray *rhs0[2];
  rhs0[0] = A_m;
  rhs0[1] = Setup;
  
  if (mexCallMATLAB(2, lhs0, 2, rhs0, "ilu"))
    {
      tmp << " In GMRES, the incomplet LU decomposition (ilu) ahs failed.\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mxArray *L1 = lhs0[0];
  mxArray *U1 = lhs0[1];
  double flags = 2;
  mxArray *z;
  if (steady_state)  /*Octave BicStab algorihtm involves a 0 division in case of a preconditionner equal to the LU decomposition of A matrix*/
    {
      mxArray *res = mult_SAT_B(Sparse_transpose(A_m), x0_m);
      double *resid = mxGetPr(res);
      double *b = mxGetPr(b_m);
      for (unsigned int i = 0; i < n; i++)
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
      for (unsigned int i = 0; i < n; i++)
        phat[i] = x0[i] + phat[i];

      /*Check the solution*/
      res = mult_SAT_B(Sparse_transpose(A_m), z);
      resid = mxGetPr(res);
      double cum_abs = 0;
      for (unsigned int i = 0; i < n; i++)
        {
          resid[i] = b[i] - resid[i];
          cum_abs += fabs(resid[i]);
        }
      //mexPrintf("cum_abs=%g\n", cum_abs);
      if (cum_abs > 1e-7)
        flags = 2;
      else
        flags = 0;
      mxDestroyArray(res);
    }
  //else
  if (flags == 2)
    {
      /*[za,flag1] = bicgstab(g1a,b,1e-6,Blck_size*periods,L1,U1);*/
      mxArray *rhs[7];
      rhs[0] = A_m;
      rhs[1] = b_m;
      rhs[2] = mxCreateDoubleScalar(1e-6);
      rhs[3] = mxCreateDoubleScalar(n);
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
  /*mexPrintf("z");
    mexCallMATLAB(0, NULL, 1, &z, "disp");*/
  mxDestroyArray(rhs0[1]);

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
        for (unsigned int i = 0; i < n; i++)
          {
            int eq = index_vara[i+Size*y_kmin];
            double yy = -(res[i] + y[eq]);
            direction[eq] = yy;
            y[eq] += slowc * yy;
          }
      else
        for (unsigned int i = 0; i < n; i++)
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
SparseMatrix::Singular_display(int block, int Size, bool steady_state, it_code_type it_code)
{
  bool zero_solution;
  Simple_Init(Size, IM_i, zero_solution);
  NonZeroElem *first;
  mxArray *rhs[1];
  rhs[0] = mxCreateDoubleMatrix(Size, Size, mxREAL);
  double *pind;
  pind = mxGetPr(rhs[0]);
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
                    else
                      if (rr != 1)
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
SparseMatrix::Solve_ByteCode_Sparse_GaussianElimination(int Size, int blck, bool steady_state, int it_)
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
      l = 0; N_max = 0;
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
      //#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
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
SparseMatrix::Solve_ByteCode_Symbolic_Sparse_GaussianElimination(int Size, bool symbolic, int Block_number)
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

  for (int t = 0; t < periods; t++)
    {
      if (record && symbolic)
        {
          if (save_op)
            {
              mxFree(save_op);
              save_op = NULL;
            }
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
          if (symbolic)
            {
              if (record)
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
                }
              nop += 2;
            }
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
          NonZeroElem **bb;
          bb = (NonZeroElem **) mxMalloc(nb_var*sizeof(first));
          for (int j = 0; j < nb_var; j++)
            {
              bb[j] = first;
              first = first->NZE_R_N;
            }

          for (int j = 0; j < nb_var; j++)
            {
              first = bb[j];
              u[first->u_index] /= piv;
              if (symbolic)
                {
                  if (record)
                    {
                      if (nop+j*2+1 >= nopa)
                        {
                          nopa = long (mem_increasing_factor*(double) nopa);
                          save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                        }
                      save_op_s = (t_save_op_s *) (&(save_op[nop+j*2]));
                      save_op_s->operat = IFDIV;
                      save_op_s->first = first->u_index;
                      save_op_s->lag = first->lag_index;
                    }
                }
            }
          mxFree(bb);
          nop += nb_var*2;
          u[b[pivj]] /= piv;
          if (symbolic)
            {
              if (record)
                {
                  if (nop+1 >= nopa)
                    {
                      nopa = long (mem_increasing_factor*(double) nopa);
                      save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                    }
                  save_op_s = (t_save_op_s *) (&(save_op[nop]));
                  save_op_s->operat = IFDIV;
                  save_op_s->first = b[pivj];
                  save_op_s->lag = 0;
                }
              nop += 2;
            }
          /*substract the elements on the non treated lines*/
          nb_eq = At_Col(i, &first);
          NonZeroElem *first_piva;
          int nb_var_piva = At_Row(pivj, &first_piva);

          NonZeroElem **bc;
          bc = (NonZeroElem **) mxMalloc(nb_eq*sizeof(first));
          int nb_eq_todo = 0;
          for (int j = 0; j < nb_eq && first; j++)
            {
              if (!line_done[first->r_index])
                bc[nb_eq_todo++] = first;
              first = first->NZE_C_N;
            }
          //#pragma omp parallel for num_threads(2) shared(nb_var_piva, first_piva, nopa, nop, save_op, record)
          for (int j = 0; j < nb_eq_todo; j++)
            {
              t_save_op_s *save_op_s_l;
              first = bc[j];
              int row = first->r_index;
              double first_elem = u[first->u_index];
              if (symbolic)
                {
                  if (record)
                    {
                      if (nop+1 >= nopa)
                        {
                          nopa = long (mem_increasing_factor*(double) nopa);
                          save_op = (int *) mxRealloc(save_op, nopa*sizeof(int));
                        }
                      save_op_s_l = (t_save_op_s *) (&(save_op[nop]));
                      save_op_s_l->operat = IFLD;
                      save_op_s_l->first = first->u_index;
                      save_op_s_l->lag = abs(first->lag_index);
                    }
                  nop += 2;
                }

              int nb_var_piv = nb_var_piva;
              NonZeroElem *first_piv = first_piva;
              NonZeroElem *first_sub;
              int nb_var_sub = At_Row(row, &first_sub);
              int l_sub = 0;
              int l_piv = 0;
              int sub_c_index = first_sub->c_index;
              int piv_c_index = first_piv->c_index;
              int tmp_lag = first_sub->lag_index;
              while (l_sub < nb_var_sub || l_piv < nb_var_piv)
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
                      if (symbolic)
                        {
                          if (record)
                            {
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
                            }
                          nop += 3;
                        }
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
                          Delete(first_sub->r_index, first_sub->c_index);
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
                          if (symbolic)
                            {
                              if (record)
                                {
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
                                }
                              nop += 3;
                            }
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

              if (symbolic)
                {
                  if (record)
                    {
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
                    }
                  nop += 3;
                }
            }
          mxFree(bc);
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
              if (save_opa && save_opaa)
                {
                  if (compare(save_op, save_opa, save_opaa, t, periods, nop, Size))
                    {
                      tbreak = t;
                      tbreak_g = tbreak;
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
                  save_opaa = (int *) mxMalloc(nop1*sizeof(int));
                  memcpy(save_opaa, save_opa, nop1*sizeof(int));
                }
              if (save_opa)
                {
                  mxFree(save_opa);
                  save_opa = NULL;
                }
              save_opa = (int *) mxMalloc(nop*sizeof(int));
              memcpy(save_opa, save_op, nop*sizeof(int));
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
    }
  mxFree(piv_v);
  mxFree(pivj_v);
  mxFree(pivk_v);
  mxFree(NR);
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
  for (int i = 0; i < y_size*(periods+y_kmin); i++)
    ya[i] = y[i];
  slowc_save = slowc;
  bksub(tbreak, last_period, Size, slowc_lbx);
  End_GE(Size);
}

bool
SparseMatrix::Simulate_Newton_One_Boundary(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, bool print_it, bool cvg, int &iter, bool steady_state, int stack_solve_algo, int solve_algo)
{
  int i, j;
  mxArray *b_m = NULL, *A_m = NULL, *x0_m = NULL;
  Clear_u();
  error_not_printed = true;
  bool singular_system = false;
  u_count_alloc_save = u_count_alloc;
  if (isnan(res1) || isinf(res1) || (res2 > 12*g0 && iter > 0))
    {
      if (iter == 0 || fabs(slowc_save) < 1e-8)
        {
          for (j = 0; j < y_size; j++)
            {
#ifdef DEBUG
              bool select = false;
#endif
              for (int i = 0; i < Size; i++)
                if (j == index_vara[i])
                  {
#ifdef DEBUG
                    select = true;
#endif
                    break;
                  }
#ifdef DEBUG
              if (select)
                mexPrintf("-> variable %s (%d) at time %d = %f direction = %f\n", get_variable(eEndogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
              else
                mexPrintf("   variable %s (%d) at time %d = %f direction = %f\n", get_variable(eEndogenous, j).c_str(), j+1, it_, y[j+it_*y_size], direction[j+it_*y_size]);
#endif
            }
          if (steady_state)
            {
              if (iter == 0)
                mexPrintf(" the initial values of endogenous variables are too far from the solution.\nChange them!\n");
              else
                mexPrintf(" dynare cannot improve the simulation in block %d at time %d (variable %d)\n", blck+1, it_+1, index_vara[max_res_idx]+1);
              mexEvalString("drawnow;");
              return singular_system;
            }
          else
            {
              ostringstream tmp;
              if (iter == 0)
                tmp << " in Simulate_Newton_One_Boundary, The initial values of endogenous variables are too far from the solution.\nChange them!\n";
              else
                tmp << " in Simulate_Newton_One_Boundary, Dynare cannot improve the simulation in block " << blck+1 << " at time " << it_+1 << " (variable " << index_vara[max_res_idx]+1 << "%d)\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      if (!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0)))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max(-gp0 / (2 * (res2 - g0 - gp0)), 0.1);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min(-b + sqrt(b*b - 3 * a * gp0) / (3 * a), 0.5 * slowc_save), 0.1 * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration++;
        }
      else
        {
          prev_slowc_save = slowc_save;
          slowc_save /= 1.1;
        }
      if (print_it)
        mexPrintf("Error: Simulation diverging, trying to correct it using slowc=%f\n", slowc_save);
      for (i = 0; i < y_size; i++)
        y[i+it_*y_size] = ya[i+it_*y_size] + slowc_save*direction[i+it_*y_size];
      iter--;
      return singular_system;
    }
  if (cvg)
    {
      return singular_system;
    }
  if (print_it)
    {
      //mexPrintf("solwc=%f g0=%f res2=%f glambda2=%f\n",slowc_save,g0, res2, glambda2);
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
              mexPrintf("MODEL STEADY STATE: (method=GMRES)\n");
              break;
            case 8:
              mexPrintf("MODEL STEADY STATE: (method=BiCGStab)\n");
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
    Simple_Init(Size, IM_i, zero_solution);
  else
    {
      b_m = mxCreateDoubleMatrix(Size, 1, mxREAL);
      if (!b_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_Newton_One_Boundary, can't allocate b_m vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      A_m = mxCreateSparse(Size, Size, min(int (IM_i.size()*2), Size*Size), mxREAL);
      if (!A_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_Newton_One_Boundary, can't allocate A_m matrix\n";
          throw FatalExceptionHandling(tmp.str());
        }
      x0_m = mxCreateDoubleMatrix(Size, 1, mxREAL);
      if (!x0_m)
        {
          ostringstream tmp;
          tmp << " in Simulate_Newton_One_Boundary, can't allocate x0_m vector\n";
          throw FatalExceptionHandling(tmp.str());
        }
      Init_Matlab_Sparse_Simple(Size, IM_i, A_m, b_m, zero_solution, x0_m);
    }
  if (zero_solution)
    {
      for (int i = 0; i < Size; i++)
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
        singular_system = Solve_ByteCode_Sparse_GaussianElimination(Size, blck, steady_state, it_);
      else if ((solve_algo == 7 && steady_state) || (stack_solve_algo == 2 && !steady_state))
        Solve_Matlab_GMRES(A_m, b_m, Size, slowc, blck, false, it_, steady_state, x0_m);
      else if ((solve_algo == 8 && steady_state) || (stack_solve_algo == 3 && !steady_state))
        Solve_Matlab_BiCGStab(A_m, b_m, Size, slowc, blck, false, it_, x0_m, steady_state);
      else if ((solve_algo == 6 && steady_state) || ((stack_solve_algo == 0 || stack_solve_algo == 1) && !steady_state))
        Solve_Matlab_LU_UMFPack(A_m, b_m, Size, slowc, false, it_);
    }
  return singular_system;
}

void
SparseMatrix::Simulate_Newton_Two_Boundaries(int blck, int y_size, int it_, int y_kmin, int y_kmax, int Size, int periods, bool print_it, bool cvg, int &iter, int minimal_solving_periods, int stack_solve_algo, unsigned int endo_name_length, char *P_endo_names)
{
  if (start_compare == 0)
    start_compare = y_kmin;
  u_count_alloc_save = u_count_alloc;
  clock_t t1 = clock();
  nop1 = 0;
  error_not_printed = true;
  mxArray *b_m = NULL, *A_m = NULL, *x0_m = NULL;
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
          //Error << filename << " stopped";
          throw FatalExceptionHandling(Error.str());
        }
      if (!(isnan(res1) || isinf(res1)) && !(isnan(g0) || isinf(g0)) && (stack_solve_algo == 4 || stack_solve_algo == 5))
        {
          if (try_at_iteration == 0)
            {
              prev_slowc_save = slowc_save;
              slowc_save = max(-gp0 / (2 * (res2 - g0 - gp0)), 0.1);
            }
          else
            {
              double t1 = res2 - gp0 * slowc_save - g0;
              double t2 = glambda2 - gp0 * prev_slowc_save - g0;
              double a = (1/(slowc_save * slowc_save) * t1 - 1/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              double b = (-prev_slowc_save/(slowc_save * slowc_save) * t1 + slowc_save/(prev_slowc_save * prev_slowc_save) * t2) / (slowc_save - prev_slowc_save);
              prev_slowc_save = slowc_save;
              slowc_save = max(min(-b + sqrt(b*b - 3 * a * gp0) / (3 * a), 0.5 * slowc_save), 0.1 * slowc_save);
            }
          glambda2 = res2;
          try_at_iteration++;
          if (slowc_save <= 0.1)
            {
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
              mexPrintf("MODEL SIMULATION: (method=GMRES)\n");
              break;
            case 3:
              mexPrintf("MODEL SIMULATION: (method=BiCGStab)\n");
              break;
            case 4:
              mexPrintf("MODEL SIMULATION: (method=Sparse LU & optimal path length)\n");
              break;
            case 5:
              mexPrintf("MODEL SIMULATION: (method=ByteCode own solver)\n");
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
          A_m = mxCreateSparse(periods*Size, periods*Size, IM_i.size()* periods*2, mxREAL);
          if (!A_m)
            {
              ostringstream tmp;
              tmp << " in Simulate_Newton_Two_Boundaries, can't allocate A_m matrix\n";
              throw FatalExceptionHandling(tmp.str());
            }
          Init_Matlab_Sparse(periods, y_kmin, y_kmax, Size, IM_i, A_m, b_m, x0_m);
        }

      if (stack_solve_algo == 0 || stack_solve_algo == 4)
        Solve_Matlab_LU_UMFPack(A_m, b_m, Size, slowc, true, 0);
      else if (stack_solve_algo == 1)
        Solve_Matlab_Relaxation(A_m, b_m, Size, slowc, true, 0);
      else if (stack_solve_algo == 2)
        Solve_Matlab_GMRES(A_m, b_m, Size, slowc, blck, true, 0, false, x0_m);
      else if (stack_solve_algo == 3)
        Solve_Matlab_BiCGStab(A_m, b_m, Size, slowc, blck, true, 0, x0_m, false);
      else if (stack_solve_algo == 5)
        Solve_ByteCode_Symbolic_Sparse_GaussianElimination(Size, symbolic, blck);
    }
  if (print_it)
    {
      clock_t t2 = clock();
      mexPrintf("(** %f milliseconds **)\n", 1000.0*(double (t2) - double (t1))/double (CLOCKS_PER_SEC));
      mexEvalString("drawnow;");
    }

  time00 = clock();
  if (tbreak_g == 0)
    tbreak_g = periods;
  return;
}

void
SparseMatrix::fixe_u(double **u, int u_count_int, int max_lag_plus_max_lead_plus_1)
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

