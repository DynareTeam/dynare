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

#ifndef INTERPRETER_HH_INCLUDED
#define INTERPRETER_HH_INCLUDED

#include <stack>
#include <vector>
#include <string>
#include <cmath>
#define BYTE_CODE
#include "CodeInterpreter.hh"
#include "SparseMatrix.hh"
#include "Evaluate.hh"
#ifdef LINBCG
# include "linbcg.hh"
#endif
#ifndef DEBUG_EX
# include <dynmex.h>
#else
# include "mex_interface.hh"
#endif

//#define DEBUGC

using namespace std;


class Interpreter : public dynSparseMatrix
{
private:
protected:
  void evaluate_a_block();
  int simulate_a_block();
  void print_a_block();
public:
  ~Interpreter();
  Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
              double *direction_arg, size_t y_size_arg,
              size_t nb_row_x_arg, size_t nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
              int maxit_arg_, double solve_tolf_arg, size_t size_of_direction_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg,
              string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
              bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg,
              bool steady_state_arg, bool print_it_arg
#ifdef CUDA
              , const int CUDA_device, cublasHandle_t cublas_handle_arg, cusparseHandle_t cusparse_handle_arg, cusparseMatDescr_t descr_arg
#endif
              );
  bool compute_blocks(string file_name, string bin_basename, bool evaluate, int block, int &nb_blocks);

  inline mxArray *
  get_jacob(int block_num)
  {
    return jacobian_block[block_num];
  };
  inline mxArray *
  get_jacob_exo(int block_num)
  {
    return jacobian_exo_block[block_num];
  };
  inline mxArray *
  get_jacob_exo_det(int block_num)
  {
    return jacobian_det_exo_block[block_num];
  };
  inline mxArray *
  get_jacob_other_endo(int block_num)
  {
    return jacobian_other_endo_block[block_num];
  };
  inline vector<double>
  get_residual()
  {
    return residual;
  };
  inline mxArray *
  get_Temporary_Terms()
  {
    return GlobalTemporaryTerms;
  };
};

#endif
