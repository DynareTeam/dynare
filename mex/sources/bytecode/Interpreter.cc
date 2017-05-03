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

#include <cstring>
#include <sstream>
#include <algorithm>
#include "Interpreter.hh"
#define BIG 1.0e+8;
#define SMALL 1.0e-5;
///#define DEBUG

Interpreter::~Interpreter()
{
}

Interpreter::Interpreter(double *params_arg, double *y_arg, double *ya_arg, double *x_arg, double *steady_y_arg, double *steady_x_arg,
                         double *direction_arg, size_t y_size_arg,
                         size_t nb_row_x_arg, size_t nb_row_xd_arg, int periods_arg, int y_kmin_arg, int y_kmax_arg,
                         int maxit_arg_, double solve_tolf_arg, size_t size_of_direction_arg, double slowc_arg, int y_decal_arg, double markowitz_c_arg,
                         string &filename_arg, int minimal_solving_periods_arg, int stack_solve_algo_arg, int solve_algo_arg,
                         bool global_temporary_terms_arg, bool print_arg, bool print_error_arg, mxArray *GlobalTemporaryTerms_arg,
                         bool steady_state_arg, bool print_it_arg, int col_x_arg
#ifdef CUDA
                         , const int CUDA_device_arg, cublasHandle_t cublas_handle_arg, cusparseHandle_t cusparse_handle_arg, cusparseMatDescr_t descr_arg
#endif
                         )
                         : dynSparseMatrix(y_size_arg, y_kmin_arg, y_kmax_arg, print_it_arg, steady_state_arg, periods_arg, minimal_solving_periods_arg, slowc_arg
#ifdef CUDA
                                        , CUDA_device_arg, cublas_handle_arg, cusparse_handle_arg, descr_arg
#endif
                                        )
{
  params = params_arg;
  y = y_arg;
  ya = ya_arg;
  x = x_arg;
  steady_y = steady_y_arg;
  steady_x = steady_x_arg;
  direction = direction_arg;
  //y_size = y_size_arg;
  nb_row_x = nb_row_x_arg;
  nb_row_xd = nb_row_xd_arg;
  periods = periods_arg;
  //y_kmax = y_kmax_arg;
  //y_kmin = y_kmin_arg;
  maxit_ = maxit_arg_;
  solve_tolf = solve_tolf_arg;
  size_of_direction = size_of_direction_arg;
  slowc = slowc_arg;
  slowc_save = slowc;
  y_decal = y_decal_arg;
  markowitz_c = markowitz_c_arg;
  filename = filename_arg;
  T = NULL;
  minimal_solving_periods = minimal_solving_periods_arg;
  stack_solve_algo = stack_solve_algo_arg;
  solve_algo = solve_algo_arg;
  global_temporary_terms = global_temporary_terms_arg;
  print = print_arg;
  col_x = col_x_arg;
  GlobalTemporaryTerms = GlobalTemporaryTerms_arg;
  print_error = print_error_arg;
  //steady_state = steady_state_arg;
  print_it = print_it_arg;

}

void
Interpreter::evaluate_a_block(bool initialization)
{
  it_code_type begining;

  switch (type)
    {
    case EVALUATE_FORWARD:
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state, */false);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state, */false);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case SOLVE_FORWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      test_mxMalloc(g1, __LINE__, __FILE__, __func__, size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_FORWARD_COMPLETE:
      if (initialization)
         {
           fixe_u(&u, u_count_int, u_count_int);
           Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
         }
#ifdef DEBUG
      mexPrintf("in SOLVE_FORWARD_COMPLETE r = mxMalloc(%d*sizeof(double))\n", size);
#endif
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[it_*y_size+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case EVALUATE_BACKWARD:
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
          if (block >= 0)
            for (int j = 0; j < size; j++)
              residual[j] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
          else
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = y[Block_Contain[j].Variable] - ya[Block_Contain[j].Variable];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
              if (block >= 0)
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+Block_Contain[j].Equation] = y[it_*y_size+Block_Contain[j].Variable] - ya[it_*y_size+Block_Contain[j].Variable];
            }
        }
      break;
    case SOLVE_BACKWARD_SIMPLE:
      g1 = (double *) mxMalloc(size*size*sizeof(double));
      test_mxMalloc(g1, __LINE__, __FILE__, __func__, size*size*sizeof(double));
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(g1);
      mxFree(r);
      break;
    case SOLVE_BACKWARD_COMPLETE:
      if (initialization)
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      if (steady_state)
        {
          compute_block_time(0, true, /*block_num, size, steady_state,*/ false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[j] = r[j];
        }
      else
        {
          begining = it_code;
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              Per_y_ = it_*y_size;
              compute_block_time(0, true, /*block_num, size, steady_state, */false);
              if (block < 0)
                for (int j = 0; j < size; j++)
                  residual[Per_y_+Block_Contain[j].Equation] = r[j];
              else
                for (int j = 0; j < size; j++)
                  residual[it_*size+j] = r[j];
            }
        }
      mxFree(r);
      break;
    case SOLVE_TWO_BOUNDARIES_SIMPLE:
    case SOLVE_TWO_BOUNDARIES_COMPLETE:
      if (initialization)
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, periods, y_kmin, y_kmax, true, stack_solve_algo, solve_algo);
        }
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      begining = it_code;
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          it_code = begining;
          compute_block_time(Per_u_, true, /*block_num, size, steady_state,*/ false);
          if (block < 0)
            for (int j = 0; j < size; j++)
              residual[it_*y_size+Block_Contain[j].Equation] = r[j];
          else
            for (int j = 0; j < size; j++)
              residual[it_*size+j] = r[j];
        }
      mxFree(r);
      break;
    }
}



int
Interpreter::simulate_a_block(vector_table_conditional_local_type vector_table_conditional_local)
{
  it_code_type begining;
  max_res = 0;
  max_res_idx = 0;
  bool cvg;
  double *y_save;
#ifdef DEBUG
  mexPrintf("simulate_a_block type = %d, periods=%d, y_kmin=%d, y_kmax=%d\n", type, periods, y_kmin, y_kmax);
  mexEvalString("drawnow;");
#endif
  switch (type)
    {
    case EVALUATE_FORWARD:
#ifdef DEBUG
      mexPrintf("EVALUATE_FORWARD\n");
      mexEvalString("drawnow;");
#endif
        evaluate_over_periods(true);
      break;
    case EVALUATE_BACKWARD:
#ifdef DEBUG
      mexPrintf("EVALUATE_BACKWARD\n");
      mexEvalString("drawnow;");
#endif
        evaluate_over_periods(false);
      break;
    case SOLVE_FORWARD_SIMPLE:
#ifdef DEBUG
      mexPrintf("SOLVE_FORWARD_SIMPLE size=%d\n", size);
      mexEvalString("drawnow;");
#endif
      solve_simple_over_periods(true);
      break;
    case SOLVE_BACKWARD_SIMPLE:
#ifdef DEBUG
      mexPrintf("SOLVE_BACKWARD_SIMPLE\n");
      mexEvalString("drawnow;");
#endif
      solve_simple_over_periods(false);
      break;
    case SOLVE_FORWARD_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_FORWARD_COMPLETE\n");
      mexEvalString("drawnow;");
#endif
      if (vector_table_conditional_local.size())
        {
          it_code_type curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      start_code = it_code;
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(true);

      mxFree(u);
      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      End_Solver();
      break;
    case SOLVE_BACKWARD_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_BACKWARD_COMPLETE\n");
      mexEvalString("drawnow;");
#endif
      if (vector_table_conditional_local.size())
        {
          it_code_type curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, 1, 0, 0, false, stack_solve_algo, solve_algo);
        }
      start_code = it_code;
      Per_u_ = 0;

      Simulate_Newton_One_Boundary(false);

      mxFree(index_equa);
      mxFree(index_vara);
      memset(direction, 0, size_of_direction);
      mxFree(u);
      End_Solver();
      break;
    case SOLVE_TWO_BOUNDARIES_SIMPLE:
    case SOLVE_TWO_BOUNDARIES_COMPLETE:
#ifdef DEBUG
      mexPrintf("SOLVE_TWO_BOUNDARIES\n");
      mexEvalString("drawnow;");
#endif
      if (steady_state)
        {
          mexPrintf("SOLVE_TWO_BOUNDARIES in a steady state model: impossible case\n");
          return ERROR_ON_EXIT;
        }
      if (vector_table_conditional_local.size())
        {
          it_code_type curr_it_code = it_code;
          evaluate_a_block(true);
          it_code = curr_it_code;
        }
      else
        {
          fixe_u(&u, u_count_int, u_count_int);
          Read_SparseMatrix(bin_base_name, size, periods, y_kmin, y_kmax, true, stack_solve_algo, solve_algo);
        }
      u_count = u_count_int*(periods+y_kmax+y_kmin);
      r = (double *) mxMalloc(size*sizeof(double));
      test_mxMalloc(r, __LINE__, __FILE__, __func__, size*sizeof(double));
      res = (double *) mxMalloc(size*periods*sizeof(double));
      test_mxMalloc(res, __LINE__, __FILE__, __func__, size*periods*sizeof(double));
      y_save = (double *) mxMalloc(y_size*sizeof(double)*(periods+y_kmax+y_kmin));
      test_mxMalloc(y_save, __LINE__, __FILE__, __func__, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
      start_code = it_code;
      iter = 0;
      if (!is_linear)
        {
          cvg = false;
          glambda2 = g0 = very_big;
          try_at_iteration = 0;
          int u_count_saved = u_count;
          while (!(cvg || (iter > maxit_)))
            {
              res2 = 0;
              res1 = 0;
              max_res = 0;
              max_res_idx = 0;
              memcpy(y_save, y, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
              if (vector_table_conditional_local.size())
                {
                  for (vector_table_conditional_local_type::iterator it1 = vector_table_conditional_local.begin(); it1 != vector_table_conditional_local.end() ; it1++)
                    {
                      if (it1->is_cond)
                        {
                          //mexPrintf("y[%d] = %f\n", it1->var_endo + y_kmin * size, y[it1->var_endo + y_kmin * size]);
                          y[it1->var_endo + y_kmin * size] = it1->constrained_value;
                        }

                    }
                }
              compute_complete_2b(false, &res1, &res2, &max_res, &max_res_idx);
              end_code = it_code;
              if (!(isnan(res1) || isinf(res1)))
                cvg = (max_res < solve_tolf);
              if (isnan(res1) || isinf(res1) || (stack_solve_algo == 4 && iter > 0))
                  memcpy(y, y_save, y_size*sizeof(double)*(periods+y_kmax+y_kmin));
              u_count = u_count_saved;
              int prev_iter = iter;
              Simulate_Newton_Two_Boundaries(block_num, symbol_table_endo_nbr, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names, vector_table_conditional_local);
              iter++;
              if (iter > prev_iter)
                {
                  g0 = res2;
                  gp0 = -res2;
                  try_at_iteration = 0;
                  slowc_save = slowc;
                }
            }
          if (!cvg)
            {
              ostringstream tmp;
              tmp << " in Solve two boundaries, convergence not achieved in block " << block_num+1 << ", after " << iter << " iterations\n";
              throw FatalExceptionHandling(tmp.str());
            }
        }
      else
        {
          res1 = 0;
          res2 = 0;
          max_res = 0; max_res_idx = 0;

          compute_complete_2b(false, &res1, &res2, &max_res, &max_res_idx);
          end_code = it_code;

          cvg = false;
          Simulate_Newton_Two_Boundaries(block_num, symbol_table_endo_nbr, y_kmin, y_kmax, size, periods, cvg, minimal_solving_periods, stack_solve_algo, endo_name_length, P_endo_names, vector_table_conditional_local);
          max_res = 0; max_res_idx = 0;
        }
      it_code = end_code;
      if (r)
        mxFree(r);
      if (y_save)
        mxFree(y_save);
      if (u)
        mxFree(u);
      if (index_vara)
        mxFree(index_vara);
      if (index_equa)
        mxFree(index_equa);
      if (res)
        mxFree(res);
      memset(direction, 0, size_of_direction);
      End_Solver();
      break;
    default:
      ostringstream tmp;
      tmp << " in simulate_a_block, Unknown type = " << type << "\n";
      throw FatalExceptionHandling(tmp.str());
      return ERROR_ON_EXIT;
    }
  return NO_ERROR_ON_EXIT;
}

void
Interpreter::print_a_block()
{
  it_code_type begining;
  if (block < 0)
    mexPrintf("\nBlock %d\n", block_num+1);
  else
    mexPrintf("\nBlock %d\n", block+1);
  mexPrintf("----------\n");
  if (steady_state)
    residual = vector<double>(size);
  else
    residual = vector<double>(size*(periods+y_kmin));
  bool go_on = true;
  bool space = false;
  while (go_on)
    {
      if (it_code->first == FENDBLOCK)
        {
          go_on = false;
          it_code++;
        }
      else
        {
          string s = print_expression(it_code, false, size, block_num, steady_state, Per_u_, it_, it_code, false);
          if (s == "if (evaluate)" || s == "else")
            space = false;
          if (s.length() > 0)
            {
              if (space)
                mexPrintf("  %s\n", s.c_str());
              else
                mexPrintf("%s\n", s.c_str());
              mexEvalString("drawnow;");
            }
          if (s == "if (evaluate)" || s == "else")
            space = true;
        }
    }
}

void
Interpreter::ReadCodeFile(string file_name, CodeLoad &code)
{
  if (steady_state)
    file_name += "_static";
  else
    file_name += "_dynamic";

  //First read and store in memory the code
  code_liste = code.get_op_code(file_name);
  EQN_block_number = code.get_block_number();
  if (!code_liste.size())
    {
      ostringstream tmp;
      tmp << " in compute_blocks, " << file_name.c_str() << " cannot be opened\n";
      throw FatalExceptionHandling(tmp.str());
    }
  if (block >= (int) code.get_block_number())
    {
      ostringstream tmp;
      tmp << " in compute_blocks, input argument block = " << block+1 << " is greater than the number of blocks in the model (" << code.get_block_number() << " see M_.block_structure_stat.block)\n";
      throw FatalExceptionHandling(tmp.str());
    }

}


void
Interpreter::check_for_controlled_exo_validity(FBEGINBLOCK_ *fb,vector<s_plan> sconstrained_extended_path)
{
  vector<unsigned int> exogenous = fb->get_exogenous();
  vector<int> endogenous = fb->get_endogenous();
  for (vector<s_plan>::iterator it = sconstrained_extended_path.begin(); it != sconstrained_extended_path.end(); it++)
    {
      if ((find(endogenous.begin(), endogenous.end(), it->exo_num) != endogenous.end()) &&  (find(exogenous.begin(), exogenous.end(), it->var_num) == exogenous.end()))
        {
          ostringstream tmp;
          tmp << "\n the conditional forecast involving as constrained variable " << get_variable(eEndogenous, it->exo_num) << " and as endogenized exogenous " << get_variable(eExogenous, it->var_num) << " that do not appear in block=" << Block_Count+1 << ")\n You should not use block in model options\n";
          throw FatalExceptionHandling(tmp.str());
        }
      else if ((find(endogenous.begin(), endogenous.end(), it->exo_num) != endogenous.end()) &&  (find(exogenous.begin(), exogenous.end(), it->var_num) != exogenous.end()) && ((fb->get_type() == EVALUATE_FORWARD) || (fb->get_type() != EVALUATE_BACKWARD)))
        {
          ostringstream tmp;
          tmp << "\n the conditional forecast cannot be implemented for the block=" << Block_Count+1 << ") that has to be evaluated instead to be solved\n You should not use block in model options\n";
          throw FatalExceptionHandling(tmp.str());
        }
      else if (find(previous_block_exogenous.begin(), previous_block_exogenous.end(), it->var_num) != previous_block_exogenous.end())
        {
          ostringstream tmp;
          tmp << "\n the conditional forecast involves in the block " << Block_Count+1 << " the endogenized exogenous " << get_variable(eExogenous, it->var_num) << " that appear also in a previous block\n You should not use block in model options\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  for (vector<unsigned int>::iterator it = exogenous.begin(); it != exogenous.end(); it++)
    previous_block_exogenous.push_back(*it);
}



bool
Interpreter::MainLoop(string bin_basename, CodeLoad code, bool evaluate, int block, bool last_call, bool constrained, vector<s_plan> sconstrained_extended_path, vector_table_conditional_local_type vector_table_conditional_local)
{
  int var;
  Block_Count = -1;
  bool go_on = true;
  double max_res_local = 0;
  int max_res_idx_local = 0;

  if (block < 0)
    {
      if (steady_state)
        residual = vector<double>(y_size);
      else
        residual = vector<double>(y_size*(periods+y_kmin));
    }

  while (go_on)
    {
      switch (it_code->first)
        {
        case FBEGINBLOCK:
          Block_Count++;
#ifdef DEBUG
          mexPrintf("---------------------------------------------------------\n");
          if (block < 0)
            mexPrintf("FBEGINBLOCK Block_Count=%d\n", Block_Count+1);
          else
            mexPrintf("FBEGINBLOCK block=%d\n", block+1);
#endif
          //it's a new block
          {
            FBEGINBLOCK_ *fb = (FBEGINBLOCK_ *) it_code->second;
            Block_Contain = fb->get_Block_Contain();
            it_code++;
            if (constrained)
              check_for_controlled_exo_validity(fb,sconstrained_extended_path);
            set_block(fb->get_size(), fb->get_type(), file_name, bin_basename, Block_Count, fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int(), block);
            if (print)
              print_a_block();
            else if (evaluate)
              {
#ifdef DEBUG
                mexPrintf("jacobian_block=mxCreateDoubleMatrix(%d, %d, mxREAL)\n", fb->get_size(), fb->get_nb_col_jacob());
#endif
                jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(), mxREAL));
                if (!steady_state)
                  {
#ifdef DEBUG
                    mexPrintf("allocates jacobian_exo_block( %d, %d, mxREAL)\n", fb->get_size(), fb->get_exo_size());
                    mexPrintf("(0) Allocating Jacobian\n");
#endif

                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_exo_jacob(), mxREAL));
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_det_exo_jacob(), mxREAL));
                    jacobian_other_endo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_other_endo_jacob(), mxREAL));
                  }
                if (block >= 0)
                  {
                    if (steady_state)
                      residual = vector<double>(fb->get_size());
                    else
                      residual = vector<double>(fb->get_size()*(periods+y_kmin));
                  }
                evaluate_a_block(true);
              }
            else
              {
#ifdef DEBUG
                mexPrintf("endo in Block_Count=%d, block=%d, type=%d, steady_state=%d, print_it=%d, Block_Count=%d, fb->get_is_linear()=%d, fb->get_endo_nbr()=%d, fb->get_Max_Lag()=%d, fb->get_Max_Lead()=%d, fb->get_u_count_int()=%d\n",
                          Block_Count, fb->get_size(), fb->get_type(), steady_state, print_it, Block_Count, fb->get_is_linear(), fb->get_endo_nbr(), fb->get_Max_Lag(), fb->get_Max_Lead(), fb->get_u_count_int());
#endif
                bool result;
                if (sconstrained_extended_path.size())
                  {
                    //mexPrintf("(1) Allocating Jacobian fb->get_size()=%d fb->get_nb_col_jacob()=%d\n", fb->get_size(), fb->get_nb_col_jacob());
                    jacobian_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_jacob(), mxREAL));
                    //mexPrintf("mxGetPr(jacobian_block[block_num])=%x\n",mxGetPr(jacobian_block[0]));
                    jacobian_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_exo_jacob(), mxREAL));
                    jacobian_det_exo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_det_exo_jacob(), mxREAL));
                    jacobian_other_endo_block.push_back(mxCreateDoubleMatrix(fb->get_size(), fb->get_nb_col_other_endo_jacob(), mxREAL));
                    residual = vector<double>(fb->get_size()*(periods+y_kmin));
                    result = simulate_a_block(vector_table_conditional_local);

                    mxDestroyArray(jacobian_block.back());
                    jacobian_block.pop_back();
                    mxDestroyArray(jacobian_exo_block.back());
                    jacobian_exo_block.pop_back();
                    mxDestroyArray(jacobian_det_exo_block.back());
                    jacobian_det_exo_block.pop_back();
                    mxDestroyArray(jacobian_other_endo_block.back());
                    jacobian_other_endo_block.pop_back();
                  }
                else
                  result = simulate_a_block(vector_table_conditional_local);
                //mexPrintf("OKe\n");
                if (max_res > max_res_local)
                  {
                    max_res_local = max_res;
                    max_res_idx_local = max_res_idx;
                  }
                if (result == ERROR_ON_EXIT)
                  return ERROR_ON_EXIT;
              }
            if (last_call)
              delete fb;
          }
          if (block >= 0)
            {

              go_on = false;
            }

          break;
        case FEND:
#ifdef DEBUG
          mexPrintf("FEND\n");
#endif
          go_on = false;
          it_code++;
          break;
        case FDIMT:
#ifdef DEBUG
          mexPrintf("FDIMT size=%d\n", ((FDIMT_ *) it_code->second)->get_size());
#endif
          var = ((FDIMT_ *) it_code->second)->get_size();
          if (T)
            mxFree(T);
          T = (double *) mxMalloc(var*(periods+y_kmin+y_kmax)*sizeof(double));
          test_mxMalloc(T, __LINE__, __FILE__, __func__, var*(periods+y_kmin+y_kmax)*sizeof(double));
          if (block >= 0)
            {
              it_code = code_liste.begin() + code.get_begin_block(block);
            }
          else
            it_code++;
          break;
        case FDIMST:
#ifdef DEBUG
          mexPrintf("FDIMST size=%d\n", ((FDIMST_ *) it_code->second)->get_size());
#endif
          var = ((FDIMST_ *) it_code->second)->get_size();
          if (T)
            mxFree(T);
          if (global_temporary_terms)
            {
              if (GlobalTemporaryTerms == NULL)
                {
                  mexPrintf("GlobalTemporaryTerms is NULL\n"); mexEvalString("drawnow;");
                }
              if (var != (int) mxGetNumberOfElements(GlobalTemporaryTerms))
                GlobalTemporaryTerms = mxCreateDoubleMatrix(var, 1, mxREAL);
              T = mxGetPr(GlobalTemporaryTerms);
            }
          else
            {
              T = (double *) mxMalloc(var*sizeof(double));
              test_mxMalloc(T, __LINE__, __FILE__, __func__, var*sizeof(double));
            }


          if (block >= 0)
            it_code = code_liste.begin() + code.get_begin_block(block);
          else
            it_code++;
          break;
        default:
          ostringstream tmp;
          tmp << " in compute_blocks, unknown command " << it_code->first << " (block=" << Block_Count << ")\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  max_res = max_res_local ;
  max_res_idx = max_res_idx_local;
  Close_SaveCode();
  return true;
}

string
Interpreter::elastic(string str, unsigned int len, bool left)
{
  if (str.length() > len)
    return str;
  else
    {
      int diff = len - str.length();
      if (diff % 2 == 0)
        {
          if (left)
            {
              //mexPrintf("(1) diff=%d\n",diff);
              str.insert(str.end(),diff-1,' ');
              str.insert(str.begin(),1,' ');
            }
          else
            {
              str.insert(str.end(),diff/2,' ');
              str.insert(str.begin(),diff/2,' ');
            }
        }
      else
        {
          if (left)
            {
              //mexPrintf("(2) diff=%d\n",diff);
              str.insert(str.end(),diff-1,' ');
              str.insert(str.begin(),1,' ');
            }
          else
            {
              str.insert(str.end(),ceil(diff/2),' ');
              str.insert(str.begin(),ceil(diff/2+1),' ');
            }
        }
      return str;
    }
}

bool
Interpreter::extended_path(string file_name, string bin_basename, bool evaluate, int block, int &nb_blocks, int nb_periods, vector<s_plan> sextended_path, vector<s_plan> sconstrained_extended_path, vector<string> dates, table_conditional_global_type table_conditional_global)
{
  CodeLoad code;
  ReadCodeFile(file_name, code);
  it_code = code_liste.begin();
  it_code_type Init_Code = code_liste.begin();
  /*size_t size_of_direction = y_size*(periods + y_kmax + y_kmin)*sizeof(double);
  double *y_save = (double *) mxMalloc(size_of_direction);
  double *x_save = (double *) mxMalloc((periods + y_kmax + y_kmin) * col_x *sizeof(double));*/
  int max_periods = max(periods, nb_periods);
  size_t size_of_direction = y_size*(/*max_periods*/ nb_periods + y_kmax + y_kmin)*sizeof(double);
  double *y_save = (double *) mxMalloc(size_of_direction);
  test_mxMalloc(y_save, __LINE__, __FILE__, __func__, size_of_direction);


  double *x_save = (double *) mxMalloc((/*max_periods*/nb_periods + y_kmax + y_kmin) * col_x *sizeof(double));
  test_mxMalloc(x_save, __LINE__, __FILE__, __func__, (/*max_periods*/nb_periods + y_kmax + y_kmin) * col_x *sizeof(double));

  vector_table_conditional_local_type vector_table_conditional_local;
  vector_table_conditional_local.clear();

  int endo_name_length_l = endo_name_length;
  for (int j = 0; j < col_x* nb_row_x; j++)
    {
      x_save[j] = x[j];
      x[j] = 0;
    }
  for (int j = 0; j < col_x; j++)
    for (int i = nb_row_x; i < nb_periods; i++)
      x_save[j * (nb_periods + y_kmin + y_kmax) + i] = 0;//x[(j * nb_row_x) + 1];
  for (int j = 0; j < col_x; j++)
    x[y_kmin + j * nb_row_x] = x_save[1 + y_kmin + j * nb_row_x];
  //for (int i = 0; i < (y_size*(periods + y_kmax + y_kmin)); i++)
  for (int j = 0; j < y_size; j++)
    for (int i = 0; i < nb_periods; i++)
      if (i < periods)
        y_save[j + (i + y_kmin) * y_size] = y[ j +  (i + y_kmin) * y_size];
      else
        y_save[j + (i + y_kmin) * y_size] = y[ j +  (periods + y_kmin) * y_size];

  if (endo_name_length_l < 8)
    endo_name_length_l = 8;
  bool old_print_it = print_it;
  print_it = false;
  ostringstream res1;
  res1 << std::scientific << 2.54656875434865131;
  int real_max_length = res1.str().length();
  int date_length = dates[0].length();
  int table_length = 2 + date_length + 3 + endo_name_length_l + 3 + real_max_length + 3 + 3 + 2 + 6 + 2;
  string line;
  line.insert(line.begin(),table_length,'-');
  line.insert(line.length(),"\n");
  if (old_print_it)
    {
       mexPrintf("\nExtended Path simulation:\n");
       mexPrintf("-------------------------\n");
       mexPrintf(line.c_str());
       string title = "|" + elastic("date",date_length+2, false) + "|" + elastic("variable",endo_name_length_l+2, false) + "|" + elastic("max. value",real_max_length+2, false) + "| iter. |" + elastic("cvg",5, false) + "|\n";
       mexPrintf(title.c_str());
       mexPrintf(line.c_str());
    }
  for (int t = 0; t < nb_periods; t++)
    {
      nb_blocks = 0;
      previous_block_exogenous.clear();
      if (old_print_it)
        {
          mexPrintf("|%s|",elastic(dates[t], date_length+2, false).c_str());
          mexEvalString("drawnow;");
        }
      for (vector<s_plan>::const_iterator it = sextended_path.begin(); it != sextended_path.end(); it++)
        {
          x[y_kmin + (it->exo_num - 1) * /*(nb_periods + y_kmax + y_kmin)*/ nb_row_x] = it->value[t];
        }

      it_code = Init_Code;
      vector_table_conditional_local.clear();
      if (table_conditional_global.size())
        vector_table_conditional_local = table_conditional_global[t];
      if (t < nb_periods)
        MainLoop(bin_basename, code, evaluate, block, false, true, sconstrained_extended_path, vector_table_conditional_local);
      else
        MainLoop(bin_basename, code, evaluate, block, true, true, sconstrained_extended_path, vector_table_conditional_local);
      for (int j = 0; j < y_size; j++)
        {
          y_save[j + (t + y_kmin) * y_size] = y[ j +  (y_kmin) * y_size];
          if (y_kmin > 0)
            y[j ] = y[ j +  (y_kmin) * y_size];
        }
      for (int j = 0; j < col_x; j++)
        {
          x_save[t + y_kmin + j * nb_row_x] = x[y_kmin + j * nb_row_x];
          if (t < nb_periods)
            x[y_kmin + j * nb_row_x] = x_save[t + 1 + y_kmin + j * nb_row_x];
        }

      if (old_print_it)
        {
          ostringstream res, res1;
          for (unsigned int i = 0; i < endo_name_length; i++)
            if (P_endo_names[CHAR_LENGTH*(max_res_idx+i*y_size)] != ' ')
              res << P_endo_names[CHAR_LENGTH*(max_res_idx+i*y_size)];
          res1 << std::scientific << max_res;
          mexPrintf("%s|%s| %4d  |  x  |\n",elastic(res.str(),endo_name_length_l+2, true).c_str(), elastic(res1.str(), real_max_length+2, false).c_str(), iter);
          mexPrintf(line.c_str());
          mexEvalString("drawnow;");
        }
    }
  print_it = old_print_it;
  /*for (int j = 0; j < y_size; j++)
    {
      for(int k = nb_periods; k < periods; k++)
        y_save[j + (k + y_kmin) * y_size] = y[ j +  ( k - (nb_periods-1) + y_kmin) * y_size];
    }*/
  for (int i = 0; i < (y_size*(nb_periods + y_kmax + y_kmin)); i++)
    y[i]  = y_save[i];
  for (int j = 0; j < col_x* nb_row_x; j++)
    x[j] = x_save[j];
  if (Init_Code->second)
    mxFree(Init_Code->second);
  if (y_save)
    mxFree(y_save);
  if (x_save)
    mxFree(x_save);
  nb_blocks = Block_Count+1;
  if (T && !global_temporary_terms)
    mxFree(T);
  return true;
}

bool
Interpreter::compute_blocks(string file_name, string bin_basename, bool evaluate, int block, int &nb_blocks)
{
  CodeLoad code;
  ReadCodeFile(file_name, code);

  //The big loop on intructions
  it_code = code_liste.begin();
  it_code_type Init_Code = it_code;
  vector<s_plan> s_plan_junk;
  vector_table_conditional_local_type vector_table_conditional_local_junk;

  MainLoop(bin_basename, code, evaluate, block, true, false, s_plan_junk, vector_table_conditional_local_junk);


  mxFree(Init_Code->second);
  nb_blocks = Block_Count+1;
  if (T && !global_temporary_terms)
    mxFree(T);
  return true;
}
