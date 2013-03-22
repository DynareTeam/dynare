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
#include <cstring>
#include "Interpreter.hh"
#include "ErrorHandling.hh"
#include <ctime>
#include <math.h>

#ifdef DEBUG_EX

using namespace std;
# include <sstream>




string
Get_Argument(const char *argv)
{
  string f(argv);
  return f;
}

#else

void (*prev_fn)(int);


string
Get_Argument(const mxArray *prhs)
{
  const mxArray *mxa = prhs;
  mwSize buflen = mwSize(mxGetM(mxa) * mxGetN(mxa) + 1);
  char *first_argument;
  first_argument = (char *) mxCalloc(buflen, sizeof(char));
  size_t status = mxGetString(mxa, first_argument, buflen);
  if (status != 0)
    mexWarnMsgTxt("Not enough space. The first argument is truncated.");
  string f(first_argument);
  mxFree(first_argument);
  return f;
}
#endif


//#include <windows.h> 
#include <stdio.h> 
 

#ifdef CUDA
int
GPU_Test_and_Info(cublasHandle_t *cublas_handle, cusparseHandle_t *cusparse_handle, cusparseMatDescr_t *descr)
{
  cudaDeviceProp deviceProp;
  int device_count, device, version, version_max = 0;
  cublasStatus_t cublas_status;
  cudaError_t cuda_error;
  *descr=0;

  /* ask cuda how many devices it can find */
  cudaGetDeviceCount(&device_count);
  if (device_count < 1)
    {
      /* if it couldn't find any fail out */
      ostringstream tmp;
      tmp << " Unable to find a CUDA device. Unable to implement CUDA solvers\n";
      throw FatalExceptionHandling(tmp.str());
    }
  else
    {
      mexPrintf("-----------------------------------------\n");
      for (int i = 0; i < device_count; i++)
        {
          cudaSetDevice(i);
          // Statistics about the GPU device
          cuda_error = cudaGetDeviceProperties(&deviceProp, i);
          if (cuda_error != cudaSuccess)
            {
               ostringstream tmp;
               tmp << "  bytecode cudaGetDeviceProperties failed\n";
               throw FatalExceptionHandling(tmp.str());
            }
          mexPrintf("> GPU device %d: \"%s\" has:\n   - %d Multi-Processors,\n   - %d threads per multiprocessor,\n", i, deviceProp.name, deviceProp.multiProcessorCount, deviceProp.maxThreadsPerMultiProcessor);
          mexEvalString("drawnow;");
          version = (deviceProp.major * 0x10 + deviceProp.minor);
          if (version >= version_max)
            {
              device = i;
              version_max = version;
            }
          mexPrintf("   - %4.2fMhz clock rate,\n   - %2.0fMb of memory,\n   - %d.%d compute capabilities.\n", double(deviceProp.clockRate) / (1024 * 1024), double(deviceProp.totalGlobalMem) / (1024 * 1024), deviceProp.major, deviceProp.minor);
          mexEvalString("drawnow;");
        }
    }
  mexPrintf("> Device %d selected\n", device);
  mexEvalString("drawnow;");

  cuda_error = cudaSetDevice(device);
  if (cuda_error != cudaSuccess)
    {
       ostringstream tmp;
       tmp << "  bytecode cudaSetDevice failed\n";
       throw FatalExceptionHandling(tmp.str());
    }

  if(version_max < 0x11)
    {
       ostringstream tmp;
       tmp << "  bytecode requires a minimum CUDA compute 1.1 capability\n";
       cudaDeviceReset();
       throw FatalExceptionHandling(tmp.str());
    }

  // Initialize CuBlas library
  cublas_status = cublasCreate(cublas_handle);
  if (cublas_status != CUBLAS_STATUS_SUCCESS)
    {
      ostringstream tmp;
      switch(cublas_status)
        {
          case CUBLAS_STATUS_NOT_INITIALIZED:
            tmp << " the CUBLAS initialization failed.\n";
            break;
          case CUBLAS_STATUS_ALLOC_FAILED:
            tmp << " the resources could not be allocated.\n";
            break;
          default:
            tmp << " unknown error during the initialization of cusparse library.\n";
        }
      throw FatalExceptionHandling(tmp.str());
    }

  // Initialize the CuSparse library
  cusparseStatus_t cusparse_status;
  cusparse_status = cusparseCreate(cusparse_handle);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS)
    {
      ostringstream tmp;
      switch(cusparse_status)
        {
          case CUSPARSE_STATUS_NOT_INITIALIZED:
            tmp << " the CUDA Runtime initialization failed.\n";
            break;
          case CUSPARSE_STATUS_ALLOC_FAILED:
            tmp <<  " the resources could not be allocated.\n";
            break;
          case CUSPARSE_STATUS_ARCH_MISMATCH:
            tmp <<  " the device compute capability (CC) is less than 1.1. The CC of at least 1.1 is required.\n";
            break;
          default:
            tmp << " unknown error during the initialization of cusparse library.\n";
        }
      throw FatalExceptionHandling(tmp.str());
    }

  // Create and setup matrix descriptor
  cusparse_status = cusparseCreateMatDescr(descr);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS)
    {
      ostringstream tmp;
      tmp << " Matrix descriptor initialization failed\n";
      throw FatalExceptionHandling(tmp.str());
    }
  cusparseSetMatType(*descr, CUSPARSE_MATRIX_TYPE_GENERAL);
  cusparseSetMatIndexBase(*descr, CUSPARSE_INDEX_BASE_ZERO);

  mexPrintf("> Driver version:\n");
  int cuda_version;
  cuda_error = cudaDriverGetVersion(&cuda_version);
  if (cuda_error  != cudaSuccess)
    {
      ostringstream tmp;
      tmp << " cudaGetVersion has failed\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mexPrintf("   - CUDA version %5.3f\n", double(cuda_version) / 1000);
  int cublas_version;
  cublas_status = cublasGetVersion(*cublas_handle, &cublas_version);
  if (cublas_status != CUBLAS_STATUS_SUCCESS)
    {
      ostringstream tmp;
      tmp << " cublasGetVersion has failed\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mexPrintf("   - CUBLAS version %5.3f\n", double(cublas_version) / 1000);
  int cusparse_version;
  cusparse_status = cusparseGetVersion(*cusparse_handle, &cusparse_version);
  if (cusparse_status != CUSPARSE_STATUS_SUCCESS)
    {
      ostringstream tmp;
      tmp << " cusparseGetVersion has failed\n";
      throw FatalExceptionHandling(tmp.str());
    }
  mexPrintf("   - CUSPARSE version %5.3f\n", double(cusparse_version) / 1000);
  mexPrintf("-----------------------------------------\n");
  return device;
}

void
GPU_close(cublasHandle_t cublas_handle, cusparseHandle_t cusparse_handle, cusparseMatDescr_t descr)
{
  cublasChk(cublasDestroy(cublas_handle),"in bytecode cublasDestroy failed\n");
  cusparseChk(cusparseDestroyMatDescr(descr), "in bytecode cusparseDestroyMatDescr failed\n");
  cusparseChk(cusparseDestroy(cusparse_handle),"in bytecode cusparseDestroy failed\n");
}

#endif
string
deblank(string x)
{
  for(int i = 0; i < x.length(); i++)
    if (x[i] == ' ')
      x.erase(i--, 1);
  return x;
}

void
Get_Arguments_and_global_variables(int nrhs,
#ifndef DEBUG_EX
                                   const mxArray *prhs[],
#else
                                   const char *prhs[],
#endif
                                   int &count_array_argument,
                                   double *yd[], size_t &row_y, size_t &col_y,
                                   double *xd[], size_t &row_x, size_t &col_x,
                                   double *params[],
                                   double *steady_yd[], size_t &steady_row_y, size_t &steady_col_y,
                                   unsigned int &periods,
#ifndef DEBUG_EX
                                   mxArray *block_structur[],
#endif
                                   bool &steady_state, bool &evaluate, int &block,
                                   mxArray *M_[], mxArray *oo_[], mxArray *options_[], bool &global_temporary_terms,
                                   bool &print,
                                   bool &print_error,
                                   mxArray *GlobalTemporaryTerms[],
                                   string *plan_struct_name, string *pfplan_struct_name)
{
  size_t pos;
#ifdef DEBUG_EX
  for (int i = 2; i < nrhs; i++)
#else
    for (int i = 0; i < nrhs; i++)
#endif
      {
#ifndef DEBUG_EX
        if (!mxIsChar(prhs[i]))
          {
            switch (count_array_argument)
              {
              case 0:
                *yd = mxGetPr(prhs[i]);
                row_y = mxGetM(prhs[i]);
                col_y = mxGetN(prhs[i]);
                break;
              case 1:
                *xd =  mxGetPr(prhs[i]);
                row_x = mxGetM(prhs[i]);
                col_x = mxGetN(prhs[i]);
                break;
              case 2:
                *params = mxGetPr(prhs[i]);
                break;
              case 3:
                *steady_yd = mxGetPr(prhs[i]);
                steady_row_y = mxGetM(prhs[i]);
                steady_col_y = mxGetN(prhs[i]);
                break;
              case 4:
                periods = int(mxGetScalar(prhs[i]));
                break;
              case 5:
                *block_structur = mxDuplicateArray(prhs[i]);
                break;
              case 6:
                global_temporary_terms = true;
                *GlobalTemporaryTerms = mxDuplicateArray(prhs[i]);
                break;
              default:
                mexPrintf("Unknown argument count_array_argument=%d\n",count_array_argument);
                break;
              }
            count_array_argument++;
          }
        else
#endif
          if (Get_Argument(prhs[i]) == "static")
            steady_state = true;
          else if (Get_Argument(prhs[i]) == "dynamic")
            steady_state = false;
          else if (Get_Argument(prhs[i]) == "evaluate")
            evaluate = true;
          else if (Get_Argument(prhs[i]) == "global_temporary_terms")
            global_temporary_terms = true;
          else if (Get_Argument(prhs[i]) == "print")
            print = true;
          else if (Get_Argument(prhs[i]) == "no_print_error")
            print_error = false;
          else
            {
              ;
              if ((pos = Get_Argument(prhs[i]).find("block")) != (int) string::npos)
                {
                  size_t pos1 = Get_Argument(prhs[i]).find("=", pos+5);
                  if (pos1 != (int) string::npos)
                    pos = pos1 + 1;
                  else
                    pos += 5;
                  block =  atoi(Get_Argument(prhs[i]).substr(pos, string::npos).c_str())-1;
                }
              else if ((pos = Get_Argument(prhs[i]).find("pfplan")) != (int) string::npos)
                {
                  size_t pos1 = Get_Argument(prhs[i]).find("=", pos+6);
                  if (pos1 != (int) string::npos)
                    pos = pos1 + 1;
                  else
                    pos += 6;
                  *pfplan_struct_name =  deblank(Get_Argument(prhs[i]).substr(pos, string::npos));
                }
              else if ((pos = Get_Argument(prhs[i]).find("plan")) != (int) string::npos)
                {
                  size_t pos1 = Get_Argument(prhs[i]).find("=", pos+4);
                  if (pos1 != (int) string::npos)
                    pos = pos1 + 1;
                  else
                    pos += 4;
                  *plan_struct_name =  deblank(Get_Argument(prhs[i]).substr(pos, string::npos));
                }
              else
                {
                  ostringstream tmp;
                  tmp << " in main, unknown argument : " << Get_Argument(prhs[i]) << "\n";
                  throw FatalExceptionHandling(tmp.str());
                }
            }
      }
  if (count_array_argument > 0 && count_array_argument < 5)
    {
      if (count_array_argument == 3 && steady_state)
        periods = 1;
      else
        {
          ostringstream tmp;
          tmp << " in main, missing arguments. All the following arguments have to be indicated y, x, params, it_, ys\n";
          throw FatalExceptionHandling(tmp.str());
        }
    }
  *M_ = mexGetVariable("global", "M_");
  if (M_ == NULL)
    {
      ostringstream tmp;
      tmp << " in main, global variable not found: M_\n";
      throw FatalExceptionHandling(tmp.str());
    }
  /* Gets variables and parameters from global workspace of Matlab */
  *oo_ = mexGetVariable("global", "oo_");
  if (oo_ == NULL)
    {
      ostringstream tmp;
      tmp << " in main, global variable not found: oo_\n";
      throw FatalExceptionHandling(tmp.str());
    }
  *options_ = mexGetVariable("global", "options_");
  if (options_ == NULL)
    {
      ostringstream tmp;
      tmp << " in main, global variable not found: options_\n";
      throw FatalExceptionHandling(tmp.str());
    }
}


#ifdef DEBUG_EX
int
main(int nrhs, const char *prhs[])
#else
/* The gateway routine */
  void
  mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
#endif
{
  mxArray *M_, *oo_, *options_;
  mxArray *GlobalTemporaryTerms;
#ifndef DEBUG_EX
  mxArray *block_structur = NULL;
#else
  int nlhs = 0;
  char *plhs[1];
  load_global((char *) prhs[1]);
#endif
  mxArray *plan_struct = NULL, *pfplan_struct = NULL;
  size_t i, row_y = 0, col_y = 0, row_x = 0, col_x = 0, nb_row_xd = 0;
  size_t steady_row_y, steady_col_y;
  int y_kmin = 0, y_kmax = 0, y_decal = 0;
  unsigned int periods = 1;
  double *direction;
  bool steady_state = false;
  bool evaluate = false;
  int block = -1;
  double *params = NULL;
  double *yd = NULL, *xd = NULL;
  int count_array_argument = 0;
  bool global_temporary_terms = false;
  bool print = false, print_error = true, print_it = false;
  double *steady_yd = NULL, *steady_xd = NULL;
  string plan, pfplan;

  vector<s_plan> splan, spfplan;

#ifdef CUDA
  int CUDA_device = -1;
  cublasHandle_t cublas_handle;
  cusparseHandle_t cusparse_handle;
  cusparseMatDescr_t descr;
#endif
  try
    {
      Get_Arguments_and_global_variables(nrhs, prhs, count_array_argument,
                                         &yd, row_y, col_y,
                                         &xd, row_x, col_x,
                                         &params,
                                         &steady_yd, steady_row_y, steady_col_y,
                                         periods,
#ifndef DEBUG_EX
                                         &block_structur,
#endif
                                         steady_state, evaluate, block,
                                         &M_, &oo_, &options_, global_temporary_terms,
                                         print, print_error, &GlobalTemporaryTerms,
                                         &plan, &pfplan);
    }
  catch (GeneralExceptionHandling &feh)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(feh.GetErrorMsg().c_str());
    }
  if (!count_array_argument)
    {
      int field = mxGetFieldNumber(M_, "params");
      if (field < 0)
        DYN_MEX_FUNC_ERR_MSG_TXT("params is not a field of M_");
      params = mxGetPr(mxGetFieldByNumber(M_, 0, field));
    }

  ErrorMsg emsg;
  


  if (plan.length()>0)
    {
      mxArray* plan_struct = mexGetVariable("base", plan.c_str());
      if (plan_struct == NULL)
        {
          string tmp = plan;
          tmp.insert(0,"Can't find the plan: ");
          DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
        }
      size_t n_plan = mxGetN(plan_struct);
      splan.resize(n_plan);
      for (int i = 0; i < n_plan; i++)
        {
          splan[i].var = "";
          splan[i].exo = "";
          mxArray* tmp = mxGetField(plan_struct, i, "exo");
          if (tmp)
            {
              char name [100];
              mxGetString(tmp, name, 100);
              splan[i].var = name;
              SymbolType variable_type;
              int exo_num = emsg.get_ID(name, &variable_type);
              if (variable_type == eExogenous || variable_type == eExogenousDet)
                splan[i].var_num = exo_num;
              else
                {
                  string tmp = name;
                  tmp.insert(0,"the variable '");
                  tmp.append("'  defined as var in plan is not an exogenous or a deterministic exogenous\n");
                  DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
                }
            }
          tmp = mxGetField(plan_struct, i, "var");
          if (tmp)
            {
              char name [100];
              mxGetString(tmp, name, 100);
              splan[i].exo = name;
              SymbolType variable_type;
              int exo_num = emsg.get_ID(name, &variable_type);
              if (variable_type == eEndogenous)
                splan[i].exo_num = exo_num;
              else
                {
                  string tmp = name;
                  tmp.insert(0,"the variable '");
                  tmp.append("'  defined as exo in plan is not an endogenous variable\n");
                  DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
                }
            }
          tmp = mxGetField(plan_struct, i, "per_value");
          if (tmp)
            {
              size_t num_shocks = mxGetM(tmp);
              (splan[i]).per_value.resize(num_shocks);
              double * per_value = mxGetPr(tmp);
              for (int j = 0; j < num_shocks; j++)
                (splan[i]).per_value[j] = make_pair(ceil(per_value[j]), per_value[j + num_shocks]);
            }
        }
      int i;
      for (vector<s_plan>::iterator it = splan.begin(); it != splan.end(); it++)
        {
          mexPrintf("----------------------------------------------------------------------------------------------------\n");
          mexPrintf("suprise n°%d\n", i+1);
          if (it->exo.length())
            mexPrintf(" plan fliping var=%s (%d) exo=%s (%d) for the following periods and with the following values:\n", it->var.c_str(), it->var_num, it->exo.c_str(), it->exo_num);
          else
            mexPrintf(" plan shocks on var=%s for the following periods and with the following values:\n", it->var.c_str());
          for (vector<pair<int, double> >::iterator it1 = it->per_value.begin(); it1 != it->per_value.end(); it1++)
            {
              mexPrintf("  %3d %10.5f\n",it1->first, it1->second);
            }
          i++;
        }
    }

  if (pfplan.length()>0)
    {
      pfplan_struct = mexGetVariable("base", pfplan.c_str());
      if (!pfplan_struct)
        {
          string tmp = pfplan;
          tmp.insert(0,"Can't find the pfplan: ");
          DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
        }
      size_t n_plan = mxGetN(pfplan_struct);
      spfplan.resize(n_plan);
      for (int i = 0; i < n_plan; i++)
        {
          spfplan[i].var = "";
          spfplan[i].exo = "";
          mxArray* tmp = mxGetField(pfplan_struct, i, "var");
          if (tmp)
            {
              char name [100];
              mxGetString(tmp, name, 100);
              spfplan[i].var = name;
              SymbolType variable_type;
              int exo_num = emsg.get_ID(name, &variable_type);
              if (variable_type == eExogenous || variable_type == eExogenousDet)
                splan[i].var_num = exo_num;
              else
                {
                  string tmp = name;
                  tmp.insert(0,"the variable '");
                  tmp.append("' defined as var in pfplan is not an exogenous or a deterministic exogenous\n");
                  DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
                }
            }
          tmp = mxGetField(pfplan_struct, i, "exo");
          if (tmp)
            {
              char name [100];
              mxGetString(tmp, name, 100);
              spfplan[i].exo = name;
              SymbolType variable_type;
              int exo_num = emsg.get_ID(name, &variable_type);
              if (variable_type == eEndogenous)
                spfplan[i].exo_num = exo_num;
              else
                {
                  string tmp = name;
                  tmp.insert(0,"the variable '");
                  tmp.append("' defined as exo in pfplan  is not an endogenous variable\n");
                  DYN_MEX_FUNC_ERR_MSG_TXT(tmp.c_str());
                }
            }
          tmp = mxGetField(pfplan_struct, i, "per_value");
          if (tmp)
            {
              size_t num_shocks = mxGetM(tmp);
              double * per_value = mxGetPr(tmp);
              (spfplan[i]).per_value.resize(num_shocks);
              for (int j = 0; j < num_shocks; j++)
                spfplan[i].per_value[j] = make_pair(ceil(per_value[j]), per_value[j+ num_shocks]);
            }
        }
      int i;
      for (vector<s_plan>::iterator it = spfplan.begin(); it != spfplan.end(); it++)
        {
          mexPrintf("----------------------------------------------------------------------------------------------------\n");
          mexPrintf("perfect foresight n°%d\n", i+1);
          if (it->exo.length())
            mexPrintf(" plan flipping var=%s (%d) exo=%s (%d) for the following periods and with the following values:\n", it->var.c_str(), it->var_num, it->exo.c_str(), it->exo_num);
          else
            mexPrintf(" plan shocks on var=%s (%d) for the following periods and with the following values:\n", it->var.c_str(), it->var_num);
          for (vector<pair<int, double> >::iterator it1 = it->per_value.begin(); it1 != it->per_value.end(); it1++)
            {
              mexPrintf("  %3d %10.5f\n",it1->first, it1->second);
            }
          i++;
        }
    }



  int field_steady_state = mxGetFieldNumber(oo_, "steady_state");
  if (field_steady_state < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("steady_state is not a field of oo_");
  int field_exo_steady_state = mxGetFieldNumber(oo_, "exo_steady_state");
  if (field_exo_steady_state < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("exo_steady_state is not a field of oo_");

  if (!steady_state)
    {
      int field_endo_simul = mxGetFieldNumber(oo_, "endo_simul");
      if (field_endo_simul < 0)
        DYN_MEX_FUNC_ERR_MSG_TXT("endo_simul is not a field of oo_");

      int field_exo_simul = mxGetFieldNumber(oo_, "exo_simul");
      if (field_exo_simul < 0)
        DYN_MEX_FUNC_ERR_MSG_TXT("exo_simul is not a field of oo_");

      if (!count_array_argument)
        {
          mxArray* endo_sim_arr = mxGetFieldByNumber(oo_, 0, field_endo_simul);
          yd = mxGetPr(endo_sim_arr);
          row_y = mxGetM(endo_sim_arr);
          col_y = mxGetN(endo_sim_arr);
          mxArray* exo_sim_arr = mxGetFieldByNumber(oo_, 0, field_exo_simul);
          xd = mxGetPr(exo_sim_arr);
          row_x = mxGetM(exo_sim_arr);
          col_x = mxGetN(exo_sim_arr);
          nb_row_xd = row_x;
        }
      int field = mxGetFieldNumber(M_, "maximum_lag");
      if (field >= 0)
        y_kmin = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field)))));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("maximum_lag is not a field of M_");
      field = mxGetFieldNumber(M_, "maximum_lead");
      if (field >= 0)
        y_kmax = int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field)))));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("maximum_lead is not a field of M_");
      field = mxGetFieldNumber(M_, "maximum_endo_lag");
      if (field >= 0)
        y_decal = max(0, y_kmin-int (floor(*(mxGetPr(mxGetFieldByNumber(M_, 0, field))))));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("maximum_endo_lag is not a field of M_");

      if (!count_array_argument)
        {
          int field = mxGetFieldNumber(options_, "periods");
          if (field >= 0)
            periods = int (floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, field)))));
          else
            DYN_MEX_FUNC_ERR_MSG_TXT("options_ is not a field of options_");
        }

      if (!steady_yd )
        {
          mxArray* steady_state_arr = mxGetFieldByNumber(oo_, 0, field_steady_state);
          steady_yd = mxGetPr(steady_state_arr);
          steady_row_y = mxGetM(steady_state_arr);
          steady_col_y = mxGetN(steady_state_arr);
        }
      steady_xd = mxGetPr(mxGetFieldByNumber(oo_, 0, field_exo_steady_state));
    }
  else
    {
      if (!count_array_argument)
        {
          mxArray* steady_state_arr = mxGetFieldByNumber(oo_, 0, field_steady_state);
          yd = mxGetPr(steady_state_arr);
          row_y = mxGetM(steady_state_arr);
          col_y = mxGetN(steady_state_arr);

          mxArray* exo_steady_state_arr = mxGetFieldByNumber(oo_, 0, field_exo_steady_state);
          xd = mxGetPr(exo_steady_state_arr);
          row_x = mxGetM(exo_steady_state_arr);
          col_x = mxGetN(exo_steady_state_arr);
          nb_row_xd = row_x;
        }
    }
  int field = mxGetFieldNumber(options_, "verbosity");
  int verbose = 0;
  if (field >= 0)
    verbose = int(*mxGetPr((mxGetFieldByNumber(options_, 0, field))));
  else
    DYN_MEX_FUNC_ERR_MSG_TXT("verbosity is not a field of options_");
  if (verbose)
    print_it = true;
  field = mxGetFieldNumber(options_, "maxit_");
  if (field < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("maxit_ is not a field of options_");
  int maxit_ = int (floor(*(mxGetPr(mxGetFieldByNumber(options_, 0, field)))));
  field = mxGetFieldNumber(options_, "slowc");
  if (field < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("slows is not a field of options_");
  double slowc = double (*(mxGetPr(mxGetFieldByNumber(options_, 0, field))));
  field = mxGetFieldNumber(options_, "markowitz");
  if (field < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("markowitz is not a field of options_");
  double markowitz_c = double (*(mxGetPr(mxGetFieldByNumber(options_, 0, field))));
  field = mxGetFieldNumber(options_, "minimal_solving_periods");
  if (field < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("minimal_solving_periods is not a field of options_");
  int minimal_solving_periods = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, field))));
  field = mxGetFieldNumber(options_, "stack_solve_algo");
  if (field < 0)
    DYN_MEX_FUNC_ERR_MSG_TXT("stack_solve_algo is not a field of options_");
  int stack_solve_algo = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, field))));
  int solve_algo;
  double solve_tolf;

  if (steady_state)
    {
      int field = mxGetFieldNumber(options_, "solve_algo");
      if (field >= 0)
        solve_algo = int (*(mxGetPr(mxGetFieldByNumber(options_, 0, field))));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("solve_algo is not a field of options_");
      field = mxGetFieldNumber(options_, "solve_tolf");
      if (field >= 0)
        solve_tolf = *(mxGetPr(mxGetFieldByNumber(options_, 0, field)));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("solve_tolf is not a field of options_");
    }
  else
    {
      solve_algo = stack_solve_algo;
      int field = mxGetFieldNumber(options_, "dynatol");
      mxArray *dynatol;
      if (field >= 0)
        dynatol = mxGetFieldByNumber(options_, 0, field);
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("dynatol is not a field of options_");
      field = mxGetFieldNumber(dynatol, "f");
      if (field >= 0)
        solve_tolf= *mxGetPr((mxGetFieldByNumber(dynatol, 0, field)));
      else
        DYN_MEX_FUNC_ERR_MSG_TXT("f is not a field of options_.dynatol");
    }
  field = mxGetFieldNumber(M_, "fname");
  mxArray *mxa;
  if (field >= 0)
    mxa = mxGetFieldByNumber(M_, 0, field);
  else
    DYN_MEX_FUNC_ERR_MSG_TXT("fname is not a field of M_");
  size_t buflen = mxGetM(mxa) * mxGetN(mxa) + 1;
  char *fname;
  fname = (char *) mxCalloc(buflen+1, sizeof(char));
  size_t status = mxGetString(mxa, fname, int(buflen));
  fname[buflen] = ' ';
  if (status != 0)
    mexWarnMsgTxt("Not enough space. Filename is truncated.");
  string file_name = fname;

#ifdef CUDA
  try
    {
      if (stack_solve_algo == 7 && !steady_state)
        CUDA_device = GPU_Test_and_Info(&cublas_handle, &cusparse_handle, &descr);
    }
  catch (GeneralExceptionHandling &feh)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(feh.GetErrorMsg().c_str());
    }
#else
  if (stack_solve_algo == 7 && !steady_state)
    DYN_MEX_FUNC_ERR_MSG_TXT("bytecode has not been compiled with CUDA option. Bytecode Can't use options_.stack_solve_algo=7\n");
#endif

  size_t size_of_direction = col_y*row_y*sizeof(double);
  double *y = (double *) mxMalloc(size_of_direction);
  double *ya = (double *) mxMalloc(size_of_direction);
  direction = (double *) mxMalloc(size_of_direction);
  memset(direction, 0, size_of_direction);
  double *x = (double *) mxMalloc(col_x*row_x*sizeof(double));
  #ifdef USE_OMP
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  #endif
  for (i = 0; i < row_x*col_x; i++)
    {
      x[i] = double (xd[i]);
    }

  #ifdef USE_OMP
  #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
  #endif
  for (i = 0; i < row_y*col_y; i++)
    {
      y[i]  = double (yd[i]);
      ya[i] = double (yd[i]);
    }
  size_t y_size = row_y;
  size_t nb_row_x = row_x;
  clock_t t0 = clock();
  Interpreter interprete(params, y, ya, x, steady_yd, steady_xd, direction, y_size, nb_row_x, nb_row_xd, periods, y_kmin, y_kmax, maxit_, solve_tolf, size_of_direction, slowc, y_decal,
                         markowitz_c, file_name, minimal_solving_periods, stack_solve_algo, solve_algo, global_temporary_terms, print, print_error, GlobalTemporaryTerms, steady_state,
                         print_it
#ifdef CUDA
                         , CUDA_device, cublas_handle, cusparse_handle, descr
#endif
                         );
  string f(fname);
  mxFree(fname);
  int nb_blocks = 0;
  double *pind;
  bool no_error = true;
  try
    {
      interprete.compute_blocks(f, f, evaluate, block, nb_blocks);
    }
  catch (GeneralExceptionHandling &feh)
    {
      DYN_MEX_FUNC_ERR_MSG_TXT(feh.GetErrorMsg().c_str());
    }

#ifdef CUDA
  if (stack_solve_algo == 7 && !steady_state)
    GPU_close(cublas_handle, cusparse_handle, descr);
#endif

  clock_t t1 = clock();
  if (!steady_state && !evaluate && no_error && print)
    mexPrintf("Simulation Time=%f milliseconds\n", 1000.0*(double (t1)-double (t0))/double (CLOCKS_PER_SEC));
#ifndef DEBUG_EX
  bool dont_store_a_structure = false;
  if (nlhs > 0)
    {
      plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
      pind = mxGetPr(plhs[0]);
      if (no_error)
        pind[0] = 0;
      else
        pind[0] = 1;
      if (nlhs > 1)
        {
          if (block >= 0)
            {
              if (evaluate)
                {
                  vector<double> residual = interprete.get_residual();
                  plhs[1] = mxCreateDoubleMatrix(int(residual.size()/double(col_y)), int(col_y), mxREAL);
                  pind = mxGetPr(plhs[1]);
                  for (i = 0; i < residual.size(); i++)
                    pind[i] = residual[i];
                }
              else
                {
                  plhs[1] = mxCreateDoubleMatrix(int(row_y), int(col_y), mxREAL);
                  pind = mxGetPr(plhs[1]);
                  for (i = 0; i < row_y*col_y; i++)
                    pind[i] = y[i];
                }
            }
          else
            {
              plhs[1] = mxCreateDoubleMatrix(int(row_y), int(col_y), mxREAL);
              pind = mxGetPr(plhs[1]);
              if (evaluate)
                {
                  vector<double> residual = interprete.get_residual();
                  for (i = 0; i < residual.size(); i++)
                    pind[i] = residual[i];
                }
              else
                for (i = 0; i < row_y*col_y; i++)
                  pind[i] = y[i];
            }
          if (nlhs > 2)
            {
              if (evaluate)
                {
                  int jacob_field_number = 0, jacob_exo_field_number = 0, jacob_exo_det_field_number = 0, jacob_other_endo_field_number = 0;
                  if (!block_structur)
                    {
                      const char *field_names[] = {"g1", "g1_x", "g1_xd", "g1_o"};
                      jacob_field_number = 0;
                      jacob_exo_field_number = 1;
                      jacob_exo_det_field_number = 2;
                      jacob_other_endo_field_number = 3;
                      mwSize dims[1] = {(mwSize)nb_blocks };
                      plhs[2] = mxCreateStructArray(1, dims, 4, field_names);
                    }
                  else if (!mxIsStruct(block_structur))
                    {
                      plhs[2] = interprete.get_jacob(0);
                      //mexCallMATLAB(0,NULL, 1, &plhs[2], "disp");
                      dont_store_a_structure = true;
                    }
                  else
                    {
                      plhs[2] = block_structur;
                      jacob_field_number = mxAddField(plhs[2], "g1");
                      if (jacob_field_number == -1)
                        DYN_MEX_FUNC_ERR_MSG_TXT("Fatal error in bytecode: in main, cannot add extra field jacob to the structArray\n");
                      jacob_exo_field_number = mxAddField(plhs[2], "g1_x");
                      if (jacob_exo_field_number == -1)
                        DYN_MEX_FUNC_ERR_MSG_TXT("Fatal error in bytecode: in main, cannot add extra field jacob_exo to the structArray\n");
                      jacob_exo_det_field_number = mxAddField(plhs[2], "g1_xd");
                      if (jacob_exo_det_field_number == -1)
                        DYN_MEX_FUNC_ERR_MSG_TXT("Fatal error in bytecode: in main, cannot add extra field jacob_exo_det to the structArray\n");
                      jacob_other_endo_field_number = mxAddField(plhs[2], "g1_o");
                      if (jacob_other_endo_field_number == -1)
                        DYN_MEX_FUNC_ERR_MSG_TXT("Fatal error in bytecode: in main, cannot add extra field jacob_other_endo to the structArray\n");
                    }
                  if (!dont_store_a_structure)
                    {
                      for (int i = 0; i < nb_blocks; i++)
                        {
                          mxSetFieldByNumber(plhs[2], i, jacob_field_number, interprete.get_jacob(i));
                          if (!steady_state)
                            {
                              mxSetFieldByNumber(plhs[2], i, jacob_exo_field_number, interprete.get_jacob_exo(i));
                              mxSetFieldByNumber(plhs[2], i, jacob_exo_det_field_number, interprete.get_jacob_exo_det(i));
                              mxSetFieldByNumber(plhs[2], i, jacob_other_endo_field_number, interprete.get_jacob_other_endo(i));
                            }
                        }
                    }
                }
              else
                {
                  plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
                  pind = mxGetPr(plhs[0]);
                  pind[0] = NAN;
                }
              if (nlhs > 3)
                {
                  plhs[3] = mxCreateDoubleMatrix(int(row_y), int(col_y), mxREAL);
                  pind = mxGetPr(plhs[3]);
                  for (i = 0; i < row_y*col_y; i++)
                    pind[i] = y[i];
                  if (nlhs > 4)
                    {
                      mxArray *GlobalTemporaryTerms = interprete.get_Temporary_Terms();
                      size_t nb_temp_terms = mxGetM(GlobalTemporaryTerms);
                      plhs[4] = mxCreateDoubleMatrix(int(nb_temp_terms), 1, mxREAL);
                      pind = mxGetPr(plhs[4]);
                      double *tt = mxGetPr(GlobalTemporaryTerms);
                      for (i = 0; i < nb_temp_terms; i++)
                        pind[i] = tt[i];
                    }
                }

            }
        }
    }
#else
  Free_global();
#endif
  if (x)
    mxFree(x);
  if (y)
    mxFree(y);
  if (ya)
    mxFree(ya);
  if (direction)
    mxFree(direction);
#ifdef _MSC_VER_
  /*fFreeResult =*/ FreeLibrary(hinstLib);
#endif
  return;
}
