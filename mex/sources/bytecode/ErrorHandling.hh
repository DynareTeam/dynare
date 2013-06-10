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

#ifndef ERROR_HANDLING
#define ERROR_HANDLING

#include <cstring>
#include <iostream>
#include <sstream>
#include <map>
#define BYTE_CODE
#include "CodeInterpreter.hh"
#ifdef DEBUG_EX
# include <math.h>
# include "mex_interface.hh"
#endif

#ifdef OCTAVE_MEX_FILE
# define CHAR_LENGTH 1
#else
# define CHAR_LENGTH 2
#endif

#ifdef _MSC_VER
#include <limits>
#define M_E 2.71828182845904523536
#define M_LOG2E 1.44269504088896340736
#define M_LOG10E 0.434294481903251827651
#define M_LN2 0.693147180559945309417
#define M_LN10 2.30258509299404568402
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define M_PI_4 0.785398163397448309616
#define M_1_PI 0.318309886183790671538
#define M_2_PI 0.636619772367581343076
#define M_1_SQRTPI 0.564189583547756286948
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2 1.41421356237309504880
#define M_SQRT_2 0.707106781186547524401
#define NAN numeric_limits<double>::quiet_NaN()

#define isnan(x) _isnan(x)
#define isinf(x) (!_finite(x))
#define fpu_error(x) (isinf(x) || isnan(x))

#define finite(x) _finite(x)

class MSVCpp_missings
{
  public:
  inline double
  asinh(double x) const
    {
      if(x==0.0)
        return 0.0;
      double ax = abs(x);
      return log(x+ax*sqrt(1.+1./(ax*ax)));
    }

  inline double
  acosh(double x) const
    {
      if(x==0.0)
        return 0.0;
      double ax = abs(x);
      return log(x+ax*sqrt(1.-1./(ax*ax)));
    }

  inline double
  atanh(double x) const
    {
      return log((1+x)/(1-x))/2;
    }

  inline double
  erf(double x) const
    {
      const double a1 = -1.26551223,   a2 = 1.00002368,
                   a3 =  0.37409196,   a4 = 0.09678418,
                   a5 = -0.18628806,   a6 = 0.27886807,
                   a7 = -1.13520398,   a8 = 1.48851587,
                   a9 = -0.82215223,  a10 = 0.17087277;
     double v = 1;
     double z = abs(x);
     if (z <= 0)
       return v;
     double t = 1 / (1 + 0.5 * z);
     v = t*exp((-z*z) +a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*(a9+t*a10)))))))));
     if (x < 0)
       v = 2 - v;
     return 1 - v;
    }

  inline double
  nearbyint(double x) const
    {
      return floor(x + 0.5);
    }

  inline double
  fmax(double x, double y) const
    {
      if (x > y)
        return x;
      else
        return y;
   }

  inline double
  fmin(double x, double y) const
    {
      if (x < y)
        return x;
      else
        return y;
    }

};
#endif


//#define DEBUG
using namespace std;

const int NO_ERROR_ON_EXIT = 0;
const int ERROR_ON_EXIT = 1;


typedef vector<pair<Tags, void * > > code_liste_type;
typedef code_liste_type::const_iterator it_code_type;


class GeneralExceptionHandling
{
  string ErrorMsg;
public:
#ifdef _MSC_VER_
  ~GeneralExceptionHandling()
  {
    FreeLibrary(hinstLib);
  };
#endif
  GeneralExceptionHandling(string ErrorMsg_arg) : ErrorMsg(ErrorMsg_arg)
  {
  };
  inline string
  GetErrorMsg()
  {
    return ErrorMsg;
  }
  inline void
  completeErrorMsg(string ErrorMsg_arg)
  {
    ErrorMsg += ErrorMsg_arg;
  }
};

class FloatingPointExceptionHandling : public GeneralExceptionHandling
{
public:
  FloatingPointExceptionHandling(string value) : GeneralExceptionHandling(string("Floating point error in bytecode: " + value))
  {
  };
};

class LogExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  LogExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log(X)"),
                                           value(value_arg)
  {
    ostringstream tmp;
    tmp << " with X=" << value << "\n";
    completeErrorMsg(tmp.str());
  };
};

class Log10ExceptionHandling : public FloatingPointExceptionHandling
{
  double value;
public:
  Log10ExceptionHandling(double value_arg) : FloatingPointExceptionHandling("log10(X)"),
                                             value(value_arg)
  {
    ostringstream tmp;
    tmp << " with X=" << value << "\n";
    completeErrorMsg(tmp.str());
  };
};

class DivideExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  DivideExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("a/X"),
                                                                  value1(value1_arg),
                                                                  value2(value2_arg)
  {
    ostringstream tmp;
    tmp << " with X=" << value2 << "\n";
    completeErrorMsg(tmp.str());
  };
};

class PowExceptionHandling : public FloatingPointExceptionHandling
{
  double value1, value2;
public:
  PowExceptionHandling(double value1_arg, double value2_arg) : FloatingPointExceptionHandling("X^a"),
                                                               value1(value1_arg),
                                                               value2(value2_arg)
  {
    ostringstream tmp;
    if (abs(value1) > 1e-10 )
      tmp << " with X=" << value1 << "\n";
    else
      tmp << " with X=" << value1 << " and a=" << value2 << "\n";
    completeErrorMsg(tmp.str());
  };
};

class UserExceptionHandling : public GeneralExceptionHandling
{
  double value;
public:
  UserExceptionHandling() : GeneralExceptionHandling("Fatal error in bytecode:")
  {
    completeErrorMsg(" User break\n");
  };
};

class FatalExceptionHandling : public GeneralExceptionHandling
{
public:
  FatalExceptionHandling(string ErrorMsg_arg) : GeneralExceptionHandling("Fatal error in bytecode:")
  {
    completeErrorMsg(ErrorMsg_arg);
  };
  FatalExceptionHandling() : GeneralExceptionHandling("")
  {
  };
};

struct s_plan
{
  string var, exo;
  int var_num, exo_num;
  vector<pair<int, double> > per_value;
};

#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#else
#include <octave/oct.h>
#include <octave/unwind-prot.h>
#endif

#ifdef _MSC_VER
class ErrorMsg : public MSVCpp_missings
#else
class ErrorMsg
#endif
{
private:
  bool is_load_variable_list;

public:
  double *y, *ya;
  int y_size;
  double *T;
  int nb_row_xd, nb_row_x;
  int y_kmin, y_kmax, periods;
  double *x, *params;
  double *u;
  double *steady_y, *steady_x;
  double *g2, *g1, *r, *res;
  vector<s_plan> splan, spfplan;
  vector<mxArray *> jacobian_block, jacobian_other_endo_block, jacobian_exo_block, jacobian_det_exo_block;
  map<unsigned int, double> TEF;
  map<pair<unsigned int, unsigned int>, double > TEFD;
  map<pair<unsigned int, pair<unsigned int, unsigned int> >, double > TEFDD;

  ExpressionType EQN_type;
  it_code_type it_code_expr;
  /*unsigned int*/size_t nb_endo, nb_exo, nb_param;
  char *P_endo_names, *P_exo_names, *P_param_names;
  size_t/*unsigned int*/ endo_name_length, exo_name_length, param_name_length;
  unsigned int EQN_equation, EQN_block, EQN_block_number;
  unsigned int EQN_dvar1, EQN_dvar2, EQN_dvar3;
  vector<pair<string, pair<SymbolType, unsigned int> > > Variable_list;

  inline
  ErrorMsg()
  {
    mxArray *M_ = mexGetVariable("global", "M_");
    nb_endo = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
    endo_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
    P_endo_names = (char *) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "endo_names")));
    nb_exo = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
    exo_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
    P_exo_names = (char *) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "exo_names")));
    nb_param = mxGetM(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
    param_name_length = mxGetN(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
    P_param_names = (char *) mxGetPr(mxGetFieldByNumber(M_, 0, mxGetFieldNumber(M_, "param_names")));
    is_load_variable_list = false;
  }

  inline string
  add_underscore_to_fpe(const string &str)
  {
    string temp;
    int pos1 = -1, pos2 = -1;
    string tmp_n(str.length(), ' ');
    for (unsigned int i = 0; i < str.length(); i++)
      {
        if (str[i] != '$' && str[i] != '£')
          temp += str[i];
        else
          {
            if (str[i] == '$')
              pos1 = int(temp.length());
            else
              pos2 = int(temp.length());
            if (pos1 >= 0 && pos2 >= 0)
              {
                tmp_n.erase(pos1, pos2-pos1+1);
                tmp_n.insert(pos1, pos2-pos1, '~');
                pos1 = pos2 = -1;
              }
          }
      }
    temp += "\n" + tmp_n;
    return temp;
  }

  inline void
  load_variable_list()
  {
    ostringstream res;
    for (unsigned int variable_num = 0; variable_num < (unsigned int)nb_endo; variable_num++)
      {
        for (unsigned int i = 0; i < endo_name_length; i++)
          if (P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)] != ' ')
            res << P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)];
        Variable_list.push_back(make_pair(res.str(), make_pair(eEndogenous, variable_num)));
      }
    for (unsigned int variable_num = 0; variable_num < (unsigned int)nb_exo; variable_num++)
      {
        for (unsigned int i = 0; i < exo_name_length; i++)
          if (P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)] != ' ')
            res << P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)];
        Variable_list.push_back(make_pair(res.str(), make_pair(eExogenous, variable_num)));
      }
  }

  inline int
  get_ID(const string variable_name, SymbolType *variable_type)
  {
    if (!is_load_variable_list)
      {
        load_variable_list();
        is_load_variable_list = true;
      }
    size_t n = Variable_list.size();
    int i = 0;
    bool notfound = true;
    while (notfound && i < n)
      {
        if (variable_name == Variable_list[i].first)
          {
            notfound = false;
            *variable_type = Variable_list[i].second.first;
            return Variable_list[i].second.second;
          }
        i++;
      }
    return(-1);
  }

  inline string
  get_variable(const SymbolType variable_type, const unsigned int variable_num) const
  {
    ostringstream res;
    switch (variable_type)
      {
      case eEndogenous:
        if (variable_num < nb_endo)
          {
            for (unsigned int i = 0; i < endo_name_length; i++)
              if (P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)] != ' ')
                res << P_endo_names[CHAR_LENGTH*(variable_num+i*nb_endo)];
          }
        else
          mexPrintf("=> Unknown endogenous variable n° %d", variable_num);
        break;
      case eExogenous:
      case eExogenousDet:
        if (variable_num < nb_exo)
          {
            for (unsigned int i = 0; i < exo_name_length; i++)
              if (P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)] != ' ')
                res << P_exo_names[CHAR_LENGTH*(variable_num+i*nb_exo)];
          }
        else
          mexPrintf("=> Unknown exogenous variable n° %d", variable_num);
        break;
      case eParameter:
        if (variable_num < nb_param)
          {
            for (unsigned int i = 0; i < param_name_length; i++)
              if (P_param_names[CHAR_LENGTH*(variable_num+i*nb_param)] != ' ')
                res << P_param_names[CHAR_LENGTH*(variable_num+i*nb_param)];
          }
        else
          mexPrintf("=> Unknown parameter n° %d", variable_num);
        break;
      default:
        break;
      }
    return (res.str());
  }

  inline string
  error_location(bool evaluate, bool steady_state, int size, int block_num, int it_, int Per_u_)
  {
    stringstream Error_loc("");
    if (!steady_state)
      switch (EQN_type)
        {
        case TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "temporary term " << EQN_equation+1 << " at time " << it_;
          break;
        case ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1 << " at time " << it_;
          else
            Error_loc << "equation " << EQN_equation+1 << " at time " << it_;
          break;
        case FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        case FirstOtherEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to other endogenous variable "  << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to other endogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        case FirstExoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        case FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        case FirstParamDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to parameter "  << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to parameter " << get_variable(eEndogenous, EQN_dvar1) << " at time " << it_;
          break;
        default:
          return ("???");
        }
    else
      switch (EQN_type)
        {
        case TemporaryTerm:
          if (EQN_block_number > 1)
            Error_loc << "temporary term " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "temporary term " << EQN_equation+1;
          break;
        case ModelEquation:
          if (EQN_block_number > 1)
            Error_loc << "equation " << EQN_equation+1 << " in block " << EQN_block+1;
          else
            Error_loc << "equation " << EQN_equation+1;
          break;
        case FirstEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to endogenous variable "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to endogenous variable " << get_variable(eEndogenous, EQN_dvar1);
          break;
        case FirstOtherEndoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to other endogenous variable "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to other endogenous variable " << get_variable(eEndogenous, EQN_dvar1);
          break;
        case FirstExoDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to exogenous variable "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to exogenous variable " << get_variable(eEndogenous, EQN_dvar1);
          break;
        case FirstExodetDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to deterministic exogenous variable "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to deterministic exogenous variable " << get_variable(eEndogenous, EQN_dvar1);
          break;
        case FirstParamDerivative:
          if (EQN_block_number > 1)
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " in block " << EQN_block+1 << " with respect to parameter "  << get_variable(eEndogenous, EQN_dvar1);
          else
            Error_loc << "first order derivative of equation " << EQN_equation+1 << " with respect to parameter " << get_variable(eEndogenous, EQN_dvar1);
          break;
        default:
          return ("???");
        }
    it_code_type it_code_ret;
    Error_loc << endl << add_underscore_to_fpe("      " + print_expression(it_code_expr, evaluate, size, block_num, steady_state, Per_u_, it_, it_code_ret, true));
    return (Error_loc.str());
  }

  inline string
  print_expression(it_code_type it_code, bool evaluate, int size, int block_num, bool steady_state, int Per_u_, int it_, it_code_type &it_code_ret, bool compute) const
  {
    int var, lag = 0, op, eq;
    stack<string> Stack;
    stack<double> Stackf;
    ostringstream tmp_out, tmp_out2;
    string v1, v2, v3;
    double v1f, v2f, v3f = 0.0;
    bool go_on = true;
    double ll;
    ExpressionType equation_type = TemporaryTerm;
    size_t found;
    double *jacob = NULL, *jacob_other_endo = NULL, *jacob_exo = NULL, *jacob_exo_det = NULL;
    external_function_type function_type = ExternalFunctionWithoutDerivative;

    if (evaluate)
      {
        jacob = mxGetPr(jacobian_block[block_num]);
        if (!steady_state)
          {
            jacob_other_endo = mxGetPr(jacobian_other_endo_block[block_num]);
            jacob_exo = mxGetPr(jacobian_exo_block[block_num]);
            jacob_exo_det = mxGetPr(jacobian_det_exo_block[block_num]);
          }
      }

    while (go_on)
      {
#ifdef OCTAVE_MEX_FILE
        OCTAVE_QUIT;
#else
	      if ( utIsInterruptPending() )
		      throw UserExceptionHandling();
#endif
        switch (it_code->first)
          {
          case FNUMEXPR:
            switch (((FNUMEXPR_ *) it_code->second)->get_expression_type())
              {
              case TemporaryTerm:
                equation_type = TemporaryTerm;
                break;
              case ModelEquation:
                equation_type = ModelEquation;
                break;
              case FirstEndoDerivative:
                equation_type = FirstEndoDerivative;
                break;
              case FirstOtherEndoDerivative:
                equation_type = FirstOtherEndoDerivative;
                break;
              case FirstExoDerivative:
                equation_type = FirstExoDerivative;
                break;
              case FirstExodetDerivative:
                equation_type = FirstExodetDerivative;
                break;
              case FirstParamDerivative:
                equation_type = FirstParamDerivative;
                break;
              case SecondEndoDerivative:
                equation_type = SecondEndoDerivative;
                break;
              case SecondExoDerivative:
                equation_type = SecondExoDerivative;
                break;
              case SecondExodetDerivative:
                equation_type = SecondExodetDerivative;
                break;
              case SecondParamDerivative:
                equation_type = SecondExodetDerivative;
                break;
              case ThirdEndoDerivative:
                equation_type = ThirdEndoDerivative;
                break;
              case ThirdExoDerivative:
                equation_type = ThirdExoDerivative;
                break;
              case ThirdExodetDerivative:
                equation_type = ThirdExodetDerivative;
                break;
              case ThirdParamDerivative:
                equation_type = ThirdExodetDerivative;
                break;
              default:
                ostringstream tmp;
                tmp << " in print_expression, expression type " << ((FNUMEXPR_ *) it_code->second)->get_expression_type() << " not implemented yet\n";
                throw FatalExceptionHandling(tmp.str());
              }
            break;
          case FLDV:
            //load a variable in the processor
            switch (((FLDV_ *) it_code->second)->get_type())
              {
              case eParameter:
                var = ((FLDV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FLDV_ Param var=%d\n", var);
                mexPrintf("get_variable(eParameter, var)=%s\n",get_variable(eParameter, var).c_str());
                mexEvalString("drawnow;");
#endif
                Stack.push(get_variable(eParameter, var));
                if (compute)
                  Stackf.push(params[var]);
                break;
              case eEndogenous:
                var = ((FLDV_ *) it_code->second)->get_pos();
                lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FLDV_ endo var=%d, lag=%d\n", var, lag);
                mexPrintf("get_variable(eEndogenous, var)=%s, compute=%d\n",get_variable(eEndogenous, var).c_str(), compute);
                mexPrintf("it_=%d, lag=%d, y_size=%d, var=%d, y=%x\n", it_, lag, y_size, var, y);
                mexEvalString("drawnow;");
#endif
                tmp_out.str("");
                if (lag > 0)
                  tmp_out << get_variable(eEndogenous, var) << "(+" << lag << ")";
                else if (lag < 0)
                  tmp_out << get_variable(eEndogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(eEndogenous, var);
                Stack.push(tmp_out.str());
                if (compute)
                  {
                    if (evaluate)
                      Stackf.push(ya[(it_+lag)*y_size+var]);
                    else
                      Stackf.push(y[(it_+lag)*y_size+var]);
                  }
                break;
              case eExogenous:
                var = ((FLDV_ *) it_code->second)->get_pos();
                lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FLDV_ exo var=%d, lag=%d", var, lag);
#endif
                tmp_out.str("");
                if (lag != 0)
                  tmp_out << get_variable(eExogenous, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(eExogenous, var);
                Stack.push(tmp_out.str());
                if (compute)
                  Stackf.push(x[it_+lag+var*nb_row_x]);
                break;
              case eExogenousDet:
                var = ((FLDV_ *) it_code->second)->get_pos();
                lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FLDV_ exo_det var=%d, lag=%d", var, lag);
#endif
                tmp_out.str("");
                if (lag != 0)
                  tmp_out << get_variable(eExogenousDet, var) << "(" << lag << ")";
                else
                  tmp_out << get_variable(eExogenousDet, var);
                Stack.push(tmp_out.str());
                if (compute)
                  Stackf.push(x[it_+lag+var*nb_row_xd]);
                break;
              case eModelLocalVariable:
                break;
              default:
                mexPrintf("FLDV: Unknown variable type\n");
              }
            break;
          case FLDSV:
          case FLDVS:
            //load a variable in the processor
            switch (((FLDSV_ *) it_code->second)->get_type())
              {
              case eParameter:
                var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FLDSV_ param var=%d", var);
#endif
                Stack.push(get_variable(eParameter, var));
                if (compute)
                  Stackf.push(params[var]);
                break;
              case eEndogenous:
                var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FLDSV_ endo var=%d", var);
#endif
                Stack.push(get_variable(eEndogenous, var));
                if (compute)
                  {
                    if (it_code->first == FLDSV)
                      {
                        if (evaluate)
                          Stackf.push(ya[var]);
                        else
                          Stackf.push(y[var]);
                      }
                    else
                      Stackf.push(steady_y[var]);
                  }
                break;
              case eExogenous:
                var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FLDSV_ exo var=%d", var);
#endif
                Stack.push(get_variable(eExogenous, var));
#ifdef DEBUG
                mexPrintf("oka var=%d, Stack.size()=%d x=%x\n", var, Stack.size(), x);
#endif
                if (compute)
                  Stackf.push(x[var]);
                break;
              case eExogenousDet:
                var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FLDSV_ exo_det var=%d", var);
#endif
                Stack.push(get_variable(eExogenousDet, var));
                if (compute)
                  Stackf.push(x[var]);
                break;
              case eModelLocalVariable:
                break;
              default:
                mexPrintf("FLDSV: Unknown variable type\n");
              }
            break;
          case FLDT:
            //load a temporary variable in the processor
            var = ((FLDT_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FLDT_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(T[var*(periods+y_kmin+y_kmax)+it_]);
            break;
          case FLDST:
            //load a temporary variable in the processor
            var = ((FLDST_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FLDST_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "T" << var+1;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(T[var]);
            break;
          case FLDU:
            //load u variable in the processor
            var = ((FLDU_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FLDU_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_)";
            Stack.push(tmp_out.str());
            var += Per_u_;
            if (compute)
              Stackf.push(u[var]);
            break;
          case FLDSU:
            //load u variable in the processor
            var = ((FLDSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FLDSU_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ")";
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(u[var]);
            break;
          case FLDR:
            var = ((FLDR_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FLDSR_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ")";
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(r[var]);
            break;
          case FLDZ:
            //load 0 in the processor
#ifdef DEBUG
            mexPrintf("FLDZ_");
#endif
            tmp_out.str("");
            tmp_out << 0;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(0.0);
            break;
          case FLDC:
            //load a numerical constant in the processor
            ll = ((FLDC_ *) it_code->second)->get_value();
            tmp_out.str("");
#ifdef DEBUG
            mexPrintf("FLDC_ ll=%f", ll);
#endif
            tmp_out << ll;
            Stack.push(tmp_out.str());
            if (compute)
              Stackf.push(ll);
            break;
          case FSTPV:
            //load a variable in the processor
            go_on = false;
            switch (((FSTPV_ *) it_code->second)->get_type())
              {
              case eParameter:
                var = ((FSTPV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FSTPV_ param var=%d", var);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eParameter, var) << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    params[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case eEndogenous:
                var = ((FSTPV_ *) it_code->second)->get_pos();
                lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FSTPV_ endo var=%d, lag=%d", var, lag);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eEndogenous, var);
                if (lag > 0)
                  tmp_out << "(+" << lag << ")";
                else if (lag < 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    y[(it_+lag)*y_size+var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case eExogenous:
                var = ((FSTPV_ *) it_code->second)->get_pos();
                lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FSTPV_ exo var=%d, lag=%d", var, lag);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eExogenous, var);
                if (lag != 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[it_+lag+var*nb_row_x]  = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case eExogenousDet:
                var = ((FSTPV_ *) it_code->second)->get_pos();
                lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
                mexPrintf("FSTPV_ exodet var=%d, lag=%d", var, lag);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eExogenousDet, var);
                if (lag != 0)
                  tmp_out << "(" << lag << ")";
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[it_+lag+var*nb_row_xd] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              default:
                mexPrintf("FSTPV: Unknown variable type\n");
              }
            break;
          case FSTPSV:
            go_on = false;
            //load a variable in the processor
            switch (((FSTPSV_ *) it_code->second)->get_type())
              {
              case eParameter:
                var = ((FSTPSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FSTPSV_ param var=%d", var);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eParameter, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    params[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case eEndogenous:
                var = ((FSTPSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FSTPSV_ endo var=%d", var);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eEndogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    y[var] = Stackf.top();
                    Stackf.pop();
                  }
                break;
              case eExogenous:
              case eExogenousDet:
                var = ((FSTPSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
                mexPrintf("FSTPSV_ exo var=%d", var);
#endif
                tmp_out2.str("");
                tmp_out2 << Stack.top();
                tmp_out.str("");
                tmp_out << get_variable(eExogenous, var);
                tmp_out << " = " << tmp_out2.str();
                Stack.pop();
                if (compute)
                  {
                    x[var]  = Stackf.top();
                    Stackf.pop();
                  }
                break;
              default:
                mexPrintf("FSTPSV: Unknown variable type\n");
              }
            break;
          case FSTPT:
            go_on = false;
            //store in a temporary variable from the processor
            var = ((FSTPT_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTPT_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                T[var*(periods+y_kmin+y_kmax)+it_] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPST:
            go_on = false;
            //store in a temporary variable from the processor
            var = ((FSTPST_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTPST_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "T" << var+1 << " = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                T[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPU:
            go_on = false;
            //store in u variable from the processor
            var = ((FSTPU_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTPU_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "u(" << var+1 << " + it_) = " << Stack.top();
            var += Per_u_;
            Stack.pop();
            if (compute)
              {
                u[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPSU:
            go_on = false;
            //store in u variable from the processor
            var = ((FSTPSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTPSU_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "u(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                u[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPR:
            go_on = false;
            //store in residual variable from the processor
            var = ((FSTPR_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTPR_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "residual(" << var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                r[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPG:
            go_on = false;
            //store in derivative (g) variable from the processor
            var = ((FSTPG_ *) it_code->second)->get_pos();
#ifdef DEBUG
            mexPrintf("FSTG_ var=%d", var);
#endif
            tmp_out.str("");
            tmp_out << "g1[" << var+1 << "] = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                g1[var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPG2:
            go_on = false;
            //store in derivative (g) variable from the processor
            eq = ((FSTPG2_ *) it_code->second)->get_row();
            var = ((FSTPG2_ *) it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FSTG2_ eq=%d var=%d", eq, var);
#endif
            tmp_out.str("");
            tmp_out << "jacob(" << eq+size*var+1 << ") = " << Stack.top();
            Stack.pop();
            if (compute)
              {
                jacob[eq + size*var] = Stackf.top();
                Stackf.pop();
              }
            break;
          case FSTPG3:
            //store in derivative (g) variable from the processor
#ifdef DEBUG
            mexPrintf("FSTPG3\n");
            mexEvalString("drawnow;");
#endif
            double r;
            unsigned int pos_col;
            go_on = false;
            if (compute)
              {
                r = Stackf.top();
                Stackf.pop();
              }
            eq = ((FSTPG3_ *) it_code->second)->get_row();
            var = ((FSTPG3_ *) it_code->second)->get_col();
            lag = ((FSTPG3_ *) it_code->second)->get_lag();
            pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
            switch (equation_type)
              {
              case FirstEndoDerivative:
#ifdef DEBUG
                mexPrintf("Endo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
#endif
                if (compute)
                  jacob[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "jacob(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case FirstOtherEndoDerivative:
                if (compute)
                  jacob_other_endo[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "jacob_other_endo(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case FirstExoDerivative:
#ifdef DEBUG
                mexPrintf("Exo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
#endif
                if (compute)
                  jacob_exo[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "jacob_exo(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              case FirstExodetDerivative:
                if (compute)
                  jacob_exo_det[eq + size*pos_col] = r;
                tmp_out.str("");
                tmp_out << "jacob_exo_det(" << eq+size*pos_col+1 << ") = " << Stack.top();
                Stack.pop();
                break;
              default:
                ostringstream tmp;
                tmp << " in compute_block_time, variable " << EQN_type << " not used yet\n";
                //throw FatalExceptionHandling(tmp.str());
                mexPrintf("%s", tmp.str().c_str());
              }
#ifdef DEBUG
            tmp_out << "=>";
            mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
            tmp_out.str("");
#endif
            break;
          case FBINARY:
            op = ((FBINARY_ *) it_code->second)->get_op_type();
            v2 = Stack.top();
            Stack.pop();
            v1 = Stack.top();
            Stack.pop();
            if (compute)
              {
                v2f = Stackf.top();
                Stackf.pop();
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (op)
              {
              case oPlus:
#ifdef DEBUG
                mexPrintf("+");
#endif
                if (compute)
                  Stackf.push(v1f + v2f);
                tmp_out.str("");
                tmp_out << v1 << " + " << v2;
                Stack.push(tmp_out.str());
                break;
              case oMinus:
#ifdef DEBUG
                mexPrintf("-");
#endif
                if (compute)
                  Stackf.push(v1f - v2f);
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " - ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oTimes:
#ifdef DEBUG
                mexPrintf("*");
#endif
                if (compute)
                  Stackf.push(v1f * v2f);
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " * ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oDivide:
#ifdef DEBUG
                mexPrintf("/");
#endif
                if (compute)
                  {
                    r = v1f / v2f;
                    Stackf.push(r);
                  }
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (compute)
                  {
                    if (isinf(r))
                      tmp_out << "$";
                    tmp_out << " / ";
                    if (isinf(r))
                      tmp_out << "£";
                  }
                else
                  tmp_out << " / ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oLess:
#ifdef DEBUG
                mexPrintf("<");
#endif
                if (compute)
                  Stackf.push(double (v1f < v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " < ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oGreater:
#ifdef DEBUG
                mexPrintf(">");
#endif
                if (compute)
                  Stackf.push(double (v1f > v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " > ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oLessEqual:
#ifdef DEBUG
                mexPrintf("<=");
#endif
                if (compute)
                  Stackf.push(double (v1f <= v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " <= ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oGreaterEqual:
#ifdef DEBUG
                mexPrintf(">=");
#endif
                if (compute)
                  Stackf.push(double (v1f >= v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " >= ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oEqualEqual:
#ifdef DEBUG
                mexPrintf("==");
#endif
                if (compute)
                  Stackf.push(double (v1f == v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " == ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oDifferent:
#ifdef DEBUG
                mexPrintf("!=");
#endif
                if (compute)
                  Stackf.push(double (v1f != v2f));
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                tmp_out << " != ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oPower:
#ifdef DEBUG
                mexPrintf("^");
#endif
                if (compute)
                  {
                    r = pow(v1f, v2f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                found = v1.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v1;
                if (found != string::npos)
                  tmp_out << ")";
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$ ^ £";
                    else
                      tmp_out << " ^ ";
                  }
                else
                  tmp_out << " ^ ";
                found = v2.find(" ");
                if (found != string::npos)
                  tmp_out << "(";
                tmp_out << v2;
                if (found != string::npos)
                  tmp_out << ")";
                Stack.push(tmp_out.str());
                break;
              case oPowerDeriv:
                {
                  v3 = Stack.top();
                  Stack.pop();
                  if (compute)
                    {
                      int derivOrder = int(nearbyint(Stackf.top()));
                      Stackf.pop();
                      if (fabs(v1f) < NEAR_ZERO && v2f > 0
                          && derivOrder > v2f
                          && fabs(v2f-nearbyint(v2f)) < NEAR_ZERO)
                        {
                          r = 0.0;
                          Stackf.push(r);
                        }
                      else
                        {
                          double dxp = pow(v1f, v2f-derivOrder);
                          for (int i = 0; i < derivOrder; i++)
                            dxp *= v2f--;
                          Stackf.push(dxp);
                          r = dxp;
                        }
                    }
                  tmp_out.str("");
                  if (compute)
                    {
                      if (isnan(r))
                        tmp_out << "$ PowerDeriv £";
                      else
                        tmp_out << "PowerDeriv";
                    }
                  else
                    tmp_out << "PowerDeriv";
                  tmp_out << "(" << v1 << ", " << v2 << ", " << v3 << ")";
                  Stack.push(tmp_out.str());
                }
#ifdef DEBUG
                tmp_out << " |PowerDeriv(" << v1 << ", " << v2 << v3 << ")|";
#endif
                break;
              case oMax:
#ifdef DEBUG
                mexPrintf("max");
#endif
                if (compute)
                  Stackf.push(max(v1f, v2f));
                tmp_out.str("");
                tmp_out << "max(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case oMin:
#ifdef DEBUG
                mexPrintf("min");
#endif
                if (compute)
                  Stackf.push(min(v1f, v2f));
                tmp_out.str("");
                tmp_out << "min(" << v1 << ", " << v2 << ")";
                Stack.push(tmp_out.str());
                break;
              case oEqual:
              default:
                mexPrintf("Error unknown Unary operator=%d\n", op); mexEvalString("drawnow;");
                ;
              }
            break;
          case FUNARY:
            op = ((FUNARY_ *) it_code->second)->get_op_type();
            v1 = Stack.top();
            Stack.pop();
            if (compute)
              {
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (op)
              {
              case oUminus:
                if (compute)
                  Stackf.push(-v1f);
                tmp_out.str("");
                tmp_out << " -" << v1;
                Stack.push(tmp_out.str());
                break;
              case oExp:
                if (compute)
                  Stackf.push(exp(v1f));
                tmp_out.str("");
                tmp_out << "exp(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oLog:
                if (compute)
                  {
                    r = log(v1f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$log£(" << v1 << ")";
                    else
                      tmp_out << "log(" << v1 << ")";
                  }
                else
                  tmp_out << "log(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oLog10:
                if (compute)
                  {
                    r = log10(v1f);
                    Stackf.push(r);
                  }
                tmp_out.str("");
                if (compute)
                  {
                    if (isnan(r))
                      tmp_out << "$log10£(" << v1 << ")";
                    else
                      tmp_out << "log10(" << v1 << ")";
                  }
                else
                  tmp_out << "log10(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oCos:
                if (compute)
                  Stackf.push(cos(v1f));
                tmp_out.str("");
                tmp_out << "cos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oSin:
                if (compute)
                  Stackf.push(sin(v1f));
                tmp_out.str("");
                tmp_out << "sin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oTan:
                if (compute)
                  Stackf.push(tan(v1f));
                tmp_out.str("");
                tmp_out << "tan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAcos:
                if (compute)
                  Stackf.push(acos(v1f));
                tmp_out.str("");
                tmp_out << "acos(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAsin:
                if (compute)
                  Stackf.push(asin(v1f));
                tmp_out.str("");
                tmp_out << "asin(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAtan:
                if (compute)
                  Stackf.push(atan(v1f));
                tmp_out.str("");
                tmp_out << "atan(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oCosh:
                if (compute)
                  Stackf.push(cosh(v1f));
                tmp_out.str("");
                tmp_out << "cosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oSinh:
                if (compute)
                  Stackf.push(sinh(v1f));
                tmp_out.str("");
                tmp_out << "sinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oTanh:
                if (compute)
                  Stackf.push(tanh(v1f));
                tmp_out.str("");
                tmp_out << "tanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAcosh:
                if (compute)
                  Stackf.push(acosh(v1f));
                tmp_out.str("");
                tmp_out << "acosh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAsinh:
                if (compute)
                  Stackf.push(asinh(v1f));
                tmp_out.str("");
                tmp_out << "asinh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oAtanh:
                if (compute)
                  Stackf.push(atanh(v1f));
                tmp_out.str("");
                tmp_out << "atanh(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oSqrt:
                if (compute)
                  Stackf.push(sqrt(v1f));
                tmp_out.str("");
                tmp_out << "sqrt(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              case oErf:
                if (compute)
                  Stackf.push(erf(v1f));
                tmp_out.str("");
                tmp_out << "erf(" << v1 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown Binary operator=%d\n", op); mexEvalString("drawnow;");
                ;
              }
            break;
          case FTRINARY:
            op = ((FTRINARY_ *) it_code->second)->get_op_type();
            v3 = Stack.top();
            Stack.pop();
            v2 = Stack.top();
            Stack.pop();
            v1 = Stack.top();
            Stack.pop();
            if (compute)
              {
                v3f = Stackf.top();
                Stackf.pop();
                v2f = Stackf.top();
                Stackf.pop();
                v1f = Stackf.top();
                Stackf.pop();
              }
            switch (op)
              {
              case oNormcdf:
                if (compute)
                  Stackf.push(0.5*(1+erf((v1f-v2f)/v3f/M_SQRT2)));
                tmp_out.str("");
                tmp_out << "normcdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              case oNormpdf:
                if (compute)
                  Stackf.push(1/(v3f*sqrt(2*M_PI)*exp(pow((v1f-v2f)/v3f, 2)/2)));
                tmp_out.str("");
                tmp_out << "normpdf(" << v1 << ", " << v2 << ", " << v3 << ")";
                Stack.push(tmp_out.str());
                break;
              default:
                mexPrintf("Error unknown Trinary operator=%d\n", op); mexEvalString("drawnow;");
              }
            break;
          case FCALL:
            {
#ifdef DEBUG
              mexPrintf("------------------------------\n");
              mexPrintf("CALL "); mexEvalString("drawnow;");
#endif
              FCALL_ *fc = (FCALL_ *) it_code->second;
              string function_name = fc->get_function_name();
#ifdef DEBUG
              mexPrintf("function_name=%s ", function_name.c_str()); mexEvalString("drawnow;");
#endif
              unsigned int nb_input_arguments = fc->get_nb_input_arguments();
#ifdef DEBUG
              mexPrintf("nb_input_arguments=%d ", nb_input_arguments); mexEvalString("drawnow;");
#endif
              unsigned int nb_output_arguments = fc->get_nb_output_arguments();
#ifdef DEBUG
              mexPrintf("nb_output_arguments=%d\n", nb_output_arguments); mexEvalString("drawnow;");
#endif

              mxArray *output_arguments[3];
              string arg_func_name = fc->get_arg_func_name();
#ifdef DEBUG
              mexPrintf("arg_func_name.length() = %d\n", arg_func_name.length());
              mexPrintf("arg_func_name.c_str() = %s\n", arg_func_name.c_str());
#endif
              unsigned int nb_add_input_arguments = fc->get_nb_add_input_arguments();
              function_type = fc->get_function_type();
#ifdef DEBUG
              mexPrintf("function_type=%d ExternalFunctionWithoutDerivative=%d\n", function_type, ExternalFunctionWithoutDerivative);
              mexEvalString("drawnow;");
#endif
              mxArray **input_arguments;
              switch (function_type)
                {
                case ExternalFunctionWithoutDerivative:
                case ExternalFunctionWithFirstDerivative:
                case ExternalFunctionWithFirstandSecondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
#ifdef DEBUG
                        mexPrintf("Stack.size()=%d\n", Stack.size());
                        mexEvalString("drawnow;");
#endif
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[nb_input_arguments - i - 1] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
#ifndef _MSC_VER
                    string ss[nb_input_arguments];
#else
                    vector<string> ss(nb_input_arguments);
#endif
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionNumericalFirstDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = (mxArray **) mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                        mxArray *vv = mxCreateString(arg_func_name.c_str());
                        input_arguments[0] = vv;
                        vv = mxCreateDoubleScalar(fc->get_row());
                        input_arguments[1] = vv;
                        vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                        for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                          {
                            double rr = Stackf.top();
#ifdef DEBUG
                            mexPrintf("i=%d rr = %f Stack.size()=%d\n", i, rr, Stack.size());
#endif
                            mxSetCell(vv, nb_add_input_arguments - (i+1), mxCreateDoubleScalar(rr));
                            Stackf.pop();
                          }
                        input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                        mexCallMATLAB(0, NULL, 1, &input_arguments[0], "disp");
                        mexCallMATLAB(0, NULL, 1, &input_arguments[1], "disp");
                        mexCallMATLAB(0, NULL, 1, &input_arguments[2], "celldisp");
                        mexPrintf("OK\n");
                        mexEvalString("drawnow;");
#endif
                        nb_input_arguments = 3;
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
#ifdef DEBUG
                        mexPrintf("*rr=%f\n", *rr);
#endif
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", {";
#ifndef _MSC_VER
                    string ss[nb_add_input_arguments];
#else
                    vector<string> ss(nb_input_arguments);
#endif
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionFirstDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[(nb_input_arguments - 1) - i] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
#ifndef _MSC_VER
                    string ss[nb_input_arguments];
#else
                    vector<string> ss(nb_input_arguments);
#endif
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionNumericalSecondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = (mxArray **) mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                        mxArray *vv = mxCreateString(arg_func_name.c_str());
                        input_arguments[0] = vv;
                        vv = mxCreateDoubleScalar(fc->get_row());
                        input_arguments[1] = vv;
                        vv = mxCreateDoubleScalar(fc->get_col());
                        input_arguments[2] = vv;
                        vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                        for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                          {
                            double rr = Stackf.top();
#ifdef DEBUG
                            mexPrintf("i=%d rr = %f\n", i, rr);
#endif
                            mxSetCell(vv, (nb_add_input_arguments - 1) - i, mxCreateDoubleScalar(rr));
                            Stackf.pop();
                          }
                        input_arguments[nb_input_arguments+nb_add_input_arguments] = vv;
#ifdef DEBUG
                        mexCallMATLAB(0, NULL, 1, &input_arguments[0], "disp");
                        mexCallMATLAB(0, NULL, 1, &input_arguments[1], "disp");
                        mexCallMATLAB(0, NULL, 1, &input_arguments[2], "celldisp");
                        mexPrintf("OK\n");
                        mexEvalString("drawnow;");
#endif
                        nb_input_arguments = 3;
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                        double *rr = mxGetPr(output_arguments[0]);
                        Stackf.push(*rr);
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
                    tmp_out << arg_func_name.c_str() << ", " << fc->get_row() << ", " << fc->get_col() << ", {";
#ifndef _MSC_VER
                    string ss[nb_add_input_arguments];
#else
                    vector<string> ss(nb_input_arguments);
#endif
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        ss[nb_add_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_add_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << "})";
                    Stack.push(tmp_out.str());
                  }
                  break;
                case ExternalFunctionSecondDerivative:
                  {
                    if (compute)
                      {
                        input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                        for (unsigned int i = 0; i < nb_input_arguments; i++)
                          {
                            mxArray *vv = mxCreateDoubleScalar(Stackf.top());
                            input_arguments[i] = vv;
                            Stackf.pop();
                          }
                        mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str());
                      }
                    tmp_out.str("");
                    tmp_out << function_name << "(";
#ifndef _MSC_VER
                    string ss[nb_input_arguments];
#else
                    vector<string> ss(nb_input_arguments);
#endif
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        ss[nb_input_arguments-i-1] = Stack.top();
                        Stack.pop();
                      }
                    for (unsigned int i = 0; i < nb_input_arguments; i++)
                      {
                        tmp_out << ss[i];
                        if (i < nb_input_arguments - 1)
                          tmp_out << ", ";
                      }
                    tmp_out << ")";
                    Stack.push(tmp_out.str());
                  }
                  break;
                }
#ifdef DEBUG
              mexPrintf("end CALL\n");
              mexEvalString("drawnow;");
#endif
              break;
            }
          case FSTPTEF:
            go_on = false;
            var = ((FSTPTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
            mexPrintf("FSTPTEF\n");
            mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
#endif
            if (compute)
              {
#ifdef DEBUG
                double rr = Stackf.top();
                mexPrintf("FSTP TEF(var-1)=%f done\n", rr);
                mexEvalString("drawnow;");
#endif
                Stackf.pop();
              }
            tmp_out.str("");
            switch (function_type)
              {
              case ExternalFunctionWithoutDerivative:
                tmp_out << "TEF(" << var << ") = " << Stack.top();
                break;
              case ExternalFunctionWithFirstDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << ") ]= " << Stack.top();
                break;
              case ExternalFunctionWithFirstandSecondDerivative:
                tmp_out << "[TEF(" << var << "), TEFD(" << var << "), TEFDD(" << var << ") ]= " << Stack.top();
                break;
              default:
                break;
              }
            Stack.pop();
#ifdef DEBUG
            mexPrintf("end FSTPEF\n");
            mexEvalString("drawnow;");
#endif
            break;
          case FLDTEF:
            var = ((FLDTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
            mexPrintf("FLDTEF\n");
            mexPrintf("var=%d Stack.size()=%d\n", var, Stackf.size());
            {
              map<unsigned int, double>::const_iterator it = TEF.find(var-1);
              mexPrintf("FLD TEF[var-1]=%f done\n", it->second);
            }
            mexEvalString("drawnow;");
#endif
            if (compute)
              {
                map<unsigned int, double>::const_iterator it = TEF.find(var-1);
                Stackf.push(it->second);
              }
            tmp_out.str("");
            tmp_out << "TEF(" << var << ")";
            Stack.push(tmp_out.str());
#ifdef DEBUG
            mexPrintf("end FLDTEF\n");
            mexEvalString("drawnow;");
#endif

            break;
          case FSTPTEFD:
            {
              go_on = false;
              unsigned int indx = ((FSTPTEFD_ *) it_code->second)->get_indx();
              unsigned int row = ((FSTPTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
              mexPrintf("FSTPTEFD\n");
              mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
              if (compute)
                {
#ifdef DEBUG
                  double rr = Stackf.top();
                  mexPrintf("FSTP TEFD[make_pair(indx, row)]=%f done\n", rr);
                  mexEvalString("drawnow;");
#endif
                  Stackf.pop();
                }
              tmp_out.str("");
              if (function_type == ExternalFunctionNumericalFirstDerivative)
                tmp_out << "TEFD(" << indx << ", " << row << ") = " << Stack.top();
              else if (function_type == ExternalFunctionFirstDerivative)
                tmp_out << "TEFD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }
            break;
          case FLDTEFD:
            {
              unsigned int indx = ((FLDTEFD_ *) it_code->second)->get_indx();
              unsigned int row = ((FLDTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
              mexPrintf("FLDTEFD\n");
              mexPrintf("indx=%d row=%d Stack.size()=%d\n", indx, row, Stack.size());
              map<pair<unsigned int, unsigned int>, double>::const_iterator it = TEFD.find(make_pair(indx, row-1));
              mexPrintf("FLD TEFD[make_pair(indx, row)]=%f done\n", it->second);
              mexEvalString("drawnow;");
#endif
              if (compute)
                {
                  map<pair<unsigned int, unsigned int>, double>::const_iterator it = TEFD.find(make_pair(indx, row-1));
                  Stackf.push(it->second);
                }
              tmp_out.str("");
              tmp_out << "TEFD(" << indx << ", " << row << ")";
              Stack.push(tmp_out.str());
            }
            break;
          case FSTPTEFDD:
            {
              go_on = false;
              unsigned int indx = ((FSTPTEFDD_ *) it_code->second)->get_indx();
              unsigned int row = ((FSTPTEFDD_ *) it_code->second)->get_row();
              unsigned int col = ((FSTPTEFDD_ *) it_code->second)->get_col();
#ifdef DEBUG
              mexPrintf("FSTPTEFD\n");
              mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
              if (compute)
                {
#ifdef DEBUG
                  double rr = Stackf.top();
                  mexPrintf("rr=%f\n", rr);
                  map<pair<unsigned int, pair<unsigned int, unsigned int> >, double>::const_iterator it = TEFDD.find(make_pair(indx, make_pair(row-1, col-1)));
                  mexPrintf("FSTP TEFDD[make_pair(indx, make_pair(row, col))]=%f done\n", it->second);
                  mexEvalString("drawnow;");
#endif
                  Stackf.pop();
                }
              tmp_out.str("");
              if (function_type == ExternalFunctionNumericalSecondDerivative)
                tmp_out << "TEFDD(" << indx << ", " << row << ", " << col << ") = " << Stack.top();
              else if (function_type == ExternalFunctionSecondDerivative)
                tmp_out << "TEFDD(" << indx << ") = " << Stack.top();
              Stack.pop();
            }

            break;
          case FLDTEFDD:
            {
              unsigned int indx = ((FLDTEFDD_ *) it_code->second)->get_indx();
              unsigned int row = ((FLDTEFDD_ *) it_code->second)->get_row();
              unsigned int col = ((FSTPTEFDD_ *) it_code->second)->get_col();
#ifdef DEBUG
              mexPrintf("FLDTEFD\n");
              mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
              map<pair<unsigned int, pair<unsigned int, unsigned int> >, double>::const_iterator it = TEFDD.find(make_pair(indx, make_pair(row-1, col-1)));
              mexPrintf("FLD TEFD[make_pair(indx, make_pair(row, col))]=%f done\n", it->second);
              mexEvalString("drawnow;");
#endif
              if (compute)
                {
                  map<pair<unsigned int, pair<unsigned int, unsigned int> >, double>::const_iterator it = TEFDD.find(make_pair(indx, make_pair(row-1, col-1)));
                  Stackf.push(it->second);
                }
              tmp_out.str("");
              tmp_out << "TEFDD(" << indx << ", " << row << ", " << col << ")";
              Stack.push(tmp_out.str());
            }
            break;
          case FJMPIFEVAL:
            tmp_out.str("");
            tmp_out << "if (~evaluate)";
            go_on = false;
            break;
          case FJMP:
            tmp_out.str("");
            tmp_out << "else";
            go_on = false;
            break;
          case FCUML:
            if (compute)
              {
                v1f = Stackf.top();
                Stackf.pop();
                v2f = Stackf.top();
                Stackf.pop();
                Stackf.push(v1f+v2f);
              }
            v1 = Stack.top();
            Stack.pop();
            v2 = Stack.top();
            Stack.pop();
            tmp_out.str("");
            tmp_out << v1 << " + " << v2;
            Stack.push(tmp_out.str());
            break;
          case FENDBLOCK:
          case FENDEQU:
            go_on = false;
            break;
          case FOK:
            break;
          default:
            ostringstream tmp;
            mexPrintf("Error it_code->first=%d unknown\n", it_code->first); mexEvalString("drawnow;");
            tmp << " in print_expression, unknown opcode " << it_code->first << "!! FENDEQU=" << FENDEQU << "\n";
            throw FatalExceptionHandling(tmp.str());
          }
        it_code++;
      }
#ifdef DEBUG
    mexPrintf("print_expression end tmp_out.str().c_str()=%s\n", tmp_out.str().c_str()); mexEvalString("drawnow;");
#endif
    it_code_ret = it_code;
    return (tmp_out.str());
  }

};

#endif
