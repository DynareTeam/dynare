#include <cstring>
#include <sstream>
#include <math.h>
#include "Evaluate.hh"

#ifdef MATLAB_MEX_FILE
extern "C" bool utIsInterruptPending();
#else
#include <octave/oct.h>
#include <octave/unwind-prot.h>
#endif

Evaluate::Evaluate()
{
  symbol_table_endo_nbr = 0;
  Block_List_Max_Lag = 0;
  Block_List_Max_Lead = 0;
  u_count_int = 0;
  block = -1;
}

Evaluate::Evaluate(const int y_size_arg, const int y_kmin_arg, const int y_kmax_arg, const bool print_it_arg, const bool steady_state_arg, const int periods_arg, const int minimal_solving_periods_arg, const double slowc_arg):
print_it(print_it_arg),  minimal_solving_periods(minimal_solving_periods_arg)
{
  symbol_table_endo_nbr = 0;
  Block_List_Max_Lag = 0;
  Block_List_Max_Lead = 0;
  u_count_int = 0;
  block = -1;
  y_size = y_size_arg;
  y_kmin = y_kmin_arg;
  y_kmax  = y_kmax_arg;
  periods = periods_arg;
  steady_state = steady_state_arg;
  slowc = slowc_arg;
}

double
Evaluate::pow1(double a, double b)
{
  double r = pow_(a, b);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = 0.0000000000000000000000001;
      if (print_error)
        throw PowExceptionHandling(a, b);
    }
  return r;
}

double
Evaluate::divide(double a, double b)
{
  double r = a / b;
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = 1e70;
      if (print_error)
        throw DivideExceptionHandling(a, b);
    }
  return r;
}

double
Evaluate::log1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = -1e70;
      if (print_error)
        throw LogExceptionHandling(a);
    }
  return r;
}

double
Evaluate::log10_1(double a)
{
  double r = log(a);
  if (isnan(r) || isinf(r))
    {
      res1 = NAN;
      r = -1e70;
      if (print_error)
        throw Log10ExceptionHandling(a);
    }
  return r;
}


void
Evaluate::compute_block_time(const int Per_u_, const bool evaluate, /*const int block_num, const int size, const bool steady_state,*/ const bool no_derivative)
{
  int var = 0, lag = 0, op;
  unsigned int eq, pos_col;
  ostringstream tmp_out;
  double v1, v2, v3;
  bool go_on = true;
  double ll;
  double rr;
  double *jacob = NULL, *jacob_other_endo = NULL, *jacob_exo = NULL, *jacob_exo_det = NULL;
  EQN_block = block_num;
  stack<double> Stack;
  external_function_type function_type = ExternalFunctionWithoutDerivative;

#ifdef DEBUG
  mexPrintf("compute_block_time\n");
#endif
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
#ifdef OCTAVE_MEX_FILE
  OCTAVE_QUIT;
#else
	if ( utIsInterruptPending() )
		throw UserExceptionHandling();
#endif
  
  while (go_on)
    {
      switch (it_code->first)
        {
        case FNUMEXPR:
#ifdef DEBUG
          mexPrintf("FNUMEXPR\n");
#endif
          it_code_expr = it_code;
          switch (((FNUMEXPR_ *) it_code->second)->get_expression_type())
            {
            case TemporaryTerm:
#ifdef DEBUG
              mexPrintf("TemporaryTerm\n");
#endif
              EQN_type = TemporaryTerm;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
#ifdef DEBUG
              mexPrintf("EQN_equation=%d\n", EQN_equation); mexEvalString("drawnow;");
#endif
              break;
            case ModelEquation:
#ifdef DEBUG
              mexPrintf("ModelEquation\n");
#endif
              EQN_type = ModelEquation;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              break;
            case FirstEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstEndoDerivative\n");
#endif
              EQN_type = FirstEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstOtherEndoDerivative:
#ifdef DEBUG
              mexPrintf("FirstOtherEndoDerivative\n");
#endif
              EQN_type = FirstOtherEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExoDerivative:
#ifdef DEBUG
              mexPrintf("FirstExoDerivative\n");
#endif
              EQN_type = FirstExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstExodetDerivative:
#ifdef DEBUG
              mexPrintf("FirstExodetDerivative\n");
#endif
              EQN_type = FirstExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              break;
            case FirstParamDerivative:
#ifdef DEBUG
              mexPrintf("FirstParamDerivative\n");
#endif
              EQN_type = FirstParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              break;
            case SecondEndoDerivative:
#ifdef DEBUG
              mexPrintf("SecondEndoDerivative\n");
#endif
              EQN_type = SecondEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExoDerivative:
#ifdef DEBUG
              mexPrintf("SecondExoDerivative\n");
#endif
              EQN_type = SecondExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondExodetDerivative:
#ifdef DEBUG
              mexPrintf("SecondExodetDerivative\n");
#endif
              EQN_type = SecondExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              break;
            case SecondParamDerivative:
#ifdef DEBUG
              mexPrintf("SecondParamDerivative\n");
#endif
              EQN_type = SecondParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              break;
            case ThirdEndoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdEndoDerivative\n");
#endif
              EQN_type = ThirdEndoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExoDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExoDerivative\n");
#endif
              EQN_type = ThirdExoDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdExodetDerivative:
#ifdef DEBUG
              mexPrintf("ThirdExodetDerivative\n");
#endif
              EQN_type = ThirdExodetDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_lag1 = ((FNUMEXPR_ *) it_code->second)->get_lag1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_lag2 = ((FNUMEXPR_ *) it_code->second)->get_lag2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              EQN_lag3 = ((FNUMEXPR_ *) it_code->second)->get_lag3();
              break;
            case ThirdParamDerivative:
#ifdef DEBUG
              mexPrintf("ThirdParamDerivative\n");
#endif
              EQN_type = ThirdParamDerivative;
              EQN_equation = ((FNUMEXPR_ *) it_code->second)->get_equation();
              EQN_dvar1 = ((FNUMEXPR_ *) it_code->second)->get_dvariable1();
              EQN_dvar2 = ((FNUMEXPR_ *) it_code->second)->get_dvariable2();
              EQN_dvar3 = ((FNUMEXPR_ *) it_code->second)->get_dvariable3();
              break;
            }
          break;
        case FLDV:
          //load a variable in the processor
          switch (((FLDV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDV Param[var=%d]\n", var);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV y[var=%d, lag=%d, it_=%d], y_size=%d evaluate=%d\n", var, lag, it_, y_size, evaluate);
#endif
              if (evaluate)
                Stack.push(ya[(it_+lag)*y_size+var]);
              else
                Stack.push(y[(it_+lag)*y_size+var]);
#ifdef DEBUG
              tmp_out << " y[" << it_+lag << ", " << var << "](" << y[(it_+lag)*y_size+var] << ")";
#endif
              break;
            case eExogenous:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
#ifdef DEBUG
              mexPrintf("FLDV x[var=%d, lag=%d, it_=%d], nb_row_x=%d evaluate=%d\n", var, lag, it_, nb_row_x, evaluate);
              tmp_out << " x[" << it_+lag << ", " << var << "](" << x[it_+lag+var*nb_row_x] << ")";
#endif
              Stack.push(x[it_+lag+var*nb_row_x]);
              break;
            case eExogenousDet:
              var = ((FLDV_ *) it_code->second)->get_pos();
              lag = ((FLDV_ *) it_code->second)->get_lead_lag();
              Stack.push(x[it_+lag+var*nb_row_xd]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDV: Unknown variable type\n");
            }
          break;
        case FLDSV:
          //load a variable in the processor
          switch (((FLDSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV Param[var=%d]=%f\n", var, params[var]);
              tmp_out << " params[" << var << "](" << params[var] << ")";
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV y[var=%d]=%f\n", var, ya[var]);
              tmp_out << " y[" << var << "](" << y[var] << ")";
#endif
              if (evaluate)
                Stack.push(ya[var]);
              else
                Stack.push(y[var]);
              break;
            case eExogenous:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV x[var=%d]\n", var);
              tmp_out << " x[" << var << "](" << x[var] << ")";
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDSV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDSV xd[var=%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDSV a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDSV: Unknown variable type\n");
            }
          break;
        case FLDVS:
          //load a variable in the processor
          switch (((FLDVS_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("params[%d]\n", var);
#endif
              Stack.push(params[var]);
              break;
            case eEndogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS steady_y[%d]\n", var);
#endif
              Stack.push(steady_y[var]);
              break;
            case eExogenous:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS x[%d] \n", var);
#endif
              Stack.push(x[var]);
              break;
            case eExogenousDet:
              var = ((FLDVS_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FLDVS xd[%d]\n", var);
#endif
              Stack.push(x[var]);
              break;
            case eModelLocalVariable:
#ifdef DEBUG
              mexPrintf("FLDVS a local variable in Block %d Stack.size()=%d", block_num, Stack.size());
              mexPrintf(" value=%f\n", Stack.top());
#endif
              break;
            default:
              mexPrintf("FLDVS: Unknown variable type\n");
            }
          break;
        case FLDT:
          //load a temporary variable in the processor
          var = ((FLDT_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("T[it_=%d var=%d, y_kmin=%d, y_kmax=%d == %d]=>%f\n", it_, var, y_kmin, y_kmax, var*(periods+y_kmin+y_kmax)+it_, var);
          tmp_out << " T[" << it_ << ", " << var << "](" << T[var*(periods+y_kmin+y_kmax)+it_] << ")";
#endif
          Stack.push(T[var*(periods+y_kmin+y_kmax)+it_]);
          break;
        case FLDST:
          //load a temporary variable in the processor
          var = ((FLDST_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDST T[%d]", var);
#endif
          Stack.push(T[var]);
#ifdef DEBUG
          mexPrintf("=%f\n", T[var]);
          tmp_out << " T[" << var << "](" << T[var] << ")";
#endif
          break;
        case FLDU:
          //load u variable in the processor
          var = ((FLDU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FLDU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDSU:
          //load u variable in the processor
          var = ((FLDSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDSU u[%d]\n", var);
          tmp_out << " u[" << var << "](" << u[var] << ")";
#endif
          Stack.push(u[var]);
          break;
        case FLDR:
          //load u variable in the processor
          var = ((FLDR_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("FLDR r[%d]\n", var);
#endif
          Stack.push(r[var]);
          break;
        case FLDZ:
          //load 0 in the processor
#ifdef DEBUG
          mexPrintf("FLDZ\n");
#endif
          Stack.push(0.0);
#ifdef DEBUG
          tmp_out << " 0";
#endif
          break;
        case FLDC:
          //load a numerical constant in the processor
          ll = ((FLDC_ *) it_code->second)->get_value();
#ifdef DEBUG
          mexPrintf("FLDC = %f\n", ll);
          tmp_out << " " << ll;
#endif

          Stack.push(ll);
          break;
        case FSTPV:
          //load a variable in the processor
          switch (((FSTPV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPV_ *) it_code->second)->get_pos();
#ifdef DEBUG
              mexPrintf("FSTPV params[%d]\n", var);
#endif
              params[var] = Stack.top();
              Stack.pop();
              break;
            case eEndogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              y[(it_+lag)*y_size+var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d, %d](%f)=%s\n", it_+lag, var, y[(it_+lag)*y_size+var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case eExogenous:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_x]  = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_x], tmp_out.str().c_str());
              tmp_out.str("");
#endif

              Stack.pop();
              break;
            case eExogenousDet:
              var = ((FSTPV_ *) it_code->second)->get_pos();
              lag = ((FSTPV_ *) it_code->second)->get_lead_lag();
              x[it_+lag+var*nb_row_xd] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[it_+lag+var*nb_row_xd], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPV: Unknown variable type\n");
            }
          break;
        case FSTPSV:
          //load a variable in the processor
          switch (((FSTPSV_ *) it_code->second)->get_type())
            {
            case eParameter:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              params[var] = Stack.top();
              Stack.pop();
              break;
            case eEndogenous:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              y[var] = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" y[%d](%f)=%s\n", var, y[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            case eExogenous:
            case eExogenousDet:
              var = ((FSTPSV_ *) it_code->second)->get_pos();
              x[var]  = Stack.top();
#ifdef DEBUG
              tmp_out << "=>";
              mexPrintf(" x[%d, %d](%f)=%s\n", it_+lag, var, x[var], tmp_out.str().c_str());
              tmp_out.str("");
#endif
              Stack.pop();
              break;
            default:
              mexPrintf("FSTPSV: Unknown variable type\n");
            }
          break;
        case FSTPT:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPT\n");
#endif
          var = ((FSTPT_ *) it_code->second)->get_pos();
          T[var*(periods+y_kmin+y_kmax)+it_] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d, %d](%f)=%s\n", it_, var, T[var*(periods+y_kmin+y_kmax)+it_], tmp_out.str().c_str());
          tmp_out.str("");
#endif

          Stack.pop();
          break;
        case FSTPST:
          //store in a temporary variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPST\n");
#endif
          var = ((FSTPST_ *) it_code->second)->get_pos();
#ifdef DEBUG
          mexPrintf("var=%d\n", var);
#endif
          T[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" T[%d](%f)=%s\n", var, T[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPU:
          //store in u variable from the processor
          var = ((FSTPU_ *) it_code->second)->get_pos();
          var += Per_u_;
#ifdef DEBUG
          mexPrintf("FSTPU\n");
          mexPrintf("var=%d\n", var);
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPSU:
          //store in u variable from the processor
          var = ((FSTPSU_ *) it_code->second)->get_pos();
#ifdef DEBUG
          /*if (var >= u_count_alloc || var < 0)
            mexPrintf("Erreur var=%d\n", var);*/
#endif
          u[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" u[%d](%f)=%s\n", var, u[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPR:
          //store in residual variable from the processor
          var = ((FSTPR_ *) it_code->second)->get_pos();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("FSTPR r[%d]", var);
          tmp_out.str("");
#endif
          r[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf("(%f)=%s\n", r[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;
        case FSTPG:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG\n");
          mexEvalString("drawnow;");
#endif
          var = ((FSTPG_ *) it_code->second)->get_pos();
          g1[var] = Stack.top();
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case FSTPG2:
          //store in the jacobian matrix
          rr = Stack.top();
          if (EQN_type != FirstEndoDerivative)
            {
              ostringstream tmp;
              tmp << " in compute_block_time, impossible case " << EQN_type << " not implement in static jacobian\n";
              throw FatalExceptionHandling(tmp.str());
            }
          eq = ((FSTPG2_ *) it_code->second)->get_row();
          var = ((FSTPG2_ *) it_code->second)->get_col();
#ifdef DEBUG
          mexPrintf("FSTPG2 eq=%d, var=%d\n", eq, var);
          mexEvalString("drawnow;");
#endif
          jacob[eq + size*var] = rr;
          break;
        case FSTPG3:
          //store in derivative (g) variable from the processor
#ifdef DEBUG
          mexPrintf("FSTPG3 Evaluate=%d\n", evaluate);
          mexEvalString("drawnow;");
          if (!evaluate)
            {
              mexPrintf("impossible case!! \n");
              mexEvalString("drawnow;pause;");
            }

#endif
          rr = Stack.top();
          switch (EQN_type)
            {
            case FirstEndoDerivative:
              eq = ((FSTPG3_ *) it_code->second)->get_row();
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Endo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;pause;");
#endif
              jacob[eq + size*pos_col] = rr;
              break;
            case FirstOtherEndoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("other_endo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;pause;");
#endif
              jacob_other_endo[eq + size*pos_col] = rr;
              break;
            case FirstExoDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;pause;");
#endif
              jacob_exo[eq + size*pos_col] = rr;
              break;
            case FirstExodetDerivative:
              //eq = ((FSTPG3_ *) it_code->second)->get_row();
              eq = EQN_equation;
              var = ((FSTPG3_ *) it_code->second)->get_col();
              lag = ((FSTPG3_ *) it_code->second)->get_lag();
              pos_col = ((FSTPG3_ *) it_code->second)->get_col_pos();
#ifdef DEBUG
              mexPrintf("Exo det eq=%d, pos_col=%d, size=%d\n", eq, pos_col, size);
              mexEvalString("drawnow;pause;");
#endif

              jacob_exo_det[eq + size*pos_col] = rr;
              break;
            default:
              ostringstream tmp;
              tmp << " in compute_block_time, variable " << EQN_type << " not used yet\n";
              throw FatalExceptionHandling(tmp.str());
            }
#ifdef DEBUG
          tmp_out << "=>";
          mexPrintf(" g1[%d](%f)=%s\n", var, g1[var], tmp_out.str().c_str());
          tmp_out.str("");
#endif
          Stack.pop();
          break;

        case FBINARY:
          op = ((FBINARY_ *) it_code->second)->get_op_type();
#ifdef DEBUG
          mexPrintf("FBINARY, op=%d\n", op);
#endif
          v2 = Stack.top();
          Stack.pop();
          v1 = Stack.top();
          Stack.pop();
          switch (op)
            {
            case oPlus:
              Stack.push(v1 + v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "+" << v2 << "|";
#endif
              break;
            case oMinus:
              Stack.push(v1 - v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "-" << v2 << "|";
#endif
              break;
            case oTimes:
              Stack.push(v1 * v2);
#ifdef DEBUG
              tmp_out << " |" << v1 << "*" << v2 << "|";
#endif
              break;
            case oDivide:
              double tmp;
#ifdef DEBUG
              mexPrintf("v1=%f / v2=%f\n", v1, v2);
#endif
              try
                {
                  tmp = divide(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |" << v1 << "/" << v2 << "|";
#endif
              break;
            case oLess:
              Stack.push(double (v1 < v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "<" << v2 << "|";
#endif
              break;
            case oGreater:
              Stack.push(double (v1 > v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">" << v2 << "|";
#endif
              break;
            case oLessEqual:
              Stack.push(double (v1 <= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "<=" << v2 << "|";
#endif
              break;
            case oGreaterEqual:
              Stack.push(double (v1 >= v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << ">=" << v2 << "|";
#endif
              break;
            case oEqualEqual:
              Stack.push(double (v1 == v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "==" << v2 << "|";
#endif
              break;
            case oDifferent:
              Stack.push(double (v1 != v2));
#ifdef DEBUG
              tmp_out << " |" << v1 << "!=" << v2 << "|";
#endif
              break;
            case oPower:
#ifdef DEBUG
              mexPrintf("pow\n");
#endif
              try
                {
                  tmp = pow1(v1, v2);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);

#ifdef DEBUG
              tmp_out << " |" << v1 << "^" << v2 << "|";
#endif
              break;
            case oPowerDeriv:
              {
                int derivOrder = int(nearbyint(Stack.top()));
                Stack.pop();
                try
                  {
                    if (fabs(v1) < NEAR_ZERO && v2 > 0
                        && derivOrder > v2
                        && fabs(v2-nearbyint(v2)) < NEAR_ZERO)
                      Stack.push(0.0);
                    else
                      {
                        double dxp = pow1(v1, v2-derivOrder);
                        for (int i = 0; i < derivOrder; i++)
                          dxp *= v2--;
                        Stack.push(dxp);
                      }
                  }
                catch (FloatingPointExceptionHandling &fpeh)
                  {
                    mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                    go_on = false;
                  }
              }

#ifdef DEBUG
              tmp_out << " |PowerDeriv(" << v1 << ", " << v2 << ")|";
#endif
              break;
            case oMax:
              Stack.push(max(v1, v2));
#ifdef DEBUG
              tmp_out << " |max(" << v1 << "," << v2 << ")|";
#endif
              break;
            case oMin:
              Stack.push(min(v1, v2));
#ifdef DEBUG
              tmp_out << " |min(" << v1 << "," << v2 << ")|";
#endif
              break;
            case oEqual:
              // Nothing to do
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown binary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
            }
          break;
        case FUNARY:
          op = ((FUNARY_ *) it_code->second)->get_op_type();
          v1 = Stack.top();
          Stack.pop();
#ifdef DEBUG
          mexPrintf("FUNARY, op=%d\n", op);
#endif
          switch (op)
            {
            case oUminus:
              Stack.push(-v1);
#ifdef DEBUG
              tmp_out << " |-(" << v1 << ")|";
#endif

              break;
            case oExp:
              Stack.push(exp(v1));
#ifdef DEBUG
              tmp_out << " |exp(" << v1 << ")|";
#endif
              break;
            case oLog:
              double tmp;
              try
                {
                  tmp = log1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
              //if (isnan(res1))

#ifdef DEBUG
              tmp_out << " |log(" << v1 << ")|";
#endif
              break;
            case oLog10:
              try
                {
                  tmp = log10_1(v1);
                }
              catch (FloatingPointExceptionHandling &fpeh)
                {
                  mexPrintf("%s      %s\n", fpeh.GetErrorMsg().c_str(), error_location(evaluate, steady_state, size, block_num, it_, Per_u_).c_str());
                  go_on = false;
                }
              Stack.push(tmp);
#ifdef DEBUG
              tmp_out << " |log10(" << v1 << ")|";
#endif
              break;
            case oCos:
              Stack.push(cos(v1));
#ifdef DEBUG
              tmp_out << " |cos(" << v1 << ")|";
#endif
              break;
            case oSin:
              Stack.push(sin(v1));
#ifdef DEBUG
              tmp_out << " |sin(" << v1 << ")|";
#endif
              break;
            case oTan:
              Stack.push(tan(v1));
#ifdef DEBUG
              tmp_out << " |tan(" << v1 << ")|";
#endif
              break;
            case oAcos:
              Stack.push(acos(v1));
#ifdef DEBUG
              tmp_out << " |acos(" << v1 << ")|";
#endif
              break;
            case oAsin:
              Stack.push(asin(v1));
#ifdef DEBUG
              tmp_out << " |asin(" << v1 << ")|";
#endif
              break;
            case oAtan:
              Stack.push(atan(v1));
#ifdef DEBUG
              tmp_out << " |atan(" << v1 << ")|";
#endif
              break;
            case oCosh:
              Stack.push(cosh(v1));
#ifdef DEBUG
              tmp_out << " |cosh(" << v1 << ")|";
#endif
              break;
            case oSinh:
              Stack.push(sinh(v1));
#ifdef DEBUG
              tmp_out << " |sinh(" << v1 << ")|";
#endif
              break;
            case oTanh:
              Stack.push(tanh(v1));
#ifdef DEBUG
              tmp_out << " |tanh(" << v1 << ")|";
#endif
              break;
            case oAcosh:
              Stack.push(acosh(v1));
#ifdef DEBUG
              tmp_out << " |acosh(" << v1 << ")|";
#endif
              break;
            case oAsinh:
              Stack.push(asinh(v1));
#ifdef DEBUG
              tmp_out << " |asinh(" << v1 << ")|";
#endif
              break;
            case oAtanh:
              Stack.push(atanh(v1));
#ifdef DEBUG
              tmp_out << " |atanh(" << v1 << ")|";
#endif
              break;
            case oSqrt:
              Stack.push(sqrt(v1));
#ifdef DEBUG
              tmp_out << " |sqrt(" << v1 << ")|";
#endif
              break;
            case oErf:
              Stack.push(erf(v1));
#ifdef DEBUG
              tmp_out << " |erf(" << v1 << ")|";

#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown unary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
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
          switch (op)
            {
            case oNormcdf:
              Stack.push(0.5*(1+erf((v1-v2)/v3/M_SQRT2)));
#ifdef DEBUG
              tmp_out << " |normcdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            case oNormpdf:
              Stack.push(1/(v3*sqrt(2*M_PI)*exp(pow((v1-v2)/v3, 2)/2)));
#ifdef DEBUG
              tmp_out << " |normpdf(" << v1 << ", " << v2 << ", " << v3 << ")|";
#endif
              break;
            default:
              {
                mexPrintf("Error\n");
                ostringstream tmp;
                tmp << " in compute_block_time, unknown trinary operator " << op << "\n";
                throw FatalExceptionHandling(tmp.str());
              }
            }
          break;

        case FPUSH:
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
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
#ifdef DEBUG
                  mexPrintf("Stack.size()=%d\n", Stack.size());
                  mexEvalString("drawnow;");
#endif
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[nb_input_arguments - i - 1] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    {
                      ostringstream tmp;
                      tmp << " external function: " << function_name << " not found";
                      throw FatalExceptionHandling(tmp.str());
                    }
 
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                  if (function_type == ExternalFunctionWithFirstDerivative || function_type == ExternalFunctionWithFirstandSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD1 = mxGetPr(output_arguments[1]);
                      size_t rows = mxGetN(output_arguments[1]);
                      for (unsigned int i = 0; i < rows; i++)
                        TEFD[make_pair(indx, i)] = FD1[i];
                    }
                  if (function_type == ExternalFunctionWithFirstandSecondDerivative)
                    {
                      unsigned int indx = fc->get_indx();
                      double *FD2 = mxGetPr(output_arguments[2]);
                      size_t rows = mxGetM(output_arguments[2]);
                      size_t cols = mxGetN(output_arguments[2]);
                      unsigned int k = 0;
                      for (unsigned int j = 0; j < cols; j++)
                        for (unsigned int i = 0; i < rows; i++)
                          TEFDD[make_pair(indx, make_pair(i, j))] = FD2[k++];
                    }
                }
                break;
              case ExternalFunctionNumericalFirstDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc((nb_input_arguments+1+nb_add_input_arguments) * sizeof(mxArray *));
                  mxArray *vv = mxCreateString(arg_func_name.c_str());
                  input_arguments[0] = vv;
                  vv = mxCreateDoubleScalar(fc->get_row());
                  input_arguments[1] = vv;
                  vv = mxCreateCellMatrix(1, nb_add_input_arguments);
                  for (unsigned int i = 0; i < nb_add_input_arguments; i++)
                    {
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f Stack.size()=%d\n", i, rr, Stack.size());
#endif
                      mxSetCell(vv, nb_add_input_arguments - (i+1), mxCreateDoubleScalar(rr));
                      Stack.pop();
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
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    {
                      ostringstream tmp;
                      tmp << " external function: " << function_name << " not found";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  double *rr = mxGetPr(output_arguments[0]);
#ifdef DEBUG
                  mexPrintf("*rr=%f\n", *rr);
#endif
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionFirstDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[(nb_input_arguments - 1) - i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    {
                      ostringstream tmp;
                      tmp << " external function: " << function_name << " not found";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  unsigned int indx = fc->get_indx();
                  double *FD1 = mxGetPr(output_arguments[0]);
                  //mexPrint
                  size_t rows = mxGetN(output_arguments[0]);
                  for (unsigned int i = 0; i < rows; i++)
                    TEFD[make_pair(indx, i)] = FD1[i];
                }
                break;
              case ExternalFunctionNumericalSecondDerivative:
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
                      double rr = Stack.top();
#ifdef DEBUG
                      mexPrintf("i=%d rr = %f\n", i, rr);
#endif
                      mxSetCell(vv, (nb_add_input_arguments - 1) - i, mxCreateDoubleScalar(rr));
                      Stack.pop();
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
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    {
                      ostringstream tmp;
                      tmp << " external function: " << function_name << " not found";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  double *rr = mxGetPr(output_arguments[0]);
                  Stack.push(*rr);
                }
                break;
              case ExternalFunctionSecondDerivative:
                {
                  input_arguments = (mxArray **) mxMalloc(nb_input_arguments * sizeof(mxArray *));
                  for (unsigned int i = 0; i < nb_input_arguments; i++)
                    {
                      mxArray *vv = mxCreateDoubleScalar(Stack.top());
                      input_arguments[i] = vv;
                      Stack.pop();
                    }
                  if (mexCallMATLAB(nb_output_arguments, output_arguments, nb_input_arguments, input_arguments, function_name.c_str()))
                    {
                      ostringstream tmp;
                      tmp << " external function: " << function_name << " not found";
                      throw FatalExceptionHandling(tmp.str());
                    }
                  unsigned int indx = fc->get_indx();
                  double *FD2 = mxGetPr(output_arguments[2]);
                  size_t rows = mxGetM(output_arguments[0]);
                  size_t cols = mxGetN(output_arguments[0]);
                  unsigned int k = 0;
                  for (unsigned int j = 0; j < cols; j++)
                    for (unsigned int i = 0; i < rows; i++)
                      TEFDD[make_pair(indx, make_pair(i, j))] = FD2[k++];
                }
                break;
              }
          }
          break;
        case FSTPTEF:
          var = ((FSTPTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FSTPTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
#endif
          TEF[var-1] = Stack.top();
#ifdef DEBUG
          mexPrintf("FSTP TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.pop();
          break;
        case FLDTEF:
          var = ((FLDTEF_ *) it_code->second)->get_number();
#ifdef DEBUG
          mexPrintf("FLDTEF\n");
          mexPrintf("var=%d Stack.size()=%d\n", var, Stack.size());
          mexPrintf("FLD TEF[var-1]=%f done\n", TEF[var-1]);
          mexEvalString("drawnow;");
#endif
          Stack.push(TEF[var-1]);
          break;
        case FSTPTEFD:
          {
            unsigned int indx = ((FSTPTEFD_ *) it_code->second)->get_indx();
            unsigned int row = ((FSTPTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionNumericalFirstDerivative)
              {
                TEFD[make_pair(indx, row-1)] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFD[make_pair(indx, row)]=%f done\n", TEFD[make_pair(indx, row-1)]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
          }

          break;
        case FLDTEFD:
          {
            unsigned int indx = ((FLDTEFD_ *) it_code->second)->get_indx();
            unsigned int row = ((FLDTEFD_ *) it_code->second)->get_row();
#ifdef DEBUG
            mexPrintf("FLDTEFD\n");
            mexPrintf("indx=%d row=%d Stack.size()=%d\n", indx, row, Stack.size());
            mexPrintf("FLD TEFD[make_pair(indx, row)]=%f done\n", TEFD[make_pair(indx, row-1)]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFD[make_pair(indx, row-1)]);
          }
          break;
        case FSTPTEFDD:
          {
            unsigned int indx = ((FSTPTEFDD_ *) it_code->second)->get_indx();
            unsigned int row = ((FSTPTEFDD_ *) it_code->second)->get_row();
            unsigned int col = ((FSTPTEFDD_ *) it_code->second)->get_col();
#ifdef DEBUG
            mexPrintf("FSTPTEFD\n");
            mexPrintf("indx=%d Stack.size()=%d\n", indx, Stack.size());
#endif
            if (function_type == ExternalFunctionNumericalSecondDerivative)
              {
                TEFDD[make_pair(indx, make_pair(row-1, col-1))] = Stack.top();
#ifdef DEBUG
                mexPrintf("FSTP TEFDD[make_pair(indx, make_pair(row, col))]=%f done\n", TEFDD[make_pair(indx, make_pair(row, col))]);
                mexEvalString("drawnow;");
#endif
                Stack.pop();
              }
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
            mexPrintf("FLD TEFD[make_pair(indx, make_pair(row, col))]=%f done\n", TEFDD[make_pair(indx, make_pair(row, col))]);
            mexEvalString("drawnow;");
#endif
            Stack.push(TEFDD[make_pair(indx, make_pair(row-1, col-1))]);
          }
          break;
        case FCUML:
          v1 = Stack.top();
          Stack.pop();
          v2 = Stack.top();
          Stack.pop();
          Stack.push(v1+v2);
          break;
        case FENDBLOCK:
          //it's the block end
#ifdef DEBUG
          mexPrintf("FENDBLOCK\n");
#endif
          go_on = false;
          break;
        case FBEGINBLOCK:
          mexPrintf("Impossible case in Bytecode\n");
          break;
        case FENDEQU:
          if (no_derivative)
            go_on = false;
          break;
        case FJMPIFEVAL:
          if (evaluate)
            {
#ifdef DEBUG
              mexPrintf("FJMPIFEVAL length=%d\n", ((FJMPIFEVAL_ *) it_code->second)->get_pos());
              mexEvalString("drawnow;");
#endif
              it_code += ((FJMPIFEVAL_ *) it_code->second)->get_pos() /* - 1*/;
            }
          break;
        case FJMP:
#ifdef DEBUG
          mexPrintf("FJMP length=%d\n", ((FJMP_ *) it_code->second)->get_pos());
          mexEvalString("drawnow;");
#endif
          it_code += ((FJMP_ *) it_code->second)->get_pos() /*- 1 */;
          break;
        case FOK:
          op = ((FOK_ *) it_code->second)->get_arg();
          if (Stack.size() > 0)
            {
              ostringstream tmp;
              tmp << " in compute_block_time, stack not empty\n";
              throw FatalExceptionHandling(tmp.str());
            }
          break;
        default:
          ostringstream tmp;
          tmp << " in compute_block_time, unknown opcode " << it_code->first << "\n";
          throw FatalExceptionHandling(tmp.str());
        }
      it_code++;
    }
#ifdef DEBUG
  mexPrintf("==> end of compute_block_time Block = %d\n", block_num);
  mexEvalString("drawnow;");
#endif
}



void
Evaluate::evaluate_over_periods(const bool forward)
{
  if (steady_state)
    compute_block_time(0, false, false);
  else
    {
      it_code_type begining = it_code;
      if (forward)
        {
          for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
            {
              it_code = begining;
              compute_block_time(0, false, false);
            }
        }
      else
        {
          for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
            {
              it_code = begining;
              compute_block_time(0, false, false);
            }
        }
    }
}

void
Evaluate::solve_simple_one_periods()
{
  bool cvg = false;
  int iter = 0;
  double ya ;
  double slowc_save = slowc;
  res1 = 0;
  while (!(cvg || (iter > maxit_)))
    {
      it_code = start_code;
      Per_y_ = it_*y_size;
      ya = y[Block_Contain[0].Variable + Per_y_];
      compute_block_time(0, false, false);
      if (!finite(res1))
        {
          res1 = NAN;
          while ((isinf(res1) || isnan(res1)) && (slowc > 1e-9) )
            {
              it_code = start_code;
              compute_block_time(0, false, false);
              if (!finite(res1))
                {
                  slowc /= 1.5;
                  mexPrintf("Reducing the path length in Newton step slowc=%f\n", slowc);
                  y[Block_Contain[0].Variable + Per_y_] = ya - slowc * divide(r[0], g1[0]);
                }
            }
        }
      double rr;
      rr = r[0];
      cvg = (fabs(rr) < solve_tolf);
      //mexPrintf("g1=%x, g1[0]=%f, type=%d, block_num=%d, it_=%d y=%x\n", g1, g1[0], type, block_num, it_, y);
      if (cvg)
        continue;
      try
        {
          y[Block_Contain[0].Variable + Per_y_] += - slowc * divide(rr, g1[0]);
        }
      catch (FloatingPointExceptionHandling &fpeh)
        {
          mexPrintf("%s      \n", fpeh.GetErrorMsg().c_str());
          mexPrintf("      Singularity in block %d", block_num+1);
        }
      iter++;
    }
  slowc = slowc_save;
  if (!cvg)
    {
      ostringstream tmp;
      tmp << " in Solve Forward simple, convergence not achieved in block " << block_num+1 << ", after " << iter << " iterations\n";
      throw FatalExceptionHandling(tmp.str());
    }
}


void
Evaluate::solve_simple_over_periods(const bool forward)
{
  g1 = (double *) mxMalloc(sizeof(double));
  r = (double *) mxMalloc(sizeof(double));
  start_code = it_code;
  if (steady_state)
    {
      it_ = 0;
      solve_simple_one_periods();
    }
  else
    {
      if (forward)
        for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
          solve_simple_one_periods();
      else
        for (it_ = periods+y_kmin-1; it_ >= y_kmin; it_--)
          solve_simple_one_periods();
    }
  mxFree(g1);
  mxFree(r);
}

void
Evaluate::set_block(const int size_arg, const int type_arg, string file_name_arg, string bin_base_name_arg, const int block_num_arg,
          const bool is_linear_arg, const int symbol_table_endo_nbr_arg, const int Block_List_Max_Lag_arg, const int Block_List_Max_Lead_arg, const int u_count_int_arg, const int block_arg)
{
  size = size_arg;
  type = type_arg;
  file_name = file_name_arg;
  bin_base_name = bin_base_name_arg;
  block_num = block_num_arg;
  is_linear = is_linear_arg;
  symbol_table_endo_nbr = symbol_table_endo_nbr_arg;
  Block_List_Max_Lag = Block_List_Max_Lag_arg;
  Block_List_Max_Lead = Block_List_Max_Lead_arg;
  u_count_int = u_count_int_arg;
  block = block_arg;
}

void
Evaluate::evaluate_complete(const bool no_derivatives)
{
  it_code = start_code;
  compute_block_time(0, false, no_derivatives);
}


void
Evaluate::compute_complete_2b(const bool no_derivatives, double *_res1, double *_res2, double *_max_res, int *_max_res_idx)
{
  bool result;
  res1 = 0;
  *_res1 = 0;
  *_res2 = 0;
  *_max_res = 0;
  for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
    {
      Per_u_ = (it_-y_kmin)*u_count_int;
      Per_y_ = it_*y_size;
      it_code = start_code;
      int shift = (it_-y_kmin) * size;
      compute_block_time(Per_u_, false, no_derivatives);
      if (!(isnan(res1) || isinf(res1)))
        {
#ifdef USE_OMP
          if (size > 10)
            {
              double res1_ = 0;
              double res2_ = 0;
              double max_res_ = 0;
              int max_res_idx_ = 0;
              #ifdef USE_OMP
              #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:res1_, res2_) shared(max_res_, max_res_idx_)
              #endif
              for (int i = 0; i < size; i++)
                {
                  double rr = r[i];
                  res[i+shift] = rr;
                  if (max_res_ < fabs(rr))
                    {
                      max_res_ = fabs(rr);
                      max_res_idx_ = i;
                    }
                  res2_ += rr*rr;
                  res1_ += fabs(rr);
                }
              *_res1 += res1_;
              *_res2 += res2_;
              if (max_res_ > *_max_res)
                {
                  *_max_res = max_res_;
                  *_max_res_idx = max_res_idx_;
                }
            }
          else
#endif
            {
              for (int i = 0; i < size; i++)
                {
                  double rr;
                  rr = r[i];
                  res[i+shift] = rr;
                  if (max_res < fabs(rr))
                    {
                      *_max_res = fabs(rr);
                      *_max_res_idx = i;
                    }
                  *_res2 += rr*rr;
                  *_res1 += fabs(rr);
                }
            }
        }
      else
        return;
    }
  return;
}


bool
Evaluate::compute_complete(const bool no_derivatives, double &_res1, double &_res2, double &_max_res, int &_max_res_idx)
{
  bool result;
  res1 = 0;
  it_code = start_code;
  compute_block_time(0, false, no_derivatives);
  if (!(isnan(res1) || isinf(res1)))
    {
#ifdef USE_OMP
      if (size > 10)
        {
          double res1_ = 0;
          double res2_ = 0;
          double max_res_ = 0;
          int max_res_idx_ = 0;
          #ifdef USE_OMP
          #pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS"))) reduction(+:res1_, res2_) shared(max_res_, max_res_idx_)
          #endif
          for (int i = 0; i < size; i++)
            {
              double rr = r[i];
              if (max_res_ < fabs(rr))
                {
                  #ifdef USE_OMP
                  #pragma omp critical
                  #endif
                    {
                      max_res_ = fabs(rr);
                      max_res_idx_ = i;
                    }
                }
              res2_ += rr*rr;
              res1_ += fabs(rr);
            }
          _res1 = res1_;
          _res2 = res2_;
          _max_res = max_res_;
          _max_res_idx = max_res_idx_;
        }
      else
#endif
        {
          _res1 = 0;
          _res2 = 0;
          _max_res = 0;
          for (int i = 0; i < size; i++)
            {
              double rr;
              rr = r[i];
              if (max_res < fabs(rr))
                {
                  _max_res = fabs(rr);
                  _max_res_idx = i;
                }
              _res2 += rr*rr;
              _res1 += fabs(rr);
            }
        }
      result = true;
    }
  else
    result = false;
  return result;
}


bool
Evaluate::compute_complete(double lambda, double *crit)
{
  double res1_ = 0, res2_ = 0, max_res_ = 0;
  //double res1 = 0, res2, max_res;
  int max_res_idx_ = 0;
  if (steady_state)
    {
      it_ = 0;
#ifdef USE_OMP
#pragma omp parallel for num_threads(atoi(getenv("DYNARE_NUM_THREADS")))
#endif
      for (int i = 0; i < size; i++)
        {
          int eq = index_vara[i];
          y[eq] = ya[eq] + lambda * direction[eq];
        }
      Per_u_ = 0;
      Per_y_ = 0;
      if (compute_complete(true, res1, res2, max_res, max_res_idx))
        {
          res2_ = res2;
          /*res1_ = res1;
          if (max_res > max_res_)
            {
              max_res = max_res_;
              max_res_idx = max_res_idx_;
            }*/
        }
      else
        return false;
    }
  else
    {
      for (int it = y_kmin; it < periods+y_kmin; it++)
        {
          for (int i = 0; i < size; i++)
            {
              int eq = index_vara[i];
              y[eq+it*y_size] = ya[eq+it*y_size] + lambda * direction[eq+it*y_size];
            }
        }
      for (it_ = y_kmin; it_ < periods+y_kmin; it_++)
        {
          Per_u_ = (it_-y_kmin)*u_count_int;
          Per_y_ = it_*y_size;
          if (compute_complete(true, res1, res2, max_res, max_res_idx))
            {
              res2_ += res2;
              res1_ += res1;
              if (max_res > max_res_)
                {
                  max_res = max_res_;
                  max_res_idx = max_res_idx_;
                }
            }
          else
            return false;
        }
    }
    mexPrintf("  lambda=%e, res2=%e\n", lambda, res2_);
  *crit = res2_/2;
  return true;
}
