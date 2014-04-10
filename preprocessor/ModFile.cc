/*
 * Copyright (C) 2006-2014 Dynare Team
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

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <typeinfo>
#include <cassert>
#ifndef _WIN32
# include <unistd.h>
#endif

#include "ModFile.hh"
#include "ConfigFile.hh"
#include "ComputingTasks.hh"

ModFile::ModFile(WarningConsolidation &warnings_arg)
  : expressions_tree(symbol_table, num_constants, external_functions_table),
    dynamic_model(symbol_table, num_constants, external_functions_table),
    trend_dynamic_model(symbol_table, num_constants, external_functions_table),
    ramsey_FOC_equations_dynamic_model(symbol_table, num_constants, external_functions_table),
    static_model(symbol_table, num_constants, external_functions_table),
    steady_state_model(symbol_table, num_constants, external_functions_table, static_model),
    linear(false), block(false), byte_code(false), use_dll(false), no_static(false), 
    differentiate_forward_vars(false),
    nonstationary_variables(false), ramsey_model_orig_eqn_nbr(0),
    warnings(warnings_arg)
{
}

ModFile::~ModFile()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    delete (*it);
}

void
ModFile::evalAllExpressions(bool warn_uninit)
{
  cout << "Evaluating expressions...";

  // Loop over all statements, and fill global eval context if relevant
  for (vector<Statement *>::const_iterator it = statements.begin(); it != statements.end(); it++)
    {
      InitParamStatement *ips = dynamic_cast<InitParamStatement *>(*it);
      if (ips)
        ips->fillEvalContext(global_eval_context);

      InitOrEndValStatement *ies = dynamic_cast<InitOrEndValStatement *>(*it);
      if (ies)
        ies->fillEvalContext(global_eval_context);

      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass)
        lpass->fillEvalContext(global_eval_context);
    }

  // Evaluate model local variables
  dynamic_model.fillEvalContext(global_eval_context);

  cout << "done" << endl;

  // Check if some symbols are not initialized, and give them a zero value then
  for (int id = 0; id <= symbol_table.maxID(); id++)
    {
      SymbolType type = symbol_table.getType(id);
      if ((type == eEndogenous || type == eExogenous || type == eExogenousDet
           || type == eParameter || type == eModelLocalVariable)
          && global_eval_context.find(id) == global_eval_context.end())
        {
          if (warn_uninit)
            warnings << "WARNING: Can't find a numeric initial value for "
                     << symbol_table.getName(id) << ", using zero" << endl;
          global_eval_context[id] = 0;
        }
    }
}

void
ModFile::addStatement(Statement *st)
{
  statements.push_back(st);
}

void
ModFile::addStatementAtFront(Statement *st)
{
  statements.insert(statements.begin(), st);
}

void
ModFile::checkPass()
{
  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->checkPass(mod_file_struct, warnings);

  // Check the steady state block
  steady_state_model.checkPass(mod_file_struct.ramsey_model_present, warnings);

  // If order option has not been set, default to 2
  if (!mod_file_struct.order_option)
    mod_file_struct.order_option = 2;

  bool stochastic_statement_present = mod_file_struct.stoch_simul_present
    || mod_file_struct.estimation_present
    || mod_file_struct.osr_present
    || mod_file_struct.ramsey_policy_present
    || mod_file_struct.discretionary_policy_present
    || mod_file_struct.calib_smoother_present;

  // Allow empty model only when doing a standalone BVAR estimation
  if (dynamic_model.equation_number() == 0
      && (mod_file_struct.check_present
          || mod_file_struct.perfect_foresight_solver_present
          || stochastic_statement_present))
    {
      cerr << "ERROR: At least one model equation must be declared!" << endl;
      exit(EXIT_FAILURE);
    }

  if (((mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present) 
       && !mod_file_struct.planner_objective_present)
      || (!(mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present)
	  && mod_file_struct.planner_objective_present))
    {
      cerr << "ERROR: A planner_objective statement must be used with a ramsey_model, a ramsey_policy or a discretionary_policy statement and vice versa." << endl;
      exit(EXIT_FAILURE);
    }

  if ((mod_file_struct.osr_present && (!mod_file_struct.osr_params_present || !mod_file_struct.optim_weights_present))
      || ((!mod_file_struct.osr_present || !mod_file_struct.osr_params_present) && mod_file_struct.optim_weights_present)
      || ((!mod_file_struct.osr_present || !mod_file_struct.optim_weights_present) && mod_file_struct.osr_params_present))
    {
      cerr << "ERROR: The osr statement must be used with osr_params and optim_weights." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.perfect_foresight_solver_present && stochastic_statement_present)
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver,simul} and one of {stoch_simul, estimation, osr, ramsey_policy, discretionary_policy}. This is not possible: one cannot mix perfect foresight context with stochastic context in the same file." << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.k_order_solver && byte_code)
    {
      cerr << "ERROR: 'k_order_solver' (which is implicit if order >= 3), is not yet compatible with 'bytecode'." << endl;
      exit(EXIT_FAILURE);
    }

  if (use_dll && (block || byte_code))
    {
      cerr << "ERROR: In 'model' block, 'use_dll' option is not compatible with 'block' or 'bytecode'" << endl;
      exit(EXIT_FAILURE);
    }

  if (block || byte_code)
    if (dynamic_model.isModelLocalVariableUsed())
      {
        cerr << "ERROR: In 'model' block, 'block' or 'bytecode' options are not yet compatible with pound expressions" << endl;
        exit(EXIT_FAILURE);
      }

  if ((stochastic_statement_present || mod_file_struct.check_present || mod_file_struct.steady_present) && no_static)
    {
      cerr << "ERROR: no_static option is incompatible with stoch_simul, estimation, osr, ramsey_policy, discretionary_policy, steady and check commands" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.dsge_var_estimated)
    if (!mod_file_struct.dsge_prior_weight_in_estimated_params)
      {
        cerr << "ERROR: When estimating a DSGE-VAR model and estimating the weight of the prior, dsge_prior_weight must "
             << "be referenced in the estimated_params block." << endl;
        exit(EXIT_FAILURE);
      }

  if (symbol_table.exists("dsge_prior_weight"))
    {
      if (symbol_table.getType("dsge_prior_weight") != eParameter)
        {
          cerr << "ERROR: dsge_prior_weight may only be used as a parameter." << endl;
          exit(EXIT_FAILURE);
        }
      else
        warnings << "WARNING: When estimating a DSGE-Var, declaring dsge_prior_weight as a "
                 << "parameter is deprecated. The preferred method is to do this via "
                 << "the dsge_var option in the estimation statement." << endl;

      if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
        {
          cerr << "ERROR: dsge_prior_weight can either be declared as a parameter (deprecated) or via the dsge_var option "
               << "to the estimation statement (preferred), but not both." << endl;
          exit(EXIT_FAILURE);
        }

      if (!mod_file_struct.dsge_prior_weight_initialized && !mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: If dsge_prior_weight is declared as a parameter, it must either be initialized or placed in the "
               << "estimated_params block." << endl;
          exit(EXIT_FAILURE);
        }

      if (mod_file_struct.dsge_prior_weight_initialized && mod_file_struct.dsge_prior_weight_in_estimated_params)
        {
          cerr << "ERROR: dsge_prior_weight cannot be both initalized and estimated." << endl;
          exit(EXIT_FAILURE);
        }
    }

  if (mod_file_struct.dsge_prior_weight_in_estimated_params)
    if (!mod_file_struct.dsge_var_estimated && !mod_file_struct.dsge_var_calibrated.empty())
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, the prior weight cannot be calibrated "
             << "via the dsge_var option in the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }
    else if (!mod_file_struct.dsge_var_estimated && !symbol_table.exists("dsge_prior_weight"))
      {
        cerr << "ERROR: If dsge_prior_weight is in the estimated_params block, it must either be declared as a parameter "
             << "(deprecated) or the dsge_var option must be passed to the estimation statement (preferred)." << endl;
        exit(EXIT_FAILURE);
      }

  if (dynamic_model.staticOnlyEquationsNbr() != dynamic_model.dynamicOnlyEquationsNbr())
    {
      cerr << "ERROR: the number of equations marked [static] must be equal to the number of equations marked [dynamic]" << endl;
      exit(EXIT_FAILURE);
    }
  
  if (dynamic_model.staticOnlyEquationsNbr() > 0 &&
      (mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present))
    {
      cerr << "ERROR: marking equations as [static] or [dynamic] is not possible with ramsey_model, ramsey_policy or discretionary_policy" << endl;
      exit(EXIT_FAILURE);
    }

  if (stochastic_statement_present &&
      (dynamic_model.isUnaryOpUsed(oSign)
       || dynamic_model.isUnaryOpUsed(oAbs)
       || dynamic_model.isBinaryOpUsed(oMax)
       || dynamic_model.isBinaryOpUsed(oMin)
       || dynamic_model.isBinaryOpUsed(oGreater)
       || dynamic_model.isBinaryOpUsed(oLess)
       || dynamic_model.isBinaryOpUsed(oGreaterEqual)
       || dynamic_model.isBinaryOpUsed(oLessEqual)
       || dynamic_model.isBinaryOpUsed(oEqualEqual)
       || dynamic_model.isBinaryOpUsed(oDifferent)))
    warnings << "WARNING: you are using a function (max, min, abs, sign) or an operator (<, >, <=, >=, ==, !=) which is unsuitable for a stochastic context; see the reference manual, section about \"Expressions\", for more details." << endl;

  // Test if some estimated parameters are used within the values of shocks
  // statements (see issue #469)
  set<int> parameters_intersect;
  set_intersection(mod_file_struct.parameters_within_shocks_values.begin(),
                   mod_file_struct.parameters_within_shocks_values.end(),
                   mod_file_struct.estimated_parameters.begin(),
                   mod_file_struct.estimated_parameters.end(),
                   inserter(parameters_intersect, parameters_intersect.begin()));
  if (parameters_intersect.size() > 0)
    {
      cerr << "ERROR: some estimated parameters (";
      for (set<int>::const_iterator it = parameters_intersect.begin();
           it != parameters_intersect.end(); )
        {
          cerr << symbol_table.getName(*it);
          if (++it != parameters_intersect.end())
            cerr << ", ";
        }
      cerr << ") also appear in the expressions defining the variance/covariance matrix of shocks; this is not allowed." << endl;
      exit(EXIT_FAILURE);
    }

  // Check if some exogenous is not used in the model block
  set<int> unusedExo = dynamic_model.findUnusedExogenous();
  if (unusedExo.size() > 1)
    {
      warnings << "WARNING: some exogenous (";
      for (set<int>::const_iterator it = unusedExo.begin();
           it != unusedExo.end(); )
        {
          warnings << symbol_table.getName(*it);
          if (++it != unusedExo.end())
            warnings << ", ";
        }
      warnings << ") are declared but not used in the model. This may lead to crashes or unexpected behaviour." << endl;
    }
}

void
ModFile::transformPass(bool nostrict)
{
  if (nostrict)
    {
      set<int> unusedEndogs = dynamic_model.findUnusedEndogenous();
      for (set<int>::iterator it = unusedEndogs.begin(); it != unusedEndogs.end(); it++)
        {
          symbol_table.changeType(*it, eUnusedEndogenous);
          warnings << "WARNING: '" << symbol_table.getName(*it)
                   << "' not used in model block, removed by nostrict command-line option" << endl;
        }
    }

  if (symbol_table.predeterminedNbr() > 0)
    dynamic_model.transformPredeterminedVariables();

  // Create auxiliary vars for Expectation operator
  dynamic_model.substituteExpectation(mod_file_struct.partial_information);

  if (nonstationary_variables)
    {
      dynamic_model.detrendEquations();
      dynamic_model.cloneDynamic(trend_dynamic_model);
      dynamic_model.removeTrendVariableFromEquations();
    }

  if (mod_file_struct.ramsey_model_present)
    {
      StaticModel *planner_objective = NULL;
      for (vector<Statement *>::iterator it = statements.begin(); it != statements.end(); it++)
        {
          PlannerObjectiveStatement *pos = dynamic_cast<PlannerObjectiveStatement *>(*it);
          if (pos != NULL)
            planner_objective = pos->getPlannerObjective();
        }
      assert(planner_objective != NULL);
      ramsey_model_orig_eqn_nbr = dynamic_model.equation_number();

      /*
        clone the model then clone the new equations back to the original because
        we have to call computeDerivIDs (in computeRamseyPolicyFOCs and computingPass)
       */
      dynamic_model.cloneDynamic(ramsey_FOC_equations_dynamic_model);
      ramsey_FOC_equations_dynamic_model.computeRamseyPolicyFOCs(*planner_objective);
      ramsey_FOC_equations_dynamic_model.replaceMyEquations(dynamic_model);
    }

  if (mod_file_struct.stoch_simul_present
      || mod_file_struct.estimation_present
      || mod_file_struct.osr_present
      || mod_file_struct.ramsey_policy_present
      || mod_file_struct.discretionary_policy_present
      || mod_file_struct.calib_smoother_present)
    {
      // In stochastic models, create auxiliary vars for leads and lags greater than 2, on both endos and exos
      dynamic_model.substituteEndoLeadGreaterThanTwo(false);
      dynamic_model.substituteExoLead(false);
      dynamic_model.substituteEndoLagGreaterThanTwo(false);
      dynamic_model.substituteExoLag(false);
    }
  else
    {
      // In deterministic models, create auxiliary vars for leads and lags endogenous greater than 2, only on endos (useless on exos)
      dynamic_model.substituteEndoLeadGreaterThanTwo(true);
      dynamic_model.substituteEndoLagGreaterThanTwo(true);
    }

  if (differentiate_forward_vars)
    dynamic_model.differentiateForwardVars(differentiate_forward_vars_subset);

  if (mod_file_struct.dsge_var_estimated || !mod_file_struct.dsge_var_calibrated.empty())
    try
      {
        int sid = symbol_table.addSymbol("dsge_prior_weight", eParameter);
        if (!mod_file_struct.dsge_var_calibrated.empty())
          addStatementAtFront(new InitParamStatement(sid,
                                                     expressions_tree.AddNonNegativeConstant(mod_file_struct.dsge_var_calibrated),
                                                     symbol_table));
      }
    catch (SymbolTable::AlreadyDeclaredException &e)
      {
        cerr << "ERROR: dsge_prior_weight should not be declared as a model variable / parameter "
             << "when the dsge_var option is passed to the estimation statement." << endl;
        exit(EXIT_FAILURE);
      }

  // Freeze the symbol table
  symbol_table.freeze();

  /*
    Enforce the same number of equations and endogenous, except in three cases:
    - ramsey_model, ramsey_policy or discretionary_policy is used
    - a BVAR command is used and there is no equation (standalone BVAR estimation)
    - nostrict option is passed and there are more endogs than equations (dealt with before freeze)
  */
  if (!(mod_file_struct.ramsey_model_present || mod_file_struct.discretionary_policy_present)
      && !(mod_file_struct.bvar_present && dynamic_model.equation_number() == 0)
      && (dynamic_model.equation_number() != symbol_table.endo_nbr()))
    {
      cerr << "ERROR: There are " << dynamic_model.equation_number() << " equations but " << symbol_table.endo_nbr() << " endogenous variables!" << endl;
      exit(EXIT_FAILURE);
    }

  if (symbol_table.exo_det_nbr() > 0 && mod_file_struct.perfect_foresight_solver_present)
    {
      cerr << "ERROR: A .mod file cannot contain both one of {perfect_foresight_solver,simul}  and varexo_det declaration (all exogenous variables are deterministic in this case)" << endl;
      exit(EXIT_FAILURE);
    }

  if (mod_file_struct.ramsey_policy_present && symbol_table.exo_det_nbr() > 0)
    {
      cerr << "ERROR: ramsey_policy is incompatible with deterministic exogenous variables" << endl;
      exit(EXIT_FAILURE);
    }

  if (!mod_file_struct.ramsey_model_present)
    cout << "Found " << dynamic_model.equation_number() << " equation(s)." << endl;
  else
    {
      cout << "Found " << ramsey_model_orig_eqn_nbr  << " equation(s)." << endl;
      cout << "Found " << dynamic_model.equation_number() << " FOC equation(s) for Ramsey Problem." << endl;
    }

  if (symbol_table.exists("dsge_prior_weight"))
    if (mod_file_struct.bayesian_irf_present)
      {
        if (symbol_table.exo_nbr() != symbol_table.observedVariablesNbr())
          {
            cerr << "ERROR: When estimating a DSGE-Var and the bayesian_irf option is passed to the estimation "
                 << "statement, the number of shocks must equal the number of observed variables." << endl;
            exit(EXIT_FAILURE);
          }
      }
    else
      if (symbol_table.exo_nbr() < symbol_table.observedVariablesNbr())
        {
          cerr << "ERROR: When estimating a DSGE-Var, the number of shocks must be "
               << "greater than or equal to the number of observed variables." << endl;
          exit(EXIT_FAILURE);
        }
}

void
ModFile::computingPass(bool no_tmp_terms, FileOutputType output)
{
  // Mod file may have no equation (for example in a standalone BVAR estimation)
  if (dynamic_model.equation_number() > 0)
    {
      if (nonstationary_variables)
        trend_dynamic_model.runTrendTest(global_eval_context);

      // Compute static model and its derivatives
      dynamic_model.toStatic(static_model);
      if (!no_static)
        {
          if (mod_file_struct.stoch_simul_present
              || mod_file_struct.estimation_present || mod_file_struct.osr_present
              || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
              || mod_file_struct.calib_smoother_present)
            static_model.set_cutoff_to_zero();

          const bool static_hessian = mod_file_struct.identification_present
            || mod_file_struct.estimation_analytic_derivation;
          const bool paramsDerivatives = mod_file_struct.identification_present
            || mod_file_struct.estimation_analytic_derivation;
          static_model.computingPass(global_eval_context, no_tmp_terms, static_hessian,
                                     false, paramsDerivatives, block, byte_code);
        }
      // Set things to compute for dynamic model
      if (mod_file_struct.perfect_foresight_solver_present || mod_file_struct.check_present
          || mod_file_struct.stoch_simul_present
          || mod_file_struct.estimation_present || mod_file_struct.osr_present
          || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
          || mod_file_struct.calib_smoother_present)
        {
          if (mod_file_struct.perfect_foresight_solver_present)
            dynamic_model.computingPass(true, false, false, false, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
          else
            {
              if (mod_file_struct.stoch_simul_present
                  || mod_file_struct.estimation_present || mod_file_struct.osr_present
                  || mod_file_struct.ramsey_model_present || mod_file_struct.identification_present
                  || mod_file_struct.calib_smoother_present)
                dynamic_model.set_cutoff_to_zero();
              if (mod_file_struct.order_option < 1 || mod_file_struct.order_option > 3)
                {
                  cerr << "ERROR: Incorrect order option..." << endl;
                  exit(EXIT_FAILURE);
                }
              bool hessian = mod_file_struct.order_option >= 2 
		|| mod_file_struct.identification_present 
		|| mod_file_struct.estimation_analytic_derivation
		|| output == second 
		|| output == third;
              bool thirdDerivatives = mod_file_struct.order_option == 3 
		|| mod_file_struct.estimation_analytic_derivation
		|| output == third;
              bool paramsDerivatives = mod_file_struct.identification_present || mod_file_struct.estimation_analytic_derivation;
              dynamic_model.computingPass(true, hessian, thirdDerivatives, paramsDerivatives, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
            }
        }
      else // No computing task requested, compute derivatives up to 2nd order by default
        dynamic_model.computingPass(true, true, false, false, global_eval_context, no_tmp_terms, block, use_dll, byte_code);
    }

  for (vector<Statement *>::iterator it = statements.begin();
       it != statements.end(); it++)
    (*it)->computingPass();
}

void
ModFile::writeOutputFiles(const string &basename, bool clear_all, bool no_log, bool no_warn, bool console, bool nograph, bool nointeractive, const ConfigFile &config_file
#if defined(_WIN32) || defined(__CYGWIN32__)
                          , bool cygwin, bool msvc
#endif
                          ) const
{
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : main Dynare file " << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl;

  if (no_warn)
    mOutputFile << "warning off" << endl; // This will be executed *after* function warning_config()

  if (clear_all)
    mOutputFile << "clear all" << endl;

  mOutputFile << "tic;" << endl
              << "global M_ oo_ options_ ys0_ ex0_ estimation_info" << endl
              << "options_ = [];" << endl
              << "M_.fname = '" << basename << "';" << endl
              << "%" << endl
              << "% Some global variables initialization" << endl
              << "%" << endl;
  config_file.writeHooks(mOutputFile);
  mOutputFile << "global_initialization;" << endl
              << "diary off;" << endl;
  if (!no_log)
    mOutputFile << "diary('" << basename << ".log');" << endl;

  if (console)
    mOutputFile << "options_.console_mode = 1;" << endl
                << "options_.nodisplay = 1;" << endl;
  if (nograph)
    mOutputFile << "options_.nograph = 1;" << endl;

  if (nointeractive)
    mOutputFile << "options_.nointeractive = 1;" << endl;
    
  cout << "Processing outputs ...";

  symbol_table.writeOutput(mOutputFile);

  // Initialize M_.Sigma_e, M_.Correlation_matrix, M_.H, and M_.Correlation_matrix_ME
  mOutputFile << "M_.Sigma_e = zeros(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl
              << "M_.Correlation_matrix = eye(" << symbol_table.exo_nbr() << ", "
              << symbol_table.exo_nbr() << ");" << endl;

  if (mod_file_struct.calibrated_measurement_errors)
    mOutputFile << "M_.H = zeros(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl
                << "M_.Correlation_matrix_ME = eye(" << symbol_table.observedVariablesNbr() << ", "
                << symbol_table.observedVariablesNbr() << ");" << endl;
  else
    mOutputFile << "M_.H = 0;" << endl
                << "M_.Correlation_matrix_ME = 1;" << endl;

  // May be later modified by a shocks block
  mOutputFile << "M_.sigma_e_is_diagonal = 1;" << endl;

  // Initialize M_.det_shocks
  mOutputFile << "M_.det_shocks = [];" << endl;

  if (linear == 1)
    mOutputFile << "options_.linear = 1;" << endl;

  mOutputFile << "options_.block=" << block << ";" << endl
              << "options_.bytecode=" << byte_code << ";" << endl
              << "options_.use_dll=" << use_dll << ";" << endl;

  if (parallel_local_files.size() > 0)
    {
      mOutputFile << "options_.parallel_info.local_files = {" << endl;
      for (size_t i = 0; i < parallel_local_files.size(); i++)
        {
          size_t j = parallel_local_files[i].find_last_of("/\\");
          if (j == string::npos)
            mOutputFile << "'', '" << parallel_local_files[i] << "';" << endl;
          else
            mOutputFile << "'" << parallel_local_files[i].substr(0, j+1) << "', '"
                        << parallel_local_files[i].substr(j+1, string::npos) << "';" << endl;
        }
      mOutputFile << "};" << endl;
    }

  config_file.writeCluster(mOutputFile);

  if (byte_code)
    mOutputFile << "if exist('bytecode') ~= 3" << endl
                << "  error('DYNARE: Can''t find bytecode DLL. Please compile it or remove the ''bytecode'' option.')" << endl
                << "end" << endl;

  // Erase possible remnants of previous runs
  unlink((basename + "_dynamic.m").c_str());
  unlink((basename + "_dynamic.cod").c_str());
  unlink((basename + "_dynamic.bin").c_str());

  unlink((basename + "_static.m").c_str());
  unlink((basename + "_static.cod").c_str());
  unlink((basename + "_static.bin").c_str());

  unlink((basename + "_steadystate2.m").c_str());
  unlink((basename + "_set_auxiliary_variables.m").c_str());

  if (!use_dll)
    {
      mOutputFile << "erase_compiled_function('" + basename + "_static');" << endl;
      mOutputFile << "erase_compiled_function('" + basename + "_dynamic');" << endl;
    }

#if defined(_WIN32) || defined(__CYGWIN32__)
  // If using USE_DLL with MSVC, check that the user didn't use a function not supported by MSVC (because MSVC doesn't comply with C99 standard)
  if (use_dll && msvc)
    {
      if (dynamic_model.isUnaryOpUsed(oAcosh))
        {
          cerr << "ERROR: acosh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAsinh))
        {
          cerr << "ERROR: asinh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isUnaryOpUsed(oAtanh))
        {
          cerr << "ERROR: atanh() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
      if (dynamic_model.isTrinaryOpUsed(oNormcdf))
        {
          cerr << "ERROR: normcdf() function is not supported with USE_DLL option and MSVC compiler; use Cygwin compiler instead." << endl;
          exit(EXIT_FAILURE);
        }
    }
#endif

  // Compile the dynamic MEX file for use_dll option
  if (use_dll)
    {
      mOutputFile << "if ~exist('OCTAVE_VERSION')" << endl;
      // Some mex commands are enclosed in an eval(), because otherwise it will make Octave fail
#if defined(_WIN32) || defined(__CYGWIN32__)
      if (msvc)
        // MATLAB/Windows + Microsoft Visual C++
        mOutputFile << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Dynamic\" " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                    << "    eval('mex -O LINKFLAGS=\"$LINKFLAGS /export:Static\" " << basename << "_static.c "<< basename << "_static_mex.c')" << endl;
      else if (cygwin)
        // MATLAB/Windows + Cygwin g++
        mOutputFile << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Dynamic >> mex.def\" " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                    << "    eval('mex -O PRELINK_CMDS1=\"echo EXPORTS > mex.def & echo mexFunction >> mex.def & echo Static >> mex.def\" " << basename << "_static.c "<< basename << "_static_mex.c')" << endl;
      else
        mOutputFile << "    error('When using the USE_DLL option, you must give either ''cygwin'' or ''msvc'' option to the ''dynare'' command')" << endl;
#else
# ifdef __linux__
      // MATLAB/Linux
      mOutputFile << "    if matlab_ver_less_than('8.3')" << endl
                  << "        eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "        eval('mex -O LDFLAGS=''-pthread -shared -Wl,--no-undefined'' " << basename << "_static.c "<< basename << "_static_mex.c')" << endl
                  << "    else" << endl
                  << "        eval('mex -O LINKEXPORT='''' " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "        eval('mex -O LINKEXPORT='''' " << basename << "_static.c "<< basename << "_static_mex.c')" << endl
                  << "    end" << endl;
# else // MacOS
      // MATLAB/MacOS
      mOutputFile << "    if matlab_ver_less_than('8.3')" << endl
                  << "        if matlab_ver_less_than('8.1')" << endl
                  << "            eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "            eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_static.c " << basename << "_static_mex.c')" << endl
                  << "        else" << endl
                  << "            eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$MW_SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "            eval('mex -O LDFLAGS=''-Wl,-twolevel_namespace -undefined error -arch \\$ARCHS -Wl,-syslibroot,\\$MW_SDKROOT -mmacosx-version-min=\\$MACOSX_DEPLOYMENT_TARGET -bundle'' "
                  << basename << "_static.c " << basename << "_static_mex.c')" << endl
                  << "        end" << endl
                  << "    else" << endl
                  << "        eval('mex -O LINKEXPORT='''' " << basename << "_dynamic.c " << basename << "_dynamic_mex.c')" << endl
                  << "        eval('mex -O LINKEXPORT='''' " << basename << "_static.c "<< basename << "_static_mex.c')" << endl
                  << "    end" << endl;
# endif
#endif
      mOutputFile << "else" << endl // Octave
                  << "    mex " << basename << "_dynamic.c " << basename << "_dynamic_mex.c" << endl
                  << "    mex " << basename << "_static.c " << basename << "_static_mex.c" << endl
                  << "end" << endl;
    }

  // Add path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "addpath " << basename << ";" << endl;

  if (mod_file_struct.ramsey_model_present)
    mOutputFile << "M_.orig_eq_nbr = " << ramsey_model_orig_eqn_nbr << ";" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeOutput(mOutputFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      if (!no_static)
        static_model.writeOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
    {
      (*it)->writeOutput(mOutputFile, basename);

      /* Special treatment for initval block: insert initial values for the
         auxiliary variables and initialize exo det */
      InitValStatement *ivs = dynamic_cast<InitValStatement *>(*it);
      if (ivs != NULL)
        {
          static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
          ivs->writeOutputPostInit(mOutputFile);
        }

      // Special treatment for endval block: insert initial values for the auxiliary variables
      EndValStatement *evs = dynamic_cast<EndValStatement *>(*it);
      if (evs != NULL)
        static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);

      // Special treatment for load params and steady state statement: insert initial values for the auxiliary variables
      LoadParamsAndSteadyStateStatement *lpass = dynamic_cast<LoadParamsAndSteadyStateStatement *>(*it);
      if (lpass && !no_static)
        static_model.writeAuxVarInitval(mOutputFile, oMatlabOutsideModel);
    }

  // Remove path for block option with M-files
  if (block && !byte_code)
    mOutputFile << "rmpath " << basename << ";" << endl;

  mOutputFile << "save('" << basename << "_results.mat', 'oo_', 'M_', 'options_');" << endl
              << "if exist('estim_params_', 'var') == 1" << endl
              << "  save('" << basename << "_results.mat', 'estim_params_', '-append');" << endl << "end" << endl
              << "if exist('bayestopt_', 'var') == 1" << endl
              << "  save('" << basename << "_results.mat', 'bayestopt_', '-append');" << endl << "end" << endl
              << "if exist('dataset_', 'var') == 1" << endl
              << "  save('" << basename << "_results.mat', 'dataset_', '-append');" << endl << "end" << endl
              << "if exist('estimation_info', 'var') == 1" << endl
              << "  save('" << basename << "_results.mat', 'estimation_info', '-append');" << endl << "end" << endl;

  config_file.writeEndParallel(mOutputFile);

  mOutputFile << endl << endl
	      << "disp(['Total computing time : ' dynsec2hms(toc) ]);" << endl;

  if (!no_warn)
    {
      if (warnings.countWarnings() > 0)
        mOutputFile << "disp('Note: " << warnings.countWarnings() << " warning(s) encountered in the preprocessor')" << endl;

      mOutputFile << "if ~isempty(lastwarn)" << endl
                  << "  disp('Note: warning(s) encountered in MATLAB/Octave code')" << endl
                  << "end" << endl;
    }
  

  if (!no_log)
    mOutputFile << "diary off" << endl;

  mOutputFile.close();

  // Create static and dynamic files
  if (dynamic_model.equation_number() > 0)
    {
      if (!no_static)
        {
          static_model.writeStaticFile(basename, block, byte_code, use_dll);
          static_model.writeParamsDerivativesFile(basename);
        }

      dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);
      dynamic_model.writeParamsDerivativesFile(basename);
    }

  // Create steady state file
  steady_state_model.writeSteadyStateFile(basename, mod_file_struct.ramsey_model_present);

  cout << "done" << endl;
}

void
ModFile::writeModelC(const string &basename, bool cuda) const
{
  string filename = basename + ".c";

  ofstream mDriverCFile;
  mDriverCFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDriverCFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  mDriverCFile << "/*" << endl
               << " * " << filename << " : Driver file for Dynare C code" << endl
               << " *" << endl
               << " * Warning : this file is generated automatically by Dynare" << endl
               << " *           from model file (.mod)" << endl
               << " */" << endl
               << endl
               << "#include \"dynare_driver.h\"" << endl
               << endl
               << "struct" << endl
               << "{" << endl;

  // Write basic info
  symbol_table.writeCOutput(mDriverCFile);

  mDriverCFile << endl << "params.resize(param_nbr);" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeCOutput(mDriverCFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      //      if (!no_static)
      //        static_model.writeCOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
      (*it)->writeCOutput(mDriverCFile, basename);

  mDriverCFile << "} DynareInfo;" << endl;
  mDriverCFile.close();

  // Write informational m file
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : informational m file" << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl
              << "disp('The following C file was successfully created:');" << endl
              << "ls preprocessorOutput.c" << endl << endl;
  mOutputFile.close();
}

void
ModFile::writeModelCC(const string &basename, bool cuda) const
{
  string filename = basename + ".cc";

  ofstream mDriverCFile;
  mDriverCFile.open(filename.c_str(), ios::out | ios::binary);
  if (!mDriverCFile.is_open())
    {
      cerr << "Error: Can't open file " << filename << " for writing" << endl;
      exit(EXIT_FAILURE);
    }

  mDriverCFile << "/*" << endl
               << " * " << filename << " : Driver file for Dynare C++ code" << endl
               << " *" << endl
               << " * Warning : this file is generated automatically by Dynare" << endl
               << " *           from model file (.mod)" << endl
               << " */" << endl
               << endl
               << "#include \"dynare_cpp_driver.hh\"" << endl
               << endl
               << "DynareInfo::DynareInfo(void)" << endl
               << "{" << endl;

  // Write basic info
  symbol_table.writeCCOutput(mDriverCFile);

  mDriverCFile << endl << "params.resize(param_nbr);" << endl;

  if (dynamic_model.equation_number() > 0)
    {
      dynamic_model.writeCOutput(mDriverCFile, basename, block, byte_code, use_dll, mod_file_struct.order_option, mod_file_struct.estimation_present);
      //      if (!no_static)
      //        static_model.writeCOutput(mOutputFile, block);
    }

  // Print statements
  for (vector<Statement *>::const_iterator it = statements.begin();
       it != statements.end(); it++)
      (*it)->writeCOutput(mDriverCFile, basename);

  mDriverCFile << "};" << endl;
  mDriverCFile.close();

  // Write informational m file
  ofstream mOutputFile;

  if (basename.size())
    {
      string fname(basename);
      fname += ".m";
      mOutputFile.open(fname.c_str(), ios::out | ios::binary);
      if (!mOutputFile.is_open())
        {
          cerr << "ERROR: Can't open file " << fname
               << " for writing" << endl;
          exit(EXIT_FAILURE);
        }
    }
  else
    {
      cerr << "ERROR: Missing file name" << endl;
      exit(EXIT_FAILURE);
    }

  mOutputFile << "%" << endl
              << "% Status : informational m file" << endl
              << "%" << endl
              << "% Warning : this file is generated automatically by Dynare" << endl
              << "%           from model file (.mod)" << endl << endl
              << "disp('The following C++ file was successfully created:');" << endl
              << "ls preprocessorOutput.cc" << endl << endl;
  mOutputFile.close();
}

void
ModFile::writeExternalFiles(const string &basename, FileOutputType output, bool cuda) const
{
  writeModelC(basename, cuda);
  steady_state_model.writeSteadyStateFileC(basename, mod_file_struct.ramsey_model_present, cuda);

  dynamic_model.writeDynamicFile(basename, block, byte_code, use_dll, mod_file_struct.order_option);

  if (!no_static)
    static_model.writeStaticFile(basename, false, false, true);


  //  static_model.writeStaticCFile(basename, block, byte_code, use_dll);
  //  static_model.writeParamsDerivativesFileC(basename, cuda);
  //  static_model.writeAuxVarInitvalC(mOutputFile, oMatlabOutsideModel, cuda);

  // dynamic_model.writeResidualsC(basename, cuda);
  // dynamic_model.writeParamsDerivativesFileC(basename, cuda);
  dynamic_model.writeFirstDerivativesC(basename, cuda);
  
  if (output == second)
    dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  else if (output == third)
    {
        dynamic_model.writeSecondDerivativesC_csr(basename, cuda);
  	dynamic_model.writeThirdDerivativesC_csr(basename, cuda);
    }
}
