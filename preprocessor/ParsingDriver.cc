/*
 * Copyright (C) 2003-2014 Dynare Team
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
#include <fstream>
#include <iostream>
#include <cassert>
#include <sstream>
#include <cmath>

#include "ParsingDriver.hh"
#include "Statement.hh"
#include "ExprNode.hh"
#include "WarningConsolidation.hh"

bool
ParsingDriver::symbol_exists_and_is_not_modfile_local_or_external_function(const char *s)
{
  if (!mod_file->symbol_table.exists(s))
    return false;

  SymbolType type = mod_file->symbol_table.getType(s);

  return (type != eModFileLocalVariable && type != eExternalFunction);
}

void
ParsingDriver::check_symbol_existence(const string &name)
{
  if (!mod_file->symbol_table.exists(name))
    error("Unknown symbol: " + name);
}

void
ParsingDriver::check_symbol_is_parameter(string *name)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  if (mod_file->symbol_table.getType(symb_id) != eParameter)
    error(*name + " is not a parameter");
}

void
ParsingDriver::set_current_data_tree(DataTree *data_tree_arg)
{
  data_tree = data_tree_arg;
  model_tree = dynamic_cast<ModelTree *>(data_tree_arg);
  dynamic_model = dynamic_cast<DynamicModel *>(data_tree_arg);
}

void
ParsingDriver::reset_data_tree()
{
  set_current_data_tree(&mod_file->expressions_tree);
}

void
ParsingDriver::reset_current_external_function_options()
{
  current_external_function_options.nargs = eExtFunSetDefaultNargs;
  current_external_function_options.firstDerivSymbID = eExtFunNotSet;
  current_external_function_options.secondDerivSymbID = eExtFunNotSet;
  current_external_function_id = eExtFunNotSet;
}

ModFile *
ParsingDriver::parse(istream &in, bool debug)
{
  mod_file = new ModFile(warnings);

  symbol_list.clear();

  reset_data_tree();
  estim_params.init(*data_tree);
  reset_current_external_function_options();

  lexer = new DynareFlex(&in);
  lexer->set_debug(debug);

  Dynare::parser parser(*this);
  parser.set_debug_level(debug);
  parser.parse();

  delete lexer;

  return mod_file;
}

void
ParsingDriver::error(const Dynare::parser::location_type &l, const string &m)
{
  cerr << "ERROR: " << *l.begin.filename << ": line " << l.begin.line;
  if (l.begin.line == l.end.line)
    if (l.begin.column == l.end.column - 1)
      cerr << ", col " << l.begin.column;
    else
      cerr << ", cols " << l.begin.column << "-" << l.end.column - 1;
  else
    cerr << ", col " << l.begin.column << " -"
         << " line " << l.end.line << ", col " << l.end.column - 1;
  cerr << ": " << m << endl;
  exit(EXIT_FAILURE);
}

void
ParsingDriver::error(const string &m)
{
  error(location, m);
}

void
ParsingDriver::warning(const string &m)
{
  warnings << "WARNING: " << location << ": " << m << endl;
}

void
ParsingDriver::declare_symbol(const string *name, SymbolType type, const string *tex_name, const string *long_name)
{
  try
    {
      if (tex_name == NULL && long_name == NULL)
        mod_file->symbol_table.addSymbol(*name, type);
      else
        if (tex_name == NULL)
          mod_file->symbol_table.addSymbol(*name, type, "", *long_name);
        else if (long_name == NULL)
          mod_file->symbol_table.addSymbol(*name, type, *tex_name, "");
        else
          mod_file->symbol_table.addSymbol(*name, type, *tex_name, *long_name);
    }
  catch (SymbolTable::AlreadyDeclaredException &e)
    {
      if (e.same_type)
        warning("Symbol " + *name + " declared twice.");
      else
        error("Symbol " + *name + " declared twice with different types!");
    }
}

void
ParsingDriver::declare_endogenous(string *name, string *tex_name, string *long_name)
{
  declare_symbol(name, eEndogenous, tex_name, long_name);
  delete name;
  if (tex_name != NULL)
    delete tex_name;
  if (long_name != NULL)
    delete long_name;
}

void
ParsingDriver::declare_exogenous(string *name, string *tex_name, string *long_name)
{
  declare_symbol(name, eExogenous, tex_name, long_name);
  delete name;
  if (tex_name != NULL)
    delete tex_name;
  if (long_name != NULL)
    delete long_name;
}

void
ParsingDriver::declare_exogenous_det(string *name, string *tex_name, string *long_name)
{
  declare_symbol(name, eExogenousDet, tex_name, long_name);
  delete name;
  if (tex_name != NULL)
    delete tex_name;
  if (long_name != NULL)
    delete long_name;
}

void
ParsingDriver::declare_parameter(string *name, string *tex_name, string *long_name)
{
  declare_symbol(name, eParameter, tex_name, long_name);
  delete name;
  if (tex_name != NULL)
    delete tex_name;
  if (long_name != NULL)
    delete long_name;
}

void
ParsingDriver::declare_statement_local_variable(string *name)
{
  if (mod_file->symbol_table.exists(*name))
    error("Symbol " + *name + " cannot be assigned within a statement " +
          "while being assigned elsewhere in the modfile");
  declare_symbol(name, eStatementDeclaredVariable, NULL, NULL);
  delete name;
}

void
ParsingDriver::declare_optimal_policy_discount_factor_parameter(expr_t exprnode)
{
  string *optimalParName_declare = new string("optimal_policy_discount_factor");
  string *optimalParName_init = new string("optimal_policy_discount_factor");
  if (mod_file->symbol_table.exists(*optimalParName_declare))
    error("Symbol optimal_policy_discount_factor is needed by Dynare when using a ramsey_model, a ramsey_policy or a discretionary_policy statement");
  declare_parameter(optimalParName_declare, NULL);
  init_param(optimalParName_init, exprnode);
}

void
ParsingDriver::begin_trend()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::declare_trend_var(bool log_trend, string *name, string *tex_name)
{
  declare_symbol(name, log_trend ? eLogTrend : eTrend, tex_name, NULL);
  declared_trend_vars.push_back(mod_file->symbol_table.getID(*name));
  delete name;
  if (tex_name != NULL)
    delete tex_name;
}

void
ParsingDriver::end_trend_var(expr_t growth_factor)
{
  try
    {
      dynamic_model->addTrendVariables(declared_trend_vars, growth_factor);
    }
  catch (DataTree::TrendException &e)
    {
      error("Trend variable " + e.name + " was declared twice.");
    }
  declared_trend_vars.clear();
  reset_data_tree();
}

void
ParsingDriver::add_predetermined_variable(string *name)
{
  try
    {
      int symb_id = mod_file->symbol_table.getID(*name);
      if (mod_file->symbol_table.getType(symb_id) != eEndogenous)
        error("Predetermined variables must be endogenous variables");

      mod_file->symbol_table.markPredetermined(symb_id);
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      error("Undeclared symbol name: " + *name);
    }
  delete name;
}

void
ParsingDriver::add_equation_tags(string *key, string *value)
{
  eq_tags.push_back(make_pair(*key, *value));
  delete key;
  delete value;
}

expr_t
ParsingDriver::add_non_negative_constant(string *constant)
{
  expr_t id = data_tree->AddNonNegativeConstant(*constant);
  delete constant;
  return id;
}

expr_t
ParsingDriver::add_nan_constant()
{
  return data_tree->NaN;
}

expr_t
ParsingDriver::add_inf_constant()
{
  return data_tree->Infinity;
}

expr_t
ParsingDriver::add_model_variable(string *name)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  delete name;
  return add_model_variable(symb_id, 0);
}

expr_t
ParsingDriver::add_model_variable(int symb_id, int lag)
{
  assert(symb_id >= 0);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type == eModFileLocalVariable)
    error("Variable " + mod_file->symbol_table.getName(symb_id) + " not allowed inside model declaration. Its scope is only outside model.");

  if (type == eExternalFunction)
    error("Symbol " + mod_file->symbol_table.getName(symb_id) + " is a function name external to Dynare. It cannot be used inside model.");

  if (type == eModelLocalVariable && lag != 0)
    error("Model local variable " + mod_file->symbol_table.getName(symb_id) + " cannot be given a lead or a lag.");

  if (dynamic_cast<StaticModel *>(model_tree) != NULL && lag != 0)
    error("Leads and lags on variables are forbidden in 'planner_objective'.");

  // It makes sense to allow a lead/lag on parameters: during steady state calibration, endogenous and parameters can be swapped
  return model_tree->AddVariable(symb_id, lag);
}

expr_t
ParsingDriver::add_expression_variable(string *name)
{
  // If symbol doesn't exist, then declare it as a mod file local variable
  if (!mod_file->symbol_table.exists(*name))
    mod_file->symbol_table.addSymbol(*name, eModFileLocalVariable);

  // This check must come after the previous one!
  if (mod_file->symbol_table.getType(*name) == eModelLocalVariable)
    error("Variable " + *name + " not allowed outside model declaration. Its scope is only inside model.");

  if (mod_file->symbol_table.getType(*name) == eTrend
      || mod_file->symbol_table.getType(*name) == eLogTrend)
    error("Variable " + *name + " not allowed outside model declaration, because it is a trend variable.");

  if (mod_file->symbol_table.getType(*name) == eExternalFunction)
    error("Symbol '" + *name + "' is the name of a MATLAB/Octave function, and cannot be used as a variable.");

  int symb_id = mod_file->symbol_table.getID(*name);
  expr_t id = data_tree->AddVariable(symb_id);

  delete name;
  return id;
}

void
ParsingDriver::declare_nonstationary_var(string *name, string *tex_name, string *long_name)
{
  if (tex_name == NULL && long_name == NULL)
    declare_endogenous(new string(*name));
  else
    if (tex_name == NULL)
      declare_endogenous(new string(*name), NULL, new string(*long_name));
    else if (long_name == NULL)
      declare_endogenous(new string(*name), new string(*tex_name));
    else
      declare_endogenous(new string(*name), new string(*tex_name), new string(*long_name));

  declared_nonstationary_vars.push_back(mod_file->symbol_table.getID(*name));
  mod_file->nonstationary_variables = true;
  delete name;
  if (tex_name != NULL)
    delete tex_name;
  if (long_name != NULL)
    delete long_name;
}

void
ParsingDriver::end_nonstationary_var(bool log_deflator, expr_t deflator)
{
  try
    {
      dynamic_model->addNonstationaryVariables(declared_nonstationary_vars, log_deflator, deflator);
    }
  catch (DataTree::TrendException &e)
    {
      error("Variable " + e.name + " was listed more than once as following a trend.");
    }

  set<int> r;
  deflator->collectVariables(eEndogenous, r);
  for (set<int>::const_iterator it = r.begin(); it != r.end(); ++it)
    if (dynamic_model->isNonstationary(*it))
      error("The deflator contains a non-stationary endogenous variable. This is not allowed. Please use only stationary endogenous and/or {log_}trend_vars.");

  declared_nonstationary_vars.clear();
  reset_data_tree();
}

void
ParsingDriver::periods(string *periods)
{
  warning("periods: this command is now deprecated and may be removed in a future version of Dynare. Please use the ''periods'' option of the ''simul'' command instead.");

  int periods_val = atoi(periods->c_str());
  mod_file->addStatement(new PeriodsStatement(periods_val));
  delete periods;
}

void
ParsingDriver::dsample(string *arg1)
{
  int arg1_val = atoi(arg1->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val));
  delete arg1;
}

void
ParsingDriver::dsample(string *arg1, string *arg2)
{
  int arg1_val = atoi(arg1->c_str());
  int arg2_val = atoi(arg2->c_str());
  mod_file->addStatement(new DsampleStatement(arg1_val, arg2_val));
  delete arg1;
  delete arg2;
}

void
ParsingDriver::init_param(string *name, expr_t rhs)
{
  check_symbol_is_parameter(name);
  int symb_id = mod_file->symbol_table.getID(*name);
  mod_file->addStatement(new InitParamStatement(symb_id, rhs, mod_file->symbol_table));
  delete name;
}

void
ParsingDriver::init_val(string *name, expr_t rhs)
{
  if (nostrict)
    if (!mod_file->symbol_table.exists(*name))
      {
        warnings << "WARNING: discarding '" << *name
                 << "' not recognized in initval or endval statement" << endl;
        delete name;
        return;
      }

  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("initval/endval: " + *name + " should be an endogenous or exogenous variable");

  init_values.push_back(make_pair(symb_id, rhs));

  delete name;
}

void
ParsingDriver::initval_file(string *filename)
{
  mod_file->addStatement(new InitvalFileStatement(*filename));
  delete filename;
}

void
ParsingDriver::hist_val(string *name, string *lag, expr_t rhs)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eEndogenous
      && type != eExogenous
      && type != eExogenousDet)
    error("histval: " + *name + " should be an endogenous or exogenous variable");

  int ilag = atoi(lag->c_str());
  pair<int, int> key(symb_id, ilag);

  if (mod_file->dynamic_model.minLagForSymbol(symb_id) > ilag - 1)
    {
      ostringstream s;
      s << ilag-1;
      error("histval: variable " + *name + " does not appear in the model with the lag " + s.str() + " (see the reference manual for the timing convention in 'histval')");
    }

  if (hist_values.find(key) != hist_values.end())
    error("hist_val: (" + *name + ", " + *lag + ") declared twice");

  hist_values[key] = rhs;

  delete name;
  delete lag;
}

void
ParsingDriver::homotopy_val(string *name, expr_t val1, expr_t val2)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (type != eParameter
      && type != eExogenous
      && type != eExogenousDet)
    error("homotopy_val: " + *name + " should be a parameter or exogenous variable");

  homotopy_values.push_back(make_pair(symb_id, make_pair(val1, val2)));

  delete name;
}

void
ParsingDriver::forecast()
{
  mod_file->addStatement(new ForecastStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::use_dll()
{
  mod_file->use_dll = true;
}

void
ParsingDriver::block()
{
  mod_file->block = true;
}

void
ParsingDriver::no_static()
{
  mod_file->no_static = true;
}

void
ParsingDriver::byte_code()
{
  mod_file->byte_code = true;
}

void
ParsingDriver::differentiate_forward_vars_all()
{
  mod_file->differentiate_forward_vars = true;
}

void
ParsingDriver::differentiate_forward_vars_some()
{
  mod_file->differentiate_forward_vars = true;
  mod_file->differentiate_forward_vars_subset = symbol_list.get_symbols();
  for (vector<string>::const_iterator it = mod_file->differentiate_forward_vars_subset.begin();
       it != mod_file->differentiate_forward_vars_subset.end(); ++it)
    {
      check_symbol_existence(*it);
      if (mod_file->symbol_table.getType(*it) != eEndogenous)
        error("Symbol " + *it + " is not an endogenous");
    }
  symbol_list.clear();
}

void
ParsingDriver::cutoff(string *value)
{
  double val = atof(value->c_str());
  mod_file->dynamic_model.cutoff = val;
  mod_file->static_model.cutoff = val;
  delete value;
}

void
ParsingDriver::mfs(string *value)
{
  int val = atoi(value->c_str());
  mod_file->dynamic_model.mfs = val;
  mod_file->static_model.mfs = val;
  delete value;
}

void
ParsingDriver::end_initval(bool all_values_required)
{
  mod_file->addStatement(new InitValStatement(init_values, mod_file->symbol_table, all_values_required));
  init_values.clear();
}

void
ParsingDriver::end_endval(bool all_values_required)
{
  mod_file->addStatement(new EndValStatement(init_values, mod_file->symbol_table, all_values_required));
  init_values.clear();
}

void
ParsingDriver::end_histval()
{
  mod_file->addStatement(new HistValStatement(hist_values, mod_file->symbol_table));
  hist_values.clear();
}

void
ParsingDriver::end_homotopy()
{
  mod_file->addStatement(new HomotopyStatement(homotopy_values, mod_file->symbol_table));
  homotopy_values.clear();
}

void
ParsingDriver::begin_model()
{
  set_current_data_tree(&mod_file->dynamic_model);
}

void
ParsingDriver::end_shocks(bool overwrite)
{
  mod_file->addStatement(new ShocksStatement(overwrite, det_shocks, var_shocks, std_shocks,
                                             covar_shocks, corr_shocks, mod_file->symbol_table));
  det_shocks.clear();
  var_shocks.clear();
  std_shocks.clear();
  covar_shocks.clear();
  corr_shocks.clear();
}

void
ParsingDriver::end_mshocks(bool overwrite)
{
  mod_file->addStatement(new MShocksStatement(overwrite, det_shocks, mod_file->symbol_table));
  det_shocks.clear();
}

void
ParsingDriver::add_det_shock(string *var, bool conditional_forecast)
{
  check_symbol_existence(*var);
  int symb_id = mod_file->symbol_table.getID(*var);
  SymbolType type = mod_file->symbol_table.getType(symb_id);

  if (conditional_forecast)
    {
      if (type != eEndogenous)
        error("conditional_forecast_paths: shocks can only be applied to endogenous variables");
    }
  else
    {
      if (type != eExogenous && type != eExogenousDet)
        error("shocks: shocks can only be applied to exogenous variables");
    }

  if (det_shocks.find(symb_id) != det_shocks.end())
    error("shocks/conditional_forecast_paths: variable " + *var + " declared twice");

  if (det_shocks_periods.size() != det_shocks_values.size())
    error("shocks/conditional_forecast_paths: variable " + *var + ": number of periods is different from number of shock values");

  vector<ShocksStatement::DetShockElement> v;

  for (size_t i = 0; i < det_shocks_periods.size(); i++)
    {
      ShocksStatement::DetShockElement dse;
      dse.period1 = det_shocks_periods[i].first;
      dse.period2 = det_shocks_periods[i].second;
      dse.value = det_shocks_values[i];
      v.push_back(dse);
    }

  det_shocks[symb_id] = v;
  
  det_shocks_periods.clear();
  det_shocks_values.clear();
  delete var;
}

void
ParsingDriver::add_stderr_shock(string *var, expr_t value)
{
  check_symbol_existence(*var);
  int symb_id = mod_file->symbol_table.getID(*var);

  if (var_shocks.find(symb_id) != var_shocks.end()
      || std_shocks.find(symb_id) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  std_shocks[symb_id] = value;

  delete var;
}

void
ParsingDriver::add_var_shock(string *var, expr_t value)
{
  check_symbol_existence(*var);
  int symb_id = mod_file->symbol_table.getID(*var);

  if (var_shocks.find(symb_id) != var_shocks.end()
      || std_shocks.find(symb_id) != std_shocks.end())
    error("shocks: variance or stderr of shock on " + *var + " declared twice");

  var_shocks[symb_id] = value;

  delete var;
}

void
ParsingDriver::add_covar_shock(string *var1, string *var2, expr_t value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);
  int symb_id1 = mod_file->symbol_table.getID(*var1);
  int symb_id2 = mod_file->symbol_table.getID(*var2);

  pair<int, int> key(symb_id1, symb_id2), key_inv(symb_id2, symb_id1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");

  covar_shocks[key] = value;

  delete var1;
  delete var2;
}

void
ParsingDriver::add_correl_shock(string *var1, string *var2, expr_t value)
{
  check_symbol_existence(*var1);
  check_symbol_existence(*var2);
  int symb_id1 = mod_file->symbol_table.getID(*var1);
  int symb_id2 = mod_file->symbol_table.getID(*var2);

  pair<int, int> key(symb_id1, symb_id2), key_inv(symb_id2, symb_id1);

  if (covar_shocks.find(key) != covar_shocks.end()
      || covar_shocks.find(key_inv) != covar_shocks.end()
      || corr_shocks.find(key) != corr_shocks.end()
      || corr_shocks.find(key_inv) != corr_shocks.end())
    error("shocks: covariance or correlation shock on variable pair (" + *var1 + ", "
          + *var2 + ") declared twice");

  corr_shocks[key] = value;

  delete var1;
  delete var2;
}

void
ParsingDriver::add_period(string *p1, string *p2)
{
  int p1_val = atoi(p1->c_str());
  int p2_val = atoi(p2->c_str());
  if (p1_val > p2_val)
    error("shocks/conditional_forecast_paths: can't have first period index greater than second index in range specification");
  det_shocks_periods.push_back(make_pair(p1_val, p2_val));
  delete p1;
  delete p2;
}

void
ParsingDriver::add_period(string *p1)
{
  int p1_val = atoi(p1->c_str());
  det_shocks_periods.push_back(make_pair(p1_val, p1_val));
  delete p1;
}

void
ParsingDriver::add_value(expr_t value)
{
  det_shocks_values.push_back(value);
}

void
ParsingDriver::add_value(string *v)
{
  expr_t id;

  if (v->at(0) == '-')
    id = data_tree->AddUMinus(data_tree->AddNonNegativeConstant(v->substr(1, string::npos)));
  else
    id = data_tree->AddNonNegativeConstant(*v);

  delete v;
  det_shocks_values.push_back(id);
}

void
ParsingDriver::begin_svar_identification()
{
  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_constants_exclusion = false;
}

void
ParsingDriver::end_svar_identification()
{
  mod_file->addStatement(new SvarIdentificationStatement(svar_ident_restrictions,
                                                         svar_upper_cholesky,
                                                         svar_lower_cholesky,
							 svar_constants_exclusion,
                                                         mod_file->symbol_table));
  svar_restriction_symbols.clear();
  svar_equation_restrictions.clear();
  svar_ident_restrictions.clear();
  svar_Qi_restriction_nbr.clear();
  svar_Ri_restriction_nbr.clear();
}

void
ParsingDriver::combine_lag_and_restriction(string *lag)
{
  int current_lag = atoi(lag->c_str());

  for (SvarIdentificationStatement::svar_identification_restrictions_t::const_iterator it = svar_ident_restrictions.begin();
       it != svar_ident_restrictions.end(); it++)
    if (it->lag == current_lag)
      error("lag " + *lag + " used more than once.");

  for (map<int, vector<int> >::const_iterator it = svar_equation_restrictions.begin();
       it != svar_equation_restrictions.end(); it++)
    for (vector<int>::const_iterator it1 = it->second.begin();
	 it1 != it->second.end(); it1++)
      {
	SvarIdentificationStatement::svar_identification_restriction new_restriction;
	new_restriction.equation = it->first;
	if (current_lag > 0)
	  new_restriction.restriction_nbr = ++svar_Ri_restriction_nbr[it->first];
	else
	  new_restriction.restriction_nbr = ++svar_Qi_restriction_nbr[it->first];
	new_restriction.lag = current_lag;
	new_restriction.variable = *it1;
	new_restriction.value = data_tree->One;
	svar_ident_restrictions.push_back(new_restriction);
      } 
  //    svar_ident_exclusion_values[make_pair(current_lag, it->first)] = it->second;

  svar_upper_cholesky = false;
  svar_lower_cholesky = false;
  svar_equation_restrictions.clear();
  delete lag;
}

void
ParsingDriver::add_restriction_in_equation(string *equation)
{
  int eqn = atoi(equation->c_str());
  if (eqn < 1)
    error("equation numbers must be greater than or equal to 1.");

  if (svar_equation_restrictions.count(eqn) > 0)
    error("equation number " + *equation + " referenced more than once under a single lag.");

  svar_equation_restrictions[eqn] = svar_restriction_symbols;

  svar_restriction_symbols.clear();
  delete equation;
}

void
ParsingDriver::add_in_svar_restriction_symbols(string *tmp_var)
{
  check_symbol_existence(*tmp_var);
  int symb_id = mod_file->symbol_table.getID(*tmp_var);

  for (vector<int>::const_iterator viit = svar_restriction_symbols.begin();
       viit != svar_restriction_symbols.end(); viit++)
    if (symb_id == *viit)
      error(*tmp_var + " restriction added twice.");

  svar_restriction_symbols.push_back(symb_id);
  delete tmp_var;
}

void 
ParsingDriver::add_restriction_equation_nbr(string *eq_nbr)
{
  svar_equation_nbr = atoi(eq_nbr->c_str());
  svar_left_handside = true;
  // reinitialize restriction type that must be set from the first restriction element
  svar_restriction_type = ParsingDriver::NOT_SET;
}

void
ParsingDriver::add_restriction_equal()
{
  if (svar_left_handside)
    svar_left_handside = false;
  else
    error("svar_identification: there are more than one EQUAL sign in a restriction equation");
}

void
ParsingDriver::add_positive_restriction_element(expr_t value, string *variable, string *lag)
{
  // if the expression is not on the left handside, change its sign
  if (!svar_left_handside)
    value = add_uminus(value);
  
  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_positive_restriction_element(string *variable, string *lag)
{
  expr_t value(data_tree->One); 

  // if the expression is not on the left handside, change its sign
  if (!svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_negative_restriction_element(expr_t value, string *variable, string *lag)
{
  // if the expression is on the left handside, change its sign
  if (svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_negative_restriction_element(string *variable, string *lag)
{
  expr_t value(data_tree->One); 

  // if the expression is on the left handside, change its sign
  if (svar_left_handside)
    value = add_uminus(value);

  add_restriction_element(value, variable, lag);
}

void
ParsingDriver::add_restriction_element(expr_t value, string *variable, string *lag)
{
  check_symbol_existence(*variable);
  int symb_id = mod_file->symbol_table.getID(*variable);

  int current_lag = atoi(lag->c_str());
  if (svar_restriction_type == ParsingDriver::NOT_SET)
    {
      if (current_lag == 0)
	{
	  svar_restriction_type = ParsingDriver::Qi_TYPE;
	  ++svar_Qi_restriction_nbr[svar_equation_nbr];
	}
      else
	{
	  svar_restriction_type = ParsingDriver::Ri_TYPE;
	  ++svar_Ri_restriction_nbr[svar_equation_nbr];
	}
    }
  else
    {
      if ((svar_restriction_type == Qi_TYPE && current_lag > 0)
	  || (svar_restriction_type == Ri_TYPE && current_lag == 0))
	error("SVAR_IDENTIFICATION: a single restrictions must affect either Qi or Ri, but not both");
    }
  SvarIdentificationStatement::svar_identification_restriction new_restriction;
  new_restriction.equation = svar_equation_nbr;
  if (current_lag > 0)
    new_restriction.restriction_nbr = svar_Ri_restriction_nbr[svar_equation_nbr];
  else
    new_restriction.restriction_nbr = svar_Qi_restriction_nbr[svar_equation_nbr];
  new_restriction.lag = current_lag;
  new_restriction.variable = symb_id;
  new_restriction.value = value;

  svar_ident_restrictions.push_back(new_restriction);
}

void
ParsingDriver::check_restriction_expression_constant(expr_t value)
{
  if (value->eval(eval_context_t()) != 0)
    error("SVAR_INDENTIFICATION restrictions must be homogenous");
}

void
ParsingDriver::add_upper_cholesky()
{
  svar_upper_cholesky = true;
}

void
ParsingDriver::add_lower_cholesky()
{
  svar_lower_cholesky = true;
}

void
ParsingDriver::add_constants_exclusion()
{
  svar_constants_exclusion = true;
}

void
ParsingDriver::do_sigma_e()
{
  warning("Sigma_e: this command is now deprecated and may be removed in a future version of Dynare. Please use the ''shocks'' command instead.");

  try
    {
      mod_file->addStatement(new SigmaeStatement(sigmae_matrix));
    }
  catch (SigmaeStatement::MatrixFormException &e)
    {
      error("Sigma_e: matrix is neither upper triangular nor lower triangular");
    }
  sigmae_matrix.clear();
}

void
ParsingDriver::end_of_row()
{
  sigmae_matrix.push_back(sigmae_row);
  sigmae_row.clear();
}

void
ParsingDriver::add_to_row_const(string *v)
{
  expr_t id;

  if (v->at(0) == '-')
    id = data_tree->AddUMinus(data_tree->AddNonNegativeConstant(v->substr(1, string::npos)));
  else
    id = data_tree->AddNonNegativeConstant(*v);

  delete v;
  sigmae_row.push_back(id);
}

void
ParsingDriver::add_to_row(expr_t v)
{
  sigmae_row.push_back(v);
}

void
ParsingDriver::steady()
{
  mod_file->addStatement(new SteadyStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::option_num(const string &name_option, string *opt1, string *opt2)
{
  if (options_list.paired_num_options.find(name_option)
      != options_list.paired_num_options.end())
    error("option " + name_option + " declared twice");

  options_list.paired_num_options[name_option] = make_pair(*opt1, *opt2);
  delete opt1;
  delete opt2;
}

void
ParsingDriver::option_num(const string &name_option, string *opt)
{
  option_num(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_num(const string &name_option, const string &opt)
{
  if (options_list.num_options.find(name_option) != options_list.num_options.end())
    error("option " + name_option + " declared twice");

  options_list.num_options[name_option] = opt;
}

void
ParsingDriver::option_str(const string &name_option, string *opt)
{
  option_str(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_str(const string &name_option, const string &opt)
{
  if (options_list.string_options.find(name_option)
      != options_list.string_options.end())
    error("option " + name_option + " declared twice");

  options_list.string_options[name_option] = opt;
}

void
ParsingDriver::option_date(const string &name_option, string *opt)
{
  option_date(name_option, *opt);
  delete opt;
}

void
ParsingDriver::option_date(const string &name_option, const string &opt)
{
  if (options_list.date_options.find(name_option)
      != options_list.date_options.end())
    error("option " + name_option + " declared twice");

  options_list.date_options[name_option] = opt;
}

void
ParsingDriver::option_symbol_list(const string &name_option)
{
  if (options_list.symbol_list_options.find(name_option)
      != options_list.symbol_list_options.end())
    error("option " + name_option + " declared twice");

  if (name_option.compare("irf_shocks")==0)
    {
      vector<string> shocks = symbol_list.get_symbols();
      for (vector<string>::const_iterator it = shocks.begin();
           it != shocks.end(); it++)
        if (mod_file->symbol_table.getType(*it) != eExogenous)
          error("Variables passed to irf_shocks must be exogenous. Caused by: " + *it);
    }

  if (name_option.compare("ms.parameters")==0)
    {
      vector<string> parameters = symbol_list.get_symbols();
      for (vector<string>::const_iterator it = parameters.begin();
           it != parameters.end(); it++)
        if (mod_file->symbol_table.getType(*it) != eParameter)
          error("Variables passed to the parameters option of the markov_switching statement must be parameters. Caused by: " + *it);
    }

  options_list.symbol_list_options[name_option] = symbol_list;
  symbol_list.clear();
}

void
ParsingDriver::option_vec_int(const string &name_option, const vector<int> *opt)
{
  if (options_list.vector_int_options.find(name_option)
      != options_list.vector_int_options.end())
    error("option " + name_option + " declared twice");

  if ((*opt).empty())
    error("option " + name_option + " was passed an empty vector.");

  options_list.vector_int_options[name_option] = *opt;
  delete opt;
}

void
ParsingDriver::linear()
{
  mod_file->linear = true;
}

void
ParsingDriver::add_in_symbol_list(string *tmp_var)
{
  if (*tmp_var != ":")
    check_symbol_existence(*tmp_var);
  symbol_list.addSymbol(*tmp_var);
  delete tmp_var;
}

void
ParsingDriver::rplot()
{
  mod_file->addStatement(new RplotStatement(symbol_list, options_list));
  options_list.clear();
  symbol_list.clear();
}

void
ParsingDriver::stoch_simul()
{
  mod_file->addStatement(new StochSimulStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::simul()
{
  mod_file->addStatement(new SimulStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::model_info()
{
  mod_file->addStatement(new ModelInfoStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::check()
{
  mod_file->addStatement(new CheckStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::add_estimated_params_element()
{
  if (estim_params.name != "dsge_prior_weight")
    {
      check_symbol_existence(estim_params.name);
      SymbolType type = mod_file->symbol_table.getType(estim_params.name);
      switch (estim_params.type)
        {
        case 1:
          if (type != eEndogenous && type != eExogenous)
            error(estim_params.name + " must be an endogenous or an exogenous variable");
          break;
        case 2:
          check_symbol_is_parameter(&estim_params.name);
          break;
        case 3:
          check_symbol_existence(estim_params.name2);
          SymbolType type2 = mod_file->symbol_table.getType(estim_params.name2);
          if ((type != eEndogenous && type != eExogenous) || type != type2)
            error(estim_params.name + " and " + estim_params.name2 + " must either be both endogenous variables or both exogenous");
          break;
        }
    }
  estim_params_list.push_back(estim_params);
  estim_params.init(*data_tree);
}

void
ParsingDriver::estimated_params()
{
  mod_file->addStatement(new EstimatedParamsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_init(bool use_calibration)
{
  mod_file->addStatement(new EstimatedParamsInitStatement(estim_params_list, mod_file->symbol_table, use_calibration));
  estim_params_list.clear();
}

void
ParsingDriver::estimated_params_bounds()
{
  mod_file->addStatement(new EstimatedParamsBoundsStatement(estim_params_list, mod_file->symbol_table));
  estim_params_list.clear();
}

void
ParsingDriver::set_unit_root_vars()
{
  mod_file->addStatement(new UnitRootVarsStatement());
  warning("''unit_root_vars'' is now obsolete; use the ''diffuse_filter'' option of ''estimation'' instead");
  symbol_list.clear();
}

void
ParsingDriver::set_time(string *arg)
{
  option_date("initial_period", arg);
  mod_file->addStatement(new SetTimeStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::estimation_data()
{
  mod_file->addStatement(new EstimationDataStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::set_subsamples(string *name1, string *name2)
{
  check_symbol_existence(*name1);
  if (!name2->empty())
    check_symbol_existence(*name2);

  mod_file->addStatement(new SubsamplesStatement(*name1, *name2, subsample_declaration_map,
                                                 mod_file->symbol_table));
  subsample_declarations[make_pair(*name1, *name2)] = subsample_declaration_map;
  subsample_declaration_map.clear();
  delete name1;
  delete name2;
}

void
ParsingDriver::copy_subsamples(string *to_name1, string *to_name2, string *from_name1, string *from_name2)
{
  check_symbol_existence(*to_name1);
  check_symbol_existence(*from_name1);
  if (!to_name2->empty())
    check_symbol_existence(*to_name2);
  if (!from_name2->empty())
    check_symbol_existence(*from_name2);

  if (subsample_declarations.find(make_pair(*from_name1,*from_name2)) == subsample_declarations.end())
    {
      string err = *from_name1;
      if (!from_name2->empty())
        err.append(",").append(*from_name2);
      error(err + " does not have an associated subsample statement.");
    }

  mod_file->addStatement(new SubsamplesEqualStatement(*to_name1, *to_name2, *from_name1, *from_name2,
                                                      mod_file->symbol_table));

  subsample_declarations[make_pair(*to_name1, *to_name2)] =
    subsample_declarations[make_pair(*from_name1, *from_name2)];

  delete to_name1;
  delete to_name2;
  delete from_name1;
  delete from_name2;
}

void
ParsingDriver::check_symbol_is_statement_variable(string *name)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  if (mod_file->symbol_table.getType(symb_id) != eStatementDeclaredVariable)
    error(*name + " is not a variable assigned in a statement");
}

void
ParsingDriver::set_subsample_name_equal_to_date_range(string *name, string *date1, string *date2)
{
  if (subsample_declaration_map.find(*name) != subsample_declaration_map.end())
    error("Symbol " + *name + " may only be assigned once in a SUBSAMPLE statement");
  subsample_declaration_map[*name] = make_pair(*date1, *date2);
  delete name;
  delete date1;
  delete date2;
}

void
ParsingDriver::check_subsample_declaration_exists(string *name1, string *subsample_name)
{
  if (subsample_name->empty())
    return;

  string *str_empty = new string ("");
  check_subsample_declaration_exists(name1, str_empty, subsample_name);
  delete str_empty;
}

void
ParsingDriver::check_subsample_declaration_exists(string *name1, string *name2, string *subsample_name)
{
  if (subsample_name->empty())
    return;

  check_symbol_existence(*name1);
  if (!name2->empty())
      check_symbol_existence(*name2);

  subsample_declarations_t::const_iterator it = subsample_declarations.find(make_pair(*name1, *name2));
  if (it == subsample_declarations.end())
    {
      it = subsample_declarations.find(make_pair(*name2, *name1));
      if (it== subsample_declarations.end())
        {
          string err = *name1;
          if (!name2->empty())
            err.append(",").append(*name2);
          error("A subsample statement has not been issued for " + err);
        }
    }

  SubsamplesStatement::subsample_declaration_map_t tmp_map = it->second;
  if (tmp_map.find(*subsample_name) == tmp_map.end())
    error("The subsample name " + *subsample_name + " was not previously declared in a subsample statement.");
}


void
ParsingDriver::set_prior(string *name, string *subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new PriorStatement(*name, *subsample_name, prior_shape, prior_variance, options_list));
  options_list.clear();
  set_prior_variance();
  prior_shape = eNoShape;
  delete name;
  delete subsample_name;
}

void
ParsingDriver::set_prior_variance(expr_t variance)
{
  prior_variance = variance;
}

void
ParsingDriver::copy_prior(string *to_declaration_type, string *to_name1, string *to_name2, string *to_subsample_name,
                          string *from_declaration_type, string *from_name1, string *from_name2, string *from_subsample_name)
{
  if (strcmp(to_declaration_type->c_str(), "par") == 0)
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1);
      if (!to_name2->empty())
        check_symbol_is_endogenous_or_exogenous(to_name2);
    }

  if (strcmp(from_declaration_type->c_str(), "par") == 0)
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1);
      if (!from_name2->empty())
        check_symbol_is_endogenous_or_exogenous(from_name2);
    }

  mod_file->addStatement(new PriorEqualStatement(*to_declaration_type, *to_name1, *to_name2, *to_subsample_name,
                                                 *from_declaration_type, *from_name1, *from_name2, *from_subsample_name,
                                                 mod_file->symbol_table));

  delete to_declaration_type;
  delete to_name1;
  delete to_name2;
  delete to_subsample_name;
  delete from_declaration_type;
  delete from_name1;
  delete from_name2;
  delete from_subsample_name;
}

void
ParsingDriver::set_options(string *name, string *subsample_name)
{
  check_symbol_is_parameter(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new OptionsStatement(*name, *subsample_name, options_list));
  options_list.clear();
  delete name;
  delete subsample_name;
}

void
ParsingDriver::copy_options(string *to_declaration_type, string *to_name1, string *to_name2, string *to_subsample_name,
                            string *from_declaration_type, string *from_name1, string *from_name2, string *from_subsample_name)
{
  if (strcmp(to_declaration_type->c_str(), "par") == 0)
    check_symbol_is_parameter(to_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(to_name1);
      if (!to_name2->empty())
        check_symbol_is_endogenous_or_exogenous(to_name2);
    }

  if (strcmp(from_declaration_type->c_str(), "par") == 0)
    check_symbol_is_parameter(from_name1);
  else
    {
      check_symbol_is_endogenous_or_exogenous(from_name1);
      if (!from_name2->empty())
        check_symbol_is_endogenous_or_exogenous(from_name2);
    }

  mod_file->addStatement(new OptionsEqualStatement(*to_declaration_type, *to_name1, *to_name2, *to_subsample_name,
                                                   *from_declaration_type, *from_name1, *from_name2, *from_subsample_name,
                                                   mod_file->symbol_table));
  delete to_declaration_type;
  delete to_name1;
  delete to_name2;
  delete to_subsample_name;
  delete from_declaration_type;
  delete from_name1;
  delete from_name2;
  delete from_subsample_name;
}

void
ParsingDriver::check_symbol_is_endogenous_or_exogenous(string *name)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  switch(mod_file->symbol_table.getType(symb_id))
    {
    case eEndogenous:
    case eExogenous:
    case eExogenousDet:
      break;
    default:
      error(*name + " is neither endogenous or exogenous.");
    }
}

void
ParsingDriver::set_std_prior(string *name, string *subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new StdPriorStatement(*name, *subsample_name, prior_shape, prior_variance,
                                               options_list, mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = eNoShape;
  delete name;
  delete subsample_name;
}

void
ParsingDriver::set_std_options(string *name, string *subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name);
  check_subsample_declaration_exists(name, subsample_name);
  mod_file->addStatement(new StdOptionsStatement(*name, *subsample_name, options_list, mod_file->symbol_table));
  options_list.clear();
  delete name;
  delete subsample_name;
}

void
ParsingDriver::set_corr_prior(string *name1, string *name2, string *subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1);
  check_symbol_is_endogenous_or_exogenous(name2);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(new CorrPriorStatement(*name1, *name2, *subsample_name, prior_shape, prior_variance,
                                                options_list, mod_file->symbol_table));
  options_list.clear();
  set_prior_variance();
  prior_shape = eNoShape;
  delete name1;
  delete name2;
  delete subsample_name;
}

void
ParsingDriver::set_corr_options(string *name1, string *name2, string *subsample_name)
{
  check_symbol_is_endogenous_or_exogenous(name1);
  check_symbol_is_endogenous_or_exogenous(name2);
  check_subsample_declaration_exists(name1, name2, subsample_name);
  mod_file->addStatement(new CorrOptionsStatement(*name1, *name2, *subsample_name, options_list, mod_file->symbol_table));
  options_list.clear();
  delete name1;
  delete name2;
  delete subsample_name;
}

void
ParsingDriver::run_estimation()
{
  mod_file->addStatement(new EstimationStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::dynare_sensitivity()
{
  mod_file->addStatement(new DynareSensitivityStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::optim_options_helper(const string &name)
{
  if (options_list.string_options.find("optim_opt") == options_list.string_options.end())
    options_list.string_options["optim_opt"] = "";
  else
    options_list.string_options["optim_opt"] += ",";
  options_list.string_options["optim_opt"] += "''" + name + "'',";
}

void
ParsingDriver::optim_options_string(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += "''" + *value + "''";
  delete name;
  delete value;
}

void
ParsingDriver::optim_options_num(string *name, string *value)
{
  optim_options_helper(*name);
  options_list.string_options["optim_opt"] += *value;
  delete name;
  delete value;
}

void
ParsingDriver::check_varobs()
{
  if (mod_file->symbol_table.observedVariablesNbr() > 0)
    error("varobs: you cannot have several 'varobs' statements in the same MOD file");
}

void
ParsingDriver::add_varobs(string *name)
{
  check_symbol_existence(*name);
  int symb_id = mod_file->symbol_table.getID(*name);
  if (mod_file->symbol_table.getType(symb_id) != eEndogenous)
    error("varobs: " + *name + " is not an endogenous variable");
  mod_file->symbol_table.addObservedVariable(symb_id);
  delete name;
}

void
ParsingDriver::set_trends()
{
  mod_file->addStatement(new ObservationTrendsStatement(trend_elements, mod_file->symbol_table));
  trend_elements.clear();
}

void
ParsingDriver::set_trend_element(string *arg1, expr_t arg2)
{
  check_symbol_existence(*arg1);
  if (trend_elements.find(*arg1) != trend_elements.end())
    error("observation_trends: " + *arg1 + " declared twice");
  trend_elements[*arg1] = arg2;
  delete arg1;
}

void
ParsingDriver::set_optim_weights(string *name, expr_t value)
{
  check_symbol_existence(*name);
  if (mod_file->symbol_table.getType(*name) != eEndogenous)
    error("optim_weights: " + *name + " isn't an endogenous variable");
  if (var_weights.find(*name) != var_weights.end())
    error("optim_weights: " + *name + " declared twice");
  var_weights[*name] = value;
  delete name;
}

void
ParsingDriver::set_optim_weights(string *name1, string *name2, expr_t value)
{
  check_symbol_existence(*name1);
  if (mod_file->symbol_table.getType(*name1) != eEndogenous)
    error("optim_weights: " + *name1 + " isn't an endogenous variable");

  check_symbol_existence(*name2);
  if (mod_file->symbol_table.getType(*name2) != eEndogenous)
    error("optim_weights: " + *name2 + " isn't an endogenous variable");

  pair<string, string> covar_key(*name1, *name2);

  if (covar_weights.find(covar_key) != covar_weights.end())
    error("optim_weights: pair of variables (" + *name1 + ", " + *name2
          + ") declared twice");

  covar_weights[covar_key] = value;
  delete name1;
  delete name2;
}

void
ParsingDriver::optim_weights()
{
  mod_file->addStatement(new OptimWeightsStatement(var_weights, covar_weights, mod_file->symbol_table));
  var_weights.clear();
  covar_weights.clear();
}

void
ParsingDriver::set_osr_params()
{
  mod_file->addStatement(new OsrParamsStatement(symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::run_osr()
{
  mod_file->addStatement(new OsrStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::run_dynatype(string *filename)
{
  mod_file->addStatement(new DynaTypeStatement(symbol_list, *filename));
  symbol_list.clear();
  delete filename;
}

void
ParsingDriver::run_dynasave(string *filename)
{
  mod_file->addStatement(new DynaSaveStatement(symbol_list, *filename));
  symbol_list.clear();
  delete filename;
}

void
ParsingDriver::run_load_params_and_steady_state(string *filename)
{
  mod_file->addStatement(new LoadParamsAndSteadyStateStatement(*filename, mod_file->symbol_table, warnings));
  delete filename;
}

void
ParsingDriver::run_save_params_and_steady_state(string *filename)
{
  mod_file->addStatement(new SaveParamsAndSteadyStateStatement(*filename));
  delete filename;
}

void
ParsingDriver::run_identification()
{
  mod_file->addStatement(new IdentificationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::add_mc_filename(string *filename, string *prior)
{
  for (ModelComparisonStatement::filename_list_t::iterator it = filename_list.begin();
       it != filename_list.end(); it++)
    if ((*it).first == *filename)
      error("model_comparison: filename " + *filename + " declared twice");
  filename_list.push_back(make_pair(*filename, *prior));
  delete filename;
  delete prior;
}

void
ParsingDriver::run_model_comparison()
{
  mod_file->addStatement(new ModelComparisonStatement(filename_list, options_list));
  filename_list.clear();
  options_list.clear();
}

void
ParsingDriver::begin_planner_objective()
{
  set_current_data_tree(new StaticModel(mod_file->symbol_table, mod_file->num_constants, mod_file->external_functions_table));
}

void
ParsingDriver::end_planner_objective(expr_t expr)
{
  // Add equation corresponding to expression
  expr_t eq = model_tree->AddEqual(expr, model_tree->Zero);
  model_tree->addEquation(eq, location.begin.line);

  mod_file->addStatement(new PlannerObjectiveStatement(dynamic_cast<StaticModel *>(model_tree)));

  reset_data_tree();
}

void
ParsingDriver::ramsey_model()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new RamseyModelStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::ramsey_policy()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new RamseyPolicyStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::discretionary_policy()
{
  if (!mod_file->symbol_table.exists("optimal_policy_discount_factor"))
    declare_optimal_policy_discount_factor_parameter(data_tree->One);
  mod_file->addStatement(new DiscretionaryPolicyStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::write_latex_dynamic_model()
{
  mod_file->addStatement(new WriteLatexDynamicModelStatement(mod_file->dynamic_model));
}

void
ParsingDriver::write_latex_static_model()
{
  mod_file->addStatement(new WriteLatexStaticModelStatement(mod_file->static_model));
}

void
ParsingDriver::bvar_density(string *maxnlags)
{
  mod_file->addStatement(new BVARDensityStatement(atoi(maxnlags->c_str()), options_list));
  options_list.clear();
  delete maxnlags;
}

void
ParsingDriver::bvar_forecast(string *nlags)
{
  mod_file->addStatement(new BVARForecastStatement(atoi(nlags->c_str()), options_list));
  options_list.clear();
  delete nlags;
}

void
ParsingDriver::sbvar()
{
  mod_file->addStatement(new SBVARStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_estimation()
{
  mod_file->addStatement(new MSSBVAREstimationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_simulation()
{
  mod_file->addStatement(new MSSBVARSimulationStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_compute_mdd()
{
  mod_file->addStatement(new MSSBVARComputeMDDStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_compute_probabilities()
{
  mod_file->addStatement(new MSSBVARComputeProbabilitiesStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_irf()
{
  mod_file->addStatement(new MSSBVARIrfStatement(symbol_list,options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::ms_forecast()
{
  mod_file->addStatement(new MSSBVARForecastStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::ms_variance_decomposition()
{
  mod_file->addStatement(new MSSBVARVarianceDecompositionStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::svar()
{
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  OptionsList::vec_int_options_t::const_iterator itv;

  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");
  if (it0 == options_list.string_options.end()
      && it1 == options_list.string_options.end()
      && it2 == options_list.string_options.end())
    error("You must pass one of 'coefficients', 'variances', or 'constants'.");

  if ((it0 != options_list.string_options.end()
       && it1 != options_list.string_options.end())
      || (it1 != options_list.string_options.end()
          && it2 != options_list.string_options.end())
      || (it0 != options_list.string_options.end()
          && it2 != options_list.string_options.end()))
    error("You may only pass one of 'coefficients', 'variances', or 'constants'.");

  it0 = options_list.num_options.find("ms.chain");
  if (it0 == options_list.num_options.end())
    error("A chain option must be passed to the svar statement.");
  else if (atoi(it0->second.c_str()) <= 0)
    error("The value passed to the chain option must be greater than zero.");

  itv = options_list.vector_int_options.find("ms.equations");
  if (itv != options_list.vector_int_options.end())
    for (vector<int>::const_iterator viit = itv->second.begin(); viit != itv->second.end(); viit++)
      if (*viit <= 0)
        error("The value(s) passed to the equation option must be greater than zero.");

  mod_file->addStatement(new SvarStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::markov_switching()
{
  OptionsList::num_options_t::const_iterator it0;

  it0 = options_list.num_options.find("ms.chain");
  if (it0 == options_list.num_options.end())
    error("A chain option must be passed to the markov_switching statement.");
  else if (atoi(it0->second.c_str()) <= 0)
    error("The value passed to the chain option must be greater than zero.");

  it0 = options_list.num_options.find("ms.number_of_regimes");
  if (it0 == options_list.num_options.end())
    error("A number_of_regimes option must be passed to the markov_switching statement.");
  else if (atoi(it0->second.c_str()) <= 0)
    error("The value passed to the number_of_regimes option must be greater than zero.");

  it0 = options_list.num_options.find("ms.duration");
  if (it0 == options_list.num_options.end())
    error("A duration option must be passed to the markov_switching statement.");

  mod_file->addStatement(new MarkovSwitchingStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::shock_decomposition()
{
  mod_file->addStatement(new ShockDecompositionStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::conditional_forecast()
{
  mod_file->addStatement(new ConditionalForecastStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::plot_conditional_forecast(string *periods)
{
  int nperiods;
  if (periods == NULL)
    nperiods = -1;
  else
    {
      nperiods = atoi(periods->c_str());
      delete periods;
    }
  mod_file->addStatement(new PlotConditionalForecastStatement(nperiods, symbol_list));
  symbol_list.clear();
}

void
ParsingDriver::conditional_forecast_paths()
{
  mod_file->addStatement(new ConditionalForecastPathsStatement(det_shocks));
  det_shocks.clear();
}

void
ParsingDriver::calib_smoother()
{
  mod_file->addStatement(new CalibSmootherStatement(symbol_list, options_list));
  symbol_list.clear();
  options_list.clear();
}

void
ParsingDriver::extended_path()
{
  mod_file->addStatement(new ExtendedPathStatement(options_list));
  options_list.clear();
}

expr_t
ParsingDriver::add_model_equal(expr_t arg1, expr_t arg2)
{
  expr_t id = model_tree->AddEqual(arg1, arg2);

  // Detect if the equation is tagged [static]
  bool is_static_only = false;  
  for (vector<pair<string, string> >::const_iterator it = eq_tags.begin();
       it != eq_tags.end(); ++it)
    if (it->first == "static")
      {
        is_static_only = true;
        break;
      }

  if (is_static_only)
    {
      if (!id->isInStaticForm())
        error("An equation tagged [static] cannot contain leads, lags, expectations or STEADY_STATE operators");
      
      dynamic_model->addStaticOnlyEquation(id, location.begin.line);
    }
  else
    model_tree->addEquation(id, location.begin.line, eq_tags);

  eq_tags.clear();
  return id;
}

expr_t
ParsingDriver::add_model_equal_with_zero_rhs(expr_t arg)
{
  return add_model_equal(arg, model_tree->Zero);
}

void
ParsingDriver::declare_and_init_model_local_variable(string *name, expr_t rhs)
{
  int symb_id;
  try
    {
      symb_id = mod_file->symbol_table.addSymbol(*name, eModelLocalVariable);
    }
  catch (SymbolTable::AlreadyDeclaredException &e)
    {
      // It can have already been declared in a steady_state_model block, check that it is indeed a ModelLocalVariable
      symb_id = mod_file->symbol_table.getID(*name);
      if (mod_file->symbol_table.getType(symb_id) != eModelLocalVariable)
        error(*name + " has wrong type, you cannot use it within as left-hand side of a pound ('#') expression");
    }

  try
    {
      model_tree->AddLocalVariable(symb_id, rhs);
    }
  catch (DataTree::LocalVariableException &e)
    {
      error("Local model variable " + *name + " declared twice.");
    }
  delete name;
}

void
ParsingDriver::change_type(SymbolType new_type, vector<string *> *var_list)
{
  for (vector<string *>::iterator it = var_list->begin();
       it != var_list->end(); it++)
    {
      int id;
      try
        {
          id = mod_file->symbol_table.getID(**it);
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          error("Unknown variable " + **it);
        }

      // Check if symbol already used in a VariableNode
      if (mod_file->expressions_tree.isSymbolUsed(id)
          || mod_file->dynamic_model.isSymbolUsed(id))
        error("You cannot modify the type of symbol " + **it + " after having used it in an expression");

      mod_file->symbol_table.changeType(id, new_type);

      delete *it;
    }
  delete var_list;
}

expr_t
ParsingDriver::add_plus(expr_t arg1, expr_t arg2)
{
  return data_tree->AddPlus(arg1, arg2);
}

expr_t
ParsingDriver::add_minus(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMinus(arg1, arg2);
}

expr_t
ParsingDriver::add_uminus(expr_t arg1)
{
  return data_tree->AddUMinus(arg1);
}

expr_t
ParsingDriver::add_times(expr_t arg1, expr_t arg2)
{
  return data_tree->AddTimes(arg1, arg2);
}

expr_t
ParsingDriver::add_divide(expr_t arg1, expr_t arg2)
{
  return data_tree->AddDivide(arg1, arg2);
}

expr_t
ParsingDriver::add_less(expr_t arg1, expr_t arg2)
{
  return data_tree->AddLess(arg1, arg2);
}

expr_t
ParsingDriver::add_greater(expr_t arg1, expr_t arg2)
{
  return data_tree->AddGreater(arg1, arg2);
}

expr_t
ParsingDriver::add_less_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddLessEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_greater_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddGreaterEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_equal_equal(expr_t arg1, expr_t arg2)
{
  return data_tree->AddEqualEqual(arg1, arg2);
}

expr_t
ParsingDriver::add_different(expr_t arg1, expr_t arg2)
{
  return data_tree->AddDifferent(arg1, arg2);
}

expr_t
ParsingDriver::add_power(expr_t arg1, expr_t arg2)
{
  return data_tree->AddPower(arg1, arg2);
}

expr_t
ParsingDriver::add_expectation(string *arg1, expr_t arg2)
{
  expr_t expectationNode;
  expectationNode = data_tree->AddExpectation(atoi(arg1->c_str()), arg2);
  delete arg1;
  return expectationNode;
}

expr_t
ParsingDriver::add_exp(expr_t arg1)
{
  return data_tree->AddExp(arg1);
}

expr_t
ParsingDriver::add_log(expr_t arg1)
{
  return data_tree->AddLog(arg1);
}

expr_t
ParsingDriver::add_log10(expr_t arg1)
{
  return data_tree->AddLog10(arg1);
}

expr_t
ParsingDriver::add_cos(expr_t arg1)
{
  return data_tree->AddCos(arg1);
}

expr_t
ParsingDriver::add_sin(expr_t arg1)
{
  return data_tree->AddSin(arg1);
}

expr_t
ParsingDriver::add_tan(expr_t arg1)
{
  return data_tree->AddTan(arg1);
}

expr_t
ParsingDriver::add_acos(expr_t arg1)
{
  return data_tree->AddAcos(arg1);
}

expr_t
ParsingDriver::add_asin(expr_t arg1)
{
  return data_tree->AddAsin(arg1);
}

expr_t
ParsingDriver::add_atan(expr_t arg1)
{
  return data_tree->AddAtan(arg1);
}

expr_t
ParsingDriver::add_cosh(expr_t arg1)
{
  return data_tree->AddCosh(arg1);
}

expr_t
ParsingDriver::add_sinh(expr_t arg1)
{
  return data_tree->AddSinh(arg1);
}

expr_t
ParsingDriver::add_tanh(expr_t arg1)
{
  return data_tree->AddTanh(arg1);
}

expr_t
ParsingDriver::add_acosh(expr_t arg1)
{
  return data_tree->AddAcosh(arg1);
}

expr_t
ParsingDriver::add_asinh(expr_t arg1)
{
  return data_tree->AddAsinh(arg1);
}

expr_t
ParsingDriver::add_atanh(expr_t arg1)
{
  return data_tree->AddAtanh(arg1);
}

expr_t
ParsingDriver::add_sqrt(expr_t arg1)
{
  return data_tree->AddSqrt(arg1);
}

expr_t
ParsingDriver::add_abs(expr_t arg1)
{
  return data_tree->AddAbs(arg1);
}

expr_t
ParsingDriver::add_sign(expr_t arg1)
{
  return data_tree->AddSign(arg1);
}

expr_t
ParsingDriver::add_max(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMax(arg1, arg2);
}

expr_t
ParsingDriver::add_min(expr_t arg1, expr_t arg2)
{
  return data_tree->AddMin(arg1, arg2);
}

expr_t
ParsingDriver::add_normcdf(expr_t arg1, expr_t arg2, expr_t arg3)
{
  return data_tree->AddNormcdf(arg1, arg2, arg3);
}

expr_t
ParsingDriver::add_normcdf(expr_t arg)
{
  return add_normcdf(arg, data_tree->Zero, data_tree->One);
}

expr_t
ParsingDriver::add_normpdf(expr_t arg1, expr_t arg2, expr_t arg3)
{
  return data_tree->AddNormpdf(arg1, arg2, arg3);
}

expr_t
ParsingDriver::add_normpdf(expr_t arg)
{
  return add_normpdf(arg, data_tree->Zero, data_tree->One);
}

expr_t
ParsingDriver::add_erf(expr_t arg1)
{
  return data_tree->AddErf(arg1);
}

expr_t
ParsingDriver::add_steady_state(expr_t arg1)
{
  return data_tree->AddSteadyState(arg1);
}

void
ParsingDriver::external_function_option(const string &name_option, string *opt)
{
  external_function_option(name_option, *opt);
  delete opt;
}

void
ParsingDriver::external_function_option(const string &name_option, const string &opt)
{
  if (name_option == "name")
    {
      if (opt.empty())
        error("An argument must be passed to the 'name' option of the external_function() statement.");
      declare_symbol(&opt, eExternalFunction, NULL, NULL);
      current_external_function_id = mod_file->symbol_table.getID(opt);
    }
  else if (name_option == "first_deriv_provided")
    {
      if (opt.empty())
        current_external_function_options.firstDerivSymbID = eExtFunSetButNoNameProvided;
      else
        {
          declare_symbol(&opt, eExternalFunction, NULL, NULL);
          current_external_function_options.firstDerivSymbID = mod_file->symbol_table.getID(opt);
        }
    }
  else if (name_option == "second_deriv_provided")
    {
      if (opt.empty())
        current_external_function_options.secondDerivSymbID = eExtFunSetButNoNameProvided;
      else
        {
          declare_symbol(&opt, eExternalFunction, NULL, NULL);
          current_external_function_options.secondDerivSymbID = mod_file->symbol_table.getID(opt);
        }
    }
  else if (name_option == "nargs")
    current_external_function_options.nargs = atoi(opt.c_str());
  else
    error("Unexpected error in ParsingDriver::external_function_option(): Please inform Dynare Team.");
}

void
ParsingDriver::external_function()
{
  if (current_external_function_id == eExtFunNotSet)
    error("The 'name' option must be passed to external_function().");

  if (current_external_function_options.secondDerivSymbID >= 0
      && current_external_function_options.firstDerivSymbID  == eExtFunNotSet)
    error("If the second derivative is provided to the external_function command, the first derivative must also be provided.");

  if (current_external_function_options.secondDerivSymbID == eExtFunSetButNoNameProvided
      && current_external_function_options.firstDerivSymbID  != eExtFunSetButNoNameProvided)
    error("If the second derivative is provided in the top-level function, the first derivative must also be provided in that function.");

  mod_file->external_functions_table.addExternalFunction(current_external_function_id, current_external_function_options, true);
  reset_current_external_function_options();
}

void
ParsingDriver::push_external_function_arg_vector_onto_stack()
{
  vector<expr_t> emptyvec;
  stack_external_function_args.push(emptyvec);
}

void
ParsingDriver::add_external_function_arg(expr_t arg)
{
  stack_external_function_args.top().push_back(arg);
}

expr_t
ParsingDriver::add_model_var_or_external_function(string *function_name, bool in_model_block)
{
  expr_t nid;
  if (mod_file->symbol_table.exists(*function_name))
    {
      if (mod_file->symbol_table.getType(*function_name) != eExternalFunction)
        {
          if (!in_model_block)
            {
              if (stack_external_function_args.top().size() > 0)
                error(string("Symbol ") + *function_name + string(" cannot take arguments."));
              else
                return add_expression_variable(function_name);
            }
          else
            { // e.g. model_var(lag) => ADD MODEL VARIABLE WITH LEAD (NumConstNode)/LAG (UnaryOpNode)
              if (stack_external_function_args.top().size() != 1)
                error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., has received more than one argument)."));

              NumConstNode *numNode = dynamic_cast<NumConstNode *>(stack_external_function_args.top().front());
              UnaryOpNode *unaryNode = dynamic_cast<UnaryOpNode *>(stack_external_function_args.top().front());

              if (numNode == NULL && unaryNode == NULL)
                error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., takes an argument that is not an integer)."));

              eval_context_t ectmp;
              double model_var_arg;
              if (unaryNode == NULL)
                {
                  try
                    {
                      model_var_arg = numNode->eval(ectmp);
                    }
                  catch (ExprNode::EvalException &e)
                    {
                      error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., takes an argument that is not an integer)."));
                    }
                }
              else
                if (unaryNode->get_op_code() != oUminus)
                  error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., takes an argument that is not an integer)."));
                else
                  {
                    try
                      {
                        model_var_arg = unaryNode->eval(ectmp);
                      }
                    catch (ExprNode::EvalException &e)
                      {
                        error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., takes an argument that is not an integer)."));
                      }
                  }

              if (model_var_arg != floor(model_var_arg))
                error(string("Symbol ") + *function_name + string(" is being treated as if it were a function (i.e., takes an argument that is not an integer)."));

              nid = add_model_variable(mod_file->symbol_table.getID(*function_name), (int) model_var_arg);
              stack_external_function_args.pop();
              delete function_name;
              return nid;
            }
        }
      else
        { // e.g. this function has already been referenced (either ad hoc or through the external_function() statement
          // => check that the information matches previously declared info
          int symb_id = mod_file->symbol_table.getID(*function_name);
          assert(mod_file->external_functions_table.exists(symb_id));

          if (in_model_block)
            if (mod_file->external_functions_table.getNargs(symb_id) == eExtFunNotSet)
              error("Before using " + *function_name
                    +"() in the model block, you must first declare it via the external_function() statement");
            else if ((int) (stack_external_function_args.top().size()) != mod_file->external_functions_table.getNargs(symb_id))
              error("The number of arguments passed to " + *function_name
                    +"() does not match those of a previous call or declaration of this function.");
        }
    }
  else
    { //First time encountering this external function i.e., not previously declared or encountered
      if (in_model_block)
        error("To use an external function (" + *function_name + ") within the model block, you must first declare it via the external_function() statement.");

      declare_symbol(function_name, eExternalFunction, NULL, NULL);
      current_external_function_options.nargs = stack_external_function_args.top().size();
      mod_file->external_functions_table.addExternalFunction(mod_file->symbol_table.getID(*function_name),
                                                             current_external_function_options, in_model_block);
      reset_current_external_function_options();
    }

  //By this point, we're sure that this function exists in the External Functions Table and is not a mod var
  int symb_id = mod_file->symbol_table.getID(*function_name);
  nid = data_tree->AddExternalFunction(symb_id, stack_external_function_args.top());
  stack_external_function_args.pop();
  delete function_name;
  return nid;
}

void
ParsingDriver::add_native(const string &s)
{
  mod_file->addStatement(new NativeStatement(s));
}

void
ParsingDriver::add_native_remove_charset(const char *s, const string &token)
{
  string str = string(s);
  size_t found = str.find(token);

  assert(found != string::npos);
  str.resize(found);
  add_native(str);
}

void
ParsingDriver::begin_steady_state_model()
{
  set_current_data_tree(&mod_file->steady_state_model);
}

void
ParsingDriver::add_steady_state_model_equal(string *varname, expr_t expr)
{
  int id;
  try
    {
      id = mod_file->symbol_table.getID(*varname);
    }
  catch (SymbolTable::UnknownSymbolNameException &e)
    {
      // Unknown symbol, declare it as a ModFileLocalVariable
      id = mod_file->symbol_table.addSymbol(*varname, eModFileLocalVariable);
    }

  SymbolType type = mod_file->symbol_table.getType(id);
  if (type != eEndogenous && type != eModFileLocalVariable && type != eParameter)
    error(*varname + " has incorrect type");

  mod_file->steady_state_model.addDefinition(id, expr);

  delete varname;
}

void
ParsingDriver::add_steady_state_model_equal_multiple(expr_t expr)
{
  const vector<string> &symbs = symbol_list.get_symbols();
  vector<int> ids;

  for (size_t i = 0; i < symbs.size(); i++)
    {
      int id;
      try
        {
          id = mod_file->symbol_table.getID(symbs[i]);
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          // Unknown symbol, declare it as a ModFileLocalVariable
          id = mod_file->symbol_table.addSymbol(symbs[i], eModFileLocalVariable);
        }
      SymbolType type = mod_file->symbol_table.getType(id);
      if (type != eEndogenous && type != eModFileLocalVariable && type != eParameter)
        error(symbs[i] + " has incorrect type");
      ids.push_back(id);
    }

  mod_file->steady_state_model.addMultipleDefinitions(ids, expr);

  symbol_list.clear();
}

void
ParsingDriver::add_graph_format(const string &name)
{
  graph_formats.addSymbol(name);
}

void
ParsingDriver::process_graph_format_option()
{
  options_list.symbol_list_options["graph_format"] = graph_formats;
  graph_formats.clear();
}

void
ParsingDriver::model_diagnostics()
{
  mod_file->addStatement(new ModelDiagnosticsStatement());
}

void
ParsingDriver::add_parallel_local_file(string *filename)
{
  mod_file->parallel_local_files.push_back(*filename);
  delete filename;
}

void
ParsingDriver::add_moment_calibration_item(string *endo1, string *endo2, string *lags, vector<string *> *range)
{
  MomentCalibration::Constraint c;

  check_symbol_existence(*endo1);
  c.endo1 = mod_file->symbol_table.getID(*endo1);
  if (mod_file->symbol_table.getType(*endo1) != eEndogenous)
    error("Variable " + *endo1 + " is not an endogenous.");
  delete endo1;

  check_symbol_existence(*endo2);
  c.endo2 = mod_file->symbol_table.getID(*endo2);
  if (mod_file->symbol_table.getType(*endo2) != eEndogenous)
    error("Variable " + *endo2 + " is not an endogenous.");
  delete endo2;

  c.lags = *lags;
  delete lags;
  
  assert(range->size() == 2);
  c.lower_bound = *((*range)[0]);
  c.upper_bound = *((*range)[1]);
  delete (*range)[0];
  delete (*range)[1];
  delete range;
  
  moment_calibration_constraints.push_back(c);
}

void ParsingDriver::end_moment_calibration()
{
  mod_file->addStatement(new MomentCalibration(moment_calibration_constraints,
                                               mod_file->symbol_table));
  moment_calibration_constraints.clear();
}

void
ParsingDriver::add_irf_calibration_item(string *endo, string *periods, string *exo, vector<string *> *range)
{
  IrfCalibration::Constraint c;

  check_symbol_existence(*endo);
  c.endo = mod_file->symbol_table.getID(*endo);
  if (mod_file->symbol_table.getType(*endo) != eEndogenous)
    error("Variable " + *endo + " is not an endogenous.");
  delete endo;

  c.periods = *periods;
  delete periods;

  check_symbol_existence(*exo);
  c.exo = mod_file->symbol_table.getID(*exo);
  if (mod_file->symbol_table.getType(*exo) != eExogenous)
    error("Variable " + *endo + " is not an exogenous.");
  delete exo;
  
  assert(range->size() == 2);
  c.lower_bound = *((*range)[0]);
  c.upper_bound = *((*range)[1]);
  delete (*range)[0];
  delete (*range)[1];
  delete range;
  
  irf_calibration_constraints.push_back(c);
}

void ParsingDriver::end_irf_calibration()
{
  mod_file->addStatement(new IrfCalibration(irf_calibration_constraints,
                                            mod_file->symbol_table));
  irf_calibration_constraints.clear();
}

void
ParsingDriver::smoother2histval()
{
  mod_file->addStatement(new Smoother2histvalStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::histval_file(string *filename)
{
  mod_file->addStatement(new HistvalFileStatement(*filename));
  delete filename;
}


void
ParsingDriver::perfect_foresight_setup()
{
  mod_file->addStatement(new PerfectForesightSetupStatement(options_list));
  options_list.clear();
}

void
ParsingDriver::perfect_foresight_solver()
{
  mod_file->addStatement(new PerfectForesightSolverStatement(options_list));
  options_list.clear();
}
