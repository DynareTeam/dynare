/*
 * Copyright (C) 2003-2013 Dynare Team
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
#include <cassert>
#include <iostream>
#include <sstream>

using namespace std;

#include "ComputingTasks.hh"
#include "Statement.hh"

#include <boost/algorithm/string/trim.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>

SteadyStatement::SteadyStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SteadyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.steady_present = true;
}

void
SteadyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "steady;\n";
}

CheckStatement::CheckStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
CheckStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "oo_.dr.eigval = check(M_,options_,oo_);" << endl;
}

void
CheckStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.check_present = true;
}

ModelInfoStatement::ModelInfoStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
ModelInfoStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  //mod_file_struct.model_info_present = true;
}

void
ModelInfoStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "model_info();\n";
}

SimulStatement::SimulStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SimulStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.simul_present = true;

  // The following is necessary to allow shocks+endval+simul in a loop
  mod_file_struct.shocks_present_but_simul_not_yet = false;
}

void
SimulStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "simul();\n";
}

StochSimulStatement::StochSimulStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
StochSimulStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.stoch_simul_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, atoi(it->second.c_str()));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;

}

void
StochSimulStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = stoch_simul(var_list_);" << endl;
}

ForecastStatement::ForecastStatement(const SymbolList &symbol_list_arg,
                                     const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
ForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "info = dyn_forecast(var_list_,'simul');" << endl;
}

RamseyPolicyStatement::RamseyPolicyStatement(const SymbolList &symbol_list_arg,
                                             const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
RamseyPolicyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.ramsey_policy_present = true;

  /* Fill in option_order of mod_file_struct
     Since ramsey policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = atoi(it->second.c_str());
      if (order > 1)
        {
          cerr << "ERROR: ramsey_policy: order > 1 is not yet implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;
}

void
RamseyPolicyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "ramsey_policy(var_list_);\n";
}

DiscretionaryPolicyStatement::DiscretionaryPolicyStatement(const SymbolList &symbol_list_arg,
							   const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
DiscretionaryPolicyStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.discretionary_policy_present = true;

  if (options_list.symbol_list_options.find("instruments") == options_list.symbol_list_options.end())
    {
      cerr << "ERROR: discretionary_policy: the instruments option is required." << endl;
      exit(EXIT_FAILURE);
    }

  /* Fill in option_order of mod_file_struct
     Since discretionary policy needs one further order of derivation (for example, for 1st order
     approximation, it needs 2nd derivatives), we add 1 to the order declared by user */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = atoi(it->second.c_str());
      if (order > 1)
        {
          cerr << "ERROR: discretionary_policy: order > 1 is not yet implemented" << endl;
          exit(EXIT_FAILURE);
        }
      mod_file_struct.order_option = max(mod_file_struct.order_option, order + 1);
    }

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;
}

void
DiscretionaryPolicyStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "discretionary_policy(var_list_);\n";
}

EstimationStatement::EstimationStatement(const SymbolList &symbol_list_arg,
                                         const OptionsList &options_list_arg,
                                         const SymbolTable &symbol_table_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    {
      int order = atoi(it->second.c_str());
          
      if (order > 2)
        {
          cerr << "ERROR: order > 2 is not supported in estimation" << endl;
          exit(EXIT_FAILURE);
        }
      
      mod_file_struct.order_option = max(mod_file_struct.order_option, order);
    }
  
  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Fill in mod_file_struct.estimation_analytic_derivation
  it = options_list.num_options.find("analytic_derivation");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.estimation_analytic_derivation = true;

  it = options_list.num_options.find("dsge_var");
  if (it != options_list.num_options.end())
    // Ensure that irf_shocks & dsge_var have not both been passed
    if (options_list.symbol_list_options.find("irf_shocks") != options_list.symbol_list_options.end())
      {
        cerr << "The irf_shocks and dsge_var options may not both be passed to estimation." << endl;
        exit(EXIT_FAILURE);
      }
    else
      // Fill in mod_file_struct.dsge_var_calibrated
      mod_file_struct.dsge_var_calibrated = it->second;

  // Fill in mod_file_struct.dsge_var_estimated
  OptionsList::string_options_t::const_iterator it_str = options_list.string_options.find("dsge_var");
  if (it_str != options_list.string_options.end())
    mod_file_struct.dsge_var_estimated = true;

  // Fill in mod_file_struct.bayesian_irf_present
  it = options_list.num_options.find("bayesian_irf");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.bayesian_irf_present = true;

  it = options_list.num_options.find("dsge_varlag");
  if (it != options_list.num_options.end())
    if (mod_file_struct.dsge_var_calibrated.empty()
        && !mod_file_struct.dsge_var_estimated)
      {
        cerr << "ERROR: The estimation statement requires a dsge_var option to be passed "
             << "if the dsge_varlag option is passed." << endl;
        exit(EXIT_FAILURE);
      }

  if (!mod_file_struct.dsge_var_calibrated.empty()
      && mod_file_struct.dsge_var_estimated)
    {
      cerr << "ERROR: An estimation statement cannot take more than one dsge_var option." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.string_options.find("datafile") == options_list.string_options.end() &&
      !mod_file_struct.estimation_data_statement_present)
    {
      cerr << "ERROR: The estimation statement requires a data file to be supplied via the datafile option." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);

  // Special treatment for order option and particle filter
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it == options_list.num_options.end())
    output << "options_.order = 1;" << endl;
  else if (atoi(it->second.c_str()) == 2)
    output << "options_.particle.status = 1;" << endl;

  // Do not check for the steady state in diffuse filter mode (#400)
  it = options_list.num_options.find("diffuse_filter");
  if (it != options_list.num_options.end() && it->second == "1")
    output << "options_.steadystate.nocheck = 1;" << endl;

  symbol_list.writeOutput("var_list_", output);
  output << "dynare_estimation(var_list_);\n";
}

DynareSensitivityStatement::DynareSensitivityStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
DynareSensitivityStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("identification");
  if (it != options_list.num_options.end()
      && it->second == "1")
    mod_file_struct.identification_present = true;
}

void
DynareSensitivityStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_gsa");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("nodisplay");
  if (it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  it = options_list.num_options.find("nograph");
  if (it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  OptionsList::string_options_t::const_iterator it2 = options_list.string_options.find("graph_format");
  if (it2 != options_list.string_options.end())
    output << "options_.graph_format = '" << it2->second << "';" << endl;
  
  output << "dynare_sensitivity(options_gsa);" << endl;
}

RplotStatement::RplotStatement(const SymbolList &symbol_list_arg,
                               const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
RplotStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "rplot(var_list_);\n";
}

UnitRootVarsStatement::UnitRootVarsStatement(void)
{
}

void
UnitRootVarsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.diffuse_filter = 1;" << endl
	 << "options_.steadystate.nocheck = 1;" << endl;
}

PeriodsStatement::PeriodsStatement(int periods_arg) : periods(periods_arg)
{
}

void
PeriodsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.periods = " << periods << ";" << endl;
}

DsampleStatement::DsampleStatement(int val1_arg) : val1(val1_arg), val2(-1)
{
}

DsampleStatement::DsampleStatement(int val1_arg, int val2_arg) : val1(val1_arg), val2(val2_arg)
{
}

void
DsampleStatement::writeOutput(ostream &output, const string &basename) const
{
  if (val2 < 0)
    output << "dsample(" << val1 << ");" << endl;
  else
    output << "dsample(" << val1 << ", " << val2 << ");" << endl;
}

EstimatedParamsStatement::EstimatedParamsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                   const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  for (vector<EstimationParams>::const_iterator it = estim_params_list.begin();
       it != estim_params_list.end(); it++)
    {
      if (it->name == "dsge_prior_weight")
        mod_file_struct.dsge_prior_weight_in_estimated_params = true;

      // Handle case of degenerate beta prior
      if (it->prior == eBeta)
        try
          {
            if (it->mean->eval(eval_context_t()) == 0.5
                && it->std->eval(eval_context_t()) == 0.5)
              {
                cerr << "ERROR: The prior density is not defined for the beta distribution when the mean = standard deviation = 0.5." << endl;
                exit(EXIT_FAILURE);
              }
          }
        catch (ExprNode::EvalException &e)
          {
            // We don't have enough information to compute the numerical value, skip the test
          }
    }

  // Check that no parameter/endogenous is declared twice in the block
  set<string> already_declared;
  set<pair<string, string> > already_declared_corr;
  for (vector<EstimationParams>::const_iterator it = estim_params_list.begin();
       it != estim_params_list.end(); it++)
    {
      if (it->type == 3) // Correlation
        {
          // Use lexical ordering for the pair of symbols
          pair<string, string> x = it->name < it->name2 ? make_pair(it->name, it->name2) : make_pair(it->name2, it->name);

          if (already_declared_corr.find(x) == already_declared_corr.end())
            already_declared_corr.insert(x);
          else
            {
              cerr << "ERROR: in `estimated_params' block, the correlation between " << it->name << " and " << it->name2 << " is declared twice." << endl;
              exit(EXIT_FAILURE);
            }
        }
      else
        {
          if (already_declared.find(it->name) == already_declared.end())
            already_declared.insert(it->name);
          else
            {
              cerr << "ERROR: in `estimated_params' block, the symbol " << it->name << " is declared twice." << endl;
              exit(EXIT_FAILURE);
            }
        }
    }
}

void
EstimatedParamsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "global estim_params_" << endl
         << "estim_params_.var_exo = [];" << endl
         << "estim_params_.var_endo = [];" << endl
         << "estim_params_.corrx = [];" << endl
         << "estim_params_.corrn = [];" << endl
         << "estim_params_.param_vals = [];" << endl;

  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      switch (it->type)
        {
        case 1:
          if (symb_type == eExogenous)
            output << "estim_params_.var_exo = [estim_params_.var_exo; ";
          else if (symb_type == eEndogenous)
            output << "estim_params_.var_endo = [estim_params_.var_endo; ";
          output << symb_id;
          break;
        case 2:
          output << "estim_params_.param_vals = [estim_params_.param_vals; "
                 << symb_id;
          break;
        case 3:
          if (symb_type == eExogenous)
            output << "estim_params_.corrx = [estim_params_.corrx; ";
          else if (symb_type == eEndogenous)
            output << "estim_params_.corrn = [estim_params_.corrn; ";
          output << symb_id << " " << symbol_table.getTypeSpecificID(it->name2)+1;
          break;
        }
      output << ", ";
      it->init_val->writeOutput(output);
      output << ", ";
      it->low_bound->writeOutput(output);
      output << ", ";
      it->up_bound->writeOutput(output);
      output << ", "
             << it->prior << ", ";
      it->mean->writeOutput(output);
      output << ", ";
      it->std->writeOutput(output);
      output << ", ";
      it->p3->writeOutput(output);
      output << ", ";
      it->p4->writeOutput(output);
      output << ", ";
      it->jscale->writeOutput(output);
      output << " ];" << endl;
    }
}

EstimatedParamsInitStatement::EstimatedParamsInitStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                           const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsInitStatement::writeOutput(ostream &output, const string &basename) const
{
  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      if (it->type < 3)
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_exo(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.var_endo(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;
              output << "estim_params_.param_vals(tmp1,2) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << " & estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ") | "
                     <<             "(estim_params_.corrx(:,2)==" << symb_id << " & estim_params_.corrx(:,1)==" << symbol_table.getTypeSpecificID(it->name2)+1 << "));" << endl;
              output << "estim_params_.corrx(tmp1,3) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << " & estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ") | "
                     <<             "(estim_params_.corrn(:,2)==" << symb_id << " & estim_params_.corrn(:,1)==" << symbol_table.getTypeSpecificID(it->name2)+1 << "));" << endl;
              output << "estim_params_.corrn(tmp1,3) = ";
              it->init_val->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

EstimatedParamsBoundsStatement::EstimatedParamsBoundsStatement(const vector<EstimationParams> &estim_params_list_arg,
                                                               const SymbolTable &symbol_table_arg) :
  estim_params_list(estim_params_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
EstimatedParamsBoundsStatement::writeOutput(ostream &output, const string &basename) const
{
  vector<EstimationParams>::const_iterator it;

  for (it = estim_params_list.begin(); it != estim_params_list.end(); it++)
    {
      int symb_id = symbol_table.getTypeSpecificID(it->name) + 1;
      SymbolType symb_type = symbol_table.getType(it->name);

      if (it->type < 3)
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find(estim_params_.var_exo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_exo(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_exo(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find(estim_params_.var_endo(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.var_endo(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.var_endo(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eParameter)
            {
              output << "tmp1 = find(estim_params_.param_vals(:,1)==" << symb_id << ");" << endl;

              output << "estim_params_.param_vals(tmp1,3) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.param_vals(tmp1,4) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
      else
        {
          if (symb_type == eExogenous)
            {
              output << "tmp1 = find((estim_params_.corrx(:,1)==" << symb_id << " & estim_params_.corrx(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ") | "
                     <<             "(estim_params_.corrx(:,2)==" << symb_id << " & estim_params_.corrx(:,1)==" << symbol_table.getTypeSpecificID(it->name2)+1 << "));" << endl;

              output << "estim_params_.corrx(tmp1,4) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrx(tmp1,5) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
          else if (symb_type == eEndogenous)
            {
              output << "tmp1 = find((estim_params_.corrn(:,1)==" << symb_id << " & estim_params_.corrn(:,2)==" << symbol_table.getTypeSpecificID(it->name2)+1 << ") | "
                     <<             "(estim_params_.corrn(:,2)==" << symb_id << " & estim_params_.corrn(:,1)==" << symbol_table.getTypeSpecificID(it->name2)+1 << "));" << endl;

              output << "estim_params_.corrn(tmp1,4) = ";
              it->low_bound->writeOutput(output);
              output << ";" << endl;

              output << "estim_params_.corrn(tmp1,5) = ";
              it->up_bound->writeOutput(output);
              output << ";" << endl;
            }
        }
    }
}

ObservationTrendsStatement::ObservationTrendsStatement(const trend_elements_t &trend_elements_arg,
                                                       const SymbolTable &symbol_table_arg) :
  trend_elements(trend_elements_arg),
  symbol_table(symbol_table_arg)
{
}

void
ObservationTrendsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_.trend_coeff_ = {};" << endl;

  trend_elements_t::const_iterator it;

  for (it = trend_elements.begin(); it != trend_elements.end(); it++)
    {
      SymbolType type = symbol_table.getType(it->first);
      if (type == eEndogenous)
        {
          output << "tmp1 = strmatch('" << it->first << "',options_.varobs,'exact');\n";
          output << "options_.trend_coeffs{tmp1} = '";
          it->second->writeOutput(output);
          output << "';" << endl;
        }
      else
        cout << "Error : Non-variable symbol used in TREND_COEFF: " << it->first << endl;
    }
}

OsrParamsStatement::OsrParamsStatement(const SymbolList &symbol_list_arg) :
  symbol_list(symbol_list_arg)
{
}

void
OsrParamsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.osr_params_present = true;
}

void
OsrParamsStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("osr_params_", output);
}

OsrStatement::OsrStatement(const SymbolList &symbol_list_arg,
                           const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
OsrStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.osr_present = true;

  // Fill in option_order of mod_file_struct
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("order");
  if (it != options_list.num_options.end())
    mod_file_struct.order_option = max(mod_file_struct.order_option, atoi(it->second.c_str()));

  // Fill in mod_file_struct.partial_information
  it = options_list.num_options.find("partial_information");
  if (it != options_list.num_options.end() && it->second == "1")
    mod_file_struct.partial_information = true;

  // Option k_order_solver (implicit when order >= 3)
  it = options_list.num_options.find("k_order_solver");
  if ((it != options_list.num_options.end() && it->second == "1")
      || mod_file_struct.order_option >= 3)
    mod_file_struct.k_order_solver = true;
}

void
OsrStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "osr(var_list_,osr_params_,obj_var_,optim_weights_);\n";
}

OptimWeightsStatement::OptimWeightsStatement(const var_weights_t &var_weights_arg,
                                             const covar_weights_t &covar_weights_arg,
                                             const SymbolTable &symbol_table_arg) :
  var_weights(var_weights_arg),
  covar_weights(covar_weights_arg),
  symbol_table(symbol_table_arg)
{
}

void
OptimWeightsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.optim_weights_present = true;
}

void
OptimWeightsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "%" << endl
         << "% OPTIM_WEIGHTS" << endl
         << "%" << endl
         << "optim_weights_ = sparse(M_.endo_nbr,M_.endo_nbr);" << endl
         << "obj_var_ = [];" << endl << endl;

  for (var_weights_t::const_iterator it = var_weights.begin();
       it != var_weights.end(); it++)
    {
      const string &name = it->first;
      const expr_t value = it->second;
      int id = symbol_table.getTypeSpecificID(name) + 1;
      output << "optim_weights_(" << id << "," << id << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id << "];\n";
    }

  for (covar_weights_t::const_iterator it = covar_weights.begin();
       it != covar_weights.end(); it++)
    {
      const string &name1 = it->first.first;
      const string &name2 = it->first.second;
      const expr_t value = it->second;
      int id1 = symbol_table.getTypeSpecificID(name1) + 1;
      int id2 = symbol_table.getTypeSpecificID(name2) + 1;
      output << "optim_weights_(" << id1 << "," << id2 << ") = ";
      value->writeOutput(output);
      output << ";" << endl;
      output << "obj_var_ = [obj_var_; " << id1 << "; " << id2 << "];\n";
    }
}

DynaSaveStatement::DynaSaveStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg)
{
}

void
DynaSaveStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynasave('" << filename
         << "',var_list_);" << endl;
}

DynaTypeStatement::DynaTypeStatement(const SymbolList &symbol_list_arg,
                                     const string &filename_arg) :
  symbol_list(symbol_list_arg),
  filename(filename_arg)
{
}

void
DynaTypeStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  output << "dynatype('" << filename
         << "',var_list_);" << endl;
}

ModelComparisonStatement::ModelComparisonStatement(const filename_list_t &filename_list_arg,
                                                   const OptionsList &options_list_arg) :
  filename_list(filename_list_arg),
  options_list(options_list_arg)
{
}

void
ModelComparisonStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);

  output << "ModelNames_ = {};" << endl;
  output << "ModelPriors_ = [];" << endl;

  for (filename_list_t::const_iterator it = filename_list.begin();
       it != filename_list.end(); it++)
    {
      output << "ModelNames_ = { ModelNames_{:} '" << (*it).first << "'};" << endl;
      output << "ModelPriors_ = [ ModelPriors_ ; " << (*it).second << "];" << endl;
    }
  output << "model_comparison(ModelNames_,ModelPriors_,oo_,options_,M_.fname);" << endl;
}

PlannerObjectiveStatement::PlannerObjectiveStatement(StaticModel *model_tree_arg) :
  model_tree(model_tree_arg)
{
}

PlannerObjectiveStatement::~PlannerObjectiveStatement()
{
  delete model_tree;
}

void
PlannerObjectiveStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  assert(model_tree->equation_number() == 1);
  mod_file_struct.planner_objective_present = true;
}

StaticModel *
PlannerObjectiveStatement::getPlannerObjective() const
{
  return model_tree;
}

void
PlannerObjectiveStatement::computingPass()
{
  model_tree->computingPass(eval_context_t(), false, true, false, false, false);
}

void
PlannerObjectiveStatement::writeOutput(ostream &output, const string &basename) const
{
  model_tree->writeStaticFile(basename + "_objective", false, false, false);
}

BVARDensityStatement::BVARDensityStatement(int maxnlags_arg, const OptionsList &options_list_arg) :
  maxnlags(maxnlags_arg),
  options_list(options_list_arg)
{
}

void
BVARDensityStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
BVARDensityStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "bvar_density(" << maxnlags << ");" << endl;
}

BVARForecastStatement::BVARForecastStatement(int nlags_arg, const OptionsList &options_list_arg) :
  nlags(nlags_arg),
  options_list(options_list_arg)
{
}

void
BVARForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
BVARForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "bvar_forecast(" << nlags << ");" << endl;
}

SBVARStatement::SBVARStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SBVARStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
SBVARStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  output << "sbvar(M_,options_);" << endl;
}

MSSBVAREstimationStatement::MSSBVAREstimationStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVAREstimationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.create_init") == options_list.num_options.end())
    if (options_list.string_options.find("datafile") == options_list.string_options.end() ||
        options_list.num_options.find("ms.initial_year") == options_list.num_options.end())
      {
        cerr << "ERROR: If you do not pass no_create_init to ms_estimation, "
             << "you must pass the datafile and initial_year options." << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVAREstimationStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl
         << "options_.datafile = '';" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_estimation(M_, options_, oo_);" << endl;
}

MSSBVARSimulationStatement::MSSBVARSimulationStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARSimulationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARSimulationStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);

  // Redeclare drop option if necessary
  OptionsList::num_options_t::const_iterator mh_replic_it = options_list.num_options.find("ms.mh_replic");
  OptionsList::num_options_t::const_iterator thinning_factor_it = options_list.num_options.find("ms.thinning_factor");
  OptionsList::num_options_t::const_iterator drop_it = options_list.num_options.find("ms.drop");
  if (mh_replic_it != options_list.num_options.end() || thinning_factor_it != options_list.num_options.end())
    if (drop_it == options_list.num_options.end())
      output << "options_.ms.drop = 0.1*options_.ms.mh_replic*options_.ms.thinning_factor;" << endl;

  output << "[options_, oo_] = ms_simulation(M_, options_, oo_);" << endl;
}

MSSBVARComputeMDDStatement::MSSBVARComputeMDDStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARComputeMDDStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;
}

void
MSSBVARComputeMDDStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_mdd(M_, options_, oo_);" << endl;
}

MSSBVARComputeProbabilitiesStatement::MSSBVARComputeProbabilitiesStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARComputeProbabilitiesStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.real_time_smoothed_probabilities") != options_list.num_options.end())
    if (options_list.num_options.find("ms.filtered_probabilities") != options_list.num_options.end())
      {
        cerr << "ERROR: You may only pass one of real_time_smoothed "
             << "and filtered_probabilities to ms_compute_probabilities." << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVARComputeProbabilitiesStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_compute_probabilities(M_, options_, oo_);" << endl;
}

MSSBVARIrfStatement::MSSBVARIrfStatement(const SymbolList &symbol_list_arg,
					 const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
MSSBVARIrfStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = false;
  bool regimes_present = false;
  bool filtered_probabilities_present = false;

  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("ms.regimes");
  if (it != options_list.num_options.end())
    regimes_present = true;

  it = options_list.num_options.find("ms.regime");
  if (it != options_list.num_options.end())
    regime_present = true;

  it = options_list.num_options.find("ms.filtered_probabilities");
  if (it != options_list.num_options.end())
    filtered_probabilities_present = true;

  if ((filtered_probabilities_present && regime_present) ||
      (filtered_probabilities_present && regimes_present) ||
      (regimes_present && regime_present))
      {
        cerr << "ERROR: You may only pass one of regime, regimes and "
             << "filtered_probabilities to ms_irf" << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVARIrfStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  symbol_list.writeOutput("var_list_", output);
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_irf(var_list_,M_, options_, oo_);" << endl;
}

MSSBVARForecastStatement::MSSBVARForecastStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARForecastStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  if (options_list.num_options.find("ms.regimes") != options_list.num_options.end())
    if (options_list.num_options.find("ms.regime") != options_list.num_options.end())
      {
        cerr << "ERROR: You may only pass one of regime and regimes to ms_forecast" << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVARForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_forecast(M_, options_, oo_);" << endl;
}

MSSBVARVarianceDecompositionStatement::MSSBVARVarianceDecompositionStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
MSSBVARVarianceDecompositionStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.bvar_present = true;

  bool regime_present = false;
  bool regimes_present = false;
  bool filtered_probabilities_present = false;

  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("ms.regimes");
  if (it != options_list.num_options.end())
    regimes_present = true;

  it = options_list.num_options.find("ms.regime");
  if (it != options_list.num_options.end())
    regime_present = true;

  it = options_list.num_options.find("ms.filtered_probabilities");
  if (it != options_list.num_options.end())
    filtered_probabilities_present = true;

  if ((filtered_probabilities_present && regime_present) ||
      (filtered_probabilities_present && regimes_present) ||
      (regimes_present && regime_present))
      {
        cerr << "ERROR: You may only pass one of regime, regimes and "
             << "filtered_probabilities to ms_variance_decomposition" << endl;
        exit(EXIT_FAILURE);
      }
}

void
MSSBVARVarianceDecompositionStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "options_ = initialize_ms_sbvar_options(M_, options_);" << endl;
  options_list.writeOutput(output);
  output << "[options_, oo_] = ms_variance_decomposition(M_, options_, oo_);" << endl;
}

IdentificationStatement::IdentificationStatement(const OptionsList &options_list_arg)
{
  options_list = options_list_arg;
  if (options_list.num_options.find("max_dim_cova_group") != options_list.num_options.end())
    if (atoi(options_list.num_options["max_dim_cova_group"].c_str()) == 0)
      {
        cerr << "ERROR: The max_dim_cova_group option to identification only accepts integers > 0." << endl;
        exit(EXIT_FAILURE);
      }
}

void
IdentificationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.identification_present = true;
}

void
IdentificationStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_ident");

  /* Ensure that nograph, nodisplay and graph_format are also set in top-level
     options_.
     \todo factorize this code between identification and dynare_sensitivity,
     and provide a generic mechanism for this situation (maybe using regexps) */
  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("nodisplay");
  if (it != options_list.num_options.end())
    output << "options_.nodisplay = " << it->second << ";" << endl;
  it = options_list.num_options.find("nograph");
  if (it != options_list.num_options.end())
    output << "options_.nograph = " << it->second << ";" << endl;
  OptionsList::string_options_t::const_iterator it2 = options_list.string_options.find("graph_format");
  if (it2 != options_list.string_options.end())
    output << "options_.graph_format = '" << it2->second << "';" << endl;

  output << "dynare_identification(options_ident);" << endl;
}

WriteLatexDynamicModelStatement::WriteLatexDynamicModelStatement(const DynamicModel &dynamic_model_arg) :
  dynamic_model(dynamic_model_arg)
{
}

void
WriteLatexDynamicModelStatement::writeOutput(ostream &output, const string &basename) const
{
  dynamic_model.writeLatexFile(basename);
}

WriteLatexStaticModelStatement::WriteLatexStaticModelStatement(const StaticModel &static_model_arg) :
  static_model(static_model_arg)
{
}

void
WriteLatexStaticModelStatement::writeOutput(ostream &output, const string &basename) const
{
  static_model.writeLatexFile(basename);
}

ShockDecompositionStatement::ShockDecompositionStatement(const SymbolList &symbol_list_arg,
                                                         const OptionsList &options_list_arg) :
  symbol_list(symbol_list_arg),
  options_list(options_list_arg)
{
}

void
ShockDecompositionStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "oo_ = shock_decomposition(M_,oo_,options_,var_list_);\n";
}

ConditionalForecastStatement::ConditionalForecastStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
ConditionalForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_cond_fcst_");
  output << "imcforecast(constrained_paths_, constrained_vars_, options_cond_fcst_, constrained_perfect_foresight_);" << endl;
}

PlotConditionalForecastStatement::PlotConditionalForecastStatement(int periods_arg, const SymbolList &symbol_list_arg) :
  periods(periods_arg),
  symbol_list(symbol_list_arg)
{
}

void
PlotConditionalForecastStatement::writeOutput(ostream &output, const string &basename) const
{
  symbol_list.writeOutput("var_list_", output);
  if (periods == -1)
    output << "plot_icforecast(var_list_,[],options_);" << endl;
  else
    output << "plot_icforecast(var_list_, " << periods << ",options_);" << endl;
}

SvarIdentificationStatement::SvarIdentificationStatement(const svar_identification_restrictions_t &restrictions_arg,
                                                         const bool &upper_cholesky_present_arg,
                                                         const bool &lower_cholesky_present_arg,
                                                         const bool &constants_exclusion_present_arg,
                                                         const SymbolTable &symbol_table_arg) :
  restrictions(restrictions_arg),
  upper_cholesky_present(upper_cholesky_present_arg),
  lower_cholesky_present(lower_cholesky_present_arg),
  constants_exclusion_present(constants_exclusion_present_arg),
  symbol_table(symbol_table_arg)
{
}

int
SvarIdentificationStatement::getMaxLag() const
{
  int max_lag = 0;
  for (svar_identification_restrictions_t::const_iterator it = restrictions.begin(); it != restrictions.end(); it++)
    if (it->lag > max_lag)
      max_lag = it->lag;

  return max_lag;
}

void
SvarIdentificationStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (!mod_file_struct.svar_identification_present)
    mod_file_struct.svar_identification_present = true;
  else
    {
      cerr << "ERROR: You may only have one svar_identification block in your .mod file." << endl;
      exit(EXIT_FAILURE);
    }

  if (upper_cholesky_present && lower_cholesky_present)
    {
      cerr << "ERROR: Within the svar_identification statement, you may only have one of "
           << "upper_cholesky and lower_cholesky." << endl;
      exit(EXIT_FAILURE);
    }
}

void
SvarIdentificationStatement::writeOutput(ostream &output, const string &basename) const
{
  assert(!(upper_cholesky_present && lower_cholesky_present));
  output << "%" << endl
         << "% SVAR IDENTIFICATION" << endl
         << "%" << endl;

  if (upper_cholesky_present)
    output << "options_.ms.upper_cholesky=1;" << endl;

  if (lower_cholesky_present)
    output << "options_.ms.lower_cholesky=1;" << endl;

  if (constants_exclusion_present)
    output << "options_.ms.constants_exclusion=1;" << endl;

  if (!upper_cholesky_present && !lower_cholesky_present)
    {
      int n = symbol_table.endo_nbr();
      int m = 1; // this is the constant, not the shocks
      int r = getMaxLag();
      int k = r*n+m;

      if (k < 1)
        {
          cerr << "ERROR: lag = " << r
               << ", number of endogenous variables = " << n
               << ", number of exogenous variables = " << m
               << ". If this is not a logical error in the specification"
               << " of the .mod file, please report it to the Dynare Team." << endl;
          exit(EXIT_FAILURE);
        }
      if (n < 1)
        {
          cerr << "ERROR: Number of endogenous variables = " << n << "< 1. If this is not a logical "
               << "error in the specification of the .mod file, please report it to the Dynare Team." << endl;
          exit(EXIT_FAILURE);
        }
      output << "options_.ms.Qi = cell(" << n << ",1);" << endl;
      output << "options_.ms.Ri = cell(" << n << ",1);" << endl;

      for (svar_identification_restrictions_t::const_iterator it = restrictions.begin(); it != restrictions.end(); it++)
        {
          assert(it->lag >= 0);
	  if (it->lag == 0)
            output << "options_.ms.Qi{" << it->equation << "}(" << it->restriction_nbr << ", " << it->variable + 1 << ") = ";
	  else
	    {
	      int col = (it->lag-1)*n+it->variable+1;
	      if (col > k)
                {
                  cerr << "ERROR: lag =" << it->lag << ", num endog vars = " << n << "current endog var index = " << it->variable << ". Index "
                       << "out of bounds. If the above does not represent a logical error, please report this to the Dyanre Team." << endl;
                  exit(EXIT_FAILURE);
                }
	      output << "options_.ms.Ri{" << it->equation << "}(" << it->restriction_nbr << ", " << col << ") = ";
	    }
          it->value->writeOutput(output);
          output << ";" << endl;
        }
    }
}

MarkovSwitchingStatement::MarkovSwitchingStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
  OptionsList::num_options_t::const_iterator it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      OptionsList::num_options_t::const_iterator it_num_regimes =
        options_list.num_options.find("ms.number_of_regimes");
      assert(it_num_regimes !=  options_list.num_options.end());
      int num_regimes = lexical_cast< int >(it_num_regimes->second);

      vector<string> tokenizedRestrictions;
      split(tokenizedRestrictions, it_num->second, is_any_of("["), token_compress_on);
      for (vector<string>::iterator it = tokenizedRestrictions.begin();
            it != tokenizedRestrictions.end(); it++ )
        if (it->size() > 0)
          {
            vector<string> restriction;
            split(restriction, *it, is_any_of("], "));
            for (vector<string>::iterator it1 = restriction.begin();
                 it1 != restriction.end(); )
              if (it1->empty())
                restriction.erase(it1);
              else
                it1++;

            if (restriction.size() != 3)
              {
                cerr << "ERROR: restrictions in the subsample statement must be specified in the form "
                     << "[current_period_regime, next_period_regime, transition_probability]" << endl;
                exit(EXIT_FAILURE);
              }

            try
              {
                int from_regime = lexical_cast< int >(restriction[0]);
                int to_regime = lexical_cast< int >(restriction[1]);
                if (from_regime > num_regimes || to_regime > num_regimes)
                  {
                    cerr << "ERROR: the regimes specified in the restrictions option must be "
                         << "<= the number of regimes specified in the number_of_regimes option" << endl;
                    exit(EXIT_FAILURE);
                  }

                if (restriction_map.find(make_pair(from_regime, to_regime)) !=
                    restriction_map.end())
                  {
                    cerr << "ERROR: two restrictions were given for: " << from_regime << ", "
                         << to_regime << endl;
                    exit(EXIT_FAILURE);
                  }

                double transition_probability = lexical_cast< double >(restriction[2]);
                if (transition_probability > 1.0)
                  {
                    cerr << "ERROR: the transition probability, " << transition_probability
                         << " must be less than 1" << endl;
                    exit(EXIT_FAILURE);
                  }
                restriction_map[make_pair(from_regime, to_regime)] = transition_probability;
              }
            catch (const bad_lexical_cast &)
              {
                cerr << "ERROR: The first two arguments for a restriction must be integers "
                     << "specifying the regime and the last must be a double specifying the "
                     << "transition probability. You wrote [" << *it << endl;
                exit(EXIT_FAILURE);
              }
          }
    }
}

void
MarkovSwitchingStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  OptionsList::num_options_t::const_iterator itChain = options_list.num_options.find("ms.chain");
  assert(itChain != options_list.num_options.end());
  int chainNumber = atoi(itChain->second.c_str());
  if (++mod_file_struct.last_markov_switching_chain != chainNumber)
    {
      cerr << "ERROR: The markov_switching chain option takes consecutive integers "
           << "beginning at 1." << endl;
      exit(EXIT_FAILURE);
    }

  OptionsList::num_options_t::const_iterator it_num = options_list.num_options.find("ms.restrictions");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      OptionsList::num_options_t::const_iterator it_num_regimes =
        options_list.num_options.find("ms.number_of_regimes");
      assert(it_num_regimes != options_list.num_options.end());
      int num_regimes = lexical_cast< int >(it_num_regimes->second);
      vector<double> col_trans_prob_sum (num_regimes, 0);
      vector<double> row_trans_prob_sum (num_regimes, 0);
      vector<bool> all_restrictions_in_row (num_regimes, true);
      vector<bool> all_restrictions_in_col (num_regimes, true);
      for (int row=0; row<num_regimes; row++)
        for (int col=0; col<num_regimes; col++)
          if (restriction_map.find(make_pair(row+1, col+1)) != restriction_map.end())
            {
              row_trans_prob_sum[row] += restriction_map[make_pair(row+1, col+1)];
              col_trans_prob_sum[col] += restriction_map[make_pair(row+1, col+1)];
            }
          else
            {
              all_restrictions_in_row[row] = false;
              all_restrictions_in_col[col] = false;
            }

      for (int i=0; i<num_regimes; i++)
        {
          if (all_restrictions_in_row[i])
          {
            if (row_trans_prob_sum[i] != 1.0)
              {
                cerr << "ERROR: When all transitions probabilities are specified for a certain "
                     << "regime, they must sum to 1" << endl;
                exit(EXIT_FAILURE);
              }
          }
        else
          if (row_trans_prob_sum[i] >= 1.0)
            {
              cerr << "ERROR: When transition probabilites are not specified for every regime, "
                   << "their sum must be < 1" << endl;
              exit(EXIT_FAILURE);
            }

        if (all_restrictions_in_col[i])
          {
            if (col_trans_prob_sum[i] != 1.0)
              {
                cerr << "ERROR: When all transitions probabilities are specified for a certain "
                     << "regime, they must sum to 1" << endl;
                exit(EXIT_FAILURE);
              }
          }
        else
          if (col_trans_prob_sum[i] >= 1.0)
            {
              cerr << "ERROR: When transition probabilites are not specified for every regime, "
                   << "their sum must be < 1" << endl;
              exit(EXIT_FAILURE);
            }
        }
    }
}

void
MarkovSwitchingStatement::writeOutput(ostream &output, const string &basename) const
{
  bool isDurationAVec = true;
  string infStr("Inf");
  OptionsList::num_options_t::const_iterator itChain, itNOR, itDuration;
  map<pair<int, int>, double >::const_iterator itR;

  itChain = options_list.num_options.find("ms.chain");
  assert(itChain != options_list.num_options.end());

  itDuration = options_list.num_options.find("ms.duration");
  assert(itDuration != options_list.num_options.end());
  if (atof(itDuration->second.c_str()) || infStr.compare(itDuration->second) == 0)
    isDurationAVec = false;
  output << "options_.ms.duration = " << itDuration->second << ";" << endl;

  itNOR = options_list.num_options.find("ms.number_of_regimes");
  assert(itNOR != options_list.num_options.end());
  for (int i = 0; i < atoi(itNOR->second.c_str()); i++)
    {
      output << "options_.ms.ms_chain(" << itChain->second << ").regime("
             << i+1 << ").duration = options_.ms.duration";
      if (isDurationAVec)
        output << "(" << i+1 << ")";
      output << ";" << endl;
    }

  int restrictions_index = 0;
  for (itR=restriction_map.begin(); itR != restriction_map.end(); itR++)
    output << "options_.ms.ms_chain(" << itChain->second << ").restrictions("
           << ++restrictions_index << ") = {[" << itR->first.first << ", "
           << itR->first.second << ", " << itR->second << "]};" << endl;
}

SvarStatement::SvarStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SvarStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");
  assert((it0 != options_list.string_options.end()
          && it1 == options_list.string_options.end()
          && it2 == options_list.string_options.end()) ||
         (it0 == options_list.string_options.end()
          && it1 != options_list.string_options.end()
          && it2 == options_list.string_options.end()) ||
         (it0 == options_list.string_options.end()
          && it1 == options_list.string_options.end()
          && it2 != options_list.string_options.end()));
}

void
SvarStatement::writeOutput(ostream &output, const string &basename) const
{
  OptionsList::num_options_t::const_iterator it0, it1, it2;
  OptionsList::vec_int_options_t::const_iterator itv;

  it0 = options_list.num_options.find("ms.chain");
  assert(it0 != options_list.num_options.end());
  output << "options_.ms.ms_chain(" << it0->second << ")";

  it0 = options_list.string_options.find("ms.coefficients");
  it1 = options_list.string_options.find("ms.variances");
  it2 = options_list.string_options.find("ms.constants");

  if (it0 != options_list.string_options.end())
    output << "." << it0->second;
  else if (it1 != options_list.string_options.end())
    output << "." << it1->second;
  else
    output << "." << it2->second;

  itv = options_list.vector_int_options.find("ms.equations");
  output << ".equations = ";
  if (itv != options_list.vector_int_options.end())
    {
      assert(itv->second.size() >= 1);
      if (itv->second.size() > 1)
        {
          output << "[";
          for (vector<int>::const_iterator viit = itv->second.begin();
               viit != itv->second.end(); viit++)
            output << *viit << ";";
          output << "];" << endl;
        }
      else
        output << itv->second.front() << ";" << endl;
    }
  else
    output << "'ALL';" << endl;
}

SetTimeStatement::SetTimeStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
SetTimeStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
}

EstimationDataStatement::EstimationDataStatement(const OptionsList &options_list_arg) :
  options_list(options_list_arg)
{
}

void
EstimationDataStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.estimation_data_statement_present = true;

  OptionsList::num_options_t::const_iterator it = options_list.num_options.find("nobs");
  if (it != options_list.num_options.end())
    if (atoi(it->second.c_str()) <= 0)
      {
        cerr << "ERROR: The nobs option of the data statement only accepts positive integers." << endl;
        exit(EXIT_FAILURE);
      }

  if (options_list.string_options.find("file") == options_list.string_options.end())
    {
      cerr << "ERROR: The file option must be passed to the data statement." << endl;
      exit(EXIT_FAILURE);
    }
}

void
EstimationDataStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output, "options_.dataset");
  if (options_list.date_options.find("first_obs") == options_list.date_options.end())
    output << "options_.dataset.firstobs = options_.initial_period;" << endl;
}

SubsamplesStatement::SubsamplesStatement(const string &name1_arg,
                                         const string &name2_arg,
                                         const subsample_declaration_map_t subsample_declaration_map_arg,
                                         const SymbolTable &symbol_table_arg) :
  name1(name1_arg),
  name2(name2_arg),
  subsample_declaration_map(subsample_declaration_map_arg),
  symbol_table(symbol_table_arg)
{
}

void
SubsamplesStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

void
SubsamplesStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "subsamples_indx = get_new_or_existing_ei_index('subsamples_index', '"
         << name1 << "','" << name2 << "');" << endl
         << "estimation_info.subsamples_index(subsamples_indx) = {'" << name1;
  if (!name2.empty())
    output << ":" << name2;
  output << "'};" << endl
         << "estimation_info.subsamples(subsamples_indx).range = {};" << endl
         << "estimation_info.subsamples(subsamples_indx).range_index = {};" << endl;

  int map_indx = 1;
  for (subsample_declaration_map_t::const_iterator it = subsample_declaration_map.begin();
       it != subsample_declaration_map.end(); it++, map_indx++)
    output << "estimation_info.subsamples(subsamples_indx).range_index(" << map_indx << ") = {'"
           << it->first << "'};" << endl
           << "estimation_info.subsamples(subsamples_indx).range(" << map_indx << ").date1 = "
           << it->second.first << endl
           << "estimation_info.subsamples(subsamples_indx).range(" << map_indx << ").date2 = "
           << it->second.second << endl;

  // Initialize associated subsample substructures in estimation_info
  const SymbolType symb_type = symbol_table.getType(name1);
  string lhs_field;
  if (symb_type == eParameter)
    lhs_field = "parameter";
  else if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field;

  if (!name2.empty())
    output << "_corr";
  output << "_prior_index', '"
         << name1 << "', '";
  if (!name2.empty())
    output << name2;
  output << "');" << endl;

  lhs_field = "estimation_info." + lhs_field;
  if (!name2.empty())
    lhs_field += "_corr";
  output << lhs_field << "_prior_index(eifind) = {'" << name1;
  if (!name2.empty())
    output << ":" << name2;
  output << "'};" << endl;

  output << lhs_field << "(eifind).subsample_prior = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).subsample_prior(1:" << subsample_declaration_map.size()
         << ") = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).range_index = estimation_info.subsamples(subsamples_indx).range_index;"
         << endl;
}

SubsamplesEqualStatement::SubsamplesEqualStatement(const string &to_name1_arg,
                                                   const string &to_name2_arg,
                                                   const string &from_name1_arg,
                                                   const string &from_name2_arg,
                                                   const SymbolTable &symbol_table_arg) :
  to_name1(to_name1_arg),
  to_name2(to_name2_arg),
  from_name1(from_name1_arg),
  from_name2(from_name2_arg),
  symbol_table(symbol_table_arg)
{
}

void
SubsamplesEqualStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "subsamples_to_indx = get_new_or_existing_ei_index('subsamples_index', '"
         << to_name1 << "','" << to_name2 << "');" << endl
         << "estimation_info.subsamples_index(subsamples_to_indx) = {'" << to_name1;
  if (!to_name2.empty())
    output << ":" << to_name2;
  output << "'};" << endl
         << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');"
         << endl
         << "estimation_info.subsamples(subsamples_to_indx) = estimation_info.subsamples(subsamples_from_indx);"
         << endl;

  // Initialize associated subsample substructures in estimation_info
  const SymbolType symb_type = symbol_table.getType(to_name1);
  string lhs_field;
  if (symb_type == eParameter)
    lhs_field = "parameter";
  else if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field;

  if (!to_name2.empty())
    output << "_corr";
  output << "_prior_index', '"
         << to_name1 << "', '";
  if (!to_name2.empty())
    output << to_name2;
  output << "');" << endl;

  lhs_field = "estimation_info." + lhs_field;
  if (!to_name2.empty())
    lhs_field += "_corr";
  output << lhs_field << "_prior_index(eifind) = {'" << to_name1;
  if (!to_name2.empty())
    output << ":" << to_name2;
  output << "'};" << endl;

  output << lhs_field << "(eifind).subsample_prior = estimation_info.empty_prior;" << endl
         << lhs_field << "(eifind).subsample_prior(1:size(estimation_info.subsamples(subsamples_to_indx).range_index,2)) = estimation_info.empty_prior;"
         << endl
         << lhs_field << "(eifind).range_index = estimation_info.subsamples(subsamples_to_indx).range_index;"
         << endl;
}

BasicPriorStatement::~BasicPriorStatement()
{
}

BasicPriorStatement::BasicPriorStatement(const string &name_arg,
                                         const string &subsample_name_arg,
                                         const PriorDistributions &prior_shape_arg,
                                         const expr_t &variance_arg,
                                         const OptionsList &options_list_arg) :
  name(name_arg),
  subsample_name(subsample_name_arg),
  prior_shape(prior_shape_arg),
  variance(variance_arg),
  options_list(options_list_arg)
{
}

void
BasicPriorStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (prior_shape == eNoShape)
    {
      cerr << "ERROR: You must pass the shape option to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  if (options_list.num_options.find("mean") == options_list.num_options.end() &&
      options_list.num_options.find("mode") == options_list.num_options.end())
    {
      cerr << "ERROR: You must pass at least one of mean and mode to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  OptionsList::num_options_t::const_iterator it_stdev = options_list.num_options.find("stdev");
  if ((it_stdev == options_list.num_options.end() && variance == NULL) ||
      (it_stdev != options_list.num_options.end() && variance != NULL))
    {
      cerr << "ERROR: You must pass exactly one of stdev and variance to the prior statement." << endl;
      exit(EXIT_FAILURE);
    }

  OptionsList::num_options_t::const_iterator it_num = options_list.num_options.find("domain");
  if (it_num != options_list.num_options.end())
    {
      using namespace boost;
      vector<string> tokenizedDomain;
      split(tokenizedDomain, it_num->second, is_any_of("[ ]"), token_compress_on);
      if (tokenizedDomain.size() != 4)
        {
          cerr << "ERROR: You must pass exactly two values to the domain option." << endl;
          exit(EXIT_FAILURE);
        }
    }
}

void
BasicPriorStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicPriorStatement::writeCommonOutput(ostream &output, const string &lhs_field) const
{
  output << lhs_field << " = estimation_info.empty_prior;" << endl;

  writeCommonOutputHelper(output, "domain", lhs_field);
  writeCommonOutputHelper(output, "interval", lhs_field);
  writeCommonOutputHelper(output, "mean", lhs_field);
  writeCommonOutputHelper(output, "median", lhs_field);
  writeCommonOutputHelper(output, "mode", lhs_field);

  assert(prior_shape != eNoShape);
  output << lhs_field << ".shape = " << prior_shape << ";" << endl;

  writeCommonOutputHelper(output, "shift", lhs_field);
  writeCommonOutputHelper(output, "stdev", lhs_field);
  writeCommonOutputHelper(output, "truncate", lhs_field);

  if (variance)
    {
      output << lhs_field << ".variance = ";
      variance->writeOutput(output);
      output << ";" << endl;
    }
}

void
BasicPriorStatement::writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  OptionsList::num_options_t::const_iterator itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << lhs_field << "." << field << " = "<< itn->second << ";" << endl;
}

void
BasicPriorStatement::writePriorOutput(ostream &output, string &lhs_field, const string &name2) const
{
  if (subsample_name.empty())
    lhs_field += ".prior(1)";
  else
    {
      output << "subsamples_indx = get_existing_subsamples_indx('" << name << "','" << name2 << "');" << endl
             << "eisind = get_subsamples_range_indx(subsamples_indx, '" << subsample_name << "');" << endl;
      lhs_field += ".subsample_prior(eisind)";
    }
  writeCommonOutput(output, lhs_field);
}

PriorStatement::PriorStatement(const string &name_arg,
                               const string &subsample_name_arg,
                               const PriorDistributions &prior_shape_arg,
                               const expr_t &variance_arg,
                               const OptionsList &options_list_arg) :
  BasicPriorStatement(name_arg, subsample_name_arg, prior_shape_arg, variance_arg, options_list_arg)
{
}

void
PriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field = "estimation_info.parameter(eifind)";
  output << "eifind = get_new_or_existing_ei_index('parameter_prior_index', '"
         << name << "', '');" << endl
         << "estimation_info.parameter_prior_index(eifind) = {'" << name << "'};" << endl;
  writePriorOutput(output, lhs_field, "");
}

StdPriorStatement::StdPriorStatement(const string &name_arg,
                                     const string &subsample_name_arg,
                                     const PriorDistributions &prior_shape_arg,
                                     const expr_t &variance_arg,
                                     const OptionsList &options_list_arg,
                                     const SymbolTable &symbol_table_arg ) :
  BasicPriorStatement(name_arg, subsample_name_arg, prior_shape_arg, variance_arg, options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
StdPriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_prior_index', '"
         << name << "', '');" << endl
         << "estimation_info." << lhs_field << "_prior_index(eifind) = {'" << name << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "(eifind)";
  writePriorOutput(output, lhs_field, "");
}

CorrPriorStatement::CorrPriorStatement(const string &name_arg1, const string &name_arg2,
                                       const string &subsample_name_arg,
                                       const PriorDistributions &prior_shape_arg,
                                       const expr_t &variance_arg,
                                       const OptionsList &options_list_arg,
                                       const SymbolTable &symbol_table_arg ) :
  BasicPriorStatement(name_arg1, subsample_name_arg, prior_shape_arg, variance_arg, options_list_arg),
  name1(name_arg2),
  symbol_table(symbol_table_arg)
{
}

void
CorrPriorStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  BasicPriorStatement::checkPass(mod_file_struct, warnings);
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).prior statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
}

void
CorrPriorStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_corr_prior_index', '"
         << name << "', '" << name1 << "');" << endl
         << "estimation_info." << lhs_field << "_corr_prior_index(eifind) = {'"
         << name << ":" << name1 << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "_corr(eifind)";
  writePriorOutput(output, lhs_field, name1);
}

PriorEqualStatement::PriorEqualStatement(const string &to_declaration_type_arg,
                                         const string &to_name1_arg,
                                         const string &to_name2_arg,
                                         const string &to_subsample_name_arg,
                                         const string &from_declaration_type_arg,
                                         const string &from_name1_arg,
                                         const string &from_name2_arg,
                                         const string &from_subsample_name_arg,
                                         const SymbolTable &symbol_table_arg) :
  to_declaration_type(to_declaration_type_arg),
  to_name1(to_name1_arg),
  to_name2(to_name2_arg),
  to_subsample_name(to_subsample_name_arg),
  from_declaration_type(from_declaration_type_arg),
  from_name1(from_name1_arg),
  from_name2(from_name2_arg),
  from_subsample_name(from_subsample_name_arg),
  symbol_table(symbol_table_arg)
{
}

void
PriorEqualStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if ((to_declaration_type != "par" && to_declaration_type != "std" && to_declaration_type != "corr") ||
      (from_declaration_type != "par" && from_declaration_type != "std" && from_declaration_type != "corr"))
    {
      cerr << "Internal Dynare Error" << endl;
      exit(EXIT_FAILURE);
    }
}

void
PriorEqualStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
PriorEqualStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field, rhs_field;

  if (to_declaration_type == "par")
    lhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(to_name1), lhs_field);

  if (from_declaration_type == "par")
    rhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(from_name1), rhs_field);


  if (to_declaration_type == "corr")
    lhs_field += "_corr";

  if (from_declaration_type == "corr")
    rhs_field += "_corr";

  output << "ei_to_ind = get_new_or_existing_ei_index('" << lhs_field << "_prior_index', '"
         << to_name1 << "', '" << to_name2<< "');" << endl
         << "ei_from_ind = get_new_or_existing_ei_index('" << rhs_field << "_prior_index', '"
         << from_name1 << "', '" << from_name2<< "');" << endl
         << "estimation_info." << lhs_field << "_prior_index(ei_to_ind) = {'" << to_name1;

  if (to_declaration_type == "corr")
    output << ":" << to_name2;
  output << "'};" << endl;

  if (to_declaration_type == "par")
    lhs_field = "parameter";

  if (from_declaration_type == "par")
    rhs_field = "parameter";

  lhs_field = "estimation_info." + lhs_field + "(ei_to_ind)";
  rhs_field = "estimation_info." + rhs_field + "(ei_from_ind)";

  if (to_subsample_name.empty())
    lhs_field += ".prior";
  else
    {
      output << "subsamples_to_indx = get_existing_subsamples_indx('" << to_name1 << "','" << to_name2 << "');" << endl
             << "ei_to_ss_ind = get_subsamples_range_indx(subsamples_to_indx, '" << to_subsample_name << "');" << endl;
      lhs_field += ".subsample_prior(ei_to_ss_ind)";
    }

  if (from_subsample_name.empty())
    rhs_field += ".prior";
  else
    {
      output << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');" << endl
             << "ei_from_ss_ind = get_subsamples_range_indx(subsamples_from_indx, '" << from_subsample_name << "');" << endl;
      rhs_field += ".subsample_prior(ei_from_ss_ind)";
    }

  output << lhs_field << " = " << rhs_field << ";" << endl;
}

BasicOptionsStatement::~BasicOptionsStatement()
{
}

BasicOptionsStatement::BasicOptionsStatement(const string &name_arg,
                                             const string &subsample_name_arg,
                                             const OptionsList &options_list_arg) :
  name(name_arg),
  subsample_name(subsample_name_arg),
  options_list(options_list_arg)
{
}

void
BasicOptionsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
}

void
BasicOptionsStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
BasicOptionsStatement::writeCommonOutput(ostream &output, const string &lhs_field) const
{
  output << lhs_field << " = estimation_info.empty_options;" << endl;

  writeCommonOutputHelper(output, "bounds", lhs_field);
  writeCommonOutputHelper(output, "init", lhs_field);
  writeCommonOutputHelper(output, "jscale", lhs_field);
}

void
BasicOptionsStatement::writeCommonOutputHelper(ostream &output, const string &field, const string &lhs_field) const
{
  OptionsList::num_options_t::const_iterator itn = options_list.num_options.find(field);
  if (itn != options_list.num_options.end())
    output << lhs_field << "." << field << " = " << itn->second << ";" << endl;
}

void
BasicOptionsStatement::writeOptionsOutput(ostream &output, string &lhs_field, const string &name2) const
{
  if (subsample_name.empty())
    lhs_field += ".options(1)";
  else
    {
      output << "subsamples_indx = get_existing_subsamples_indx('" << name << "','" << name2 << "');" << endl
             << "eisind = get_subsamples_range_indx(subsamples_indx, '" << subsample_name << "');" << endl;
      lhs_field += ".subsample_options(eisind)";
    }
  writeCommonOutput(output, lhs_field);
}

OptionsStatement::OptionsStatement(const string &name_arg,
                                   const string &subsample_name_arg,
                                   const OptionsList &options_list_arg) :
  BasicOptionsStatement(name_arg, subsample_name_arg, options_list_arg)
{
}

void
OptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field = "estimation_info.parameter(eifind)";
  output << "eifind = get_new_or_existing_ei_index('parameter_options_index', '"
         << name << "', '');" << endl
         << "estimation_info.parameter_options_index(eifind) = {'" << name << "'};" << endl;
  writeOptionsOutput(output, lhs_field, "");
}

StdOptionsStatement::StdOptionsStatement(const string &name_arg,
                                         const string &subsample_name_arg,
                                         const OptionsList &options_list_arg,
                                         const SymbolTable &symbol_table_arg ) :
  BasicOptionsStatement(name_arg, subsample_name_arg, options_list_arg),
  symbol_table(symbol_table_arg)
{
}

void
StdOptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);
  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_options_index', '"
         << name << "', '');" << endl
         << "estimation_info." << lhs_field << "_options_index(eifind) = {'" << name << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "(eifind)";
  writeOptionsOutput(output, lhs_field, "");
}

CorrOptionsStatement::CorrOptionsStatement(const string &name_arg1, const string &name_arg2,
                                           const string &subsample_name_arg,
                                           const OptionsList &options_list_arg,
                                           const SymbolTable &symbol_table_arg ) :
  BasicOptionsStatement(name_arg1, subsample_name_arg, options_list_arg),
  name1(name_arg2),
  symbol_table(symbol_table_arg)
{
}

void
CorrOptionsStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (symbol_table.getType(name) != symbol_table.getType(name1))
    {
      cerr << "ERROR: In the corr(A,B).options statement, A and B must be of the same type. "
           << "In your case, " << name << " and " << name1 << " are of different "
           << "types." << endl;
      exit(EXIT_FAILURE);
    }
}

void
CorrOptionsStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field;
  get_base_name(symbol_table.getType(name), lhs_field);

  output << "eifind = get_new_or_existing_ei_index('" << lhs_field << "_corr_options_index', '"
         << name << "', '" << name1 << "');" << endl
         << "estimation_info." << lhs_field << "_corr_options_index(eifind) = {'"
         << name << ":" << name1 << "'};" << endl;

  lhs_field = "estimation_info." + lhs_field + "_corr(eifind)";
  writeOptionsOutput(output, lhs_field, name1);
}

OptionsEqualStatement::OptionsEqualStatement(const string &to_declaration_type_arg,
                                             const string &to_name1_arg,
                                             const string &to_name2_arg,
                                             const string &to_subsample_name_arg,
                                             const string &from_declaration_type_arg,
                                             const string &from_name1_arg,
                                             const string &from_name2_arg,
                                             const string &from_subsample_name_arg,
                                             const SymbolTable &symbol_table_arg) :
  to_declaration_type(to_declaration_type_arg),
  to_name1(to_name1_arg),
  to_name2(to_name2_arg),
  to_subsample_name(to_subsample_name_arg),
  from_declaration_type(from_declaration_type_arg),
  from_name1(from_name1_arg),
  from_name2(from_name2_arg),
  from_subsample_name(from_subsample_name_arg),
  symbol_table(symbol_table_arg)
{
}

void
OptionsEqualStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if ((to_declaration_type != "par" && to_declaration_type != "std" && to_declaration_type != "corr") ||
      (from_declaration_type != "par" && from_declaration_type != "std" && from_declaration_type != "corr"))
    {
      cerr << "Internal Dynare Error" << endl;
      exit(EXIT_FAILURE);
    }
}

void
OptionsEqualStatement::get_base_name(const SymbolType symb_type, string &lhs_field) const
{
  if (symb_type == eExogenous || symb_type == eExogenousDet)
    lhs_field = "structural_innovation";
  else
    lhs_field = "measurement_error";
}

void
OptionsEqualStatement::writeOutput(ostream &output, const string &basename) const
{
  string lhs_field, rhs_field;

  if (to_declaration_type == "par")
    lhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(to_name1), lhs_field);

  if (from_declaration_type == "par")
    rhs_field = "parameter";
  else
    get_base_name(symbol_table.getType(from_name1), rhs_field);


  if (to_declaration_type == "corr")
    lhs_field += "_corr";

  if (from_declaration_type == "corr")
    rhs_field += "_corr";

  output << "ei_to_ind = get_new_or_existing_ei_index('" << lhs_field << "_options_index', '"
         << to_name1 << "', '" << to_name2<< "');" << endl
         << "ei_from_ind = get_new_or_existing_ei_index('" << rhs_field << "_options_index', '"
         << from_name1 << "', '" << from_name2<< "');" << endl
         << "estimation_info." << lhs_field << "_options_index(ei_to_ind) = {'" << to_name1;

  if (to_declaration_type == "corr")
    output << ":" << to_name2;
  output << "'};" << endl;

  if (to_declaration_type == "par")
    lhs_field = "parameter";

  if (from_declaration_type == "par")
    rhs_field = "parameter";

  lhs_field = "estimation_info." + lhs_field + "(ei_to_ind)";
  rhs_field = "estimation_info." + rhs_field + "(ei_from_ind)";

  if (to_subsample_name.empty())
    lhs_field += ".options";
  else
    {
      output << "subsamples_to_indx = get_existing_subsamples_indx('" << to_name1 << "','" << to_name2 << "');" << endl
             << "ei_to_ss_ind = get_subsamples_range_indx(subsamples_to_indx, '" << to_subsample_name << "');" << endl;
      lhs_field += ".subsample_options(ei_to_ss_ind)";
    }

  if (from_subsample_name.empty())
    rhs_field += ".options";
  else
    {
      output << "subsamples_from_indx = get_existing_subsamples_indx('" << from_name1 << "','" << from_name2 << "');" << endl
             << "ei_from_ss_ind = get_subsamples_range_indx(subsamples_from_indx, '" << from_subsample_name << "');" << endl;
      rhs_field += ".subsample_options(ei_from_ss_ind)";
    }

  output << lhs_field << " = " << rhs_field << ";" << endl;
}

CalibSmootherStatement::CalibSmootherStatement(const SymbolList &symbol_list_arg,
                                               const OptionsList &options_list_arg)
  : symbol_list(symbol_list_arg), options_list(options_list_arg)
{
}

void
CalibSmootherStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.calib_smoother_present = true;
}

void
CalibSmootherStatement::writeOutput(ostream &output, const string &basename) const
{
  options_list.writeOutput(output);
  symbol_list.writeOutput("var_list_", output);
  output << "options_.smoother = 1;" << endl;
  output << "options_.order = 1;" << endl;
  output << "evaluate_smoother('calibration',var_list_);" << endl;
}

ExtendedPathStatement::ExtendedPathStatement(const OptionsList &options_list_arg)
  : options_list(options_list_arg)
{
}

void
ExtendedPathStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.extended_path_present = true;

  if (options_list.num_options.find("periods") == options_list.num_options.end())
    {
      cerr << "ERROR: the 'periods' option of 'extended_path' is mandatory" << endl;
      exit(EXIT_FAILURE);
    }
}

void
ExtendedPathStatement::writeOutput(ostream &output, const string &basename) const
{
  // Beware: options do not have the same name in the interface and in the M code...

  for (OptionsList::num_options_t::const_iterator it = options_list.num_options.begin();
       it != options_list.num_options.end(); ++it)
    if (it->first != string("periods"))
      output << "options_." << it->first << " = " << it->second << ";" << endl;

  output << "oo_.endo_simul = [ oo_.steady_state, extended_path([], " << options_list.num_options.find("periods")->second
         << ") ];" << endl
         << "oo_.exo_simul = oo_.ep.shocks;" << endl;
}

ModelDiagnosticsStatement::ModelDiagnosticsStatement()
{
}

void
ModelDiagnosticsStatement::writeOutput(ostream &output, const string &basename) const
{
  output << "model_diagnostics(M_,options_,oo_);" << endl;
}
