/*
 * Copyright (C) 2003-2015 Dynare Team
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

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>

#include "NumericalInitialization.hh"

InitParamStatement::InitParamStatement(int symb_id_arg,
                                       const expr_t param_value_arg,
                                       const SymbolTable &symbol_table_arg) :
  symb_id(symb_id_arg),
  param_value(param_value_arg),
  symbol_table(symbol_table_arg)
{
}

void
InitParamStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  if (symbol_table.getName(symb_id) == "dsge_prior_weight")
    mod_file_struct.dsge_prior_weight_initialized = true;
}

void
InitParamStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << symbol_table.getName(symb_id) << " = M_.params( " << id << " );\n";
}

void
InitParamStatement::writeM_Output(ostream &output, const string &basename, bool minimal_workspace)
{
  int id = symbol_table.getTypeSpecificID(symb_id) + 1;
  output << "M_.params( " << id << " ) = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  if (!minimal_workspace)
    output << symbol_table.getName(symb_id) << " = M_.params( " << id << " );" << endl;
}

void
InitParamStatement::writeCOutput(ostream &output, const string &basename)
{
  int id = symbol_table.getTypeSpecificID(symb_id);
  output << "params[ " << id << " ] = ";
  param_value->writeOutput(output);
  output << ";" << endl;
  output << "double " << symbol_table.getName(symb_id) << " = params[ " << id << " ];" << endl;
}

void
InitParamStatement::fillEvalContext(eval_context_t &eval_context) const
{
  try
    {
      eval_context[symb_id] = param_value->eval(eval_context);
    }
  catch (ExprNode::EvalException &e)
    {
      // Do nothing
    }
}

InitOrEndValStatement::InitOrEndValStatement(const init_values_t &init_values_arg,
                                             const SymbolTable &symbol_table_arg,
                                             const bool &all_values_required_arg) :
  init_values(init_values_arg),
  symbol_table(symbol_table_arg),
  all_values_required(all_values_required_arg)
{
}

void
InitOrEndValStatement::fillEvalContext(eval_context_t &eval_context) const
{
  for (init_values_t::const_iterator it = init_values.begin();
       it != init_values.end(); it++)
    {
      try
        {
          eval_context[it->first] = (it->second)->eval(eval_context);
        }
      catch (ExprNode::EvalException &e)
        {
          // Do nothing
        }
    }
}

set<int>
InitOrEndValStatement::getUninitializedVariables(SymbolType type)
{
  set<int> unused;
  if (!all_values_required)
    return unused;

  if (type == eEndogenous)
    unused = symbol_table.getEndogenous();
  else if (type == eExogenous)
    unused = symbol_table.getExogenous();
  else
    {
      cerr << "ERROR: Shouldn't arrive here." << endl;
      exit(EXIT_FAILURE);
    }

  set<int>::iterator sit;
  for (init_values_t::const_iterator it = init_values.begin();
       it != init_values.end(); it++)
    {
      sit = unused.find(it->first);
      if (sit != unused.end())
        unused.erase(sit);
    }
  return unused;
}

void
InitOrEndValStatement::writeInitValues(ostream &output) const
{
  for (init_values_t::const_iterator it = init_values.begin();
       it != init_values.end(); it++)
    {
      const int symb_id = it->first;
      const expr_t expression = it->second;

      SymbolType type = symbol_table.getType(symb_id);
      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      if (type == eEndogenous)
        output << "oo_.steady_state";
      else if (type == eExogenous)
        output << "oo_.exo_steady_state";
      else if (type == eExogenousDet)
        output << "oo_.exo_det_steady_state";

      output << "( " << tsid << " ) = ";
      expression->writeOutput(output);
      output << ";" << endl;
    }
}

InitValStatement::InitValStatement(const init_values_t &init_values_arg,
                                   const SymbolTable &symbol_table_arg,
                                   const bool &all_values_required_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg, all_values_required_arg)
{
}

void
InitValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  set<int> exogs = getUninitializedVariables(eExogenous);
  set<int> endogs = getUninitializedVariables(eEndogenous);

  if (endogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following endogenous variables in initval:";
      for (set<int>::const_iterator it = endogs.begin(); it != endogs.end(); it++)
        cerr << " " << symbol_table.getName(*it);
      cerr << endl;
    }

  if (exogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following exogenous variables in initval:";
      for (set<int>::const_iterator it = exogs.begin(); it != exogs.end(); it++)
        cerr << " " << symbol_table.getName(*it);
      cerr << endl;
    }

  if (endogs.size() > 0 || exogs.size() > 0)
    exit(EXIT_FAILURE);
}

void
InitValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% INITVAL instructions" << endl
         << "%" << endl;
  // Writing initval block to set initial values for variables
  output << "options_.initval_file = 0;" << endl;

  writeInitValues(output);
}

void
InitValStatement::writeOutputPostInit(ostream &output) const
{
  output << "if M_.exo_nbr > 0;" << endl
         << "\too_.exo_simul = [ones(M_.maximum_lag,1)*oo_.exo_steady_state'];" << endl
         <<"end;" << endl
         << "if M_.exo_det_nbr > 0;" << endl
         << "\too_.exo_det_simul = [ones(M_.maximum_lag,1)*oo_.exo_det_steady_state'];" << endl
         <<"end;" << endl;
}

EndValStatement::EndValStatement(const init_values_t &init_values_arg,
                                 const SymbolTable &symbol_table_arg,
                                 const bool &all_values_required_arg) :
  InitOrEndValStatement(init_values_arg, symbol_table_arg, all_values_required_arg)
{
}

void
EndValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  set<int> exogs = getUninitializedVariables(eExogenous);
  set<int> endogs = getUninitializedVariables(eEndogenous);

  if (endogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following endogenous variables in endval:";
      for (set<int>::const_iterator it = endogs.begin(); it != endogs.end(); it++)
        cerr << " " << symbol_table.getName(*it);
      cerr << endl;
    }

  if (exogs.size() > 0)
    {
      cerr << "ERROR: You have not set the following exogenous variables in endval:";
      for (set<int>::const_iterator it = exogs.begin(); it != exogs.end(); it++)
        cerr << " " << symbol_table.getName(*it);
      cerr << endl;
    }

  if (endogs.size() > 0 || exogs.size() > 0)
    exit(EXIT_FAILURE);
}

void
EndValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% ENDVAL instructions" << endl
         << "%" << endl;
  // Writing endval block to set terminal values for variables
  output << "ys0_= oo_.steady_state;" << endl
         << "ex0_ = oo_.exo_steady_state;" << endl;

  writeInitValues(output);
}

HistValStatement::HistValStatement(const hist_values_t &hist_values_arg,
                                   const SymbolTable &symbol_table_arg) :
  hist_values(hist_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HistValStatement::checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings)
{
  mod_file_struct.histval_present = true;
}

void
HistValStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (hist_values_t::const_iterator it = hist_values.begin();
       it != hist_values.end(); it++)
    {
      int symb_id = it->first.first;
      int lag = it->first.second;
      const expr_t expression = it->second;

      SymbolType type = symbol_table.getType(symb_id);

      // For a lag greater than 1 on endo, or for any exo, lookup for auxiliary variable
      if ((type == eEndogenous && lag < 0) || type == eExogenous)
        {
          try
            {
              // This function call must remain the 1st statement in this block
              symb_id = symbol_table.searchAuxiliaryVars(symb_id, lag);
              lag = 0;
              type = eEndogenous;
            }
          catch (SymbolTable::SearchFailedException &e)
            {
              if (type == eEndogenous)
                {
                  cerr << "HISTVAL: internal error of Dynare, please contact the developers";
                  exit(EXIT_FAILURE);
                }
              // We don't fail for exogenous, because they are not replaced by
              // auxiliary variables in deterministic mode.
            }
        }

      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      if (type == eExogenous)
        {
          output << "oo_.exo_simul( M_.maximum_lag + " << lag << ", " << tsid << " ) = ";
          expression->writeOutput(output);
          output << ";" << endl;
        }
      else if (type != eExogenousDet)
        {
          output << "oo_.exo_det_simul( M_.maximum_lag + " << lag  << ", " << tsid << " ) = ";
          expression->writeOutput(output);
          output << ";" << endl;
        }
    }
}

void
HistValStatement::writeM_Output(ostream &output, const string &basename, bool minimal_workspace)
{
  output << "%" << endl
         << "% HISTVAL instructions" << endl
         << "%" << endl
         << "M_.endo_histval = zeros(M_.endo_nbr,M_.maximum_endo_lag);" << endl;

  for (hist_values_t::const_iterator it = hist_values.begin();
       it != hist_values.end(); it++)
    {
      int symb_id = it->first.first;
      int lag = it->first.second;
      const expr_t expression = it->second;

      SymbolType type = symbol_table.getType(symb_id);

      // For a lag greater than 1 on endo, or for any exo, lookup for auxiliary variable
      if ((type == eEndogenous && lag < 0) || type == eExogenous)
        {
          try
            {
              // This function call must remain the 1st statement in this block
              symb_id = symbol_table.searchAuxiliaryVars(symb_id, lag);
              lag = 0;
              type = eEndogenous;
            }
          catch (SymbolTable::SearchFailedException &e)
            {
              if (type == eEndogenous)
                {
                  cerr << "HISTVAL: internal error of Dynare, please contact the developers";
                  exit(EXIT_FAILURE);
                }
              // We don't fail for exogenous, because they are not replaced by
              // auxiliary variables in deterministic mode.
            }
        }

      int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      if (type == eEndogenous)
        {
          output << "M_.endo_histval( " << tsid << ", M_.maximum_endo_lag + " << lag << ") = ";
          expression->writeOutput(output);
          output << ";" << endl;
        }
    }
}

InitvalFileStatement::InitvalFileStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
InitvalFileStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% INITVAL_FILE statement" << endl
         << "%" << endl
         << "options_.initval_file = 1;" << endl
         << "initvalf('" << filename << "');" << endl;
}

HistvalFileStatement::HistvalFileStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
HistvalFileStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "histvalf('" << filename << "');" << endl;
}

HomotopyStatement::HomotopyStatement(const homotopy_values_t &homotopy_values_arg,
                                     const SymbolTable &symbol_table_arg) :
  homotopy_values(homotopy_values_arg),
  symbol_table(symbol_table_arg)
{
}

void
HomotopyStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "%" << endl
         << "% HOMOTOPY_SETUP instructions" << endl
         << "%" << endl
         << "options_.homotopy_values = [];" << endl;

  for (homotopy_values_t::const_iterator it = homotopy_values.begin();
       it != homotopy_values.end(); it++)
    {
      const int &symb_id = it->first;
      const expr_t expression1 = it->second.first;
      const expr_t expression2 = it->second.second;

      const SymbolType type = symbol_table.getType(symb_id);
      const int tsid = symbol_table.getTypeSpecificID(symb_id) + 1;

      output << "options_.homotopy_values = vertcat(options_.homotopy_values, [ " << type << ", " << tsid << ", ";
      if (expression1 != NULL)
        expression1->writeOutput(output);
      else
        output << "NaN";
      output << ", ";
      expression2->writeOutput(output);
      output << "]);" << endl;
    }
}

SaveParamsAndSteadyStateStatement::SaveParamsAndSteadyStateStatement(const string &filename_arg) :
  filename(filename_arg)
{
}

void
SaveParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  output << "save_params_and_steady_state('" << filename << "');" << endl;
}

LoadParamsAndSteadyStateStatement::LoadParamsAndSteadyStateStatement(const string &filename,
                                                                     const SymbolTable &symbol_table_arg,
                                                                     WarningConsolidation &warnings) :
  symbol_table(symbol_table_arg)
{
  cout << "Reading " << filename << "." << endl;

  ifstream f;
  f.open(filename.c_str(), ios::in);
  if (f.fail())
    {
      cerr << "ERROR: Can't open " << filename << endl;
      exit(EXIT_FAILURE);
    }

  while (true)
    {
      string symb_name, value;
      f >> symb_name >> value;
      if (f.eof())
        break;

      try
        {
          int symb_id = symbol_table.getID(symb_name);
          content[symb_id] = value;
        }
      catch (SymbolTable::UnknownSymbolNameException &e)
        {
          warnings << "WARNING: Unknown symbol " << symb_name << " in " << filename << endl;
        }
    }
  f.close();
}

void
LoadParamsAndSteadyStateStatement::writeOutput(ostream &output, const string &basename, bool minimal_workspace) const
{
  for (map<int, string>::const_iterator it = content.begin();
       it != content.end(); it++)
    {
      switch (symbol_table.getType(it->first))
        {
        case eParameter:
          break;
        case eEndogenous:
          output << "oo_.steady_state";
          break;
        case eExogenous:
          output << "oo_.exo_steady_state";
          break;
        case eExogenousDet:
          output << "oo_.exo_det_steady_state";
          break;
        default:
          cerr << "ERROR: Unsupported variable type for " << symbol_table.getName(it->first) << " in load_params_and_steady_state" << endl;
          exit(EXIT_FAILURE);
        }
      int tsid = symbol_table.getTypeSpecificID(it->first) + 1;

      switch (symbol_table.getType(it->first))
        {
        case eParameter:
          break;
        default:
          output << "(" << tsid << ") = " << it->second << ";" << endl;
        }
    }
}

void
LoadParamsAndSteadyStateStatement::writeM_Output(ostream &output, const string &basename, bool minimal_workspace)
{
  for (map<int, string>::const_iterator it = content.begin();
       it != content.end(); it++)
    {
      switch (symbol_table.getType(it->first))
        {
        case eParameter:
          output << "M_.params";
          break;
        case eEndogenous:
        case eExogenous:
        case eExogenousDet:
          break;
        default:
          cerr << "ERROR: Unsupported variable type for " << symbol_table.getName(it->first) << " in load_params_and_steady_state" << endl;
          exit(EXIT_FAILURE);
        }
      int tsid = symbol_table.getTypeSpecificID(it->first) + 1;

      switch (symbol_table.getType(it->first))
        {
        case eParameter:
          output << "(" << tsid << ") = " << it->second << ";" << endl;
          break;
        default:
          break;
        }
    }
}

void
LoadParamsAndSteadyStateStatement::fillEvalContext(eval_context_t &eval_context) const
{
  for (map<int, string>::const_iterator it = content.begin();
       it != content.end(); it++)
    eval_context[it->first] = atof(it->second.c_str());
}
