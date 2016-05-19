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

#ifndef _NUMERICALINITIALIZATION_HH
#define _NUMERICALINITIALIZATION_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "SymbolTable.hh"
#include "ExprNode.hh"
#include "Statement.hh"

class InitParamStatement : public Statement
{
private:
  const int symb_id;
  const expr_t param_value;
  const SymbolTable &symbol_table;
public:
  InitParamStatement(int symb_id_arg, const expr_t param_value_arg,
                     const SymbolTable &symbol_table_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
  virtual void writeJuliaOutput(ostream &output, const string &basename);
  virtual void writeCOutput(ostream &output, const string &basename);
  //! Fill eval context with parameter value
  void fillEvalContext(eval_context_t &eval_context) const;
};

class InitOrEndValStatement : public Statement
{
public:
  /*!
    We use a vector instead of a map, since the order of declaration matters:
    an initialization can depend on a previously initialized variable inside the block
  */
  typedef vector<pair<int, expr_t> > init_values_t;
protected:
  const init_values_t init_values;
  const SymbolTable &symbol_table;
  const bool all_values_required;
public:
  InitOrEndValStatement(const init_values_t &init_values_arg,
                        const SymbolTable &symbol_table_arg,
                        const bool &all_values_required_arg);
  //! Return set of unused variables by type
  set<int> getUninitializedVariables(SymbolType type);
  //! Fill eval context with variables values
  void fillEvalContext(eval_context_t &eval_context) const;
protected:
  void writeInitValues(ostream &output) const;
};

class InitValStatement : public InitOrEndValStatement
{
public:
  InitValStatement(const init_values_t &init_values_arg,
                   const SymbolTable &symbol_table_arg,
                   const bool &all_values_required_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
  //! Writes initializations for oo_.exo_simul and oo_.exo_det_simul
  void writeOutputPostInit(ostream &output) const;
};

class EndValStatement : public InitOrEndValStatement
{
public:
  EndValStatement(const init_values_t &init_values_arg,
                  const SymbolTable &symbol_table_arg,
                  const bool &all_values_required_arg);
  //! Workaround for trac ticket #35
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class HistValStatement : public Statement
{
public:
  /*!
    Contrary to Initval and Endval, we use a map, since it is impossible to reuse
    a given initialization value in a second initialization inside the block.
    Maps pairs (symbol_id, lag) to expr_t
  */
  typedef map<pair<int, int>, expr_t> hist_values_t;
private:
  const hist_values_t hist_values;
  const SymbolTable &symbol_table;
public:
  HistValStatement(const hist_values_t &hist_values_arg,
                   const SymbolTable &symbol_table_arg);
  //! Workaround for trac ticket #157
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class InitvalFileStatement : public Statement
{
private:
  const string filename;
public:
  InitvalFileStatement(const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class HistvalFileStatement : public Statement
{
private:
  const string filename;
public:
  HistvalFileStatement(const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class HomotopyStatement : public Statement
{
public:
  //! Stores the declarations of homotopy_setup
  /*! Order matter so we use a vector. First expr_t can be NULL if no initial value given. */
  typedef vector<pair<int, pair<expr_t, expr_t> > > homotopy_values_t;
private:
  const homotopy_values_t homotopy_values;
  const SymbolTable &symbol_table;
public:
  HomotopyStatement(const homotopy_values_t &homotopy_values_arg,
                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class SaveParamsAndSteadyStateStatement : public Statement
{
private:
  const string filename;
public:
  SaveParamsAndSteadyStateStatement(const string &filename_arg);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
};

class LoadParamsAndSteadyStateStatement : public Statement
{
private:
  const SymbolTable &symbol_table;
  //! Content of the file
  /*! Maps symbol ID to numeric value (stored as string) */
  map<int, string> content;
public:
  LoadParamsAndSteadyStateStatement(const string &filename,
                                    const SymbolTable &symbol_table_arg,
                                    WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
  //! Fill eval context with parameters/variables values
  void fillEvalContext(eval_context_t &eval_context) const;
};

#endif
