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

#ifndef _SHOCKS_HH
#define _SHOCKS_HH

using namespace std;

#include <string>
#include <vector>
#include <map>

#include "Statement.hh"
#include "SymbolTable.hh"
#include "ExprNode.hh"

class AbstractShocksStatement : public Statement
{
public:
  struct DetShockElement
  {
    int period1;
    int period2;
    expr_t value;
  };
  //The boolean element indicates if the shock is a surprise (false) or a perfect foresight (true) shock.
  //This boolean is used only in case of conditional forecast with extended path method (simulation_type = deterministic).
  typedef map<int, vector<DetShockElement> > det_shocks_t;
protected:
  //! Is this statement a "mshocks" statement ? (instead of a "shocks" statement)
  const bool mshocks;
  //! Does this "shocks" statement replace the previous ones?
  const bool overwrite;
  const det_shocks_t det_shocks;
  const SymbolTable &symbol_table;
  void writeDetShocks(ostream &output) const;

  AbstractShocksStatement(bool mshocks_arg, bool overwrite_arg,
                          const det_shocks_t &det_shocks_arg,
                          const SymbolTable &symbol_table_arg);
};

class ShocksStatement : public AbstractShocksStatement
{
public:
  typedef map<int, expr_t> var_and_std_shocks_t;
  typedef map<pair<int, int>, expr_t> covar_and_corr_shocks_t;
private:
  const var_and_std_shocks_t var_shocks, std_shocks;
  const covar_and_corr_shocks_t covar_shocks, corr_shocks;
  void writeVarOrStdShock(ostream &output, var_and_std_shocks_t::const_iterator &it, bool stddev) const;
  void writeVarAndStdShocks(ostream &output) const;
  void writeCovarOrCorrShock(ostream &output, covar_and_corr_shocks_t::const_iterator &it, bool corr) const;
  void writeCovarAndCorrShocks(ostream &output) const;
  bool has_calibrated_measurement_errors() const;
public:
  ShocksStatement(bool overwrite_arg,
                  const det_shocks_t &det_shocks_arg,
                  const var_and_std_shocks_t &var_shocks_arg,
                  const var_and_std_shocks_t &std_shocks_arg,
                  const covar_and_corr_shocks_t &covar_shocks_arg,
                  const covar_and_corr_shocks_t &corr_shocks_arg,
                  const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
};

class MShocksStatement : public AbstractShocksStatement
{
public:
  MShocksStatement(bool overwrite_arg,
                   const det_shocks_t &det_shocks_arg,
                   const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class ConditionalForecastPathsStatement : public Statement
{
private:
  const AbstractShocksStatement::det_shocks_t paths;
  int path_length;
public:
  ConditionalForecastPathsStatement(const AbstractShocksStatement::det_shocks_t &paths_arg);
  virtual void checkPass(ModFileStructure &mod_file_struct, WarningConsolidation &warnings);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class MomentCalibration : public Statement
{
public:
  struct Constraint
  {
    int endo1, endo2;
    string lags;
    string lower_bound, upper_bound;
  };
  typedef vector<Constraint> constraints_t;
private:
  constraints_t constraints;
  const SymbolTable &symbol_table;
public:
  MomentCalibration(const constraints_t &constraints_arg,
                    const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

class IrfCalibration : public Statement
{
public:
  struct Constraint
  {
    int endo;
    int exo;
    string periods, lower_bound, upper_bound;
  };
  typedef vector<Constraint> constraints_t;
private:
  constraints_t constraints;
  const SymbolTable &symbol_table;
public:
  IrfCalibration(const constraints_t &constraints_arg,
                 const SymbolTable &symbol_table_arg);
  virtual void writeOutput(ostream &output, const string &basename) const;
};

#endif
