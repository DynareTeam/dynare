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

#ifndef _MOD_FILE_HH
#define _MOD_FILE_HH

using namespace std;

#include <ostream>
#include <ctime>

#include "SymbolTable.hh"
#include "NumericalConstants.hh"
#include "NumericalInitialization.hh"
#include "StaticModel.hh"
#include "DynamicModel.hh"
#include "SteadyStateModel.hh"
#include "Statement.hh"
#include "ExternalFunctionsTable.hh"
#include "ConfigFile.hh"
#include "WarningConsolidation.hh"
#include "FileOutputType.hh"

//! The abstract representation of a "mod" file
class ModFile
{
public:
  ModFile(WarningConsolidation &warnings_arg);
  ~ModFile();
  //! Symbol table
  SymbolTable symbol_table;
  //! External Functions table
  ExternalFunctionsTable external_functions_table;
  //! Numerical constants table
  NumericalConstants num_constants;
  //! Expressions outside model block
  DataTree expressions_tree;
  //! Dynamic model, as declared in the "model" block
  DynamicModel dynamic_model;
  //! A copy of Dynamic model, for testing trends declared by user
  DynamicModel trend_dynamic_model;
  //! A model in which to create the FOC for the ramsey problem
  DynamicModel ramsey_FOC_equations_dynamic_model;
  //! Static model, as derived from the "model" block when leads and lags have been removed
  StaticModel static_model;
  //! Static model, as declared in the "steady_state_model" block if present
  SteadyStateModel steady_state_model;
  //! Option linear
  bool linear;

  //! Is the model block decomposed?
  bool block;

  //! Is the model stored in bytecode format (byte_code=true) or in a M-file (byte_code=false)
  bool byte_code;

  //! Is the model stored in a MEX file ? (option "use_dll" of "model")
  bool use_dll;

  //! Is the static model have to computed (no_static=false) or not (no_static=true). Option of 'model'
  bool no_static;

  //! Is the 'differentiate_forward_vars' option used?
  bool differentiate_forward_vars;

  /*! If the 'differentiate_forward_vars' option is used, contains the set of
      endogenous with respect to which to do the transformation;
      if empty, means that the transformation must be applied to all endos
      with a lead */
  vector<string> differentiate_forward_vars_subset;

  //! Are nonstationary variables present ?
  bool nonstationary_variables;

  //! Global evaluation context
  /*! Filled using initval blocks and parameters initializations */
  eval_context_t global_eval_context;

  //! Stores the original number of equations in the model_block
  int ramsey_model_orig_eqn_nbr;

  //! Stores the list of extra files to be transefered during a parallel run
  /*! (i.e. option parallel_local_files of model block) */
  vector<string> parallel_local_files;

private:
  //! List of statements
  vector<Statement *> statements;
  //! Structure of the mod file
  ModFileStructure mod_file_struct;
  //! Warnings Encountered
  WarningConsolidation &warnings;

public:
  //! Add a statement
  void addStatement(Statement *st);
  //! Add a statement at the front of the statements vector
  void addStatementAtFront(Statement *st);
  //! Evaluate all the statements
  /*! \param warn_uninit Should a warning be displayed for uninitialized endogenous/exogenous/parameters ? */
  void evalAllExpressions(bool warn_uninit);
  //! Do some checking and fills mod_file_struct
  /*! \todo add check for number of equations and endogenous if ramsey_policy is present */
  void checkPass();
  //! Perform some transformations on the model (creation of auxiliary vars and equations)
  void transformPass(bool nostrict);
  //! Execute computations
  /*! \param no_tmp_terms if true, no temporary terms will be computed in the static and dynamic files */
  void computingPass(bool no_tmp_terms, FileOutputType output);
  //! Writes Matlab/Octave output files
  /*!
    \param basename The base name used for writing output files. Should be the name of the mod file without its extension
    \param clear_all Should a "clear all" instruction be written to output ?
    \param console Are we in console mode ?
    \param nograph Should we build the figures?
    \param nointeractive Should Dynare request user input?
    \param cygwin Should the MEX command of use_dll be adapted for Cygwin?
    \param msvc Should the MEX command of use_dll be adapted for MSVC?
  */
  void writeOutputFiles(const string &basename, bool clear_all, bool no_log, bool no_warn, bool console, bool nograph, bool nointeractive, const ConfigFile &config_file
#if defined(_WIN32) || defined(__CYGWIN32__)
                        , bool cygwin, bool msvc
#endif
                        ) const;
  //! Writes C output files only => No further Matlab processing
  void writeCOutputFiles(const string &basename) const;
  void writeModelC(const string &basename, bool cuda) const;
  //! Writes Cpp output files only => No further Matlab processing
  void writeCCOutputFiles(const string &basename) const;
  void writeModelCC(const string &basename, bool cuda) const;
  void writeExternalFiles(const string &basename, FileOutputType output, bool cuda) const;
};

#endif // ! MOD_FILE_HH
