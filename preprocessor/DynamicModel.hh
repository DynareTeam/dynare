/*
 * Copyright (C) 2003-2012 Dynare Team
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

#ifndef _DYNAMICMODEL_HH
#define _DYNAMICMODEL_HH

using namespace std;
#define ZERO_BAND 1e-8

#include <fstream>

#include "StaticModel.hh"

//! Stores a dynamic model
class DynamicModel : public ModelTree
{
private:
  typedef map<pair<int, int>, int> deriv_id_table_t;
  //! Maps a pair (symbol_id, lag) to a deriv ID
  deriv_id_table_t deriv_id_table;
  //! Maps a deriv ID to a pair (symbol_id, lag)
  vector<pair<int, int> > inv_deriv_id_table;

  //! Maps a deriv_id to the column index of the dynamic Jacobian
  /*! Contains only endogenous, exogenous and exogenous deterministic */
  map<int, int> dyn_jacobian_cols_table;

  //! Maximum lag and lead over all types of variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_lag, max_lead;
  //! Maximum lag and lead over endogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_endo_lag, max_endo_lead;
  //! Maximum lag and lead over exogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_exo_lag, max_exo_lead;
  //! Maximum lag and lead over deterministic exogenous variables (positive values)
  /*! Set by computeDerivIDs() */
  int max_exo_det_lag, max_exo_det_lead;

  //! Number of columns of dynamic jacobian
  /*! Set by computeDerivID()s and computeDynJacobianCols() */
  int dynJacobianColsNbr;
  //! Temporary terms for block decomposed models
  vector< vector<temporary_terms_t> > v_temporary_terms;

  vector<temporary_terms_inuse_t> v_temporary_terms_inuse;

  //! Store the derivatives or the chainrule derivatives:map<pair< equation, pair< variable, lead_lag >, expr_t>
  typedef map< pair< int, pair< int, int> >, expr_t> first_chain_rule_derivatives_t;
  first_chain_rule_derivatives_t first_chain_rule_derivatives;

  //! Writes dynamic model file (Matlab version)
  void writeDynamicMFile(const string &dynamic_basename) const;
  //! Writes dynamic model file (C version)
  /*! \todo add third derivatives handling */
  void writeDynamicCFile(const string &dynamic_basename, const int order) const;
  //! Writes dynamic model file when SparseDLL option is on
  void writeSparseDynamicMFile(const string &dynamic_basename, const string &basename) const;
  //! Writes the dynamic model equations and its derivatives
  /*! \todo add third derivatives handling in C output */
  void writeDynamicModel(ostream &DynamicOutput, bool use_dll) const;
  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(const string &dynamic_basename) const;
  //! Writes the code of the Block reordred structure of the model in virtual machine bytecode
  void writeModelEquationsCode_Block(string &file_name, const string &bin_basename, const map_idx_t &map_idx) const;
  //! Writes the code of the model in virtual machine bytecode
  void writeModelEquationsCode(string &file_name, const string &bin_basename, const map_idx_t &map_idx) const;

  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  //void evaluateJacobian(const eval_context_t &eval_context, jacob_map *j_m, bool dynamic);

  //! return a map on the block jacobian
  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> get_Derivatives(int block);
  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  void computeChainRuleJacobian(blocks_derivatives_t &blocks_derivatives);

  string reform(string name) const;
  map_idx_t map_idx;

  //! sorts the temporary terms in the blocks order
  void computeTemporaryTermsOrdered();

  //! creates a mapping from the index of temporary terms to a natural index
  void computeTemporaryTermsMapping();
  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, int lag, const map_idx_t &map_idx) const;
  //! Write chain rule derivative code of an equation w.r. to a variable
  void compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int var, int lag, const map_idx_t &map_idx) const;

  //! Get the type corresponding to a derivation ID
  virtual SymbolType getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the lag corresponding to a derivation ID
  virtual int getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the symbol ID corresponding to a derivation ID
  virtual int getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Compute the column indices of the dynamic Jacobian
  void computeDynJacobianCols(bool jacobianExo);
  //! Computes derivatives of the Jacobian w.r. to trend vars and tests that they are equal to zero
  void testTrendDerivativesEqualToZero(const eval_context_t &eval_context);
  //! Collect only the first derivatives
  map<pair<int, pair<int, int> >, expr_t> collect_first_order_derivatives_endogenous();

  //! Allocates the derivation IDs for all dynamic variables of the model
  /*! Also computes max_{endo,exo}_{lead_lag}, and initializes dynJacobianColsNbr to the number of dynamic endos */
  void computeDerivIDs();

  //! Write chain rule derivative of a recursive equation w.r. to a variable
  void writeChainRuleDerivative(ostream &output, int eq, int var, int lag, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;

  //! Collecte the derivatives w.r. to endogenous of the block, to endogenous of previouys blocks and to exogenous
  void collect_block_first_order_derivatives();

  //! Collecte the informations about exogenous, deterministic exogenous and endogenous from the previous block for each block
  void collectBlockVariables();

  //! Factorized code for substitutions of leads/lags
  /*! \param[in] type determines which type of variables is concerned */
  void substituteLeadLagInternal(aux_var_t type, bool deterministic_model);

private:
  //! Indicate if the temporary terms are computed for the overall model (true) or not (false). Default value true
  bool global_temporary_terms;

  //! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a expr_t on the new normalized equation
  equation_type_and_normalized_equation_t equation_type_and_normalized_equation;

  //! for each block contains pair< Simulation_Type, pair < Block_Size, Recursive_part_Size > >
  block_type_firstequation_size_mfs_t block_type_firstequation_size_mfs;

  //! for all blocks derivatives description
  blocks_derivatives_t blocks_derivatives;

  //! The jacobian without the elements below the cutoff
  dynamic_jacob_map_t dynamic_jacobian;

  //! Vector indicating if the block is linear in endogenous variable (true) or not (false)
  vector<bool> blocks_linear;

  //! Map the derivatives for a block pair<lag, make_pair(make_pair(eq, var)), expr_t>
  typedef map<pair< int, pair<int, int> >, expr_t> derivative_t;
  //! Vector of derivative for each blocks
  vector<derivative_t> derivative_endo, derivative_other_endo, derivative_exo, derivative_exo_det;

  //!List for each block and for each lag-lead all the other endogenous variables and exogenous variables
  typedef set<int> var_t;
  typedef map<int, var_t> lag_var_t;
  vector<lag_var_t> other_endo_block, exo_block, exo_det_block;

  //!List for each block the exogenous variables
  vector<pair<var_t, int> > block_var_exo;

  map< int, map<int, int> > block_exo_index, block_det_exo_index, block_other_endo_index;

  //! for each block described the number of static, forward, backward and mixed variables in the block
  /*! pair< pair<static, forward>, pair<backward,mixed> > */
  vector<pair< pair<int, int>, pair<int, int> > > block_col_type;

  //! List for each variable its block number and its maximum lag and lead inside the block
  vector<pair<int, pair<int, int> > > variable_block_lead_lag;
  //! List for each equation its block number
  vector<int> equation_block;

  //!Maximum lead and lag for each block on endogenous of the block, endogenous of the previous blocks, exogenous and deterministic exogenous
  vector<pair<int, int> > endo_max_leadlag_block, other_endo_max_leadlag_block, exo_max_leadlag_block, exo_det_max_leadlag_block, max_leadlag_block;

public:
  DynamicModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg, ExternalFunctionsTable &external_functions_table_argx);
  //! Adds a variable node
  /*! This implementation allows for non-zero lag */
  virtual VariableNode *AddVariable(int symb_id, int lag = 0);
  
  //! Execute computations (variable sorting + derivation)
  /*!
    \param jacobianExo whether derivatives w.r. to exo and exo_det should be in the Jacobian (derivatives w.r. to endo are always computed)
    \param hessian whether 2nd derivatives w.r. to exo, exo_det and endo should be computed (implies jacobianExo = true)
    \param thirdDerivatives whether 3rd derivatives w.r. to endo/exo/exo_det should be computed (implies jacobianExo = true)
    \param paramsDerivatives whether 2nd derivatives w.r. to a pair (endo/exo/exo_det, parameter) should be computed (implies jacobianExo = true)
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the dynamic files
  */
  void computingPass(bool jacobianExo, bool hessian, bool thirdDerivatives, bool paramsDerivatives,
                     const eval_context_t &eval_context, bool no_tmp_terms, bool block, bool use_dll, bool bytecode);
  //! Writes model initialization and lead/lag incidence matrix to output
  void writeOutput(ostream &output, const string &basename, bool block, bool byte_code, bool use_dll, int order, bool estimation_present) const;

  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File_Block(const string &dynamic_basename, const string &bin_basename,
                                   const int &num, int &u_count_int, bool &file_open, bool is_two_boundaries) const;
  //! Writes dynamic model file
  void writeDynamicFile(const string &basename, bool block, bool bytecode, bool use_dll, int order) const;
  //! Writes file containing parameters derivatives
  void writeParamsDerivativesFile(const string &basename) const;
  //! Converts to static model (only the equations)
  /*! It assumes that the static model given in argument has just been allocated */
  void toStatic(StaticModel &static_model) const;

  //! Copies a dynamic model (only the equations)
  /*! It assumes that the dynamic model given in argument has just been allocated */
  void cloneDynamic(DynamicModel &dynamic_model) const;

  //! Replaces model equations with derivatives of Lagrangian w.r.t. endogenous
  void computeRamseyPolicyFOCs(const StaticModel &static_model);
  //! Replaces the model equations in dynamic_model with those in this model
  void replaceMyEquations(DynamicModel &dynamic_model) const;

  //! Writes LaTeX file with the equations of the dynamic model
  void writeLatexFile(const string &basename) const;

  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
  virtual int getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException);
  virtual void addAllParamDerivId(set<int> &deriv_id_set);
  void print_trend_vars();

  //! Returns true indicating that this is a dynamic model
  virtual bool
  isDynamic() const
  {
    return true;
  };

  //! Drive test of detrended equations
  void runTrendTest(const eval_context_t &eval_context);

  //! Transforms the model by removing all leads greater or equal than 2 on endos
  /*! Note that this can create new lags on endos and exos */
  void substituteEndoLeadGreaterThanTwo(bool deterministic_model);

  //! Transforms the model by removing all lags greater or equal than 2 on endos
  void substituteEndoLagGreaterThanTwo(bool deterministic_model);

  //! Transforms the model by removing all leads on exos
  /*! Note that this can create new lags on endos and exos */
  void substituteExoLead(bool deterministic_model);

  //! Transforms the model by removing all lags on exos
  void substituteExoLag(bool deterministic_model);

  //! Transforms the model by removing all oExpectation
  void substituteExpectation(bool partial_information_model);

  //! Transform the model by add auxiliary variables for log and pow expressions
  void substituteLogPow(void);

  //! Transforms the model by decreasing the lead/lag of predetermined variables in model equations by one
  void transformPredeterminedVariables();

  //! Transforms the model by removing trends specified by the user
  void detrendEquations();

  //! Transforms the model by replacing trend variables with a 1
  void removeTrendVariableFromEquations();

  //! Fills eval context with values of model local variables and auxiliary variables
  void fillEvalContext(eval_context_t &eval_context) const;

  //! Return the number of blocks
  virtual unsigned int
  getNbBlocks() const
  {
    return (block_type_firstequation_size_mfs.size());
  };
  //! Determine the simulation type of each block
  virtual BlockSimulationType
  getBlockSimulationType(int block_number) const
  {
    return (block_type_firstequation_size_mfs[block_number].first.first);
  };
  //! Return the first equation number of a block
  virtual unsigned int
  getBlockFirstEquation(int block_number) const
  {
    return (block_type_firstequation_size_mfs[block_number].first.second);
  };
  //! Return the size of the block block_number
  virtual unsigned int
  getBlockSize(int block_number) const
  {
    return (block_type_firstequation_size_mfs[block_number].second.first);
  };
  //! Return the number of exogenous variable in the block block_number
  virtual unsigned int
  getBlockExoSize(int block_number) const
  {
    return (block_var_exo[block_number].first.size());
  };
  //! Return the number of colums in the jacobian matrix for exogenous variable in the block block_number
  virtual unsigned int
  getBlockExoColSize(int block_number) const
  {
    return (block_var_exo[block_number].second);
  };
  //! Return the number of feedback variable of the block block_number
  virtual unsigned int
  getBlockMfs(int block_number) const
  {
    return (block_type_firstequation_size_mfs[block_number].second.second);
  };
  //! Return the maximum lag in a block
  virtual unsigned int
  getBlockMaxLag(int block_number) const
  {
    return (block_lag_lead[block_number].first);
  };
  //! Return the maximum lead in a block
  virtual unsigned int
  getBlockMaxLead(int block_number) const
  {
    return (block_lag_lead[block_number].second);
  };
  //! Return the type of equation (equation_number) belonging to the block block_number
  virtual EquationType
  getBlockEquationType(int block_number, int equation_number) const
  {
    return (equation_type_and_normalized_equation[equation_reordered[block_type_firstequation_size_mfs[block_number].first.second+equation_number]].first);
  };
  //! Return true if the equation has been normalized
  virtual bool
  isBlockEquationRenormalized(int block_number, int equation_number) const
  {
    return (equation_type_and_normalized_equation[equation_reordered[block_type_firstequation_size_mfs[block_number].first.second+equation_number]].first == E_EVALUATE_S);
  };
  //! Return the expr_t of the equation equation_number belonging to the block block_number
  virtual expr_t
  getBlockEquationExpr(int block_number, int equation_number) const
  {
    return (equations[equation_reordered[block_type_firstequation_size_mfs[block_number].first.second+equation_number]]);
  };
  //! Return the expr_t of the renormalized equation equation_number belonging to the block block_number
  virtual expr_t
  getBlockEquationRenormalizedExpr(int block_number, int equation_number) const
  {
    return (equation_type_and_normalized_equation[equation_reordered[block_type_firstequation_size_mfs[block_number].first.second+equation_number]].second);
  };
  //! Return the original number of equation equation_number belonging to the block block_number
  virtual int
  getBlockEquationID(int block_number, int equation_number) const
  {
    return (equation_reordered[block_type_firstequation_size_mfs[block_number].first.second+equation_number]);
  };
  //! Return the original number of variable variable_number belonging to the block block_number
  virtual int
  getBlockVariableID(int block_number, int variable_number) const
  {
    return (variable_reordered[block_type_firstequation_size_mfs[block_number].first.second+variable_number]);
  };
  //! Return the original number of the exogenous variable varexo_number belonging to the block block_number
  virtual int
  getBlockVariableExoID(int block_number, int variable_number) const
  {
    map<int, var_t>::const_iterator it = exo_block[block_number].find(variable_number);
    return (it->first);
  };
  //! Return the position of equation_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialEquationID(int block_number, int equation_number) const
  {
    return ((int) inv_equation_reordered[equation_number] - (int) block_type_firstequation_size_mfs[block_number].first.second);
  };
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialVariableID(int block_number, int variable_number) const
  {
    return ((int) inv_variable_reordered[variable_number] - (int) block_type_firstequation_size_mfs[block_number].first.second);
  };
  //! Return the block number containing the endogenous variable variable_number
  int
  getBlockVariableID(int variable_number) const
  {
    return (variable_block_lead_lag[variable_number].first);
  };
  //! Return the position of the exogenous variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialExogenousID(int block_number, int variable_number) const
  {
    map< int, map<int, int> >::const_iterator it = block_exo_index.find(block_number);
    if (it != block_exo_index.end())
      {
        map<int, int>::const_iterator it1 = it->second.find(variable_number);
        if (it1 != it->second.end())
          return it1->second;
        else
          return -1;
      }
    else
      return (-1);
  };
  //! Return the position of the deterministic exogenous variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialDetExogenousID(int block_number, int variable_number) const
  {
    map< int, map<int, int> >::const_iterator it = block_det_exo_index.find(block_number);
    if (it != block_det_exo_index.end())
      {
        map<int, int>::const_iterator it1 = it->second.find(variable_number);
        if (it1 != it->second.end())
          return it1->second;
        else
          return -1;
      }
    else
      return (-1);
  };
  //! Return the position of the other endogenous variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialOtherEndogenousID(int block_number, int variable_number) const
  {
    map< int, map<int, int> >::const_iterator it = block_other_endo_index.find(block_number);
    if (it != block_other_endo_index.end())
      {
        map<int, int>::const_iterator it1 = it->second.find(variable_number);
        if (it1 != it->second.end())
          return it1->second;
        else
          return -1;
      }
    else
      return (-1);
  };
  bool isModelLocalVariableUsed() const;
};

#endif
