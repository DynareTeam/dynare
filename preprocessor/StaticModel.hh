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

#ifndef _STATIC_MODEL_HH
#define _STATIC_MODEL_HH

using namespace std;

#include <fstream>

#include "ModelTree.hh"

//! Stores a static model, as derived from the "model" block when leads and lags have been removed
class StaticModel : public ModelTree
{
private:
  //! Temporary terms for the file containing parameters dervicatives
  temporary_terms_t params_derivs_temporary_terms;

  //! global temporary terms for block decomposed models
  vector<vector<temporary_terms_t> > v_temporary_terms;

  //! local temporary terms for block decomposed models
  vector<vector<temporary_terms_t> > v_temporary_terms_local;

  vector<temporary_terms_inuse_t> v_temporary_terms_inuse;

  typedef map< pair< int, pair< int, int> >, expr_t> first_chain_rule_derivatives_t;
  first_chain_rule_derivatives_t first_chain_rule_derivatives;

  //! Writes static model file (standard Matlab version)
  void writeStaticMFile(const string &static_basename) const;

  //! Writes static model file (C version)
  void writeStaticCFile(const string &func_name) const;

  //! Writes the static model equations and its derivatives
  void writeStaticModel(ostream &StaticOutput, bool use_dll) const;

  //! Writes the static function calling the block to solve (Matlab version)
  void writeStaticBlockMFSFile(const string &basename) const;

  //! Writes the Block reordred structure of the model in M output
  void writeModelEquationsOrdered_M(const string &dynamic_basename) const;

  //! Writes the code of the Block reordred structure of the model in virtual machine bytecode
  void writeModelEquationsCode_Block(const string file_name, const string bin_basename, map_idx_t map_idx, vector<map_idx_t> map_idx2) const;

  //! Writes the code of the model in virtual machine bytecode
  void writeModelEquationsCode(const string file_name, const string bin_basename, map_idx_t map_idx) const;

  //! Computes jacobian and prepares for equation normalization
  /*! Using values from initval/endval blocks and parameter initializations:
    - computes the jacobian for the model w.r. to contemporaneous variables
    - removes edges of the incidence matrix when derivative w.r. to the corresponding variable is too close to zero (below the cutoff)
  */
  void evaluateJacobian(const eval_context_t &eval_context, jacob_map_t *j_m, bool dynamic);

  map_idx_t map_idx;

  vector<map_idx_t> map_idx2;

  //! sorts the temporary terms in the blocks order
  void computeTemporaryTermsOrdered();
  //! creates a mapping from the index of temporary terms to a natural index
  void computeTemporaryTermsMapping(temporary_terms_t &temporary_terms, map_idx_t &map_idx);

  //! Write derivative code of an equation w.r. to a variable
  void compileDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int symb_id, map_idx_t &map_idx, temporary_terms_t temporary_terms) const;
  //! Write chain rule derivative code of an equation w.r. to a variable
  void compileChainRuleDerivative(ofstream &code_file, unsigned int &instruction_number, int eq, int var, int lag, map_idx_t &map_idx, temporary_terms_t temporary_terms) const;

  //! Get the type corresponding to a derivation ID
  virtual SymbolType getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the lag corresponding to a derivation ID
  virtual int getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Get the symbol ID corresponding to a derivation ID
  virtual int getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException);
  //! Compute the column indices of the static Jacobian
  void computeStatJacobianCols();
  //! return a map on the block jacobian
  map<pair<pair<int, pair<int, int> >, pair<int, int> >, int> get_Derivatives(int block);
  //! Computes chain rule derivatives of the Jacobian w.r. to endogenous variables
  void computeChainRuleJacobian(blocks_derivatives_t &blocks_derivatives);
  //! Collect only the first derivatives
  map<pair<int, pair<int, int> >, expr_t> collect_first_order_derivatives_endogenous();

  //! Write chain rule derivative of a recursive equation w.r. to a variable
  void writeChainRuleDerivative(ostream &output, int eq, int var, int lag, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;

  //! Collecte the derivatives w.r. to endogenous of the block, to endogenous of previouys blocks and to exogenous
  void collect_block_first_order_derivatives();

protected:
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

  //!List for each block and for each lag-leag all the other endogenous variables and exogenous variables
  typedef set<int> var_t;
  typedef map<int, var_t> lag_var_t;
  vector<lag_var_t> other_endo_block, exo_block, exo_det_block;

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
  StaticModel(SymbolTable &symbol_table_arg, NumericalConstants &num_constants, ExternalFunctionsTable &external_functions_table_arg);

  //! Writes information on block decomposition when relevant
  void writeM_Output(ostream &output, bool block) const;

  //! Execute computations (variable sorting + derivation)
  /*!
    \param eval_context evaluation context for normalization
    \param no_tmp_terms if true, no temporary terms will be computed in the static files
    \param hessian whether 2nd derivatives w.r. to exo, exo_det and endo should be computed
    \param paramsDerivatives whether 2nd derivatives w.r. to a pair (endo/exo/exo_det, parameter) should be computed
  */
  void computingPass(const eval_context_t &eval_context, bool no_tmp_terms, bool hessian, bool thirdDerivatices, bool paramsDerivatives, bool block, bool bytecode);

  //! Adds informations for simulation in a binary file for a block decomposed model
  void Write_Inf_To_Bin_File_Block(const string &static_basename, const string &bin_basename, const int &num,
                                   int &u_count_int, bool &file_open) const;

  //! Writes static model file
  void writeStaticFile(const string &basename, bool block, bool bytecode, bool use_dll) const;

  //! Writes file containing static parameters derivatives
  void writeParamsDerivativesFile(const string &basename) const;

  //! Writes LaTeX file with the equations of the static model
  void writeLatexFile(const string &basename) const;

  //! Writes initializations in oo_.steady_state or steady state file for the auxiliary variables
  void writeAuxVarInitval(ostream &output, ExprNodeOutputType output_type) const;

  //! Writes definition of the auxiliary variables in a M file
  void writeAuxVarRecursiveDefinitions(const string &basename) const;

  virtual int getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException);
  virtual void addAllParamDerivId(set<int> &deriv_id_set);

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
  virtual unsigned int getBlockExoSize(int block_number) const
  {
    return 0;
  };
  //! Return the number of colums in the jacobian matrix for exogenous variable in the block block_number
  virtual unsigned int getBlockExoColSize(int block_number) const
  {
    return 0;
  }
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
    return 0;
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
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialExogenousID(int block_number, int variable_number) const
  {
    return -1;
  };
  //! Return the position of the deterministic exogenous variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialDetExogenousID(int block_number, int variable_number) const
  {
    return -1;
  };
  //! Return the position of the other endogenous variable_number in the block number belonging to the block block_number
  virtual int
  getBlockInitialOtherEndogenousID(int block_number, int variable_number) const
  {
    return -1;
  };
};

#endif
