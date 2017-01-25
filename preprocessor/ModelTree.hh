/*
 * Copyright (C) 2003-2017 Dynare Team
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

#ifndef _MODELTREE_HH
#define _MODELTREE_HH

using namespace std;

#include <string>
#include <vector>
#include <deque>
#include <map>
#include <ostream>

#include "DataTree.hh"
#include "ExtendedPreprocessorTypes.hh"

//! Vector describing equations: BlockSimulationType, if BlockSimulationType == EVALUATE_s then a expr_t on the new normalized equation
typedef vector<pair<EquationType, expr_t > > equation_type_and_normalized_equation_t;

//! Vector describing variables: max_lag in the block, max_lead in the block
typedef vector<pair< int, int> > lag_lead_vector_t;

//! for each block contains pair< pair<Simulation_Type, first_equation>, pair < Block_Size, Recursive_part_Size > >
typedef vector<pair< pair< BlockSimulationType, int>, pair<int, int> > > block_type_firstequation_size_mfs_t;

//! for a block contains derivatives pair< pair<block_equation_number, block_variable_number> , pair<lead_lag, expr_t> >
typedef vector< pair<pair<int, int>, pair< int, expr_t > > > block_derivatives_equation_variable_laglead_nodeid_t;

//! for all blocks derivatives description
typedef vector<block_derivatives_equation_variable_laglead_nodeid_t> blocks_derivatives_t;

//! Shared code for static and dynamic models
class ModelTree : public DataTree
{
  friend class DynamicModel;
  friend class StaticModel;
protected:
  //! Stores declared and generated auxiliary equations
  vector<BinaryOpNode *> equations;

  //! Stores line numbers of declared equations; -1 means undefined
  vector<int> equations_lineno;

  //! Only stores generated auxiliary equations, in an order meaningful for evaluation
  deque<BinaryOpNode *> aux_equations;

  //! Stores equation tags
  vector<pair<int, pair<string, string> > > equation_tags;

  //! Number of non-zero derivatives
  int NNZDerivatives[3];

  typedef map<pair<int, int>, expr_t> first_derivatives_t;
  //! First order derivatives
  /*! First index is equation number, second is variable w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Variable indices are those of the getDerivID() method.
  */
  first_derivatives_t first_derivatives;

  typedef map<pair<int, pair<int, int> >, expr_t> second_derivatives_t;
  //! Second order derivatives
  /*! First index is equation number, second and third are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only second order derivatives where var1 >= var2 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  second_derivatives_t second_derivatives;

  typedef map<pair<int, pair<int, pair<int, int> > >, expr_t> third_derivatives_t;
  //! Third order derivatives
  /*! First index is equation number, second, third and fourth are variables w.r. to which is computed the derivative.
    Only non-null derivatives are stored in the map.
    Contains only third order derivatives where var1 >= var2 >= var3 (for obvious symmetry reasons).
    Variable indices are those of the getDerivID() method.
  */
  third_derivatives_t third_derivatives;

  //! Derivatives of the residuals w.r. to parameters
  /*! First index is equation number, second is parameter.
    Only non-null derivatives are stored in the map.
    Parameter indices are those of the getDerivID() method.
  */
  first_derivatives_t residuals_params_derivatives;

  //! Cross reference information
  map<int, ExprNode::EquationInfo> xrefs;
  map<int, set<int> > xref_param;
  map<int, set<int> > xref_endo;
  map<int, set<int> > xref_exo;
  map<int, set<int> > xref_exo_det;

  //! Second derivatives of the residuals w.r. to parameters
  /*! First index is equation number, second and third indeces are parameters.
    Only non-null derivatives are stored in the map.
    Parameter indices are those of the getDerivID() method.
  */
  second_derivatives_t residuals_params_second_derivatives;

  //! Derivatives of the jacobian w.r. to parameters
  /*! First index is equation number, second is endo/exo/exo_det variable, and third is parameter.
    Only non-null derivatives are stored in the map.
    Variable and parameter indices are those of the getDerivID() method.
  */
  second_derivatives_t jacobian_params_derivatives;

  //! Second derivatives of the jacobian w.r. to parameters
  /*! First index is equation number, second is endo/exo/exo_det variable, and third and fourth are parameters.
    Only non-null derivatives are stored in the map.
    Variable and parameter indices are those of the getDerivID() method.
  */
  third_derivatives_t jacobian_params_second_derivatives;

  //! Derivatives of the hessian w.r. to parameters
  /*! First index is equation number, first and second are endo/exo/exo_det variable, and third is parameter.
    Only non-null derivatives are stored in the map.
    Variable and parameter indices are those of the getDerivID() method.
  */
  third_derivatives_t hessian_params_derivatives;


  //! Temporary terms for the static/dynamic file (those which will be noted Txxxx)
  temporary_terms_t temporary_terms;
  temporary_terms_t temporary_terms_res;
  temporary_terms_t temporary_terms_g1;
  temporary_terms_t temporary_terms_g2;
  temporary_terms_t temporary_terms_g3;

  //! Temporary terms for the file containing parameters derivatives
  temporary_terms_t params_derivs_temporary_terms;
  temporary_terms_t params_derivs_temporary_terms_res;
  temporary_terms_t params_derivs_temporary_terms_g1;
  temporary_terms_t params_derivs_temporary_terms_res2;
  temporary_terms_t params_derivs_temporary_terms_g12;
  temporary_terms_t params_derivs_temporary_terms_g2;


  //! Trend variables and their growth factors
  map<int, expr_t> trend_symbols_map;

  //! for all trends; the boolean is true if this is a log-trend, false otherwise
  typedef map<int, pair<bool, expr_t> > nonstationary_symbols_map_t;

  //! Nonstationary variables and their deflators
  nonstationary_symbols_map_t nonstationary_symbols_map;

  //! vector of block reordered variables and equations
  vector<int> equation_reordered, variable_reordered, inv_equation_reordered, inv_variable_reordered;

  //! the file containing the model and the derivatives code
  ofstream code_file;

  //! Computes 1st derivatives
  /*! \param vars the derivation IDs w.r. to which compute the derivatives */
  void computeJacobian(const set<int> &vars);
  //! Computes 2nd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 1st derivatives */
  void computeHessian(const set<int> &vars);
  //! Computes 3rd derivatives
  /*! \param vars the derivation IDs w.r. to which derive the 2nd derivatives */
  void computeThirdDerivatives(const set<int> &vars);
  //! Computes derivatives of the Jacobian and Hessian w.r. to parameters
  void computeParamsDerivatives(int paramsDerivsOrder);
  //! Write derivative of an equation w.r. to a variable
  void writeDerivative(ostream &output, int eq, int symb_id, int lag, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;
  //! Computes temporary terms (for all equations and derivatives)
  void computeTemporaryTerms(bool is_matlab);
  //! Computes temporary terms for the file containing parameters derivatives
  void computeParamsDerivativesTemporaryTerms();
//! Writes temporary terms
  void writeTemporaryTerms(const temporary_terms_t &tt, const temporary_terms_t &ttm1, ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const;
  //! Compiles temporary terms
  void compileTemporaryTerms(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, map_idx_t map_idx, bool dynamic, bool steady_dynamic) const;
  //! Adds informations for simulation in a binary file
  void Write_Inf_To_Bin_File(const string &basename, int &u_count_int, bool &file_open, bool is_two_boundaries, int block_mfs) const;
  //! Fixes output when there are more than 32 nested parens, Issue #1201
  void fixNestedParenthesis(ostringstream &output, map<string, string> &tmp_paren_vars, bool &message_printed) const;
  //! Tests if string contains more than 32 nested parens, Issue #1201
  bool testNestedParenthesis(const string &str) const;
  //! Writes model local variables
  /*! No temporary term is used in the output, so that local parameters declarations can be safely put before temporary terms declaration in the output files */
  void writeModelLocalVariables(ostream &output, ExprNodeOutputType output_type, deriv_node_temp_terms_t &tef_terms) const;
  //! Writes model equations
  void writeModelEquations(ostream &output, ExprNodeOutputType output_type) const;
  //! Compiles model equations
  void compileModelEquations(ostream &code_file, unsigned int &instruction_number, const temporary_terms_t &tt, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const;

  //! Writes LaTeX model file
  void writeLatexModelFile(const string &basename, ExprNodeOutputType output_type) const;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is equation number, second index is endogenous type specific ID */
  typedef map<pair<int, int>, double> jacob_map_t;

  //! Sparse matrix of double to store the values of the Jacobian
  /*! First index is lag, second index is equation number, third index is endogenous type specific ID */
  typedef map<pair<int, pair<int, int> >, expr_t> dynamic_jacob_map_t;

  //! Normalization of equations
  /*! Maps endogenous type specific IDs to equation numbers */
  vector<int> endo2eq;

  //! number of equation in the prologue and in the epilogue
  unsigned int epilogue, prologue;

  //! for each block contains pair< max_lag, max_lead>
  lag_lead_vector_t block_lag_lead;

  //! Compute the matching between endogenous and variable using the jacobian contemporaneous_jacobian
  /*!
    \param contemporaneous_jacobian Jacobian used as an incidence matrix: all elements declared in the map (even if they are zero), are used as vertices of the incidence matrix
    \return True if a complete normalization has been achieved
  */
  bool computeNormalization(const jacob_map_t &contemporaneous_jacobian, bool verbose);

  //! Try to compute the matching between endogenous and variable using a decreasing cutoff
  /*!
    Applied to the jacobian contemporaneous_jacobian and stop when a matching is found.
    If no matching is found using a strictly positive cutoff, then a zero cutoff is applied (i.e. use a symbolic normalization); in that case, the method adds zeros in the jacobian matrices to reflect all the edges in the symbolic incidence matrix.
    If no matching is found with a zero cutoff close to zero an error message is printout.
  */
  void computeNonSingularNormalization(jacob_map_t &contemporaneous_jacobian, double cutoff, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian);

  //! Try to normalized each unnormalized equation (matched endogenous variable only on the LHS)
  void computeNormalizedEquations(multimap<int, int> &endo2eqs) const;
  //! Compute cross references
  void computeXrefs();
  //! Help computeXrefs to compute the reverse references (i.e. param->eqs, endo->eqs, etc)
  void computeRevXref(map<int, set<int> > &xrefset, const set<int> &eiref, int eqn);
  //! Write cross references
  void writeXrefs(ostream &output) const;
  //! Write reverse cross references
  void writeRevXrefs(ostream &output, const map<int, set<int> > &xrefmap, const string &type) const;
  //! Evaluate the jacobian and suppress all the elements below the cutoff
  void evaluateAndReduceJacobian(const eval_context_t &eval_context, jacob_map_t &contemporaneous_jacobian, jacob_map_t &static_jacobian, dynamic_jacob_map_t &dynamic_jacobian, double cutoff, bool verbose);
  //! Search the equations and variables belonging to the prologue and the epilogue of the model
  void computePrologueAndEpilogue(const jacob_map_t &static_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered);
  //! Determine the type of each equation of model and try to normalized the unnormalized equation using computeNormalizedEquations
  equation_type_and_normalized_equation_t equationTypeDetermination(const map<pair<int, pair<int, int> >, expr_t> &first_order_endo_derivatives, const vector<int> &Index_Var_IM, const vector<int> &Index_Equ_IM, int mfs) const;
  //! Compute the block decomposition and for a non-recusive block find the minimum feedback set
  void computeBlockDecompositionAndFeedbackVariablesForEachBlock(const jacob_map_t &static_jacobian, const dynamic_jacob_map_t &dynamic_jacobian, vector<int> &equation_reordered, vector<int> &variable_reordered, vector<pair<int, int> > &blocks, const equation_type_and_normalized_equation_t &Equation_Type, bool verbose_, bool select_feedback_variable, int mfs, vector<int> &inv_equation_reordered, vector<int> &inv_variable_reordered, lag_lead_vector_t &equation_lag_lead, lag_lead_vector_t &variable_lag_lead_t, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed) const;
  //! Reduce the number of block merging the same type equation in the prologue and the epilogue and determine the type of each block
  block_type_firstequation_size_mfs_t reduceBlocksAndTypeDetermination(const dynamic_jacob_map_t &dynamic_jacobian, vector<pair<int, int> > &blocks, const equation_type_and_normalized_equation_t &Equation_Type, const vector<int> &variable_reordered, const vector<int> &equation_reordered, vector<unsigned int> &n_static, vector<unsigned int> &n_forward, vector<unsigned int> &n_backward, vector<unsigned int> &n_mixed, vector<pair< pair<int, int>, pair<int, int> > > &block_col_type);
  //! Determine the maximum number of lead and lag for the endogenous variable in a bloc
  void getVariableLeadLagByBlock(const dynamic_jacob_map_t &dynamic_jacobian, const vector<int> &components_set, int nb_blck_sim, lag_lead_vector_t &equation_lead_lag, lag_lead_vector_t &variable_lead_lag, const vector<int> &equation_reordered, const vector<int> &variable_reordered) const;
  //! Print an abstract of the block structure of the model
  void printBlockDecomposition(const vector<pair<int, int> > &blocks) const;
  //! Determine for each block if it is linear or not
  vector<bool> BlockLinear(const blocks_derivatives_t &blocks_derivatives, const vector<int> &variable_reordered) const;

  //! Determine the simulation type of each block
  virtual BlockSimulationType getBlockSimulationType(int block_number) const = 0;
  //! Return the number of blocks
  virtual unsigned int getNbBlocks() const = 0;
  //! Return the first equation number of a block
  virtual unsigned int getBlockFirstEquation(int block_number) const = 0;
  //! Return the size of the block block_number
  virtual unsigned int getBlockSize(int block_number) const = 0;
  //! Return the number of exogenous variable in the block block_number
  virtual unsigned int getBlockExoSize(int block_number) const = 0;
  //! Return the number of colums in the jacobian matrix for exogenous variable in the block block_number
  virtual unsigned int getBlockExoColSize(int block_number) const = 0;
  //! Return the number of feedback variable of the block block_number
  virtual unsigned int getBlockMfs(int block_number) const = 0;
  //! Return the maximum lag in a block
  virtual unsigned int getBlockMaxLag(int block_number) const = 0;
  //! Return the maximum lead in a block
  virtual unsigned int getBlockMaxLead(int block_number) const = 0;
  inline void setBlockLeadLag(int block, int max_lag, int max_lead) 
    {
       block_lag_lead[block] = make_pair(max_lag, max_lead);
    };
  
  //! Return the type of equation (equation_number) belonging to the block block_number
  virtual EquationType getBlockEquationType(int block_number, int equation_number) const = 0;
  //! Return true if the equation has been normalized
  virtual bool isBlockEquationRenormalized(int block_number, int equation_number) const = 0;
  //! Return the expr_t of the equation equation_number belonging to the block block_number
  virtual expr_t getBlockEquationExpr(int block_number, int equation_number) const = 0;
  //! Return the expr_t of the renormalized equation equation_number belonging to the block block_number
  virtual expr_t getBlockEquationRenormalizedExpr(int block_number, int equation_number) const = 0;
  //! Return the original number of equation equation_number belonging to the block block_number
  virtual int getBlockEquationID(int block_number, int equation_number) const = 0;
  //! Return the original number of variable variable_number belonging to the block block_number
  virtual int getBlockVariableID(int block_number, int variable_number) const = 0;
  //! Return the original number of the exogenous variable varexo_number belonging to the block block_number
  virtual int getBlockVariableExoID(int block_number, int variable_number) const = 0;
  //! Return the position of equation_number in the block number belonging to the block block_number
  virtual int getBlockInitialEquationID(int block_number, int equation_number) const = 0;
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialVariableID(int block_number, int variable_number) const = 0;
  //! Return the position of variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialExogenousID(int block_number, int variable_number) const = 0;
  //! Return the position of the deterministic exogenous variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialDetExogenousID(int block_number, int variable_number) const = 0;
  //! Return the position of the other endogenous variable_number in the block number belonging to the block block_number
  virtual int getBlockInitialOtherEndogenousID(int block_number, int variable_number) const = 0;
  //! Initialize equation_reordered & variable_reordered
  void initializeVariablesAndEquations();
public:
  ModelTree(SymbolTable &symbol_table_arg, NumericalConstants &num_constants_arg, ExternalFunctionsTable &external_functions_table_arg);
  //! Absolute value under which a number is considered to be zero
  double cutoff;
  //! Compute the minimum feedback set
  /*!   0 : all endogenous variables are considered as feedback variables
    1 : the variables belonging to non normalized equation are considered as feedback variables
    2 : the variables belonging to a non linear equation are considered as feedback variables
    3 : the variables belonging to a non normalizable non linear equation are considered as feedback variables
    default value = 0 */
  int mfs;
  //! Declare a node as an equation of the model; also give its line number
  void addEquation(expr_t eq, int lineno);
  //! Declare a node as an equation of the model, also giving its tags
  void addEquation(expr_t eq, int lineno, vector<pair<string, string> > &eq_tags);
  //! Declare a node as an auxiliary equation of the model, adding it at the end of the list of auxiliary equations
  void addAuxEquation(expr_t eq);
  //! Returns the number of equations in the model
  int equation_number() const;
  //! Adds a trend variable with its growth factor
  void addTrendVariables(vector<int> trend_vars, expr_t growth_factor) throw (TrendException);
  //! Adds a nonstationary variables with their (common) deflator
  void addNonstationaryVariables(vector<int> nonstationary_vars, bool log_deflator, expr_t deflator) throw (TrendException);
  //! Is a given variable non-stationary?
  bool isNonstationary(int symb_id) const;
  void set_cutoff_to_zero();
  //! Helper for writing the Jacobian elements in MATLAB and C
  /*! Writes either (i+1,j+1) or [i+j*no_eq] */
  void jacobianHelper(ostream &output, int eq_nb, int col_nb, ExprNodeOutputType output_type) const;
  //! Helper for writing the sparse Hessian or third derivatives in MATLAB and C
  /*! If order=2, writes either v2(i+1,j+1) or v2[i+j*NNZDerivatives[1]]
    If order=3, writes either v3(i+1,j+1) or v3[i+j*NNZDerivatives[2]] */
  void sparseHelper(int order, ostream &output, int row_nb, int col_nb, ExprNodeOutputType output_type) const;
  inline static std::string
  c_Equation_Type(int type)
  {
    char c_Equation_Type[4][13] =
      {
        "E_UNKNOWN   ",
        "E_EVALUATE  ",
        "E_EVALUATE_S",
        "E_SOLVE     "
      };
    return (c_Equation_Type[type]);
  };

  inline static std::string
  BlockType0(BlockType type)
  {
    switch (type)
      {
      case SIMULTANS:
        return ("SIMULTANEOUS TIME SEPARABLE  ");
        break;
      case PROLOGUE:
        return ("PROLOGUE                     ");
        break;
      case EPILOGUE:
        return ("EPILOGUE                     ");
        break;
      case SIMULTAN:
        return ("SIMULTANEOUS TIME UNSEPARABLE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };

  inline static std::string
  BlockSim(int type)
  {
    switch (type)
      {
      case EVALUATE_FORWARD:
        return ("EVALUATE FORWARD             ");
        break;
      case EVALUATE_BACKWARD:
        return ("EVALUATE BACKWARD            ");
        break;
      case SOLVE_FORWARD_SIMPLE:
        return ("SOLVE FORWARD SIMPLE         ");
        break;
      case SOLVE_BACKWARD_SIMPLE:
        return ("SOLVE BACKWARD SIMPLE        ");
        break;
      case SOLVE_TWO_BOUNDARIES_SIMPLE:
        return ("SOLVE TWO BOUNDARIES SIMPLE  ");
        break;
      case SOLVE_FORWARD_COMPLETE:
        return ("SOLVE FORWARD COMPLETE       ");
        break;
      case SOLVE_BACKWARD_COMPLETE:
        return ("SOLVE BACKWARD COMPLETE      ");
        break;
      case SOLVE_TWO_BOUNDARIES_COMPLETE:
        return ("SOLVE TWO BOUNDARIES COMPLETE");
        break;
      default:
        return ("UNKNOWN                      ");
        break;
      }
  };
};

#endif
