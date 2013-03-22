/*
 * Copyright (C) 2007-2012 Dynare Team
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

#ifndef _EXPR_NODE_HH
#define _EXPR_NODE_HH

using namespace std;

#include <set>
#include <map>
#include <vector>
#include <ostream>

#include "SymbolTable.hh"
#include "CodeInterpreter.hh"
#include "ExternalFunctionsTable.hh"

class DataTree;
class VariableNode;
class BinaryOpNode;

typedef class ExprNode *expr_t;

struct Model_Block;

struct ExprNodeLess;

//! Type for set of temporary terms
/*! They are ordered by index number thanks to ExprNodeLess */
typedef set<expr_t, ExprNodeLess> temporary_terms_t;

//! set of temporary terms used in a block
typedef set<int> temporary_terms_inuse_t;

typedef map<int, int> map_idx_t;

//! Type for evaluation contexts
/*! The key is a symbol id. Lags are assumed to be null */
typedef map<int, double> eval_context_t;

//! Type for tracking first/second derivative functions that have already been written as temporary terms
typedef map<pair<int, vector<expr_t> >, int> deriv_node_temp_terms_t;

//! Possible types of output when writing ExprNode(s)
enum ExprNodeOutputType
  {
    oMatlabStaticModel,                           //!< Matlab code, static model
    oMatlabDynamicModel,                          //!< Matlab code, dynamic model
    oMatlabStaticModelSparse,                     //!< Matlab code, static block decomposed model
    oMatlabDynamicModelSparse,                    //!< Matlab code, dynamic block decomposed model
    oCDynamicModel,                               //!< C code, dynamic model
    oCStaticModel,                                //!< C code, static model
    oMatlabOutsideModel,                          //!< Matlab code, outside model block (for example in initval)
    oLatexStaticModel,                            //!< LaTeX code, static model
    oLatexDynamicModel,                           //!< LaTeX code, dynamic model
    oLatexDynamicSteadyStateOperator,             //!< LaTeX code, dynamic model, inside a steady state operator
    oMatlabDynamicSteadyStateOperator,            //!< Matlab code, dynamic model, inside a steady state operator
    oMatlabDynamicSparseSteadyStateOperator,      //!< Matlab code, dynamic block decomposed model, inside a steady state operator
    oCDynamicSteadyStateOperator,                 //!< C code, dynamic model, inside a steady state operator
    oSteadyStateFile                              //!< Matlab code, in the generated steady state file
  };

#define IS_MATLAB(output_type) ((output_type) == oMatlabStaticModel     \
                                || (output_type) == oMatlabDynamicModel \
                                || (output_type) == oMatlabOutsideModel \
                                || (output_type) == oMatlabStaticModelSparse \
                                || (output_type) == oMatlabDynamicModelSparse \
                                || (output_type) == oMatlabDynamicSteadyStateOperator \
                                || (output_type) == oMatlabDynamicSparseSteadyStateOperator \
                                || (output_type) == oSteadyStateFile)

#define IS_C(output_type) ((output_type) == oCDynamicModel || (output_type) == oCStaticModel || (output_type) == oCDynamicSteadyStateOperator)

#define IS_LATEX(output_type) ((output_type) == oLatexStaticModel       \
                               || (output_type) == oLatexDynamicModel   \
                               || (output_type) == oLatexDynamicSteadyStateOperator)

/* Equal to 1 for Matlab langage, or to 0 for C language. Not defined for LaTeX.
   In Matlab, array indexes begin at 1, while they begin at 0 in C */
#define ARRAY_SUBSCRIPT_OFFSET(output_type) ((int) IS_MATLAB(output_type))

// Left and right array subscript delimiters: '(' and ')' for Matlab, '[' and ']' for C
#define LEFT_ARRAY_SUBSCRIPT(output_type) (IS_MATLAB(output_type) ? '(' : '[')
#define RIGHT_ARRAY_SUBSCRIPT(output_type) (IS_MATLAB(output_type) ? ')' : ']')

// Left and right parentheses
#define LEFT_PAR(output_type) (IS_LATEX(output_type) ? "\\left(" : "(")
#define RIGHT_PAR(output_type) (IS_LATEX(output_type) ? "\\right)" : ")")

// Computing cost above which a node can be declared a temporary term
#define MIN_COST_MATLAB (40*90)
#define MIN_COST_C (40*4)
#define MIN_COST(is_matlab) ((is_matlab) ? MIN_COST_MATLAB : MIN_COST_C)

//! Base class for expression nodes
class ExprNode
{
  friend class DataTree;
  friend class DynamicModel;
  friend class StaticModel;
  friend class ModelTree;
  friend class ExprNodeLess;
  friend class NumConstNode;
  friend class VariableNode;
  friend class UnaryOpNode;
  friend class BinaryOpNode;
  friend class TrinaryOpNode;
  friend class ExternalFunctionNode;
private:
  //! Computes derivative w.r. to a derivation ID (but doesn't store it in derivatives map)
  /*! You shoud use getDerivative() to get the benefit of symbolic a priori and of caching */
  virtual expr_t computeDerivative(int deriv_id) = 0;

protected:
  //! Reference to the enclosing DataTree
  DataTree &datatree;

  //! Index number
  int idx;

  //! Is the data member non_null_derivatives initialized ?
  bool preparedForDerivation;

  //! Set of derivation IDs with respect to which the derivative is potentially non-null
  set<int> non_null_derivatives;

  //! Used for caching of first order derivatives (when non-null)
  map<int, expr_t> derivatives;

  //! Cost of computing current node
  /*! Nodes included in temporary_terms are considered having a null cost */
  virtual int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const;

public:
  ExprNode(DataTree &datatree_arg);
  virtual ~ExprNode();

  //! Initializes data member non_null_derivatives
  virtual void prepareForDerivation() = 0;

  //! Returns derivative w.r. to derivation ID
  /*! Uses a symbolic a priori to pre-detect null derivatives, and caches the result for other derivatives (to avoid computing it several times)
    For an equal node, returns the derivative of lhs minus rhs */
  expr_t getDerivative(int deriv_id);

  //! Computes derivatives by applying the chain rule for some variables
  /*!
    \param deriv_id The derivation ID with respect to which we are derivating
    \param recursive_variables Contains the derivation ID for which chain rules must be applied. Keys are derivation IDs, values are equations of the form x=f(y) where x is the key variable and x doesn't appear in y
  */
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables) = 0;

  //! Returns precedence of node
  /*! Equals 100 for constants, variables, unary ops, and temporary terms */
  virtual int precedence(ExprNodeOutputType output_t, const temporary_terms_t &temporary_terms) const;

  //! Fills temporary_terms set, using reference counts
  /*! A node will be marked as a temporary term if it is referenced at least two times (i.e. has at least two parents), and has a computing cost (multiplied by reference count) greater to datatree.min_cost */
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms, and specifiying the set of already written external functions
  /*!
    \param[in] output the output stream
    \param[in] output_type the type of output (MATLAB, C, LaTeX...)
    \param[in] temporary_terms the nodes that are marked as temporary terms
    \param[in,out] tef_terms the set of already written external function nodes
  */
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const = 0;

  //! Writes output of node (with no temporary terms and with "outside model" output type)
  void writeOutput(ostream &output) const;

  //! Writes output of node (with no temporary terms)
  void writeOutput(ostream &output, ExprNodeOutputType output_type) const;

  //! Writes output of node, using a Txxx notation for nodes in temporary_terms
  void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;

  //! Writes the output for an external function, ensuring that the external function is called as few times as possible using temporary terms
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;

  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;

  //! Computes the set of all variables of a given symbol type in the expression
  /*!
    Variables are stored as integer pairs of the form (symb_id, lag).
    They are added to the set given in argument.
    Note that model local variables are substituted by their expression in the computation
    (and added if type_arg = ModelLocalVariable).
  */
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const = 0;

  //! Computes the set of endogenous variables in the expression
  /*!
    Endogenous are stored as integer pairs of the form (type_specific_id, lag).
    They are added to the set given in argument.
    Note that model local variables are substituted by their expression in the computation.
  */
  virtual void collectEndogenous(set<pair<int, int> > &result) const;

  //! Computes the set of exogenous variables in the expression
  /*!
    Exogenous are stored as integer pairs of the form (type_specific_id, lag).
    They are added to the set given in argument.
    Note that model local variables are substituted by their expression in the computation.
  */
  virtual void collectExogenous(set<pair<int, int> > &result) const;

  //! Computes the set of model local variables in the expression
  /*!
    Symbol IDs of these model local variables are added to the set given in argument.
    Note that this method is called recursively on the expressions associated to the model local variables detected.
  */
  virtual void collectModelLocalVariables(set<int> &result) const;

  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const = 0;

  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;

  class EvalException

  {
  };

  class EvalExternalFunctionException : public EvalException

  {
  };

  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException) = 0;
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const = 0;
  void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic) const;
  //! Creates a static version of this node
  /*!
    This method duplicates the current node by creating a similar node from which all leads/lags have been stripped,
    adds the result in the static_datatree argument (and not in the original datatree), and returns it.
  */
  virtual expr_t toStatic(DataTree &static_datatree) const = 0;
  //! Try to normalize an equation linear in its endogenous variable
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > > &List_of_Op_RHS) const = 0;

  //! Returns the maximum lead of endogenous in this expression
  /*! Always returns a non-negative value */
  virtual int maxEndoLead() const = 0;

  //! Returns the maximum lead of exogenous in this expression
  /*! Always returns a non-negative value */
  virtual int maxExoLead() const = 0;

  //! Returns the maximum lag of endogenous in this expression
  /*! Always returns a non-negative value */
  virtual int maxEndoLag() const = 0;

  //! Returns the maximum lag of exogenous in this expression
  /*! Always returns a non-negative value */
  virtual int maxExoLag() const = 0;

  //! Returns the relative period of the most forward term in this expression
  /*! A negative value means that the expression contains only lagged variables */
  virtual int maxLead() const = 0;
  virtual void print_deflator();

  //! Returns a new expression where all the leads/lags have been shifted backwards by the same amount
  /*!
    Only acts on endogenous, exogenous, exogenous det
    \param[in] n The number of lags by which to shift
    \return The same expression except that leads/lags have been shifted backwards
  */
  virtual expr_t decreaseLeadsLags(int n) const = 0;

  //! Type for the substitution map used in the process of creating auxiliary vars for leads >= 2
  typedef map<const ExprNode *, const VariableNode *> subst_table_t;

  //! Creates auxiliary endo lead variables corresponding to this expression
  /*!
    If maximum endogenous lead >= 3, this method will also create intermediary auxiliary var, and will add the equations of the form aux1 = aux2(+1) to the substitution table.
    \pre This expression is assumed to have maximum endogenous lead >= 2
    \param[in,out] subst_table The table to which new auxiliary variables and their correspondance will be added
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
    \return The new variable node corresponding to the current expression
  */
  VariableNode *createEndoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;

  //! Creates auxiliary exo lead variables corresponding to this expression
  /*!
    If maximum exogenous lead >= 2, this method will also create intermediary auxiliary var, and will add the equations of the form aux1 = aux2(+1) to the substitution table.
    \pre This expression is assumed to have maximum exogenous lead >= 1
    \param[in,out] subst_table The table to which new auxiliary variables and their correspondance will be added
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
    \return The new variable node corresponding to the current expression
  */
  VariableNode *createExoLeadAuxiliaryVarForMyself(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;

  //! Constructs a new expression where sub-expressions with max endo lead >= 2 have been replaced by auxiliary variables
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.

    If the method detects a sub-expr which needs to be substituted, two cases are possible:
    - if this expr is in the table, then it will use the corresponding variable and return the substituted expression
    - if this expr is not in the table, then it will create an auxiliary endogenous variable, add the substitution in the table and return the substituted expression

    \return A new equivalent expression where sub-expressions with max endo lead >= 2 have been replaced by auxiliary variables
  */
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const = 0;

  //! Constructs a new expression where endo variables with max endo lag >= 2 have been replaced by auxiliary variables
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
  */
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

  //! Constructs a new expression where exogenous variables with a lead have been replaced by auxiliary variables
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
  */
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const = 0;
  //! Constructs a new expression where exogenous variables with a lag have been replaced by auxiliary variables
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
  */
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const = 0;

  //! Constructs a new expression where the expectation operator has been replaced by auxiliary variables
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
    \param[in] partial_information_model Are we substituting in a partial information model?
  */
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const = 0;

  //! Constructs a new expression where log expression are replaced by an auxiliary variable and pow argument by an exponential expression
  /*!
    \param[in,out] subst_table Map used to store expressions that have already be substituted and their corresponding variable, in order to avoid creating two auxiliary variables for the same sub-expr.
    \param[out] neweqs Equations to be added to the model to match the creation of auxiliary variables.
  */
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const = 0;

  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const = 0;

  //! Return true if the nodeID is a numerical constant equal to value and false otherwise
  /*!
    \param[in] value of the numerical constante
    \param[out] the boolean equal to true if NodeId is a constant equal to value
  */
  virtual bool isNumConstNodeEqualTo(double value) const = 0;

  //! Returns true if the expression contains one or several endogenous variable
  virtual bool containsEndogenous(void) const = 0;

  //! Return true if the nodeID is a variable withe a type equal to type_arg, a specific variable id aqual to varfiable_id and a lag equal to lag_arg and false otherwise
  /*!
    \param[in] the type (type_arg), specifique variable id (variable_id and the lag (lag_arg)
    \param[out] the boolean equal to true if NodeId is the variable
  */
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const = 0;

  //! Replaces the Trend var with datatree.One
  virtual expr_t replaceTrendVar() const = 0;

  //! Constructs a new expression where the variable indicated by symb_id has been detrended
  /*!
    \param[in] symb_id indicating the variable to be detrended
    \param[in] expr_t indicating the trend
    \param[out] expr_t the new binary op pointing to a detrended variable
  */
  virtual expr_t detrend(int symb_id, expr_t trend) const = 0;

  //! Add ExprNodes to the provided datatree
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const = 0;

  //! Move a trend variable with lag/lead to time t by dividing/multiplying by its growth factor
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const = 0;
};

//! Object used to compare two nodes (using their indexes)
struct ExprNodeLess
{
  bool
  operator()(expr_t arg1, expr_t arg2) const
  {
    return arg1->idx < arg2->idx;
  }
};

//! Numerical constant node
/*! The constant is necessarily non-negative (this is enforced at the NumericalConstants class level) */
class NumConstNode : public ExprNode
{
private:
  //! Id from numerical constants table
  const int id;
  virtual expr_t computeDerivative(int deriv_id);
public:
  NumConstNode(DataTree &datatree_arg, int id_arg);
  int
  get_id() const
  {
    return id;
  };
  virtual void prepareForDerivation();
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
  virtual void print_deflator();
};

//! Symbol or variable node
class VariableNode : public ExprNode
{
  friend class UnaryOpNode;
private:
  //! Id from the symbol table
  const int symb_id;
  const SymbolType type;
  //! A positive value is a lead, a negative is a lag
  const int lag;
  virtual expr_t computeDerivative(int deriv_id);
public:
  VariableNode(DataTree &datatree_arg, int symb_id_arg, int lag_arg);
  virtual void prepareForDerivation();
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual void print_deflator();
  SymbolType
  get_type() const
  {
    return type;
  };
  int
  get_symb_id() const
  {
    return symb_id;
  };
  int get_lag() const { return lag; };
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
};

//! Unary operator node
class UnaryOpNode : public ExprNode
{
private:
  const expr_t arg;
  //! Stores the information set. Only used for expectation operator
  const int expectation_information_set;
  //! Only used for oSteadyStateParamDeriv and oSteadyStateParam2ndDeriv
  const int param1_symb_id, param2_symb_id;
  const UnaryOpcode op_code;
  virtual expr_t computeDerivative(int deriv_id);
  virtual int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const;
  //! Returns the derivative of this node if darg is the derivative of the argument
  expr_t composeDerivatives(expr_t darg, int deriv_id);
public:
  UnaryOpNode(DataTree &datatree_arg, UnaryOpcode op_code_arg, const expr_t arg_arg, int expectation_information_set_arg, int param1_symb_id_arg, int param2_symb_id_arg);
  virtual void prepareForDerivation();
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  static double eval_opcode(UnaryOpcode op_code, double v) throw (EvalException, EvalExternalFunctionException);
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  virtual void print_deflator();
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  //! Returns operand
  expr_t
  get_arg() const
  {
    return (arg);
  };
  //! Returns op code
  UnaryOpcode
  get_op_code() const
  {
    return (op_code);
  };
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  //! Creates another UnaryOpNode with the same opcode, but with a possibly different datatree and argument
  expr_t buildSimilarUnaryOpNode(expr_t alt_arg, DataTree &alt_datatree) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
};

//! Binary operator node
class BinaryOpNode : public ExprNode
{
private:
  const expr_t arg1, arg2;
  const BinaryOpcode op_code;
  virtual expr_t computeDerivative(int deriv_id);
  virtual int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const;
  //! Returns the derivative of this node if darg1 and darg2 are the derivatives of the arguments
  expr_t composeDerivatives(expr_t darg1, expr_t darg2);
  const int powerDerivOrder;
public:
  BinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
               BinaryOpcode op_code_arg, const expr_t arg2_arg);
  BinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
               BinaryOpcode op_code_arg, const expr_t arg2_arg, int powerDerivOrder);
  virtual void prepareForDerivation();
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  static double eval_opcode(double v1, BinaryOpcode op_code, double v2, int derivOrder) throw (EvalException, EvalExternalFunctionException);
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  virtual expr_t Compute_RHS(expr_t arg1, expr_t arg2, int op, int op_type) const;
  virtual void print_deflator();
  //! Returns first operand
  expr_t
  get_arg1() const
  {
    return (arg1);
  };
  //! Returns second operand
  expr_t
  get_arg2() const
  {
    return (arg2);
  };
  //! Returns op code
  BinaryOpcode
  get_op_code() const
  {
    return (op_code);
  };
  int
  get_power_deriv_order() const
  {
    return powerDerivOrder;
  }
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  //! Creates another BinaryOpNode with the same opcode, but with a possibly different datatree and arguments
  expr_t buildSimilarBinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, DataTree &alt_datatree) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
  //! Function to write out the oPowerNode in expr_t terms as opposed to writing out the function itself
  expr_t unpackPowerDeriv() const;
  //! Returns MULT_i*(lhs-rhs) = 0, creating multiplier MULT_i
  expr_t addMultipliersToConstraints(int i);
  //! Returns the non-zero hand-side of an equation (that must have a hand side equal to zero)
  expr_t getNonZeroPartofEquation() const;
};

//! Trinary operator node
class TrinaryOpNode : public ExprNode
{
  friend class ModelTree;
private:
  const expr_t arg1, arg2, arg3;
  const TrinaryOpcode op_code;
  virtual expr_t computeDerivative(int deriv_id);
  virtual int cost(const temporary_terms_t &temporary_terms, bool is_matlab) const;
  //! Returns the derivative of this node if darg1, darg2 and darg3 are the derivatives of the arguments
  expr_t composeDerivatives(expr_t darg1, expr_t darg2, expr_t darg3);
public:
  TrinaryOpNode(DataTree &datatree_arg, const expr_t arg1_arg,
                TrinaryOpcode op_code_arg, const expr_t arg2_arg, const expr_t arg3_arg);
  virtual void prepareForDerivation();
  virtual int precedence(ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  static double eval_opcode(double v1, TrinaryOpcode op_code, double v2, double v3) throw (EvalException, EvalExternalFunctionException);
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual void print_deflator();
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  //! Creates another TrinaryOpNode with the same opcode, but with a possibly different datatree and arguments
  expr_t buildSimilarTrinaryOpNode(expr_t alt_arg1, expr_t alt_arg2, expr_t alt_arg3, DataTree &alt_datatree) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
};

//! External function node
class ExternalFunctionNode : public ExprNode
{
private:
  virtual expr_t computeDerivative(int deriv_id);
  virtual expr_t composeDerivatives(const vector<expr_t> &dargs);
protected:
  //! Thrown when trying to access an unknown entry in external_function_node_map
  class UnknownFunctionNameAndArgs
  {
  };
  const int symb_id;
  const vector<expr_t> arguments;
  //! Returns true if the given external function has been written as a temporary term
  bool alreadyWrittenAsTefTerm(int the_symb_id, deriv_node_temp_terms_t &tef_terms) const;
  //! Returns the index in the tef_terms map of this external function
  int getIndxInTefTerms(int the_symb_id, deriv_node_temp_terms_t &tef_terms) const throw (UnknownFunctionNameAndArgs);
  //! Helper function to write output arguments of any given external function
  void writeExternalFunctionArguments(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
public:
  ExternalFunctionNode(DataTree &datatree_arg, int symb_id_arg,
                       const vector<expr_t> &arguments_arg);
  virtual void prepareForDerivation();
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void collectVariables(SymbolType type_arg, set<pair<int, int> > &result) const;
  virtual void collectTemporary_terms(const temporary_terms_t &temporary_terms, temporary_terms_inuse_t &temporary_terms_inuse, int Curr_Block) const;
  virtual double eval(const eval_context_t &eval_context) const throw (EvalException, EvalExternalFunctionException);
  unsigned int compileExternalFunctionArguments(ostream &CompileCode, unsigned int &instruction_number,
                                                bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                                const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                                deriv_node_temp_terms_t &tef_terms) const;
  virtual void print_deflator();

  virtual void compile(ostream &CompileCode, unsigned int &instruction_number, bool lhs_rhs, const temporary_terms_t &temporary_terms, const map_idx_t &map_idx, bool dynamic, bool steady_dynamic, deriv_node_temp_terms_t &tef_terms) const;
  virtual expr_t toStatic(DataTree &static_datatree) const;
  virtual pair<int, expr_t> normalizeEquation(int symb_id_endo, vector<pair<int, pair<expr_t, expr_t> > >  &List_of_Op_RHS) const;
  virtual expr_t getChainRuleDerivative(int deriv_id, const map<int, expr_t> &recursive_variables);
  virtual int maxEndoLead() const;
  virtual int maxExoLead() const;
  virtual int maxEndoLag() const;
  virtual int maxExoLag() const;
  virtual int maxLead() const;
  virtual expr_t decreaseLeadsLags(int n) const;
  virtual expr_t substituteEndoLeadGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteEndoLagGreaterThanTwo(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExoLead(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool deterministic_model) const;
  virtual expr_t substituteExoLag(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs) const;
  virtual expr_t substituteExpectation(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs, bool partial_information_model) const;
  virtual expr_t substituteLogPow(subst_table_t &subst_table, vector<BinaryOpNode *> &neweqs1, vector<BinaryOpNode *> &neweqs2) const;
  virtual expr_t buildSimilarExternalFunctionNode(vector<expr_t> &alt_args, DataTree &alt_datatree) const;
  virtual expr_t decreaseLeadsLagsPredeterminedVariables() const;
  virtual bool isNumConstNodeEqualTo(double value) const;
  virtual bool containsEndogenous(void) const;
  virtual bool isVariableNodeEqualTo(SymbolType type_arg, int variable_id, int lag_arg) const;
  virtual void writePrhs(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms, const string &ending) const;
  virtual expr_t replaceTrendVar() const;
  virtual expr_t detrend(int symb_id, expr_t trend) const;
  virtual expr_t cloneDynamic(DataTree &dynamic_datatree) const;
  virtual expr_t removeTrendLeadLag(map<int, expr_t> trend_symbols_map) const;
};

class FirstDerivExternalFunctionNode : public ExternalFunctionNode
{
private:
  const int inputIndex;
  virtual expr_t composeDerivatives(const vector<expr_t> &dargs);
public:
  FirstDerivExternalFunctionNode(DataTree &datatree_arg,
                                 int top_level_symb_id_arg,
                                 const vector<expr_t> &arguments_arg,
                                 int inputIndex_arg);
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void compile(ostream &CompileCode, unsigned int &instruction_number,
                       bool lhs_rhs, const temporary_terms_t &temporary_terms,
                       const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                       deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
  virtual void compileExternalFunctionOutput(ostream &CompileCode, unsigned int &instruction_number,
                                             bool lhs_rhs, const temporary_terms_t &temporary_terms,
                                             const map_idx_t &map_idx, bool dynamic, bool steady_dynamic,
                                             deriv_node_temp_terms_t &tef_terms) const;
  virtual void print_deflator();
};

class SecondDerivExternalFunctionNode : public ExternalFunctionNode
{
private:
  const int inputIndex1;
  const int inputIndex2;
  virtual expr_t computeDerivative(int deriv_id);
public:
  SecondDerivExternalFunctionNode(DataTree &datatree_arg,
                                  int top_level_symb_id_arg,
                                  const vector<expr_t> &arguments_arg,
                                  int inputIndex1_arg,
                                  int inputIndex2_arg);
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count, temporary_terms_t &temporary_terms, bool is_matlab) const;
  virtual void computeTemporaryTerms(map<expr_t, int> &reference_count,
                                     temporary_terms_t &temporary_terms,
                                     map<expr_t, pair<int, int> > &first_occurence,
                                     int Curr_block,
                                     vector< vector<temporary_terms_t> > &v_temporary_terms,
                                     int equation) const;
  virtual void writeOutput(ostream &output, ExprNodeOutputType output_type, const temporary_terms_t &temporary_terms, deriv_node_temp_terms_t &tef_terms) const;
  virtual void writeExternalFunctionOutput(ostream &output, ExprNodeOutputType output_type,
                                           const temporary_terms_t &temporary_terms,
                                           deriv_node_temp_terms_t &tef_terms) const;
 virtual void print_deflator();
};

#endif
