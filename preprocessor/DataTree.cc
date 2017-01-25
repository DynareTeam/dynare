/*
 * Copyright (C) 2003-2016 Dynare Team
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

#include "DataTree.hh"

DataTree::DataTree(SymbolTable &symbol_table_arg,
                   NumericalConstants &num_constants_arg,
                   ExternalFunctionsTable &external_functions_table_arg) :
  symbol_table(symbol_table_arg),
  num_constants(num_constants_arg),
  external_functions_table(external_functions_table_arg),
  node_counter(0)
{
  Zero = AddNonNegativeConstant("0");
  One = AddNonNegativeConstant("1");
  Two = AddNonNegativeConstant("2");

  MinusOne = AddUMinus(One);

  NaN = AddNonNegativeConstant("NaN");
  Infinity = AddNonNegativeConstant("Inf");
  MinusInfinity = AddUMinus(Infinity);

  Pi = AddNonNegativeConstant("3.141592653589793");
}

DataTree::~DataTree()
{
  for (node_list_t::iterator it = node_list.begin(); it != node_list.end(); it++)
    delete *it;
}

expr_t
DataTree::AddNonNegativeConstant(const string &value)
{
  int id = num_constants.AddNonNegativeConstant(value);

  num_const_node_map_t::iterator it = num_const_node_map.find(id);
  if (it != num_const_node_map.end())
    return it->second;
  else
    return new NumConstNode(*this, id);
}

VariableNode *
DataTree::AddVariableInternal(int symb_id, int lag)
{
  variable_node_map_t::iterator it = variable_node_map.find(make_pair(symb_id, lag));
  if (it != variable_node_map.end())
    return it->second;
  else
    return new VariableNode(*this, symb_id, lag);
}

bool
DataTree::ParamUsedWithLeadLagInternal() const
{
  for (variable_node_map_t::const_iterator it = variable_node_map.begin();
       it != variable_node_map.end(); it++)
    if (symbol_table.getType(it->first.first) == eParameter && it->first.second != 0)
      return true;
  return false;
}

VariableNode *
DataTree::AddVariable(int symb_id, int lag)
{
  assert(lag == 0);
  return AddVariableInternal(symb_id, lag);
}

expr_t
DataTree::AddPlus(expr_t iArg1, expr_t iArg2)
{
  if (iArg1 != Zero && iArg2 != Zero)
    {
      // Simplify x+(-y) in x-y
      UnaryOpNode *uarg2 = dynamic_cast<UnaryOpNode *>(iArg2);
      if (uarg2 != NULL && uarg2->get_op_code() == oUminus)
        return AddMinus(iArg1, uarg2->get_arg());

      // To treat commutativity of "+"
      // Nodes iArg1 and iArg2 are sorted by index
      if (iArg1->idx > iArg2->idx)
        {
          expr_t tmp = iArg1;
          iArg1 = iArg2;
          iArg2 = tmp;
        }
      return AddBinaryOp(iArg1, oPlus, iArg2);
    }
  else if (iArg1 != Zero)
    return iArg1;
  else if (iArg2 != Zero)
    return iArg2;
  else
    return Zero;
}

expr_t
DataTree::AddMinus(expr_t iArg1, expr_t iArg2)
{
  if (iArg2 == Zero)
    return iArg1;

  if (iArg1 == Zero)
    return AddUMinus(iArg2);

  if (iArg1 == iArg2)
    return Zero;

  return AddBinaryOp(iArg1, oMinus, iArg2);
}

expr_t
DataTree::AddUMinus(expr_t iArg1)
{
  if (iArg1 != Zero)
    {
      // Simplify -(-x) in x
      UnaryOpNode *uarg = dynamic_cast<UnaryOpNode *>(iArg1);
      if (uarg != NULL && uarg->get_op_code() == oUminus)
        return uarg->get_arg();

      return AddUnaryOp(oUminus, iArg1);
    }
  else
    return Zero;
}

expr_t
DataTree::AddTimes(expr_t iArg1, expr_t iArg2)
{
  if (iArg1 == MinusOne)
    return AddUMinus(iArg2);
  else if (iArg2 == MinusOne)
    return AddUMinus(iArg1);
  else if (iArg1 != Zero && iArg1 != One && iArg2 != Zero && iArg2 != One)
    {
      // To treat commutativity of "*"
      // Nodes iArg1 and iArg2 are sorted by index
      if (iArg1->idx > iArg2->idx)
        {
          expr_t tmp = iArg1;
          iArg1 = iArg2;
          iArg2 = tmp;
        }
      return AddBinaryOp(iArg1, oTimes, iArg2);
    }
  else if (iArg1 != Zero && iArg1 != One && iArg2 == One)
    return iArg1;
  else if (iArg2 != Zero && iArg2 != One && iArg1 == One)
    return iArg2;
  else if (iArg2 == One && iArg1 == One)
    return One;
  else
    return Zero;
}

expr_t
DataTree::AddDivide(expr_t iArg1, expr_t iArg2) throw (DivisionByZeroException)
{
  if (iArg2 == One)
    return iArg1;

  // This test should be before the next two, otherwise 0/0 won't be rejected
  if (iArg2 == Zero)
    {
      cerr << "ERROR: Division by zero!" << endl;
      throw DivisionByZeroException();
    }

  if (iArg1 == Zero)
    return Zero;

  if (iArg1 == iArg2)
    return One;

  return AddBinaryOp(iArg1, oDivide, iArg2);
}

expr_t
DataTree::AddLess(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oLess, iArg2);
}

expr_t
DataTree::AddGreater(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oGreater, iArg2);
}

expr_t
DataTree::AddLessEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oLessEqual, iArg2);
}

expr_t
DataTree::AddGreaterEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oGreaterEqual, iArg2);
}

expr_t
DataTree::AddEqualEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oEqualEqual, iArg2);
}

expr_t
DataTree::AddDifferent(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oDifferent, iArg2);
}

expr_t
DataTree::AddPower(expr_t iArg1, expr_t iArg2)
{
  if (iArg1 != Zero && iArg2 != Zero && iArg1 != One && iArg2 != One)
    return AddBinaryOp(iArg1, oPower, iArg2);
  else if (iArg1 == One)
    return One;
  else if (iArg2 == One)
    return iArg1;
  else if (iArg2 == Zero)
    return One;
  else
    return Zero;
}

expr_t
DataTree::AddPowerDeriv(expr_t iArg1, expr_t iArg2, int powerDerivOrder)
{
  assert(powerDerivOrder > 0);
  return AddBinaryOp(iArg1, oPowerDeriv, iArg2, powerDerivOrder);
}

expr_t
DataTree::AddExp(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oExp, iArg1);
  else
    return One;
}

expr_t
DataTree::AddLog(expr_t iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    return AddUnaryOp(oLog, iArg1);
  else if (iArg1 == One)
    return Zero;
  else
    {
      cerr << "ERROR: log(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
DataTree::AddLog10(expr_t iArg1)
{
  if (iArg1 != Zero && iArg1 != One)
    return AddUnaryOp(oLog10, iArg1);
  else if (iArg1 == One)
    return Zero;
  else
    {
      cerr << "ERROR: log10(0) not defined!" << endl;
      exit(EXIT_FAILURE);
    }
}

expr_t
DataTree::AddCos(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oCos, iArg1);
  else
    return One;
}

expr_t
DataTree::AddSin(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSin, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddTan(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oTan, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAcos(expr_t iArg1)
{
  if (iArg1 != One)
    return AddUnaryOp(oAcos, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAsin(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAsin, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAtan(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAtan, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddCosh(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oCosh, iArg1);
  else
    return One;
}

expr_t
DataTree::AddSinh(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSinh, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddTanh(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oTanh, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAcosh(expr_t iArg1)
{
  if (iArg1 != One)
    return AddUnaryOp(oAcosh, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAsinh(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAsinh, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAtanh(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oAtanh, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddSqrt(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oSqrt, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddAbs(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;
  if (iArg1 == One)
    return One;
  else
    return AddUnaryOp(oAbs, iArg1);
}

expr_t
DataTree::AddSign(expr_t iArg1)
{
  if (iArg1 == Zero)
    return Zero;
  if (iArg1 == One)
    return One;
  else
    return AddUnaryOp(oSign, iArg1);
}

expr_t
DataTree::AddErf(expr_t iArg1)
{
  if (iArg1 != Zero)
    return AddUnaryOp(oErf, iArg1);
  else
    return Zero;
}

expr_t
DataTree::AddMax(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oMax, iArg2);
}

expr_t
DataTree::AddMin(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oMin, iArg2);
}

expr_t
DataTree::AddNormcdf(expr_t iArg1, expr_t iArg2, expr_t iArg3)
{
  return AddTrinaryOp(iArg1, oNormcdf, iArg2, iArg3);
}

expr_t
DataTree::AddNormpdf(expr_t iArg1, expr_t iArg2, expr_t iArg3)
{
  return AddTrinaryOp(iArg1, oNormpdf, iArg2, iArg3);
}

expr_t
DataTree::AddSteadyState(expr_t iArg1)
{
  return AddUnaryOp(oSteadyState, iArg1);
}

expr_t
DataTree::AddSteadyStateParamDeriv(expr_t iArg1, int param_symb_id)
{
  return AddUnaryOp(oSteadyStateParamDeriv, iArg1, 0, param_symb_id);
}

expr_t
DataTree::AddSteadyStateParam2ndDeriv(expr_t iArg1, int param1_symb_id, int param2_symb_id)
{
  return AddUnaryOp(oSteadyStateParam2ndDeriv, iArg1, 0, param1_symb_id, param2_symb_id);
}

expr_t
DataTree::AddExpectation(int iArg1, expr_t iArg2)
{
  return AddUnaryOp(oExpectation, iArg2, iArg1);
}

expr_t
DataTree::AddEqual(expr_t iArg1, expr_t iArg2)
{
  return AddBinaryOp(iArg1, oEqual, iArg2);
}

void
DataTree::AddLocalVariable(int symb_id, expr_t value) throw (LocalVariableException)
{
  assert(symbol_table.getType(symb_id) == eModelLocalVariable);

  // Throw an exception if symbol already declared
  map<int, expr_t>::iterator it = local_variables_table.find(symb_id);
  if (it != local_variables_table.end())
    throw LocalVariableException(symbol_table.getName(symb_id));

  local_variables_table[symb_id] = value;
}

expr_t
DataTree::AddExternalFunction(int symb_id, const vector<expr_t> &arguments)
{
  assert(symbol_table.getType(symb_id) == eExternalFunction);

  external_function_node_map_t::iterator it = external_function_node_map.find(make_pair(arguments, symb_id));
  if (it != external_function_node_map.end())
    return it->second;

  return new ExternalFunctionNode(*this, symb_id, arguments);
}

expr_t
DataTree::AddFirstDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index)
{
  assert(symbol_table.getType(top_level_symb_id) == eExternalFunction);

  first_deriv_external_function_node_map_t::iterator it
    = first_deriv_external_function_node_map.find(make_pair(make_pair(arguments, input_index),
                                                            top_level_symb_id));
  if (it != first_deriv_external_function_node_map.end())
    return it->second;

  return new FirstDerivExternalFunctionNode(*this, top_level_symb_id, arguments, input_index);
}

expr_t
DataTree::AddSecondDerivExternalFunction(int top_level_symb_id, const vector<expr_t> &arguments, int input_index1, int input_index2)
{
  assert(symbol_table.getType(top_level_symb_id) == eExternalFunction);

  second_deriv_external_function_node_map_t::iterator it
    = second_deriv_external_function_node_map.find(make_pair(make_pair(arguments,
                                                                       make_pair(input_index1, input_index2)),
                                                             top_level_symb_id));
  if (it != second_deriv_external_function_node_map.end())
    return it->second;

  return new SecondDerivExternalFunctionNode(*this, top_level_symb_id, arguments, input_index1, input_index2);
}

bool
DataTree::isSymbolUsed(int symb_id) const
{
  for (variable_node_map_t::const_iterator it = variable_node_map.begin();
       it != variable_node_map.end(); it++)
    if (it->first.first == symb_id)
      return true;

  if (local_variables_table.find(symb_id) != local_variables_table.end())
    return true;

  return false;
}

int
DataTree::getDerivID(int symb_id, int lag) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

SymbolType
DataTree::getTypeByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

int
DataTree::getLagByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

int
DataTree::getSymbIDByDerivID(int deriv_id) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

void
DataTree::addAllParamDerivId(set<int> &deriv_id_set)
{
}

int
DataTree::getDynJacobianCol(int deriv_id) const throw (UnknownDerivIDException)
{
  throw UnknownDerivIDException();
}

bool
DataTree::isUnaryOpUsed(UnaryOpcode opcode) const
{
  for (unary_op_node_map_t::const_iterator it = unary_op_node_map.begin();
       it != unary_op_node_map.end(); it++)
    if (it->first.first.second == opcode)
      return true;

  return false;
}

bool
DataTree::isBinaryOpUsed(BinaryOpcode opcode) const
{
  for (binary_op_node_map_t::const_iterator it = binary_op_node_map.begin();
       it != binary_op_node_map.end(); it++)
    if (it->first.second == opcode)
      return true;

  return false;
}

bool
DataTree::isTrinaryOpUsed(TrinaryOpcode opcode) const
{
  for (trinary_op_node_map_t::const_iterator it = trinary_op_node_map.begin();
       it != trinary_op_node_map.end(); it++)
    if (it->first.second == opcode)
      return true;

  return false;
}

bool
DataTree::isExternalFunctionUsed(int symb_id) const
{
  for (external_function_node_map_t::const_iterator it = external_function_node_map.begin();
       it != external_function_node_map.end(); it++)
    if (it->first.second == symb_id)
      return true;

  return false;
}

bool
DataTree::isFirstDerivExternalFunctionUsed(int symb_id) const
{
  for (first_deriv_external_function_node_map_t::const_iterator it = first_deriv_external_function_node_map.begin();
       it != first_deriv_external_function_node_map.end(); it++)
    if (it->first.second == symb_id)
      return true;

  return false;
}

bool
DataTree::isSecondDerivExternalFunctionUsed(int symb_id) const
{
  for (second_deriv_external_function_node_map_t::const_iterator it = second_deriv_external_function_node_map.begin();
       it != second_deriv_external_function_node_map.end(); it++)
    if (it->first.second == symb_id)
      return true;

  return false;
}

int
DataTree::minLagForSymbol(int symb_id) const
{
  int r = 0;
  for (variable_node_map_t::const_iterator it = variable_node_map.begin();
       it != variable_node_map.end(); ++it)
    if (it->first.first == symb_id && it->first.second < r)
      r = it->first.second;
  return r;
}

void
DataTree::writePowerDerivCHeader(ostream &output) const
{
  if (isBinaryOpUsed(oPowerDeriv))
    output << "double getPowerDeriv(double, double, int);" << endl;
}

void
DataTree::writePowerDeriv(ostream &output) const
{
  if (isBinaryOpUsed(oPowerDeriv))
    output << "/*" << endl
           << " * The k-th derivative of x^p" << endl
           << " */" << endl
           << "double getPowerDeriv(double x, double p, int k)" << endl
           << "{" << endl
           << "#ifdef _MSC_VER" << endl
           << "# define nearbyint(x) (fabs((x)-floor(x)) < fabs((x)-ceil(x)) ? floor(x) : ceil(x))" << endl
           << "#endif" << endl
           << "  if ( fabs(x) < " << NEAR_ZERO << " && p > 0 && k > p && fabs(p-nearbyint(p)) < " << NEAR_ZERO << " )" << endl
           << "    return 0.0;" << endl
           << "  else" << endl
           << "    {" << endl
           << "      int i = 0;" << endl
           << "      double dxp = pow(x, p-k);" << endl
           << "      for (; i<k; i++)" << endl
           << "        dxp *= p--;" << endl
           << "      return dxp;" << endl
           << "    }" << endl
           << "}" << endl;
}

void
DataTree::writeNormcdfCHeader(ostream &output) const
{
#if defined(_WIN32) || defined(__CYGWIN32__) || defined(__MINGW32__)
  if (isTrinaryOpUsed(oNormcdf))
    output << "#ifdef _MSC_VER" << endl
           << "double normcdf(double);" << endl
           << "#endif" << endl;
#endif
}

void
DataTree::writeNormcdf(ostream &output) const
{
#if defined(_WIN32) || defined(__CYGWIN32__) || defined(__MINGW32__)
  if (isTrinaryOpUsed(oNormcdf))
    output << endl
           << "#ifdef _MSC_VER" << endl
           << "/*" << endl
           << " * Define normcdf for MSVC compiler" << endl
           << " */" << endl
           << "double normcdf(double x)" << endl
           << "{" << endl
           << "#if _MSC_VER >= 1700" << endl
           << "  return 0.5 * erfc(-x * M_SQRT1_2);" << endl
           << "#else" << endl
           << "  // From http://www.johndcook.com/blog/cpp_phi" << endl
           << "  double a1 =  0.254829592;" << endl
           << "  double a2 = -0.284496736;" << endl
           << "  double a3 =  1.421413741;" << endl
           << "  double a4 = -1.453152027;" << endl
           << "  double a5 =  1.061405429;" << endl
           << "  double p  =  0.3275911;" << endl
           << "  int sign = (x < 0) ? -1 : 1;" << endl
           << "  x = fabs(x)/sqrt(2.0);" << endl
           << "  // From the Handbook of Mathematical Functions by Abramowitz and Stegun, formula 7.1.26" << endl
           << "  double t = 1.0/(1.0 + p*x);" << endl
           << "  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);" << endl
           << "  return 0.5*(1.0 + sign*y);" << endl
           << "#endif" << endl
           << "}" << endl
           << "#endif" << endl;
#endif
}
