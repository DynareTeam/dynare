/*
 * Copyright (C) 2007-2011 Dynare Team
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

#ifndef _CODEINTERPRETER_HH
#define _CODEINTERPRETER_HH
//#define DEBUGL
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstring>
#include <vector>
#ifdef LINBCG
# include "linbcg.hh"
#endif
#ifdef BYTE_CODE
# ifndef DEBUG_EX
#  include "mex.h"
# else
#  include "mex_interface.hh"
# endif
#endif

#include <stdint.h>

#define NEAR_ZERO (1e-12)

using namespace std;

/**
 * \enum Tags
 * \brief The differents flags of the bytecode
 */
enum Tags
  {
    FLDZ,         //!< Stores zero in the stack - 0 (0)
    FLDC,         //!< Stores a constant term in the stack - 1 (1)

    FDIMT,        //!< Defines the number of temporary terms - dynamic context (the period has to be indicated) - 2 (2)
    FDIMST,       //!< Defines the number of temporary terms - static context (the period hasn't to be indicated) - 3  (3)
    FLDT,         //!< Stores a temporary term in the stack - dynamic context (the period has to be indicated) - 4 (4)
    FLDST,        //!< Stores a temporary term in the stack - static context (the period hasn't to be indicated) - 5 (5)
    FSTPT,        //!< Loads a temporary term from the stack - dynamic context (the period has to be indicated) - 6 (6)
    FSTPST,       //!< Loads a temporary term from the stack - static context (the period hasn't to be indicated) - 7 (7)

    FLDU,         //!< Stores an element of the vector U in the stack - dynamic context (the period has to be indicated) - 8 (8)
    FLDSU,        //!< Stores an element of the vector U in the stack - static context (the period hasn't to be indicated) - 9 (9)
    FSTPU,        //!< Loads an element of the vector U from the stack - dynamic context (the period has to be indicated) - A (10)
    FSTPSU,       //!< Loads an element of the vector U from the stack - static context (the period hasn't to be indicated) - B (11)

    FLDV,         //!< Stores a variable (described in SymbolType) in the stack - dynamic context (the period has to be indicated) - C (12)
    FLDSV,        //!< Stores a variable (described in SymbolType) in the stack - static context (the period hasn't to be indicated) - D (13)
    FLDVS,        //!< Stores a variable (described in SymbolType) in the stack - dynamic context but inside the STEADYSTATE function (the period hasn't to be indicated) - E (14)
    FSTPV,        //!< Loads a variable (described in SymbolType) from the stack - dynamic context (the period has to be indicated) - F (15)
    FSTPSV,       //!< Loads a variable (described in SymbolType) from the stack - static context (the period hasn't to be indicated) - 10 (16)

    FLDR,         //!< Stores a residual in the stack - 11 (17)
    FSTPR,        //!< Loads a residual from the stack - 12 (18)

    FSTPG,        //!< Loads a derivative from the stack - 13 (19)
    FSTPG2,       //!< Loads a derivative matrix for static model from the stack - 14 (20)
    FSTPG3,       //!< Loads a derivative matrix for a dynamic model from the stack - 15 (21)
    FSTPG4,       //!< Loads a second order derivative matrix for a dynamic model from the stack - 16 (22)

    FUNARY,       //!< A Unary operator - 17 (23)
    FBINARY,      //!< A binary operator - 18 (24)
    FTRINARY,     //!< A trinary operator - 19 (25)

    FCUML,        //!< Cumulates the result - 1A (26)

    FJMPIFEVAL,   //!< Jump if evaluate = true - 1B (27)
    FJMP,         //!< Jump - 1C (28)

    FBEGINBLOCK,  //!< Defines the begining of a model block - 1D (29)
    FENDBLOCK,    //!< Defines the end of a model block - 1E (30)
    FENDEQU,      //!< Defines the last equation of the block. For block that has to be solved, the derivatives appear just after this flag - 1F (31)
    FEND,         //!< Defines the end of the model code - 20 (32)

    FOK,          //!< Used for debugging purpose - 21 (33)

    FNUMEXPR,     //!< Store the expression type and references - 22 (34)

    FCALL,        //!< Call an external function - 23 (35)
    FPUSH,        //!< Push a double in the stack - 24 (36)
    FPOP,         //!< Pop a double from the stack - 25 (37)
    FLDTEF,       //!< Stores the result of an external function in the stack - 26 (38)
    FSTPTEF,      //!< Loads the result of an external function from the stack- 27 (39)
    FLDTEFD,      //!< Stores the result of an external function in the stack - 28 (40)
    FSTPTEFD,     //!< Loads the result of an external function from the stack- 29 (41)
    FLDTEFDD,     //!< Stores the result of an external function in the stack - 28 (42)
    FSTPTEFDD     //!< Loads the result of an external function from the stack- 29 (43)

  };

enum BlockType
  {
    SIMULTANS,  //!< Simultaneous time separable block
    PROLOGUE,   //!< Prologue block (one equation at the beginning, later merged)
    EPILOGUE,   //!< Epilogue block (one equation at the beginning, later merged)
    SIMULTAN    //!< Simultaneous time unseparable block
  };

enum EquationType
  {
    E_UNKNOWN,              //!< Unknown equation type
    E_EVALUATE,             //!< Simple evaluation, normalized variable on left-hand side
    E_EVALUATE_S,           //!< Simple evaluation, normalize using the first order derivative
    E_SOLVE                 //!< No simple evaluation of the equation, it has to be solved
  };

enum BlockSimulationType
  {
    UNKNOWN,                      //!< Unknown simulation type
    EVALUATE_FORWARD,             //!< Simple evaluation, normalized variable on left-hand side, forward
    EVALUATE_BACKWARD,            //!< Simple evaluation, normalized variable on left-hand side, backward
    SOLVE_FORWARD_SIMPLE,         //!< Block of one equation, newton solver needed, forward
    SOLVE_BACKWARD_SIMPLE,        //!< Block of one equation, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_SIMPLE,  //!< Block of one equation, newton solver needed, forward & ackward
    SOLVE_FORWARD_COMPLETE,       //!< Block of several equations, newton solver needed, forward
    SOLVE_BACKWARD_COMPLETE,      //!< Block of several equations, newton solver needed, backward
    SOLVE_TWO_BOUNDARIES_COMPLETE //!< Block of several equations, newton solver needed, forward and backwar
  };

//! Enumeration of possible symbol types
/*! Warning: do not to change existing values for 0 to 4: the values matter for homotopy_setup command */
enum SymbolType
  {
    eEndogenous = 0,                //!< Endogenous
    eExogenous = 1,                 //!< Exogenous
    eExogenousDet = 2,              //!< Exogenous deterministic
    eParameter = 4,                 //!< Parameter
    eModelLocalVariable = 10,       //!< Local variable whose scope is model (pound expression)
    eModFileLocalVariable = 11,     //!< Local variable whose scope is mod file (model excluded)
    eExternalFunction = 12,         //!< External (user-defined) function
    eTrend = 13,                    //!< Trend variable
    eStatementDeclaredVariable = 14 //!< Local variable assigned within a Statement (see subsample statement for example)
  };

enum ExpressionType
  {
    TemporaryTerm,
    ModelEquation,
    FirstEndoDerivative,
    FirstOtherEndoDerivative,
    FirstExoDerivative,
    FirstExodetDerivative,
    FirstParamDerivative,
    SecondEndoDerivative,
    SecondExoDerivative,
    SecondExodetDerivative,
    SecondParamDerivative,
    ThirdEndoDerivative,
    ThirdExoDerivative,
    ThirdExodetDerivative,
    ThirdParamDerivative
  };

enum UnaryOpcode
  {
    oUminus,
    oExp,
    oLog,
    oLog10,
    oCos,
    oSin,
    oTan,
    oAcos,
    oAsin,
    oAtan,
    oCosh,
    oSinh,
    oTanh,
    oAcosh,
    oAsinh,
    oAtanh,
    oSqrt,
    oAbs,
    oSign,
    oSteadyState,
    oSteadyStateParamDeriv, // for the derivative of the STEADY_STATE operator w.r.t. to a parameter
    oSteadyStateParam2ndDeriv, // for the 2nd derivative of the STEADY_STATE operator w.r.t. to a parameter
    oExpectation,
    oErf
  };

enum BinaryOpcode
  {
    oPlus,
    oMinus,
    oTimes,
    oDivide,
    oPower,
    oPowerDeriv, // for the derivative of the power function (see trac ticket #78)
    oEqual,
    oMax,
    oMin,
    oLess,
    oGreater,
    oLessEqual,
    oGreaterEqual,
    oEqualEqual,
    oDifferent
  };

enum TrinaryOpcode
  {
    oNormcdf,
    oNormpdf
  };

enum external_function_type
{
  ExternalFunctionWithoutDerivative,
  ExternalFunctionWithFirstDerivative,
  ExternalFunctionWithFirstandSecondDerivative,
  ExternalFunctionNumericalFirstDerivative,
  ExternalFunctionFirstDerivative,
  ExternalFunctionNumericalSecondDerivative,
  ExternalFunctionSecondDerivative
};

enum PriorDistributions
  {
    eNoShape = 0,
    eBeta = 1,
    eGamma = 2,
    eNormal = 3,
    eInvGamma = 4,
    eInvGamma1 = 4,
    eUniform = 5,
    eInvGamma2 = 6
  };

struct Block_contain_type
{
  int Equation, Variable, Own_Derivative;
};

#pragma pack(push, 1)
class TagWithoutArgument
{
protected:
  uint8_t op_code;
public:
  inline
  TagWithoutArgument(uint8_t op_code_arg) : op_code(op_code_arg)
  {
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this));
    instruction_number++;
  };
};

template < class T1 >
class TagWithOneArgument
{
protected:
  uint8_t op_code;
  T1 arg1;
public:
  inline
  TagWithOneArgument(uint8_t op_code_arg) : op_code(op_code_arg)
  {
  };
  inline
  TagWithOneArgument(uint8_t op_code_arg, T1 arg_arg1) : op_code(op_code_arg), arg1(arg_arg1)
  {
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(TagWithOneArgument));
    instruction_number++;
  };
};

template < class T1, class T2 >
class TagWithTwoArguments
{
protected:
  uint8_t op_code;
  T1 arg1;
  T2 arg2;
public:
  inline
  TagWithTwoArguments(uint8_t op_code_arg) : op_code(op_code_arg)
  {
  };
  inline
  TagWithTwoArguments(uint8_t op_code_arg, T1 arg_arg1, T2 arg_arg2) : op_code(op_code_arg), arg1(arg_arg1), arg2(arg_arg2)
  {
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this));
    instruction_number++;
  };
};

template < class T1, class T2, class T3 >
class TagWithThreeArguments
{
protected:
  uint8_t op_code;
  T1 arg1;
  T2 arg2;
  T3 arg3;
public:
  inline
  TagWithThreeArguments(uint8_t op_code_arg) : op_code(op_code_arg)
  {
  };
  inline
  TagWithThreeArguments(uint8_t op_code_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3) : op_code(op_code_arg), arg1(arg_arg1), arg2(arg_arg2), arg3(arg_arg3)
  {
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this));
    instruction_number++;
  };
};

template < class T1, class T2, class T3, class T4 >
class TagWithFourArguments
{
protected:
  uint8_t op_code;
  T1 arg1;
  T2 arg2;
  T3 arg3;
  T4 arg4;
public:
  inline
  TagWithFourArguments(uint8_t op_code_arg) : op_code(op_code_arg)
  {
  };
  inline
  TagWithFourArguments(uint8_t op_code_arg, T1 arg_arg1, T2 arg_arg2, T3 arg_arg3, T4 arg_arg4) : op_code(op_code_arg), arg1(arg_arg1), arg2(arg_arg2), arg3(arg_arg3), arg4(arg_arg4)
  {
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(*this));
    instruction_number++;
  };
};

class FLDZ_ : public TagWithoutArgument
{
public:
  inline
  FLDZ_() : TagWithoutArgument(FLDZ)
  {
  };
};

class FEND_ : public TagWithoutArgument
{
public:
  inline
  FEND_() : TagWithoutArgument(FEND)
  {
  };
};

class FENDBLOCK_ : public TagWithoutArgument
{
public:
  inline
  FENDBLOCK_() : TagWithoutArgument(FENDBLOCK)
  {
  };
};

class FENDEQU_ : public TagWithoutArgument
{
public:
  inline
  FENDEQU_() : TagWithoutArgument(FENDEQU)
  {
  };
};

class FCUML_ : public TagWithoutArgument
{
public:
  inline
  FCUML_() : TagWithoutArgument(FCUML)
  {
  };
};

class FPUSH_ : public TagWithoutArgument
{
public:
  inline
  FPUSH_() : TagWithoutArgument(FPUSH)
  {
  };
};

class FPOP_ : public TagWithoutArgument
{
public:
  inline
  FPOP_() : TagWithoutArgument(FPOP)
  {
  };
};

class FDIMT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FDIMT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMT)
  {
  };
  inline
  FDIMT_(unsigned int size_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMT, size_arg)
  {
  };
  inline unsigned int
  get_size()
  {
    return arg1;
  };
};

class FDIMST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FDIMST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMST)
  {
  };
  inline
  FDIMST_(const unsigned int size_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FDIMST, size_arg)
  {
  };
  inline unsigned int
  get_size()
  {
    return arg1;
  };
};

class FLDC_ : public TagWithOneArgument<double>
{
public:
  inline
  FLDC_() : TagWithOneArgument<double>::TagWithOneArgument(FLDC)
  {
  };
  inline
  FLDC_(const double value_arg) : TagWithOneArgument<double>::TagWithOneArgument(FLDC, value_arg)
  {
  };
  inline double
  get_value()
  {
    return arg1;
  };
};

class FLDU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDU)
  {
  };
  inline
  FLDU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDU, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FLDSU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDSU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDSU)
  {
  };
  inline
  FLDSU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDSU, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FLDR_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDR_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDR)
  {
  };
  inline
  FLDR_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDR, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FLDT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDT)
  {
  };
  inline
  FLDT_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDT, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FLDST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDST)
  {
  };
  inline
  FLDST_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDST, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPT_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPT_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPT)
  {
  };
  inline
  FSTPT_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPT, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPST_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPST_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPST)
  {
  };
  inline
  FSTPST_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPST, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPR_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPR_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPR)
  {
  };
  inline
  FSTPR_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPR, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPU)
  {
  };
  inline
  FSTPU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPU, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPSU_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPSU_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPSU)
  {
  };
  inline
  FSTPSU_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPSU, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPG_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPG, 0)
  {
  };
  inline
  FSTPG_(const unsigned int pos_arg) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPG, pos_arg)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  };
};

class FSTPG2_ : public TagWithTwoArguments<unsigned int, unsigned int>
{
public:
  inline
  FSTPG2_() : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FSTPG2, 0, 0)
  {
  };
  inline
  FSTPG2_(const unsigned int pos_arg1, const unsigned int pos_arg2) : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FSTPG2, pos_arg1, pos_arg2)
  {
  };
  inline unsigned int
  get_row()
  {
    return arg1;
  };
  inline unsigned int
  get_col()
  {
    return arg2;
  };
};

class FSTPG3_ : public TagWithFourArguments<unsigned int, unsigned int, int, unsigned int>
{
public:
  inline
  FSTPG3_() : TagWithFourArguments<unsigned int, unsigned int, int, unsigned int>::TagWithFourArguments(FSTPG3, 0, 0, 0, 0)
  {
  };
  inline
  FSTPG3_(const unsigned int pos_arg1, const unsigned int pos_arg2, const int pos_arg3, const unsigned int pos_arg4) : TagWithFourArguments<unsigned int, unsigned int, int, unsigned int>::TagWithFourArguments(FSTPG3, pos_arg1, pos_arg2, pos_arg3, pos_arg4)
  {
  };
  inline unsigned int
  get_row()
  {
    return arg1;
  };
  inline unsigned int
  get_col()
  {
    return arg2;
  };
  inline int
  get_lag()
  {
    return arg2;
  };
  inline unsigned int
  get_col_pos()
  {
    return arg4;
  };
};

class FUNARY_ : public TagWithOneArgument<uint8_t>
{
public:
  inline
  FUNARY_() : TagWithOneArgument<uint8_t>::TagWithOneArgument(FUNARY)
  {
  };
  inline
  FUNARY_(uint8_t op_type_arg) : TagWithOneArgument<uint8_t>::TagWithOneArgument(FUNARY, op_type_arg)
  {
  };
  inline uint8_t
  get_op_type()
  {
    return arg1;
  };
};

class FBINARY_ : public TagWithOneArgument<uint8_t>
{
public:
  inline
  FBINARY_() : TagWithOneArgument<uint8_t>::TagWithOneArgument(FBINARY)
  {
  };
  inline
  FBINARY_(const int op_type_arg) : TagWithOneArgument<uint8_t>::TagWithOneArgument(FBINARY, op_type_arg)
  {
  };
  inline uint8_t
  get_op_type()
  {
    return arg1;
  };
};

class FTRINARY_ : public TagWithOneArgument<uint8_t>
{
public:
  inline
  FTRINARY_() : TagWithOneArgument<uint8_t>::TagWithOneArgument(FTRINARY)
  {
  };
  inline
  FTRINARY_(const int op_type_arg) : TagWithOneArgument<uint8_t>::TagWithOneArgument(FTRINARY, op_type_arg)
  {
  };
  inline uint8_t
  get_op_type()
  {
    return arg1;
  };
};

class FOK_ : public TagWithOneArgument<int>
{
public:
  inline
  FOK_() : TagWithOneArgument<int>::TagWithOneArgument(FOK)
  {
  };
  inline
  FOK_(const int arg_arg) : TagWithOneArgument<int>::TagWithOneArgument(FOK, arg_arg)
  {
  };
  inline int
  get_arg()
  {
    return arg1;
  };
};

class FJMPIFEVAL_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FJMPIFEVAL_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FJMPIFEVAL)
  {
  };
  inline
  FJMPIFEVAL_(unsigned int arg_pos) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FJMPIFEVAL, arg_pos)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  }
};

class FJMP_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FJMP_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FJMP)
  {
  };
  inline
  FJMP_(unsigned int arg_pos) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FJMP, arg_pos)
  {
  };
  inline unsigned int
  get_pos()
  {
    return arg1;
  }
};

class FLDTEF_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FLDTEF_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDTEF)
  {
  };
  inline
  FLDTEF_(unsigned int number) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FLDTEF, number)
  {
  };
  inline unsigned int
  get_number()
  {
    return arg1;
  }
};

class FSTPTEF_ : public TagWithOneArgument<unsigned int>
{
public:
  inline
  FSTPTEF_() : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPTEF)
  {
  };
  inline
  FSTPTEF_(unsigned int number) : TagWithOneArgument<unsigned int>::TagWithOneArgument(FSTPTEF, number)
  {
  };
  inline unsigned int
  get_number()
  {
    return arg1;
  }
};

class FLDTEFD_ : public TagWithTwoArguments<unsigned int, unsigned int>
{
public:
  inline
  FLDTEFD_() : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FLDTEFD)
  {
  };
  inline
  FLDTEFD_(unsigned int indx, unsigned int row) : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FLDTEFD, indx, row)
  {
  };
  inline unsigned int
  get_indx()
  {
    return arg1;
  };
  inline unsigned int
  get_row()
  {
    return arg2;
  };
};

class FSTPTEFD_ : public TagWithTwoArguments<unsigned int, unsigned int>
{
public:
  inline
  FSTPTEFD_() : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FSTPTEFD)
  {
  };
  inline
  FSTPTEFD_(unsigned int indx, unsigned int row) : TagWithTwoArguments<unsigned int, unsigned int>::TagWithTwoArguments(FSTPTEFD, indx, row)
  {
  };
  inline unsigned int
  get_indx()
  {
    return arg1;
  };
  inline unsigned int
  get_row()
  {
    return arg2;
  };
};

class FLDTEFDD_ : public TagWithThreeArguments<unsigned int, unsigned int, unsigned int>
{
public:
  inline
  FLDTEFDD_() : TagWithThreeArguments<unsigned int, unsigned int, unsigned int>::TagWithThreeArguments(FLDTEFDD)
  {
  };
  inline
  FLDTEFDD_(unsigned int indx, unsigned int row, unsigned int col) : TagWithThreeArguments<unsigned int, unsigned int, unsigned int>::TagWithThreeArguments(FLDTEFDD, indx, row, col)
  {
  };
  inline unsigned int
  get_indx()
  {
    return arg1;
  };
  inline unsigned int
  get_row()
  {
    return arg2;
  };
  inline unsigned int
  get_col()
  {
    return arg3;
  };
};

class FSTPTEFDD_ : public TagWithThreeArguments<unsigned int, unsigned int, unsigned int>
{
public:
  inline
  FSTPTEFDD_() : TagWithThreeArguments<unsigned int, unsigned int, unsigned int>::TagWithThreeArguments(FSTPTEFDD)
  {
  };
  inline
  FSTPTEFDD_(unsigned int indx, unsigned int row, unsigned int col) : TagWithThreeArguments<unsigned int, unsigned int, unsigned int>::TagWithThreeArguments(FSTPTEF, indx, row, col)
  {
  };
  inline unsigned int
  get_indx()
  {
    return arg1;
  };
  inline unsigned int
  get_row()
  {
    return arg2;
  };
  inline unsigned int
  get_col()
  {
    return arg3;
  };
};

class FLDVS_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline
  FLDVS_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDVS)
  {
  };
  inline
  FLDVS_(uint8_t type_arg, const unsigned int pos_arg) : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDVS, type_arg, pos_arg)
  {
  };
  inline uint8_t
  get_type()
  {
    return arg1;
  };
  inline unsigned int
  get_pos()
  {
    return arg2;
  };
};

class FLDSV_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline
  FLDSV_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDSV)
  {
  };
  inline
  FLDSV_(const uint8_t type_arg, const unsigned int pos_arg) :
    TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FLDSV, type_arg, pos_arg)
  {
  };
  inline uint8_t
  get_type()
  {
    return arg1;
  };
  inline unsigned int
  get_pos()
  {
    return arg2;
  };
};

class FSTPSV_ : public TagWithTwoArguments<uint8_t, unsigned int>
{
public:
  inline
  FSTPSV_() : TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FSTPSV)
  {
  };
  inline
  FSTPSV_(const uint8_t type_arg, const unsigned int pos_arg) :
    TagWithTwoArguments<uint8_t, unsigned int>::TagWithTwoArguments(FSTPSV, type_arg, pos_arg)
  {
  };
  inline uint8_t
  get_type()
  {
    return arg1;
  };
  inline unsigned int
  get_pos()
  {
    return arg2;
  };
};

class FLDV_ : public TagWithThreeArguments<uint8_t, unsigned int, int>
{
public:
  inline
  FLDV_() : TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV)
  {
  };
  inline
  FLDV_(const int type_arg, const unsigned int pos_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV, type_arg, pos_arg, 0)
  {
  };
  inline
  FLDV_(const int type_arg, const unsigned int pos_arg, const int lead_lag_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FLDV, type_arg, pos_arg, lead_lag_arg)
  {
  };
  inline uint8_t
  get_type()
  {
    return arg1;
  };
  inline unsigned int
  get_pos()
  {
    return arg2;
  };
  inline int
  get_lead_lag()
  {
    return arg3;
  };
};

class FSTPV_ : public TagWithThreeArguments<uint8_t, unsigned int, int>
{
public:
  inline
  FSTPV_() : TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV)
  {
  };
  inline
  FSTPV_(const int type_arg, const unsigned int pos_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV, type_arg, pos_arg, 0)
  {
  };
  inline
  FSTPV_(const int type_arg, const unsigned int pos_arg, const int lead_lag_arg) :
    TagWithThreeArguments<uint8_t, unsigned int, int>::TagWithThreeArguments(FSTPV, type_arg, pos_arg, lead_lag_arg)
  {
  };
  inline uint8_t
  get_type()
  {
    return arg1;
  };
  inline unsigned int
  get_pos()
  {
    return arg2;
  };
  inline int
  get_lead_lag()
  {
    return arg3;
  };
};

class FCALL_ : public TagWithFourArguments<unsigned int, unsigned int, string, unsigned int>
{
  string func_name;
  string arg_func_name;
  unsigned int add_input_arguments, row, col;
  external_function_type function_type;
public:
  inline
  FCALL_() : TagWithFourArguments<unsigned int, unsigned int, string, unsigned int>::TagWithFourArguments(FCALL)
  {
    arg_func_name = "";
    add_input_arguments = 0;
    row = 0;
    col = 0;
    function_type = ExternalFunctionWithoutDerivative;
  };
  inline
  FCALL_(unsigned int nb_output_arguments, unsigned int nb_input_arguments, string f_name, unsigned int indx) :
    TagWithFourArguments<unsigned int, unsigned int, string, unsigned int>::TagWithFourArguments(FCALL, nb_output_arguments, nb_input_arguments, f_name, indx)
  {
    arg_func_name = "";
    add_input_arguments = 0;
    row = 0;
    col = 0;
    function_type = ExternalFunctionWithoutDerivative;
    func_name = f_name;
  };
  inline string
  get_function_name()
  {
    //printf("get_function_name => func_name=%s\n",func_name.c_str());fflush(stdout);
    return func_name;
  };
  inline unsigned int
  get_nb_output_arguments()
  {
    return arg1;
  };
  inline unsigned int
  get_nb_input_arguments()
  {
    return arg2;
  };
  inline unsigned int
  get_indx()
  {
    return arg4;
  };
  inline void
  set_arg_func_name(string arg_arg_func_name)
  {
    arg_func_name = arg_arg_func_name;
  };
  inline string
  get_arg_func_name()
  {
    return arg_func_name;
  };
  inline void
  set_nb_add_input_arguments(unsigned int arg_add_input_arguments)
  {
    add_input_arguments = arg_add_input_arguments;
  };
  inline unsigned int
  get_nb_add_input_arguments()
  {
    return add_input_arguments;
  };
  inline void
  set_row(unsigned int arg_row)
  {
    row = arg_row;
  };
  inline unsigned int
  get_row()
  {
    return row;
  }
  inline void
  set_col(unsigned int arg_col)
  {
    col = arg_col;
  };
  inline unsigned int
  get_col()
  {
    return col;
  };
  inline void
  set_function_type(external_function_type arg_function_type)
  {
    function_type = arg_function_type;
  };
  inline external_function_type
  get_function_type()
  {
    return (function_type);
  }
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(&op_code), sizeof(op_code));
    CompileCode.write(reinterpret_cast<char *>(&arg1), sizeof(arg1));
    CompileCode.write(reinterpret_cast<char *>(&arg2), sizeof(arg2));
    CompileCode.write(reinterpret_cast<char *>(&arg4), sizeof(arg4));
    CompileCode.write(reinterpret_cast<char *>(&add_input_arguments), sizeof(add_input_arguments));
    CompileCode.write(reinterpret_cast<char *>(&row), sizeof(row));
    CompileCode.write(reinterpret_cast<char *>(&col), sizeof(col));
    CompileCode.write(reinterpret_cast<char *>(&function_type), sizeof(function_type));
    size_t size = func_name.size();
    CompileCode.write(reinterpret_cast<char *>(&size), sizeof(int));
    const char *name = func_name.c_str();
    CompileCode.write(reinterpret_cast<const char *>(name), func_name.size());
    size = arg_func_name.size();
    CompileCode.write(reinterpret_cast<char *>(&size), sizeof(int));
    name = arg_func_name.c_str();
    CompileCode.write(reinterpret_cast<const char *>(name), arg_func_name.size());
    instruction_number++;
  };
#ifdef BYTE_CODE

  inline uint8_t *
  load(uint8_t *code)
  {
    op_code = FCALL; code += sizeof(op_code);
    memcpy(&arg1, code, sizeof(arg1)); code += sizeof(arg1);
    memcpy(&arg2, code, sizeof(arg2)); code += sizeof(arg2);
    memcpy(&arg4, code, sizeof(arg4)); code += sizeof(arg4);
    memcpy(&add_input_arguments, code, sizeof(add_input_arguments)); code += sizeof(add_input_arguments);
    memcpy(&row, code, sizeof(row)); code += sizeof(row);
    memcpy(&col, code, sizeof(col)); code += sizeof(col);
    memcpy(&function_type, code, sizeof(function_type)); code += sizeof(function_type);
    int size;
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    char *name = (char *) mxMalloc((size+1)*sizeof(char));
    memcpy(name, code, size); code += size;
    name[size] = 0;
    func_name = name;
    mxFree(name);
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    name = (char *) mxMalloc((size+1)*sizeof(char));
    memcpy(name, code, size); code += size;
    name[size] = 0;
    arg_func_name = name;
    mxFree(name);
    return code;
  }
#endif
};

class FNUMEXPR_ : public TagWithOneArgument<ExpressionType>
{
private:
  unsigned int equation;
  uint16_t dvariable1, dvariable2, dvariable3;
  int8_t lag1, lag2, lag3;
public:
  inline
  FNUMEXPR_() : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR)
  {
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                               dvariable1(0), dvariable2(0), dvariable3(0), lag1(0), lag2(0), lag3(0)
  {
    equation = equation_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                                                            dvariable2(0), dvariable3(0), lag1(0), lag2(0), lag3(0)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg, int lag1_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                                                                          dvariable2(0), dvariable3(0), lag2(0), lag3(0)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
    lag1 = lag1_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg, unsigned int dvariable2_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                                                                                         dvariable3(0), lag1(0), lag2(0), lag3(0)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
    dvariable2 = dvariable2_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg, int lag1_arg, unsigned int dvariable2_arg, int lag2_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                                                                                                                     dvariable3(0), lag3(0)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
    lag1 = lag1_arg;
    dvariable2 = dvariable2_arg;
    lag2 = lag2_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg, unsigned int dvariable2_arg, unsigned int dvariable3_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type),
                                                                                                                                                                      lag1(0), lag2(0), lag3(0)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
    dvariable2 = dvariable2_arg;
    dvariable3 = dvariable3_arg;
  };
  inline
  FNUMEXPR_(const ExpressionType expression_type, unsigned int equation_arg, unsigned int dvariable1_arg, int lag1_arg, unsigned int dvariable2_arg, int lag2_arg, unsigned int dvariable3_arg, int lag3_arg) : TagWithOneArgument<ExpressionType>::TagWithOneArgument(FNUMEXPR, expression_type)
  {
    equation = equation_arg;
    dvariable1 = dvariable1_arg;
    lag1 = lag1_arg;
    dvariable2 = dvariable2_arg;
    lag2 = lag2_arg;
    dvariable3 = dvariable3_arg;
    lag3 = lag3_arg;
  };
  inline ExpressionType
  get_expression_type()
  {
    return arg1;
  }
  inline unsigned int
  get_equation()
  {
    return equation;
  };
  inline unsigned int
  get_dvariable1()
  {
    return dvariable1;
  };
  inline int
  get_lag1()
  {
    return lag1;
  };
  inline unsigned int
  get_dvariable2()
  {
    return dvariable2;
  };
  inline int
  get_lag2()
  {
    return lag2;
  };
  inline unsigned int
  get_dvariable3()
  {
    return dvariable3;
  };
  inline int
  get_lag3()
  {
    return lag3;
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(this), sizeof(FNUMEXPR_));
    instruction_number++;
  };
};

class FBEGINBLOCK_
{
private:
  uint8_t op_code;
  int size;
  uint8_t type;
  vector<int> variable;
  vector<int> equation;
  vector<unsigned int> other_endogenous;
  vector<unsigned int> exogenous;
  vector<unsigned int> det_exogenous;
  bool is_linear;
  vector<Block_contain_type> Block_Contain_;
  int endo_nbr;
  int Max_Lag;
  int Max_Lead;
  int u_count_int;
  int nb_col_jacob;
  unsigned int det_exo_size, exo_size, other_endo_size;
  unsigned int nb_col_other_endo_jacob;
public:
  inline
  FBEGINBLOCK_()
  {
    op_code = FBEGINBLOCK; size = 0; type = UNKNOWN; /*variable = NULL; equation = NULL;*/
    is_linear = false; endo_nbr = 0; Max_Lag = 0; Max_Lead = 0; u_count_int = 0; nb_col_jacob = 0;
  };
  inline
  FBEGINBLOCK_(unsigned int size_arg, BlockSimulationType type_arg, int unsigned first_element, int unsigned block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int endo_nbr_arg, int Max_Lag_arg, int Max_Lead_arg, int &u_count_int_arg, int nb_col_jacob_arg,
               unsigned int det_exo_size_arg, unsigned int exo_size_arg, unsigned int other_endo_size_arg, unsigned int nb_col_other_endo_jacob_arg,
               const vector<unsigned int> &det_exogenous_arg, const vector<unsigned int> &exogenous_arg, const vector<unsigned int> &other_endogenous_arg)
  {
    op_code = FBEGINBLOCK; size = size_arg; type = type_arg;
    variable = vector<int>(variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size));
    equation = vector<int>(equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size));
    det_exogenous = vector<unsigned int>(det_exogenous_arg);
    exogenous = vector<unsigned int>(exogenous_arg);
    other_endogenous = vector<unsigned int>(other_endogenous_arg);
    is_linear = is_linear_arg; endo_nbr = endo_nbr_arg; Max_Lag = Max_Lag_arg; Max_Lead = Max_Lead_arg; u_count_int = u_count_int_arg;
    nb_col_jacob = nb_col_jacob_arg; det_exo_size = det_exo_size_arg; exo_size = exo_size_arg; other_endo_size = other_endo_size_arg;
    nb_col_other_endo_jacob = nb_col_other_endo_jacob_arg;
  };
  inline
  FBEGINBLOCK_(unsigned int size_arg, BlockSimulationType type_arg, int unsigned first_element, int unsigned block_size,
               const vector<int> &variable_arg, const vector<int> &equation_arg,
               bool is_linear_arg, int endo_nbr_arg, int Max_Lag_arg, int Max_Lead_arg, int &u_count_int_arg, int nb_col_jacob_arg)
  {
    op_code = FBEGINBLOCK; size = size_arg; type = type_arg;
    variable = vector<int>(variable_arg.begin()+first_element, variable_arg.begin()+(first_element+block_size));
    equation = vector<int>(equation_arg.begin()+first_element, equation_arg.begin()+(first_element+block_size));
    is_linear = is_linear_arg; endo_nbr = endo_nbr_arg; Max_Lag = Max_Lag_arg; Max_Lead = Max_Lead_arg; u_count_int = u_count_int_arg;
    nb_col_jacob = nb_col_jacob_arg;
    det_exo_size = 0; exo_size = 0; other_endo_size = 0;
    nb_col_other_endo_jacob = 0;
  }
  inline unsigned int
  get_size()
  {
    return size;
  };
  inline uint8_t
  get_type()
  {
    return type;
  };
  inline bool
  get_is_linear()
  {
    return is_linear;
  };
  inline int
  get_endo_nbr()
  {
    return endo_nbr;
  };
  inline int
  get_Max_Lag()
  {
    return Max_Lag;
  };
  inline int
  get_Max_Lead()
  {
    return Max_Lead;
  };
  inline int
  get_u_count_int()
  {
    return u_count_int;
  };
  inline vector<Block_contain_type>
  get_Block_Contain()
  {
    return Block_Contain_;
  };
  inline int
  get_nb_col_jacob()
  {
    return nb_col_jacob;
  };
  inline unsigned int
  get_exo_size()
  {
    return exo_size;
  };
  inline unsigned int
  get_det_exo_size()
  {
    return det_exo_size;
  };
  inline unsigned int
  get_other_endo_size()
  {
    return other_endo_size;
  };
  inline unsigned int
  get_nb_col_other_endo_jacob()
  {
    return nb_col_other_endo_jacob;
  };
  inline void
  write(ostream &CompileCode, unsigned int &instruction_number)
  {
    CompileCode.write(reinterpret_cast<char *>(&op_code), sizeof(op_code));
    CompileCode.write(reinterpret_cast<char *>(&size), sizeof(size));
    CompileCode.write(reinterpret_cast<char *>(&type), sizeof(type));
    for (int i = 0; i < size; i++)
      {
        CompileCode.write(reinterpret_cast<char *>(&variable[i]), sizeof(variable[0]));
        CompileCode.write(reinterpret_cast<char *>(&equation[i]), sizeof(equation[0]));
      }
    if (type == SOLVE_TWO_BOUNDARIES_SIMPLE || type == SOLVE_TWO_BOUNDARIES_COMPLETE
        || type == SOLVE_BACKWARD_COMPLETE || type == SOLVE_FORWARD_COMPLETE)
      {
        CompileCode.write(reinterpret_cast<char *>(&is_linear), sizeof(is_linear));
        CompileCode.write(reinterpret_cast<char *>(&endo_nbr), sizeof(endo_nbr));
        CompileCode.write(reinterpret_cast<char *>(&Max_Lag), sizeof(Max_Lag));
        CompileCode.write(reinterpret_cast<char *>(&Max_Lead), sizeof(Max_Lead));
        CompileCode.write(reinterpret_cast<char *>(&u_count_int), sizeof(u_count_int));
      }
    CompileCode.write(reinterpret_cast<char *>(&nb_col_jacob), sizeof(nb_col_jacob));
    CompileCode.write(reinterpret_cast<char *>(&det_exo_size), sizeof(det_exo_size));
    CompileCode.write(reinterpret_cast<char *>(&exo_size), sizeof(exo_size));
    CompileCode.write(reinterpret_cast<char *>(&other_endo_size), sizeof(other_endo_size));
    CompileCode.write(reinterpret_cast<char *>(&nb_col_other_endo_jacob), sizeof(nb_col_other_endo_jacob));

    for (unsigned int i = 0; i < det_exo_size; i++)
      CompileCode.write(reinterpret_cast<char *>(&det_exogenous[i]), sizeof(det_exogenous[0]));
    for (unsigned int i = 0; i < exo_size; i++)
      CompileCode.write(reinterpret_cast<char *>(&exogenous[i]), sizeof(exogenous[0]));
    for (unsigned int i = 0; i < other_endo_size; i++)
      CompileCode.write(reinterpret_cast<char *>(&other_endogenous[i]), sizeof(other_endogenous[0]));
    instruction_number++;
  };
#ifdef BYTE_CODE

  inline uint8_t *
  load(uint8_t *code)
  {
    op_code = FBEGINBLOCK; code += sizeof(op_code);
    memcpy(&size, code, sizeof(size)); code += sizeof(size);
    memcpy(&type, code, sizeof(type)); code += sizeof(type);
    for (int i = 0; i < size; i++)
      {
        Block_contain_type bc;
        memcpy(&bc.Variable, code, sizeof(bc.Variable)); code += sizeof(bc.Variable);
        memcpy(&bc.Equation, code, sizeof(bc.Equation)); code += sizeof(bc.Equation);
        Block_Contain_.push_back(bc);
      }
    if (type == SOLVE_TWO_BOUNDARIES_SIMPLE || type == SOLVE_TWO_BOUNDARIES_COMPLETE
        || type == SOLVE_BACKWARD_COMPLETE || type == SOLVE_FORWARD_COMPLETE)
      {
        memcpy(&is_linear, code, sizeof(is_linear)); code += sizeof(is_linear);
        memcpy(&endo_nbr, code, sizeof(endo_nbr)); code += sizeof(endo_nbr);
        memcpy(&Max_Lag, code, sizeof(Max_Lag)); code += sizeof(Max_Lag);
        memcpy(&Max_Lead, code, sizeof(Max_Lead)); code += sizeof(Max_Lead);
        memcpy(&u_count_int, code, sizeof(u_count_int)); code += sizeof(u_count_int);
      }
    memcpy(&nb_col_jacob, code, sizeof(nb_col_jacob)); code += sizeof(nb_col_jacob);
    memcpy(&det_exo_size, code, sizeof(det_exo_size)); code += sizeof(det_exo_size);
    memcpy(&exo_size, code, sizeof(exo_size)); code += sizeof(exo_size);
    memcpy(&other_endo_size, code, sizeof(other_endo_size)); code += sizeof(other_endo_size);
    memcpy(&nb_col_other_endo_jacob, code, sizeof(nb_col_other_endo_jacob)); code += sizeof(nb_col_other_endo_jacob);

    for (unsigned int i = 0; i < det_exo_size; i++)
      {
        unsigned int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        det_exogenous.push_back(tmp_i);
      }
    for (unsigned int i = 0; i < exo_size; i++)
      {
        unsigned int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        exogenous.push_back(tmp_i);
      }
    for (unsigned int i = 0; i < other_endo_size; i++)
      {
        unsigned int tmp_i;
        memcpy(&tmp_i, code, sizeof(tmp_i)); code += sizeof(tmp_i);
        other_endogenous.push_back(tmp_i);
      }
    return code;
  };
#endif
};

#ifdef BYTE_CODE
typedef vector<pair<Tags, void * > > tags_liste_t;
class CodeLoad
{
private:
  uint8_t *code;
  unsigned int nb_blocks;
  vector<size_t> begin_block;
public:

  inline unsigned int
  get_block_number()
  {
    return nb_blocks;
  };

  size_t inline
  get_begin_block(int block)
  {
    return begin_block[block];
  }
  inline void *
  get_current_code()
  {
    return code;
  };
  inline tags_liste_t
  get_op_code(string file_name)
  {
    tags_liste_t tags_liste;
    ifstream CompiledCode;
    streamoff Code_Size;
    CompiledCode.open((file_name + ".cod").c_str(), std::ios::in | std::ios::binary| std::ios::ate);
    if (!CompiledCode.is_open())
      {
        return tags_liste;
      }
    Code_Size = CompiledCode.tellg();
    CompiledCode.seekg(std::ios::beg);
    code = (uint8_t *) mxMalloc(Code_Size);
    CompiledCode.seekg(0);
    CompiledCode.read(reinterpret_cast<char *>(code), Code_Size);
    CompiledCode.close();
    nb_blocks = 0;
    bool done = false;
    int instruction = 0;
    while (!done)
      {
        switch (*code)
          {
          case FLDZ:
# ifdef DEBUGL
            mexPrintf("FLDZ = %d size = %d\n", FLDZ, sizeof(FLDZ_));
# endif
            tags_liste.push_back(make_pair(FLDZ, code));
            code += sizeof(FLDZ_);
            break;
          case FEND:
# ifdef DEBUGL
            mexPrintf("FEND\n");
# endif
            tags_liste.push_back(make_pair(FEND, code));
            code += sizeof(FEND_);
            done = true;
            break;
          case FENDBLOCK:
# ifdef DEBUGL
            mexPrintf("FENDBLOCK\n");
# endif
            tags_liste.push_back(make_pair(FENDBLOCK, code));
            code += sizeof(FENDBLOCK_);
            break;
          case FENDEQU:
# ifdef DEBUGL
            mexPrintf("FENDEQU\n");
# endif
            tags_liste.push_back(make_pair(FENDEQU, code));
            code += sizeof(FENDEQU_);
            break;
          case FCUML:
# ifdef DEBUGL
            mexPrintf("FCUML\n");
# endif
            tags_liste.push_back(make_pair(FCUML, code));
            code += sizeof(FCUML_);
            break;
          case FDIMT:
# ifdef DEBUGL
            mexPrintf("FDIMT = %d size = %d\n", FDIMT, sizeof(FDIMT_));
# endif
            tags_liste.push_back(make_pair(FDIMT, code));
            code += sizeof(FDIMT_);
            break;
          case FDIMST:
# ifdef DEBUGL
            mexPrintf("FDIMST\n");
# endif
            tags_liste.push_back(make_pair(FDIMST, code));
            code += sizeof(FDIMST_);
            break;
          case FNUMEXPR:
# ifdef DEBUGL
            mexPrintf("FNUMEXPR\n");
# endif
            tags_liste.push_back(make_pair(FNUMEXPR, code));
            code += sizeof(FNUMEXPR_);
            break;
          case FLDC:
# ifdef DEBUGL
            mexPrintf("FLDC\n");
# endif
            tags_liste.push_back(make_pair(FLDC, code));
            code += sizeof(FLDC_);
            break;
          case FLDU:
# ifdef DEBUGL
            mexPrintf("FLDU\n");
# endif
            tags_liste.push_back(make_pair(FLDU, code));
            code += sizeof(FLDU_);
            break;
          case FLDSU:
# ifdef DEBUGL
            mexPrintf("FLDSU\n");
# endif
            tags_liste.push_back(make_pair(FLDSU, code));
            code += sizeof(FLDSU_);
            break;
          case FLDR:
# ifdef DEBUGL
            mexPrintf("FLDR\n");
# endif
            tags_liste.push_back(make_pair(FLDR, code));
            code += sizeof(FLDR_);
            break;
          case FLDT:
# ifdef DEBUGL
            mexPrintf("FLDT\n");
# endif
            tags_liste.push_back(make_pair(FLDT, code));
            code += sizeof(FLDT_);
            break;
          case FLDST:
# ifdef DEBUGL
            mexPrintf("FLDST\n");
# endif
            tags_liste.push_back(make_pair(FLDST, code));
            code += sizeof(FLDST_);
            break;
          case FSTPT:
# ifdef DEBUGL
            mexPrintf("FSTPT = %d size = %d\n", FSTPT, sizeof(FSTPT_));
# endif
            tags_liste.push_back(make_pair(FSTPT, code));
            code += sizeof(FSTPT_);
            break;
          case FSTPST:
# ifdef DEBUGL
            mexPrintf("FSTPST\n");
# endif
            tags_liste.push_back(make_pair(FSTPST, code));
            code += sizeof(FSTPST_);
            break;
          case FSTPR:
# ifdef DEBUGL
            mexPrintf("FSTPR\n");
# endif
            tags_liste.push_back(make_pair(FSTPR, code));
            code += sizeof(FSTPR_);
            break;
          case FSTPU:
# ifdef DEBUGL
            mexPrintf("FSTPU\n");
# endif
            tags_liste.push_back(make_pair(FSTPU, code));
            code += sizeof(FSTPU_);
            break;
          case FSTPSU:
# ifdef DEBUGL
            mexPrintf("FSTPSU\n");
# endif
            tags_liste.push_back(make_pair(FSTPSU, code));
            code += sizeof(FSTPSU_);
            break;
          case FSTPG:
# ifdef DEBUGL
            mexPrintf("FSTPG\n");
# endif
            tags_liste.push_back(make_pair(FSTPG, code));
            code += sizeof(FSTPG_);
            break;
          case FSTPG2:
# ifdef DEBUGL
            mexPrintf("FSTPG2\n");
# endif
            tags_liste.push_back(make_pair(FSTPG2, code));
            code += sizeof(FSTPG2_);
            break;
          case FSTPG3:
# ifdef DEBUGL
            mexPrintf("FSTPG3\n");
# endif
            tags_liste.push_back(make_pair(FSTPG3, code));
            code += sizeof(FSTPG3_);
            break;
          case FUNARY:
# ifdef DEBUGL
            mexPrintf("FUNARY\n");
# endif
            tags_liste.push_back(make_pair(FUNARY, code));
            code += sizeof(FUNARY_);
            break;
          case FBINARY:
# ifdef DEBUGL
            mexPrintf("FBINARY\n");
# endif
            tags_liste.push_back(make_pair(FBINARY, code));
            code += sizeof(FBINARY_);
            break;
          case FTRINARY:
# ifdef DEBUGL
            mexPrintf("FTRINARY\n");
# endif
            tags_liste.push_back(make_pair(FTRINARY, code));
            code += sizeof(FTRINARY_);
            break;
          case FOK:
# ifdef DEBUGL
            mexPrintf("FOK\n");
# endif
            tags_liste.push_back(make_pair(FOK, code));
            code += sizeof(FOK_);
            break;
          case FLDVS:
# ifdef DEBUGL
            mexPrintf("FLDVS\n");
# endif
            tags_liste.push_back(make_pair(FLDVS, code));
            code += sizeof(FLDVS_);
            break;
          case FLDSV:
# ifdef DEBUGL
            mexPrintf("FLDSV\n");
# endif
            tags_liste.push_back(make_pair(FLDSV, code));
            code += sizeof(FLDSV_);
            break;
          case FSTPSV:
# ifdef DEBUGL
            mexPrintf("FSTPSV\n");
# endif
            tags_liste.push_back(make_pair(FSTPSV, code));
            code += sizeof(FSTPSV_);
            break;
          case FLDV:
# ifdef DEBUGL
            mexPrintf("FLDV\n");
# endif
            tags_liste.push_back(make_pair(FLDV, code));
            code += sizeof(FLDV_);
            break;
          case FSTPV:
# ifdef DEBUGL
            mexPrintf("FSTPV\n");
# endif
            tags_liste.push_back(make_pair(FSTPV, code));
            code += sizeof(FSTPV_);
            break;
          case FBEGINBLOCK:
# ifdef DEBUGL
            mexPrintf("FBEGINBLOCK\n");
# endif
            {
              FBEGINBLOCK_ *fbegin_block = new FBEGINBLOCK_;

              code = fbegin_block->load(code);

              begin_block.push_back(tags_liste.size());
              tags_liste.push_back(make_pair(FBEGINBLOCK, fbegin_block));
              nb_blocks++;
            }
            break;
          case FJMPIFEVAL:
# ifdef DEBUGL
            mexPrintf("FJMPIFEVAL\n");
# endif
            tags_liste.push_back(make_pair(FJMPIFEVAL, code));
            code += sizeof(FJMPIFEVAL_);
            break;
          case FJMP:
# ifdef DEBUGL
            mexPrintf("FJMP\n");
# endif
            tags_liste.push_back(make_pair(FJMP, code));
            code += sizeof(FJMP_);
            break;
          case FCALL:
            {
# ifdef DEBUGL
              mexPrintf("FCALL\n");
# endif
              FCALL_ *fcall = new FCALL_;

              code = fcall->load(code);

              tags_liste.push_back(make_pair(FCALL, fcall));
# ifdef DEBUGL
              mexPrintf("FCALL finish\n"); mexEvalString("drawnow;");
              mexPrintf("-- *code=%d\n", *code); mexEvalString("drawnow;");
# endif
            }
            break;
          case FPUSH:
# ifdef DEBUGL
            mexPrintf("FPUSH\n");
# endif
            tags_liste.push_back(make_pair(FPUSH, code));
            code += sizeof(FPUSH_);
            break;
          case FPOP:
# ifdef DEBUGL
            mexPrintf("FPOP\n");
# endif
            tags_liste.push_back(make_pair(FPOP, code));
            code += sizeof(FPOP_);
            break;
          case FLDTEF:
# ifdef DEBUGL
            mexPrintf("FLDTEF\n");
# endif
            tags_liste.push_back(make_pair(FLDTEF, code));
            code += sizeof(FLDTEF_);
            break;
          case FSTPTEF:
# ifdef DEBUGL
            mexPrintf("FSTPTEF\n");
# endif
            tags_liste.push_back(make_pair(FSTPTEF, code));
            code += sizeof(FSTPTEF_);
            break;
          case FLDTEFD:
# ifdef DEBUGL
            mexPrintf("FLDTEFD\n");
# endif
            tags_liste.push_back(make_pair(FLDTEFD, code));
            code += sizeof(FLDTEFD_);
            break;
          case FSTPTEFD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFD\n");
# endif
            tags_liste.push_back(make_pair(FSTPTEFD, code));
            code += sizeof(FSTPTEFD_);
            break;
          case FLDTEFDD:
# ifdef DEBUGL
            mexPrintf("FLDTEFDD\n");
# endif
            tags_liste.push_back(make_pair(FLDTEFDD, code));
            code += sizeof(FLDTEFDD_);
            break;
          case FSTPTEFDD:
# ifdef DEBUGL
            mexPrintf("FSTPTEFDD\n");
# endif
            tags_liste.push_back(make_pair(FSTPTEFDD, code));
            code += sizeof(FSTPTEFDD_);
            break;
          default:
            mexPrintf("Unknown Tag value=%d code=%x\n", *code, code);
            done = true;
          }
        instruction++;
      }
    return tags_liste;
  };
};
#endif
#pragma pack(pop)
#endif

