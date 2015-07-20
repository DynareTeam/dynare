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

#ifndef _SIGMAEINITIALIZATION_HH
#define _SIGMAEINITIALIZATION_HH

#include <string>
#include <vector>

#include "ExprNode.hh"
#include "Statement.hh"

//! Stores a Sigma_e statement
class SigmaeStatement : public Statement
{
public:
  //! Matrix form (lower or upper triangular) enum
  enum matrix_form_t
    {
      eLower = 0,              //!< Lower triangular matrix
      eUpper = 1               //!< Upper triangular matrix
    };
  //! Type of a matrix row
  typedef vector<expr_t> row_t;
  //! Type of a complete matrix
  typedef vector<row_t> matrix_t;

  //! An exception indicating that a matrix is neither upper triangular nor lower triangular
  class MatrixFormException
  {
  };
private:
  //! The matrix
  const matrix_t matrix;
  //! Matrix form (lower or upper)
  const matrix_form_t matrix_form;

  //! Returns the type (upper or lower triangular) of a given matrix
  /*! Throws an exception if it is neither upper triangular nor lower triangular */
  static matrix_form_t determineMatrixForm(const matrix_t &matrix) throw (MatrixFormException);

public:
  SigmaeStatement(const matrix_t &matrix_arg) throw (MatrixFormException);
  virtual void writeOutput(ostream &output, const string &basename, bool minimal_workspace) const;
  virtual void writeM_Output(ostream &output, const string &basename, bool minimal_workspace);
};

#endif
