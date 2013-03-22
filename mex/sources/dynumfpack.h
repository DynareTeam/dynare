/*
 * Defines the prototypes for BLAS Fortran functions.
 *
 * Also defines a typedef blas_int to be used for all integers passed to BLAS
 * functions.
 *
 * When used in the context of a MATLAB MEX file, you must define MATLAB_MEX_FILE
 * and MATLAB_VERSION (for version 7.4, define it to 0x0704).
 *
 *
 * Copyright (C) 2009-2011 Dynare Team
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

#ifndef _DYNUMFPACK_H
#define _DYNUMFPACK_H

/* Starting from version 7.8, MATLAB BLAS expects ptrdiff_t arguments for integers */
#if defined(MATLAB_MEX_FILE) && MATLAB_VERSION >= 0x0708
# ifdef __cplusplus
#  include <cstddef>
# else
#  include <stddef.h>
# endif
#endif
/*
#if defined(MATLAB_MEX_FILE) && defined(_WIN32) && !defined(_MSC_VER)
# define FORTRAN_WRAPPER(x) x
#else
# define FORTRAN_WRAPPER(x) x
#endif
*/
#ifdef __cplusplus
extern "C" {
#endif


/* -------------------------------------------------------------------------- */
/* size of Info and Control arrays */
/* -------------------------------------------------------------------------- */

/* These might be larger in future versions, since there are only 3 unused
 * entries in Info, and no unused entries in Control. */

#define UMFPACK_INFO 90
#define UMFPACK_CONTROL 20
/* used in all UMFPACK_report_* routines: */
#define UMFPACK_PRL 0			/* print level */
/* returned by all routines that use Info: */
#define UMFPACK_OK (0)
#define UMFPACK_STATUS 0	/* UMFPACK_OK, or other result */


typedef long long int SuiteSparse_long;

#define umfpack_dl_defaults FORTRAN_WRAPPER(umfpack_dl_defaults)
void umfpack_dl_defaults(double Control[UMFPACK_CONTROL]);

#define umfpack_dl_symbolic FORTRAN_WRAPPER(umfpack_dl_symbolic)
SuiteSparse_long umfpack_dl_symbolic(SuiteSparse_long n_row, SuiteSparse_long n_col,
                                     const SuiteSparse_long Ap [ ], const SuiteSparse_long Ai [ ],
                                     const double Ax [ ], void **Symbolic,
                                     const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

#define umfpack_dl_numeric FORTRAN_WRAPPER(umfpack_dl_numeric)
SuiteSparse_long umfpack_dl_numeric(const SuiteSparse_long Ap [ ], const SuiteSparse_long Ai [ ],
                                    const double Ax [ ], void *Symbolic, void **Numeric,
                                    const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

#define umfpack_dl_solve FORTRAN_WRAPPER(umfpack_dl_solve)
SuiteSparse_long umfpack_dl_solve(SuiteSparse_long sys, const SuiteSparse_long Ap [ ],
                                  const SuiteSparse_long Ai [ ], const double Ax [ ],
                                  double X [ ], const double B [ ], void *Numeric,
                                  const double Control [UMFPACK_CONTROL], double Info [UMFPACK_INFO]);

#define umfpack_dl_report_info FORTRAN_WRAPPER(umfpack_dl_report_info)
void umfpack_dl_report_info(const double Control [UMFPACK_CONTROL],
                            const double Info [UMFPACK_INFO]);

#define umfpack_dl_report_status FORTRAN_WRAPPER(umfpack_dl_report_status)
void umfpack_dl_report_status(const double Control [UMFPACK_CONTROL],
                              SuiteSparse_long status);

#define umfpack_dl_free_symbolic FORTRAN_WRAPPER(umfpack_dl_free_symbolic)
void umfpack_dl_free_symbolic(void **Symbolic);

#define umfpack_dl_free_numeric FORTRAN_WRAPPER(umfpack_dl_free_numeric)
void umfpack_dl_free_numeric(void **Numeric);


#define umfpack_dl_load_symbolic  FORTRAN_WRAPPER(umfpack_dl_load_symbolic )
SuiteSparse_long umfpack_dl_load_symbolic (void **Symbolic, char *filename) ;

#define umfpack_dl_load_numeric  FORTRAN_WRAPPER(umfpack_dl_load_numeric )
SuiteSparse_long umfpack_dl_load_numeric (void **Numeric, char *filename) ;

#define umfpack_dl_save_symbolic  FORTRAN_WRAPPER(umfpack_dl_save_symbolic )
SuiteSparse_long umfpack_dl_save_symbolic (void *Symbolic, char *filename) ;

#define umfpack_dl_save_numeric  FORTRAN_WRAPPER(umfpack_dl_save_numeric )
SuiteSparse_long umfpack_dl_save_numeric (void *Numeric, char *filename) ;


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif /* DYNUMFPACK */
