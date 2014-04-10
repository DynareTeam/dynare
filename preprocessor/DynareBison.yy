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

%skeleton "lalr1.cc"
%require "2.5"
%defines

%code top {
class ParsingDriver;
}

%code requires {
#include "ExprNode.hh"
#include "CodeInterpreter.hh"
}

%name-prefix "Dynare"

%parse-param { ParsingDriver &driver }
%lex-param { ParsingDriver &driver }

%locations
%initial-action
{
  // Initialize the locations' filenames to the filename maintained by the lexer
  @$.begin.filename = @$.end.filename = &(driver.lexer->filename);
}

%debug
%error-verbose

%union
{
  string *string_val;
  expr_t node_val;
  SymbolType symbol_type_val;
  vector<string *> *vector_string_val;
  vector<int> *vector_int_val;
  PriorDistributions prior_distributions_val;
};

%code {
/* Little hack: we redefine the macro which computes the locations, because
   we need to access the location from within the parsing driver for error
   and warning messages. */
#define YYLLOC_DEFAULT(Current, Rhs, N)                         \
  do {                                                          \
    if (N)                                                      \
      {                                                         \
        (Current).begin = YYRHSLOC(Rhs, 1).begin;               \
        (Current).end   = YYRHSLOC(Rhs, N).end;                 \
      }                                                         \
    else                                                        \
      {                                                         \
        (Current).begin = (Current).end = YYRHSLOC(Rhs, 0).end;	\
      }                                                         \
    driver.location = (Current);                                \
  } while(false)

#include "ParsingDriver.hh"

/* this "connects" the bison parser in the driver to the flex scanner class
 * object. it defines the yylex() function call to pull the next token from the
 * current lexer object of the driver context. */
#undef yylex
#define yylex driver.lexer->lex
}

%token AIM_SOLVER ANALYTIC_DERIVATION AR AUTOCORR
%token BAYESIAN_IRF BETA_PDF BLOCK USE_CALIBRATION
%token BVAR_DENSITY BVAR_FORECAST
%token BVAR_PRIOR_DECAY BVAR_PRIOR_FLAT BVAR_PRIOR_LAMBDA
%token BVAR_PRIOR_MU BVAR_PRIOR_OMEGA BVAR_PRIOR_TAU BVAR_PRIOR_TRAIN
%token BVAR_REPLIC BYTECODE ALL_VALUES_REQUIRED
%token CALIB_SMOOTHER CHANGE_TYPE CHECK CONDITIONAL_FORECAST CONDITIONAL_FORECAST_PATHS CONF_SIG CONSTANT CONTROLLED_VAREXO CORR COVAR CUTOFF CYCLE_REDUCTION LOGARITHMIC_REDUCTION
%token CONSIDER_ALL_ENDOGENOUS CONSIDER_ONLY_OBSERVED
%token DATAFILE FILE DOUBLING DR_CYCLE_REDUCTION_TOL DR_LOGARITHMIC_REDUCTION_TOL DR_LOGARITHMIC_REDUCTION_MAXITER DR_ALGO DROP DSAMPLE DYNASAVE DYNATYPE CALIBRATION DIFFERENTIATE_FORWARD_VARS
%token END ENDVAL EQUAL ESTIMATION ESTIMATED_PARAMS ESTIMATED_PARAMS_BOUNDS ESTIMATED_PARAMS_INIT EXTENDED_PATH ENDOGENOUS_PRIOR
%token FILENAME FILTER_STEP_AHEAD FILTERED_VARS FIRST_OBS LAST_OBS SET_TIME
%token <string_val> FLOAT_NUMBER DATES
%token DEFAULT FIXED_POINT
%token FORECAST K_ORDER_SOLVER INSTRUMENTS SHIFT MEAN STDEV VARIANCE MODE INTERVAL SHAPE DOMAINN
%token GAMMA_PDF GRAPH GRAPH_FORMAT CONDITIONAL_VARIANCE_DECOMPOSITION NOCHECK STD
%token HISTVAL HISTVAL_FILE HOMOTOPY_SETUP HOMOTOPY_MODE HOMOTOPY_STEPS HOMOTOPY_FORCE_CONTINUE HP_FILTER HP_NGRID HYBRID
%token IDENTIFICATION INF_CONSTANT INITVAL INITVAL_FILE BOUNDS JSCALE INIT INFILE INVARS
%token <string_val> INT_NUMBER
%token INV_GAMMA_PDF INV_GAMMA1_PDF INV_GAMMA2_PDF IRF IRF_SHOCKS IRF_PLOT_THRESHOLD IRF_CALIBRATION
%token KALMAN_ALGO KALMAN_TOL SUBSAMPLES OPTIONS TOLF
%token LAPLACE LIK_ALGO LIK_INIT LINEAR LOAD_IDENT_FILES LOAD_MH_FILE LOAD_PARAMS_AND_STEADY_STATE LOGLINEAR LOGDATA LYAPUNOV
%token LYAPUNOV_FIXED_POINT_TOL LYAPUNOV_DOUBLING_TOL LYAPUNOV_SQUARE_ROOT_SOLVER_TOL LOG_DEFLATOR LOG_TREND_VAR LOG_GROWTH_FACTOR MARKOWITZ MARGINAL_DENSITY MAX MAXIT
%token MFS MH_CONF_SIG MH_DROP MH_INIT_SCALE MH_JSCALE MH_MODE MH_NBLOCKS MH_REPLIC MH_RECOVER POSTERIOR_MAX_SUBSAMPLE_DRAWS MIN MINIMAL_SOLVING_PERIODS
%token MODE_CHECK MODE_CHECK_NEIGHBOURHOOD_SIZE MODE_CHECK_SYMMETRIC_PLOTS MODE_CHECK_NUMBER_OF_POINTS MODE_COMPUTE MODE_FILE MODEL MODEL_COMPARISON MODEL_INFO MSHOCKS ABS SIGN
%token MODEL_DIAGNOSTICS MODIFIEDHARMONICMEAN MOMENTS_VARENDO DIFFUSE_FILTER SUB_DRAWS TAPER_STEPS GEWEKE_INTERVAL MCMC_JUMPING_COVARIANCE MOMENT_CALIBRATION
%token <string_val> NAME
%token NAN_CONSTANT NO_STATIC NOBS NOCONSTANT NODISPLAY NOCORR NODIAGNOSTIC NOFUNCTIONS NO_HOMOTOPY
%token NOGRAPH NOMOMENTS NOPRINT NORMAL_PDF SAVE_DRAWS
%token OBSERVATION_TRENDS OPTIM OPTIM_WEIGHTS ORDER OSR OSR_PARAMS MAX_DIM_COVA_GROUP ADVANCED OUTFILE OUTVARS OVERWRITE
%token PARALLEL_LOCAL_FILES PARAMETERS PARAMETER_SET PARTIAL_INFORMATION PERFECT_FORESIGHT PERIODS PERIOD PLANNER_OBJECTIVE PLOT_CONDITIONAL_FORECAST PLOT_PRIORS PREFILTER PRESAMPLE
%token PERFECT_FORESIGHT_SETUP PERFECT_FORESIGHT_SOLVER
%token PRINT PRIOR_MC PRIOR_TRUNC PRIOR_MODE PRIOR_MEAN POSTERIOR_MODE POSTERIOR_MEAN POSTERIOR_MEDIAN PRUNING
%token <string_val> QUOTED_STRING
%token QZ_CRITERIUM QZ_ZERO_THRESHOLD FULL DSGE_VAR DSGE_VARLAG DSGE_PRIOR_WEIGHT TRUNCATE
%token RELATIVE_IRF REPLIC SIMUL_REPLIC RPLOT SAVE_PARAMS_AND_STEADY_STATE PARAMETER_UNCERTAINTY
%token SHOCKS SHOCK_DECOMPOSITION SIGMA_E SIMUL SIMUL_ALGO SIMUL_SEED ENDOGENOUS_TERMINAL_PERIOD
%token SMOOTHER SMOOTHER2HISTVAL SQUARE_ROOT_SOLVER STACK_SOLVE_ALGO STEADY_STATE_MODEL SOLVE_ALGO SOLVER_PERIODS
%token STDERR STEADY STOCH_SIMUL SURPRISE SYLVESTER SYLVESTER_FIXED_POINT_TOL REGIMES REGIME
%token TEX RAMSEY_MODEL RAMSEY_POLICY PLANNER_DISCOUNT DISCRETIONARY_POLICY DISCRETIONARY_TOL
%token <string_val> TEX_NAME
%token UNIFORM_PDF UNIT_ROOT_VARS USE_DLL USEAUTOCORR GSA_SAMPLE_FILE USE_UNIVARIATE_FILTERS_IF_SINGULARITY_IS_DETECTED
%token VALUES VAR VAREXO VAREXO_DET VAROBS PREDETERMINED_VARIABLES
%token WRITE_LATEX_DYNAMIC_MODEL WRITE_LATEX_STATIC_MODEL
%token XLS_SHEET XLS_RANGE LONG_NAME
%left COMMA
%left EQUAL_EQUAL EXCLAMATION_EQUAL
%left LESS GREATER LESS_EQUAL GREATER_EQUAL
%left PLUS MINUS
%left TIMES DIVIDE
%left UMINUS UPLUS
%nonassoc POWER
%token EXP LOG LN LOG10 SIN COS TAN ASIN ACOS ATAN SINH COSH TANH ERF
%token ASINH ACOSH ATANH SQRT NORMCDF NORMPDF STEADY_STATE EXPECTATION
/* GSA analysis */
%token DYNARE_SENSITIVITY MORRIS STAB REDFORM PPRIOR PRIOR_RANGE PPOST ILPTAU MORRIS_NLIV
%token MORRIS_NTRA NSAM LOAD_REDFORM LOAD_RMSE LOAD_STAB ALPHA2_STAB KSSTAT LOGTRANS_REDFORM THRESHOLD_REDFORM
%token KSSTAT_REDFORM ALPHA2_REDFORM NAMENDO NAMLAGENDO NAMEXO RMSE LIK_ONLY VAR_RMSE PFILT_RMSE ISTART_RMSE
%token ALPHA_RMSE ALPHA2_RMSE
/* end of GSA analysis*/
%token FREQ INITIAL_YEAR INITIAL_SUBPERIOD FINAL_YEAR FINAL_SUBPERIOD DATA VLIST LOG_VAR PERCENT_VAR
%token VLISTLOG VLISTPER
%token RESTRICTION RESTRICTION_FNAME CROSS_RESTRICTIONS NLAGS CONTEMP_REDUCED_FORM REAL_PSEUDO_FORECAST
%token DUMMY_OBS NSTATES INDXSCALESSTATES NO_BAYESIAN_PRIOR SPECIFICATION SIMS_ZHA
%token <string_val> ALPHA BETA ABAND NINV CMS NCMS CNUM GAMMA INV_GAMMA INV_GAMMA1 INV_GAMMA2 NORMAL UNIFORM EPS PDF FIG DR NONE PRIOR PRIOR_VARIANCE HESSIAN IDENTITY_MATRIX DIRICHLET
%token GSIG2_LMDM Q_DIAG FLAT_PRIOR NCSK NSTD
%token INDXPARR INDXOVR INDXAP APBAND INDXIMF IMFBAND INDXFORE FOREBAND INDXGFOREHAT INDXGIMFHAT
%token INDXESTIMA INDXGDLS EQ_MS FILTER_COVARIANCE FILTER_DECOMPOSITION
%token EQ_CMS TLINDX TLNUMBER BANACT RESTRICTIONS
%token OUTPUT_FILE_TAG DRAWS_NBR_BURN_IN_1 DRAWS_NBR_BURN_IN_2 HORIZON
%token SBVAR TREND_VAR DEFLATOR GROWTH_FACTOR MS_IRF MS_VARIANCE_DECOMPOSITION
%token MS_ESTIMATION MS_SIMULATION MS_COMPUTE_MDD MS_COMPUTE_PROBABILITIES MS_FORECAST
%token SVAR_IDENTIFICATION EQUATION EXCLUSION LAG UPPER_CHOLESKY LOWER_CHOLESKY MONTHLY QUARTERLY
%token MARKOV_SWITCHING CHAIN DURATION NUMBER_OF_REGIMES NUMBER_OF_LAGS
%token SVAR COEFF COEFFICIENTS VARIANCES CONSTANTS EQUATIONS
%token EXTERNAL_FUNCTION EXT_FUNC_NAME EXT_FUNC_NARGS FIRST_DERIV_PROVIDED SECOND_DERIV_PROVIDED
%token SELECTED_VARIABLES_ONLY COVA_COMPUTE SIMULATION_FILE_TAG FILE_TAG
%token NO_ERROR_BANDS ERROR_BAND_PERCENTILES SHOCKS_PER_PARAMETER NO_CREATE_INIT
%token SHOCK_DRAWS FREE_PARAMETERS MEDIAN DATA_OBS_NBR NEIGHBORHOOD_WIDTH PVALUE_KS PVALUE_CORR
%token FILTERED_PROBABILITIES REAL_TIME_SMOOTHED
%token PROPOSAL_TYPE PROPOSAL_UPPER_BOUND PROPOSAL_LOWER_BOUND PROPOSAL_DRAWS USE_MEAN_CENTER
%token ADAPTIVE_MH_DRAWS THINNING_FACTOR COEFFICIENTS_PRIOR_HYPERPARAMETERS
%token CONVERGENCE_STARTING_VALUE CONVERGENCE_ENDING_VALUE CONVERGENCE_INCREMENT_VALUE
%token MAX_ITERATIONS_STARTING_VALUE MAX_ITERATIONS_INCREMENT_VALUE MAX_BLOCK_ITERATIONS
%token MAX_REPEATED_OPTIMIZATION_RUNS FUNCTION_CONVERGENCE_CRITERION
%token PARAMETER_CONVERGENCE_CRITERION NUMBER_OF_LARGE_PERTURBATIONS NUMBER_OF_SMALL_PERTURBATIONS
%token NUMBER_OF_POSTERIOR_DRAWS_AFTER_PERTURBATION MAX_NUMBER_OF_STAGES
%token RANDOM_FUNCTION_CONVERGENCE_CRITERION RANDOM_PARAMETER_CONVERGENCE_CRITERION

%type <node_val> expression expression_or_empty
%type <node_val> equation hand_side
%type <string_val> non_negative_number signed_number signed_integer date_str
%type <string_val> filename symbol vec_of_vec_value vec_value_list date_expr
%type <string_val> vec_value_1 vec_value signed_inf signed_number_w_inf
%type <string_val> range vec_value_w_inf vec_value_1_w_inf named_var
%type <string_val> integer_range signed_integer_range
%type <symbol_type_val> change_type_arg
%type <vector_string_val> change_type_var_list subsamples_eq_opt prior_eq_opt options_eq_opt calibration_range
%type <vector_int_val> vec_int_elem vec_int_1 vec_int vec_int_number
%type <prior_distributions_val> prior_pdf prior_distribution
%%

%start statement_list;

statement_list : statement
               | statement_list statement
               ;

statement : parameters
          | var
          | varexo
          | varexo_det
          | predetermined_variables
          | change_type
          | periods
          | model
          | initval
          | initval_file
          | endval
          | histval
          | init_param
          | shocks
          | mshocks
          | sigma_e
          | steady
          | check
          | simul
          | stoch_simul
          | estimation
          | estimated_params
          | estimated_params_bounds
          | estimated_params_init
          | set_time
          | data
          | prior
          | prior_eq
          | subsamples
          | subsamples_eq
          | options
          | options_eq
          | varobs
          | observation_trends
          | unit_root_vars
          | dsample
          | rplot
          | optim_weights
          | osr_params
          | osr
          | dynatype
          | dynasave
          | model_comparison
          | model_info
          | planner_objective
          | ramsey_model
          | ramsey_policy
          | discretionary_policy
          | bvar_density
          | bvar_forecast
          | sbvar
          | dynare_sensitivity
          | homotopy_setup
          | forecast
          | load_params_and_steady_state
          | save_params_and_steady_state
          | identification
          | write_latex_dynamic_model
          | write_latex_static_model
          | shock_decomposition
          | conditional_forecast
          | conditional_forecast_paths
          | plot_conditional_forecast
          | svar_identification
          | markov_switching
          | svar
          | external_function
          | steady_state_model
          | trend_var
          | log_trend_var
          | ms_estimation
          | ms_simulation
          | ms_compute_mdd
          | ms_compute_probabilities
          | ms_forecast
          | ms_irf
          | ms_variance_decomposition
          | calib_smoother
          | extended_path
          | model_diagnostics
          | moment_calibration
          | irf_calibration
          | smoother2histval
          | histval_file
          | perfect_foresight_setup
          | perfect_foresight_solver
          ;

dsample : DSAMPLE INT_NUMBER ';'
          { driver.dsample($2); }
        | DSAMPLE INT_NUMBER INT_NUMBER ';'
          { driver.dsample($2, $3); }
        ;

rplot : RPLOT symbol_list ';' { driver.rplot(); };

trend_var : TREND_VAR '(' GROWTH_FACTOR EQUAL { driver.begin_trend(); } hand_side ')' trend_var_list ';'
            { driver.end_trend_var($6); }
          ;

trend_var_list : trend_var_list symbol
                 { driver.declare_trend_var(false, $2); }
               | trend_var_list COMMA symbol
                 { driver.declare_trend_var(false, $3); }
               | symbol
                 { driver.declare_trend_var(false, $1); }
               | trend_var_list symbol TEX_NAME
                 { driver.declare_trend_var(false, $2, $3); }
               | trend_var_list COMMA symbol TEX_NAME
                 { driver.declare_trend_var(false, $3, $4); }
               | symbol TEX_NAME
                 { driver.declare_trend_var(false, $1, $2); }
               ;

log_trend_var : LOG_TREND_VAR '(' LOG_GROWTH_FACTOR EQUAL { driver.begin_trend(); } hand_side ')' log_trend_var_list ';'
                { driver.end_trend_var($6); }
              ;

log_trend_var_list : log_trend_var_list symbol
                     { driver.declare_trend_var(true, $2); }
                   | log_trend_var_list COMMA symbol
                     { driver.declare_trend_var(true, $3); }
                   | symbol
                     { driver.declare_trend_var(true, $1); }
                   | log_trend_var_list symbol TEX_NAME
                     { driver.declare_trend_var(true, $2, $3); }
                   | log_trend_var_list COMMA symbol TEX_NAME
                     { driver.declare_trend_var(true, $3, $4); }
                   | symbol TEX_NAME
                     { driver.declare_trend_var(true, $1, $2); }
                   ;

var : VAR var_list ';'
    | VAR '(' DEFLATOR EQUAL { driver.begin_trend(); } hand_side ')' nonstationary_var_list ';'
      { driver.end_nonstationary_var(false, $6); }
    | VAR '(' LOG_DEFLATOR EQUAL { driver.begin_trend(); } hand_side ')' nonstationary_var_list ';'
      { driver.end_nonstationary_var(true, $6); }
    ;

nonstationary_var_list : nonstationary_var_list symbol
                         { driver.declare_nonstationary_var($2); }
                       | nonstationary_var_list COMMA symbol
                         { driver.declare_nonstationary_var($3); }
                       | symbol
                         { driver.declare_nonstationary_var($1); }
                       |  nonstationary_var_list symbol named_var
                         { driver.declare_nonstationary_var($2, NULL, $3); }
                       | nonstationary_var_list COMMA symbol named_var
                         { driver.declare_nonstationary_var($3, NULL, $4); }
                       | symbol named_var
                         { driver.declare_nonstationary_var($1, NULL, $2); }
                       | nonstationary_var_list symbol TEX_NAME
                         { driver.declare_nonstationary_var($2, $3); }
                       | nonstationary_var_list COMMA symbol TEX_NAME
                         { driver.declare_nonstationary_var($3, $4); }
                       | symbol TEX_NAME
                         { driver.declare_nonstationary_var($1, $2); }
                       | nonstationary_var_list symbol TEX_NAME named_var
                         { driver.declare_nonstationary_var($2, $3, $4); }
                       | nonstationary_var_list COMMA symbol TEX_NAME named_var
                         { driver.declare_nonstationary_var($3, $4, $5); }
                       | symbol TEX_NAME named_var
                         { driver.declare_nonstationary_var($1, $2, $3); }
                       ;

varexo : VAREXO varexo_list ';';

varexo_det : VAREXO_DET varexo_det_list ';';

predetermined_variables : PREDETERMINED_VARIABLES predetermined_variables_list ';';

parameters : PARAMETERS parameter_list ';';

named_var : '(' LONG_NAME EQUAL QUOTED_STRING ')'
            { $$ = $4; }
          ;

var_list : var_list symbol
           { driver.declare_endogenous($2); }
         | var_list COMMA symbol
           { driver.declare_endogenous($3); }
         | symbol
           { driver.declare_endogenous($1); }
         | var_list symbol named_var
           { driver.declare_endogenous($2, NULL, $3); }
         | var_list COMMA symbol named_var
           { driver.declare_endogenous($3, NULL, $4); }
         | symbol named_var
           { driver.declare_endogenous($1, NULL, $2); }
         | var_list symbol TEX_NAME
           { driver.declare_endogenous($2, $3); }
         | var_list COMMA symbol TEX_NAME
           { driver.declare_endogenous($3, $4); }
         | symbol TEX_NAME
           { driver.declare_endogenous($1, $2); }
         | var_list symbol TEX_NAME named_var
           { driver.declare_endogenous($2, $3, $4); }
         | var_list COMMA symbol TEX_NAME named_var
           { driver.declare_endogenous($3, $4, $5); }
         | symbol TEX_NAME named_var
           { driver.declare_endogenous($1, $2, $3); }
         ;

varexo_list : varexo_list symbol
              { driver.declare_exogenous($2); }
            | varexo_list COMMA symbol
              { driver.declare_exogenous($3); }
            | symbol
              { driver.declare_exogenous($1); }
            | varexo_list symbol named_var
              { driver.declare_exogenous($2, NULL, $3); }
            | varexo_list COMMA symbol named_var
              { driver.declare_exogenous($3, NULL, $4); }
            | symbol named_var
              { driver.declare_exogenous($1, NULL, $2); }
            | varexo_list symbol TEX_NAME
              { driver.declare_exogenous($2, $3); }
            | varexo_list COMMA symbol TEX_NAME
              { driver.declare_exogenous($3, $4); }
            | symbol TEX_NAME
              { driver.declare_exogenous($1, $2); }
            | varexo_list symbol TEX_NAME named_var
              { driver.declare_exogenous($2, $3, $4); }
            | varexo_list COMMA symbol TEX_NAME named_var
              { driver.declare_exogenous($3, $4, $5); }
            | symbol TEX_NAME named_var
              { driver.declare_exogenous($1, $2, $3); }
            ;

varexo_det_list : varexo_det_list symbol
                  { driver.declare_exogenous_det($2); }
                | varexo_det_list COMMA symbol
                  { driver.declare_exogenous_det($3); }
                | symbol
                  { driver.declare_exogenous_det($1); }
                | varexo_det_list symbol named_var
                  { driver.declare_exogenous_det($2, NULL, $3); }
                | varexo_det_list COMMA symbol named_var
                  { driver.declare_exogenous_det($3, NULL, $4); }
                | symbol named_var
                  { driver.declare_exogenous_det($1, NULL, $2); }
                | varexo_det_list symbol TEX_NAME
                  { driver.declare_exogenous_det($2, $3); }
                | varexo_det_list COMMA symbol TEX_NAME
                  { driver.declare_exogenous_det($3, $4); }
                | symbol TEX_NAME
                   { driver.declare_exogenous_det($1, $2); }
                | varexo_det_list symbol TEX_NAME named_var
                  { driver.declare_exogenous_det($2, $3, $4); }
                | varexo_det_list COMMA symbol TEX_NAME named_var
                  { driver.declare_exogenous_det($3, $4, $5); }
                | symbol TEX_NAME named_var
                   { driver.declare_exogenous_det($1, $2, $3); }
                ;

parameter_list : parameter_list symbol
                 { driver.declare_parameter($2); }
               | parameter_list COMMA symbol
                 { driver.declare_parameter($3); }
               | symbol
                 { driver.declare_parameter($1); }
               | parameter_list symbol named_var
                 { driver.declare_parameter($2, NULL, $3); }
               | parameter_list COMMA symbol named_var
                 { driver.declare_parameter($3, NULL, $4); }
               | symbol named_var
                 { driver.declare_parameter($1, NULL, $2); }
               | parameter_list symbol TEX_NAME
                 { driver.declare_parameter($2, $3); }
               | parameter_list COMMA symbol TEX_NAME
                 { driver.declare_parameter($3, $4); }
               | symbol TEX_NAME
                 { driver.declare_parameter($1, $2); }
               | parameter_list symbol TEX_NAME named_var
                 { driver.declare_parameter($2, $3, $4); }
               | parameter_list COMMA symbol TEX_NAME named_var
                 { driver.declare_parameter($3, $4, $5); }
               | symbol TEX_NAME named_var
                 { driver.declare_parameter($1, $2, $3); }
               ;

predetermined_variables_list : predetermined_variables_list symbol
                               { driver.add_predetermined_variable($2); }
                             | predetermined_variables_list COMMA symbol
                               { driver.add_predetermined_variable($3); }
                             | symbol
                               { driver.add_predetermined_variable($1); }
                             ;

change_type : CHANGE_TYPE '(' change_type_arg ')' change_type_var_list ';'
              { driver.change_type($3, $5); }
            ;

change_type_arg : PARAMETERS
                  { $$ = eParameter; }
                | VAR
                  { $$ = eEndogenous; }
                | VAREXO
                  { $$ = eExogenous; }
                | VAREXO_DET
                  { $$ = eExogenousDet; }
                ;

change_type_var_list : symbol
                       { $$ = new vector<string *>(); $$->push_back($1); }
                     | change_type_var_list symbol
                       { $$ = $1; $1->push_back($2); }
                     | change_type_var_list COMMA symbol
                       { $$ = $1; $1->push_back($3); }
                     ;

periods : PERIODS INT_NUMBER ';'
          { driver.periods($2); }
        | PERIODS EQUAL INT_NUMBER ';'
          { driver.periods($3); }
        ;

init_param : symbol EQUAL expression ';' { driver.init_param($1, $3); };

expression : '(' expression ')'
             { $$ = $2;}
           | symbol
             { $$ = driver.add_expression_variable($1); }
           | non_negative_number
             { $$ = driver.add_non_negative_constant($1); }
           | expression PLUS expression
             { $$ = driver.add_plus($1, $3); }
           | expression MINUS expression
             { $$ = driver.add_minus($1, $3); }
           | expression DIVIDE expression
             { $$ = driver.add_divide($1, $3); }
           | expression TIMES expression
             { $$ = driver.add_times($1, $3); }
           | expression POWER expression
             { $$ = driver.add_power($1, $3); }
           | expression LESS expression
             { $$ = driver.add_less($1, $3); }
           | expression GREATER expression
             { $$ = driver.add_greater($1, $3); }
           | expression LESS_EQUAL expression
             { $$ = driver.add_less_equal($1, $3); }
           | expression GREATER_EQUAL expression
             { $$ = driver.add_greater_equal($1, $3); }
           | expression EQUAL_EQUAL expression
             { $$ = driver.add_equal_equal($1, $3); }
           | expression EXCLAMATION_EQUAL expression
             { $$ = driver.add_different($1, $3); }
           | MINUS expression %prec UMINUS
             { $$ = driver.add_uminus($2); }
           | PLUS expression %prec UPLUS
             { $$ = $2; }
           | EXP '(' expression ')'
             { $$ = driver.add_exp($3); }
           | LOG '(' expression ')'
             { $$ = driver.add_log($3); }
           | LN '(' expression ')'
             { $$ = driver.add_log($3); }
           | LOG10 '(' expression ')'
             { $$ = driver.add_log10($3); }
           | SIN '(' expression ')'
             { $$ = driver.add_sin($3); }
           | COS '(' expression ')'
             { $$ = driver.add_cos($3); }
           | TAN '(' expression ')'
             { $$ = driver.add_tan($3); }
           | ASIN '(' expression ')'
             { $$ = driver.add_asin($3); }
           | ACOS '(' expression ')'
             { $$ = driver.add_acos($3); }
           | ATAN '(' expression ')'
             { $$ = driver.add_atan($3); }
           | SQRT '(' expression ')'
             { $$ = driver.add_sqrt($3); }
           | ABS '(' expression ')'
             { $$ = driver.add_abs($3); }
           | SIGN '(' expression ')'
             { $$ = driver.add_sign($3); }
           | MAX '(' expression COMMA expression ')'
             { $$ = driver.add_max($3, $5); }
           | MIN '(' expression COMMA expression ')'
             { $$ = driver.add_min($3, $5); }
           | symbol { driver.push_external_function_arg_vector_onto_stack(); } '(' comma_expression ')'
             { $$ = driver.add_model_var_or_external_function($1, false); }
           | NORMCDF '(' expression COMMA expression COMMA expression ')'
             { $$ = driver.add_normcdf($3, $5, $7); }
           | NORMCDF '(' expression ')'
             { $$ = driver.add_normcdf($3); }
           | NORMPDF '(' expression COMMA expression COMMA expression ')'
             { $$ = driver.add_normpdf($3, $5, $7); }
           | NORMPDF '(' expression ')'
             { $$ = driver.add_normpdf($3); }
           | ERF '(' expression ')'
             { $$ = driver.add_erf($3); }
           | NAN_CONSTANT
             { $$ = driver.add_nan_constant(); }
           | INF_CONSTANT
             { $$ = driver.add_inf_constant(); }
           ;

comma_expression : expression
                   { driver.add_external_function_arg($1); }
                 | comma_expression COMMA expression
                   { driver.add_external_function_arg($3); }
                 ;

expression_or_empty : {$$ = driver.add_nan_constant();}
                    | expression
	            ;

initval : INITVAL ';' initval_list END ';'
          { driver.end_initval(false); }
        | INITVAL '(' ALL_VALUES_REQUIRED ')' ';' initval_list END ';'
          { driver.end_initval(true); }
        ;

initval_file : INITVAL_FILE '(' FILENAME EQUAL filename ')' ';'
               { driver.initval_file($5); }
             ;

endval : ENDVAL ';' initval_list END ';'
         { driver.end_endval(false); }
       | ENDVAL '(' ALL_VALUES_REQUIRED ')' ';' initval_list END ';'
         { driver.end_endval(true); }
       ;

initval_list : initval_list initval_elem
             | initval_elem
             ;

initval_elem : symbol EQUAL expression ';' { driver.init_val($1, $3); };

histval : HISTVAL ';' histval_list END ';' { driver.end_histval(); };

histval_list : histval_list histval_elem
             | histval_elem
             ;

histval_elem : symbol '(' signed_integer ')' EQUAL expression ';' { driver.hist_val($1, $3, $6); };

histval_file : HISTVAL_FILE '(' FILENAME EQUAL filename ')' ';'
               { driver.histval_file($5); }
             ;

model_options : BLOCK { driver.block(); }
              | o_cutoff
							| o_mfs
              | BYTECODE { driver.byte_code(); }
              | USE_DLL { driver.use_dll(); }
              | NO_STATIC { driver.no_static();}
              | DIFFERENTIATE_FORWARD_VARS { driver.differentiate_forward_vars_all(); }
              | DIFFERENTIATE_FORWARD_VARS EQUAL '(' symbol_list ')' { driver.differentiate_forward_vars_some(); }
              | o_linear
              | PARALLEL_LOCAL_FILES EQUAL '(' parallel_local_filename_list ')'
              ;

model_options_list : model_options_list COMMA model_options
                   | model_options
                   ;

model : MODEL ';' { driver.begin_model(); }
        equation_list END ';' { driver.reset_data_tree(); }
      | MODEL '(' model_options_list ')' ';' { driver.begin_model(); }
        equation_list END ';' { driver.reset_data_tree(); }
      ;

equation_list : equation_list equation
              | equation_list pound_expression
              | equation
              | pound_expression
              ;

equation : hand_side EQUAL hand_side ';'
           { $$ = driver.add_model_equal($1, $3); }
         | hand_side ';'
           { $$ = driver.add_model_equal_with_zero_rhs($1); }
         | '[' tags_list ']' hand_side EQUAL hand_side ';'
           { $$ = driver.add_model_equal($4, $6); }
         | '[' tags_list ']' hand_side ';'
           { $$ = driver.add_model_equal_with_zero_rhs($4); }
         ;

tags_list : tags_list COMMA tag_pair
          | tag_pair
          ;

tag_pair : NAME EQUAL QUOTED_STRING
           { driver.add_equation_tags($1, $3); }
         | NAME
           { driver.add_equation_tags($1, new string()); }
         ;

hand_side : '(' hand_side ')'
            { $$ = $2;}
          | symbol
            { $$ = driver.add_model_variable($1); }
          | non_negative_number
            { $$ = driver.add_non_negative_constant($1); }
          | hand_side PLUS hand_side
            { $$ = driver.add_plus($1, $3); }
          | hand_side MINUS hand_side
            { $$ = driver.add_minus($1, $3); }
          | hand_side DIVIDE hand_side
            { $$ = driver.add_divide($1, $3); }
          | hand_side TIMES hand_side
            { $$ = driver.add_times($1, $3); }
          | hand_side LESS hand_side
            { $$ = driver.add_less($1, $3); }
          | hand_side GREATER hand_side
            { $$ = driver.add_greater($1, $3); }
          | hand_side LESS_EQUAL hand_side
            { $$ = driver.add_less_equal($1, $3); }
          | hand_side GREATER_EQUAL hand_side
            { $$ = driver.add_greater_equal($1, $3); }
          | hand_side EQUAL_EQUAL hand_side
            { $$ = driver.add_equal_equal($1, $3); }
          | hand_side EXCLAMATION_EQUAL hand_side
            { $$ = driver.add_different($1, $3); }
          | hand_side POWER hand_side
            { $$ = driver.add_power($1, $3); }
          | EXPECTATION '(' signed_integer ')''(' hand_side ')'
	    { $$ = driver.add_expectation($3, $6); }
          | MINUS hand_side %prec UMINUS
            { $$ = driver.add_uminus($2); }
          | PLUS hand_side
            { $$ = $2; }
          | EXP '(' hand_side ')'
            { $$ = driver.add_exp($3); }
          | LOG '(' hand_side ')'
            { $$ = driver.add_log($3); }
          | LN '(' hand_side ')'
            { $$ = driver.add_log($3); }
          | LOG10 '(' hand_side ')'
            { $$ = driver.add_log10($3); }
          | SIN '(' hand_side ')'
            { $$ = driver.add_sin($3); }
          | COS '(' hand_side ')'
            { $$ = driver.add_cos($3); }
          | TAN '(' hand_side ')'
            { $$ = driver.add_tan($3); }
          | ASIN '(' hand_side ')'
            { $$ = driver.add_asin($3); }
          | ACOS '(' hand_side ')'
            { $$ = driver.add_acos($3); }
          | ATAN '(' hand_side ')'
            { $$ = driver.add_atan($3); }
          | SQRT '(' hand_side ')'
            { $$ = driver.add_sqrt($3); }
          | ABS '(' hand_side ')'
            { $$ = driver.add_abs($3); }
          | SIGN '(' hand_side ')'
            { $$ = driver.add_sign($3); }
          | MAX '(' hand_side COMMA hand_side ')'
            { $$ = driver.add_max($3, $5); }
          | MIN '(' hand_side COMMA hand_side ')'
            { $$ = driver.add_min($3, $5); }
          | symbol { driver.push_external_function_arg_vector_onto_stack(); } '(' comma_hand_side ')'
            { $$ = driver.add_model_var_or_external_function($1, true); }
          | NORMCDF '(' hand_side COMMA hand_side COMMA hand_side ')'
            { $$ = driver.add_normcdf($3, $5, $7); }
          | NORMCDF '(' hand_side ')'
            { $$ = driver.add_normcdf($3); }
          | NORMPDF '(' hand_side COMMA hand_side COMMA hand_side ')'
            { $$ = driver.add_normpdf($3, $5, $7); }
          | NORMPDF '(' hand_side ')'
            { $$ = driver.add_normpdf($3); }
          | ERF '(' hand_side ')'
            { $$ = driver.add_erf($3); }
          | STEADY_STATE '(' hand_side ')'
            { $$ = driver.add_steady_state($3); }
          ;

comma_hand_side : hand_side
                  { driver.add_external_function_arg($1); }
                | comma_hand_side COMMA hand_side
                  { driver.add_external_function_arg($3); }
                ;

pound_expression: '#' symbol EQUAL hand_side ';'
                  { driver.declare_and_init_model_local_variable($2, $4); };

shocks : SHOCKS ';' shock_list END ';' { driver.end_shocks(false); }
       | SHOCKS '(' OVERWRITE ')' ';' shock_list END ';' { driver.end_shocks(true); }
       ;

shock_list : shock_list shock_elem
           | shock_elem
           ;

shock_elem : det_shock_elem
           | VAR symbol ';' STDERR expression ';'
             { driver.add_stderr_shock($2, $5); }
           | VAR symbol EQUAL expression ';'
             { driver.add_var_shock($2, $4); }
           | VAR symbol COMMA symbol EQUAL expression ';'
             { driver.add_covar_shock($2, $4, $6); }
           | CORR symbol COMMA symbol EQUAL expression ';'
             { driver.add_correl_shock($2, $4, $6); }
           ;

det_shock_elem : VAR symbol ';' PERIODS period_list ';' VALUES value_list ';'
                 { driver.add_det_shock($2, false); }
               ;

svar_identification : SVAR_IDENTIFICATION {driver.begin_svar_identification();} ';' svar_identification_list END ';'
                      { driver.end_svar_identification(); }
                    ;

svar_identification_list : svar_identification_list svar_identification_elem
                         | svar_identification_elem
                         ;

svar_identification_elem : EXCLUSION LAG INT_NUMBER ';' svar_equation_list
                           { driver.combine_lag_and_restriction($3); }
                         | EXCLUSION CONSTANTS ';'
                           { driver.add_constants_exclusion(); }
                         | RESTRICTION EQUATION INT_NUMBER COMMA
			 { driver.add_restriction_equation_nbr($3);}
                           restriction_expression EQUAL
                                {driver.add_restriction_equal();}
                           restriction_expression ';'
                         | UPPER_CHOLESKY ';'
                           { driver.add_upper_cholesky(); }
                         | LOWER_CHOLESKY ';'
                           { driver.add_lower_cholesky(); }
                         ;

svar_equation_list : svar_equation_list EQUATION INT_NUMBER COMMA svar_var_list ';'
                     { driver.add_restriction_in_equation($3); }
                   | EQUATION INT_NUMBER COMMA svar_var_list ';'
                     { driver.add_restriction_in_equation($2); }
                   ;

svar_var_list : svar_var_list COMMA symbol
                { driver.add_in_svar_restriction_symbols($3); }
              | symbol
                { driver.add_in_svar_restriction_symbols($1); }
              ;

restriction_expression : expression {driver.check_restriction_expression_constant($1);}
                       | restriction_expression_1
                       ;

restriction_expression_1 : restriction_elem_expression
                         | restriction_expression_1 restriction_elem_expression
                         ;   

restriction_elem_expression : COEFF '(' symbol COMMA INT_NUMBER ')'
                                 { driver.add_positive_restriction_element($3,$5);}
                            | PLUS COEFF '(' symbol COMMA INT_NUMBER ')'
		                 { driver.add_positive_restriction_element($4,$6);}
                            | MINUS COEFF '(' symbol COMMA INT_NUMBER ')'
		                 { driver.add_negative_restriction_element($4,$6);}
                            | expression TIMES COEFF '(' symbol COMMA INT_NUMBER ')'
                                 { driver.add_positive_restriction_element($1,$5,$7);}
                            ;

markov_switching : MARKOV_SWITCHING '(' ms_options_list ')' ';'
                   { driver.markov_switching(); }
                 ;

ms_options_list : ms_options_list COMMA ms_options
                | ms_options
                ;

ms_options : o_chain
           | o_duration
           | o_restrictions
           | o_number_of_regimes
           | o_number_of_lags
           | o_parameters
           ;

svar : SVAR '(' svar_options_list ')' ';'
       { driver.svar(); }
     ;

svar_options_list : svar_options_list COMMA svar_options
                  | svar_options
                  ;

svar_options : o_coefficients
             | o_variances
             | o_equations
             | o_chain
             ;

mshocks : MSHOCKS ';' mshock_list END ';' { driver.end_mshocks(false); }
        | MSHOCKS '(' OVERWRITE ')' ';' mshock_list END ';' { driver.end_mshocks(true); }
        ;

mshock_list : mshock_list det_shock_elem
            | det_shock_elem
            ;

period_list : period_list COMMA INT_NUMBER
              { driver.add_period($3); }
            | period_list INT_NUMBER
              { driver.add_period($2); }
            | period_list COMMA INT_NUMBER ':' INT_NUMBER
              { driver.add_period($3, $5); }
            | period_list INT_NUMBER ':' INT_NUMBER
              { driver.add_period($2, $4); }
            | INT_NUMBER ':' INT_NUMBER
              { driver.add_period($1, $3); }
            | INT_NUMBER
              { driver.add_period($1); }
            ;

sigma_e : SIGMA_E EQUAL '[' triangular_matrix ']' ';' { driver.do_sigma_e(); };

value_list : value_list COMMA '(' expression ')'
             { driver.add_value($4); }
           | value_list '(' expression ')'
             { driver.add_value($3); }
           | '(' expression ')'
             { driver.add_value($2); }
           | value_list COMMA signed_number
             { driver.add_value($3); }
           | value_list signed_number
             { driver.add_value($2); }
           | signed_number
             { driver.add_value($1); }
           ;

triangular_matrix : triangular_matrix ';' triangular_row
                    { driver.end_of_row(); }
                  | triangular_row
                    { driver.end_of_row(); }
                  ;

triangular_row : triangular_row COMMA '(' expression ')'
                 { driver.add_to_row($4); }
               | triangular_row COMMA signed_number
                 { driver.add_to_row_const($3); }
               | triangular_row '(' expression ')'
                 { driver.add_to_row($3); }
               | triangular_row signed_number
                 { driver.add_to_row_const($2); }
               | '(' expression ')'
                 { driver.add_to_row($2); }
               | signed_number
                 { driver.add_to_row_const($1); }
               ;

steady : STEADY ';'
         { driver.steady(); }
       | STEADY '(' steady_options_list ')' ';'
         { driver.steady(); }
       ;

steady_options_list : steady_options_list COMMA steady_options
                    | steady_options
                    ;

steady_options : o_solve_algo
               | o_homotopy_mode
               | o_homotopy_steps
               | o_homotopy_force_continue
               | o_markowitz
               | o_steady_maxit
               | o_nocheck
               ;

check : CHECK ';'
        { driver.check(); }
      | CHECK '(' check_options_list ')' ';'
        { driver.check(); }
      ;

check_options_list : check_options_list COMMA check_options
                   | check_options
                   ;

check_options : steady_options
	      | o_qz_zero_threshold
	      ;

model_info : MODEL_INFO ';'
             { driver.model_info(); }
            | MODEL_INFO '(' model_info_options_list ')' ';'
              { driver.model_info(); }
           ;

model_info_options_list : model_info_options_list COMMA model_info_options
                   | model_info_options
                   ;
model_info_options :

perfect_foresight_setup : PERFECT_FORESIGHT_SETUP ';'
                          { driver.perfect_foresight_setup(); }
                        | PERFECT_FORESIGHT_SETUP '(' perfect_foresight_setup_options_list ')' ';'
                          { driver.perfect_foresight_setup(); }
                        ;

perfect_foresight_setup_options_list : perfect_foresight_setup_options_list COMMA perfect_foresight_setup_options
                                     | perfect_foresight_setup_options
                                     ;

perfect_foresight_setup_options : o_periods
                                | o_datafile
                                ;

perfect_foresight_solver : PERFECT_FORESIGHT_SOLVER ';'
                          { driver.perfect_foresight_solver(); }
                         | PERFECT_FORESIGHT_SOLVER '(' perfect_foresight_solver_options_list ')' ';'
                          { driver.perfect_foresight_solver(); }
                         ;

perfect_foresight_solver_options_list : perfect_foresight_solver_options_list COMMA perfect_foresight_solver_options
                                     | perfect_foresight_solver_options
                                     ;

perfect_foresight_solver_options : o_stack_solve_algo
                                 | o_markowitz
                                 | o_minimal_solving_periods
                                 | o_simul_maxit
	                         | o_endogenous_terminal_period
                                 | o_no_homotopy
                                 ;

simul : SIMUL ';'
        { driver.simul(); }
      | SIMUL '(' simul_options_list ')' ';'
        { driver.simul(); }
      ;

simul_options_list : simul_options_list COMMA simul_options
                   | simul_options
                   ;

simul_options : perfect_foresight_setup_options
              | perfect_foresight_solver_options
              ;

external_function : EXTERNAL_FUNCTION '(' external_function_options_list ')' ';'
                    { driver.external_function(); }
                  ;

external_function_options_list : external_function_options_list COMMA external_function_options
                               | external_function_options
                               ;

external_function_options : o_ext_func_name
                          | o_ext_func_nargs
                          | o_first_deriv_provided
                          | o_second_deriv_provided
                          ;

stoch_simul : STOCH_SIMUL ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL symbol_list ';'
              { driver.stoch_simul(); }
            | STOCH_SIMUL '(' stoch_simul_options_list ')' symbol_list ';'
              { driver.stoch_simul(); }
            ;

stoch_simul_options_list : stoch_simul_options_list COMMA stoch_simul_options
                         | stoch_simul_options
                         ;

stoch_simul_primary_options : o_dr_algo
                            | o_solve_algo
                            | o_simul_algo
                            | o_linear
                            | o_order
                            | o_replic
                            | o_drop
                            | o_ar
                            | o_nocorr
                            | o_nofunctions
                            | o_nomoments
                            | o_nograph
                            | o_nodisplay
                            | o_graph_format
                            | o_irf
                            | o_irf_shocks
                            | o_relative_irf
                            | o_hp_filter
                            | o_hp_ngrid
                            | o_periods
                            | o_simul
                            | o_simul_seed
                            | o_simul_replic
                            | o_qz_criterium
                            | o_qz_zero_threshold
                            | o_print
                            | o_noprint
                            | o_aim_solver
                            | o_partial_information
                            | o_conditional_variance_decomposition
                            | o_k_order_solver
                            | o_pruning
                            | o_sylvester
                            | o_sylvester_fixed_point_tol
                            | o_dr
                            | o_dr_cycle_reduction_tol
                            | o_dr_logarithmic_reduction_tol
                            | o_dr_logarithmic_reduction_maxiter
                            | o_irf_plot_threshold
                            ;

stoch_simul_options : stoch_simul_primary_options
                    | o_loglinear
                    ;

symbol_list : symbol_list symbol
               { driver.add_in_symbol_list($2); }
             | symbol_list COMMA symbol
               { driver.add_in_symbol_list($3); }
             | symbol
               { driver.add_in_symbol_list($1); }
             ;

symbol_list_ext : symbol_list
                | ':'
                  {
                    string *colon = new string(":");
                    driver.add_in_symbol_list(colon);
                  }
                ;

signed_integer : PLUS INT_NUMBER
                 { $$ = $2; }
               | MINUS INT_NUMBER
                { $2->insert(0, "-"); $$ = $2; }
               | INT_NUMBER
               ;

non_negative_number : INT_NUMBER
                    | FLOAT_NUMBER
                    ;

signed_number : PLUS non_negative_number
               { $$ = $2; }
              | MINUS non_negative_number
               { $2->insert(0, "-"); $$ = $2; }
              | non_negative_number
              ;

signed_inf : PLUS INF_CONSTANT
             { $$ = new string ("Inf"); }
           | MINUS INF_CONSTANT
             { $$ = new string ("-Inf"); }
           | INF_CONSTANT
             { $$ = new string ("Inf"); }
           ;

signed_number_w_inf : signed_inf
                    | signed_number
                    ;

estimated_params : ESTIMATED_PARAMS ';' estimated_list END ';' { driver.estimated_params(); };

estimated_list : estimated_list estimated_elem
                 { driver.add_estimated_params_element(); }
               | estimated_elem
                 { driver.add_estimated_params_element(); }
               ;

estimated_elem : estimated_elem1 COMMA estimated_elem2 ';';

estimated_elem1 : STDERR symbol
                  {
                    driver.estim_params.type = 1;
                    driver.estim_params.name = *$2;
                    delete $2;
                  }
                | symbol
                  {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = *$1;
                    delete $1;
                  }
                | CORR symbol COMMA symbol
                  {
                    driver.estim_params.type = 3;
                    driver.estim_params.name = *$2;
                    driver.estim_params.name2 = *$4;
                    delete $2;
                    delete $4;
                  }
                | DSGE_PRIOR_WEIGHT
                  {
                    driver.estim_params.type = 2;
                    driver.estim_params.name = "dsge_prior_weight";
                  }
                ;

estimated_elem2 : prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.prior = $1;
                  }
                | expression_or_empty COMMA prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.prior = $3;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA prior_pdf COMMA estimated_elem3
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.low_bound = $3;
                    driver.estim_params.up_bound = $5;
                    driver.estim_params.prior = $7;
                  }
                | expression
                  {
                    driver.estim_params.init_val = $1;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.init_val = $1;
                    driver.estim_params.low_bound = $3;
                    driver.estim_params.up_bound = $5;
                  }
                ;

estimated_elem3 : expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression_or_empty
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                    driver.estim_params.p4 = $7;
                  }
                | expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression_or_empty COMMA expression
                  {
                    driver.estim_params.mean = $1;
                    driver.estim_params.std = $3;
                    driver.estim_params.p3 = $5;
                    driver.estim_params.p4 = $7;
                    driver.estim_params.jscale = $9;
                  }
                ;

estimated_params_init : ESTIMATED_PARAMS_INIT ';' estimated_init_list END ';'
                        { driver.estimated_params_init(); }
                      | ESTIMATED_PARAMS_INIT '(' USE_CALIBRATION ')' ';' END ';'
                        { driver.estimated_params_init(true); }
                      | ESTIMATED_PARAMS_INIT '(' USE_CALIBRATION ')' ';' estimated_init_list END ';'
                        { driver.estimated_params_init(true); }
                      ;

estimated_init_list : estimated_init_list estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    | estimated_init_elem
                      { driver.add_estimated_params_element(); }
                    ;

estimated_init_elem : STDERR symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 1;
                        driver.estim_params.name = *$2;
                        driver.estim_params.init_val = $4;
                        delete $2;
                      }
                    | CORR symbol COMMA symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 3;
                        driver.estim_params.name = *$2;
                        driver.estim_params.name2 = *$4;
                        driver.estim_params.init_val = $6;
                        delete $2;
                        delete $4;
                      }
                    | symbol COMMA expression ';'
                      {
                        driver.estim_params.type = 2;
                        driver.estim_params.name = *$1;
                        driver.estim_params.init_val = $3;
                        delete $1;
                      }
                    ;

estimated_params_bounds : ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END ';'
                          { driver.estimated_params_bounds(); };

estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      | estimated_bounds_elem
                        { driver.add_estimated_params_element(); }
                      ;

estimated_bounds_elem : STDERR symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 1;
                          driver.estim_params.name = *$2;
                          driver.estim_params.low_bound = $4;
                          driver.estim_params.up_bound = $6;
                          delete $2;
                        }
                      | CORR symbol COMMA symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 3;
                          driver.estim_params.name = *$2;
                          driver.estim_params.name2 = *$4;
                          driver.estim_params.low_bound = $6;
                          driver.estim_params.up_bound = $8;
                          delete $2;
                          delete $4;
                        }
                      | symbol COMMA expression COMMA expression ';'
                        {
                          driver.estim_params.type = 2;
                          driver.estim_params.name = *$1;
                          driver.estim_params.low_bound = $3;
                          driver.estim_params.up_bound = $5;
                          delete $1;
                        }
                      ;

prior_distribution : BETA
                     { $$ = eBeta; }
                   | GAMMA
                     { $$ = eGamma; }
                   | NORMAL
                     { $$ = eNormal; }
                   | INV_GAMMA
                     { $$ = eInvGamma; }
                   | INV_GAMMA1
                     { $$ = eInvGamma1; }
                   | UNIFORM
                     { $$ = eUniform; }
                   | INV_GAMMA2
                     { $$ = eInvGamma2; }
                   | DIRICHLET
                     { $$ = eDirichlet; }
                   ;

prior_pdf : BETA_PDF
            { $$ = eBeta; }
          | GAMMA_PDF
            { $$ = eGamma; }
          | NORMAL_PDF
            { $$ = eNormal; }
          | INV_GAMMA_PDF
            { $$ = eInvGamma; }
          | INV_GAMMA1_PDF
            { $$ = eInvGamma1; }
          | UNIFORM_PDF
            { $$ = eUniform; }
          | INV_GAMMA2_PDF
            { $$ = eInvGamma2; }
          ;

date_str : DATES { $$ = $1; }

date_expr : date_str
            { $$ = $1; }
          | date_expr PLUS INT_NUMBER
            { $$ = $1; $$->append("+").append(*$3); }
          ;

set_time : SET_TIME '(' date_expr ')' ';'
           { driver.set_time($3); }
         ;

data : DATA '(' data_options_list ')'';'
       { driver.estimation_data(); }
     ;

data_options_list : data_options_list COMMA data_options
                  | data_options
                  ;

data_options : o_file
             | o_new_estimation_data_first_obs
             | o_last_obs
             | o_new_estimation_data_nobs
             | o_xls_sheet
             | o_xls_range
             ;

subsamples : subsamples_eq_opt '(' subsamples_name_list ')' ';'
             { driver.set_subsamples($1->at(0), $1->at(1)); }
           ;

subsamples_eq : subsamples_eq_opt EQUAL subsamples_eq_opt ';'
                {
                  driver.copy_subsamples($1->at(0), $1->at(1), $3->at(0), $3->at(1));
                  delete $1;
                  delete $3;
                }
              ;

subsamples_eq_opt : symbol '.' SUBSAMPLES
                    {
                      $$ = new vector<string *>();
                      $$->push_back($1);
                      $$->push_back(new string (""));
                    }
                  | STD '(' symbol ')' '.'SUBSAMPLES
                    {
                      $$ = new vector<string *>();
                      $$->push_back($3);
                      $$->push_back(new string (""));
                    }
                  | CORR '(' symbol COMMA symbol ')' '.' SUBSAMPLES
                    {
                      $$ = new vector<string *>();
                      $$->push_back($3);
                      $$->push_back($5);
                    }
                 ;

subsamples_name_list : subsamples_name_list COMMA o_subsample_name
                     | o_subsample_name
                     ;

prior : symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_prior($1, new string ("")); }
      | symbol '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_prior($1, $3); }
      | STD '(' symbol ')' '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_std_prior($3, new string ("")); }
      | STD '(' symbol ')' '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_std_prior($3, $6); }
      | CORR '(' symbol COMMA symbol ')' '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_corr_prior($3, $5, new string ("")); }
      | CORR '(' symbol COMMA symbol ')' '.' symbol '.' PRIOR { driver.set_prior_variance(); driver.prior_shape = eNoShape; } '(' prior_options_list ')' ';'
        { driver.set_corr_prior($3, $5, $8); }
      ;

prior_options_list : prior_options_list COMMA prior_options
                   | prior_options
                   ;

prior_options : o_shift
              | o_mean
              | o_median
              | o_stdev
              | o_truncate
              | o_variance
              | o_mode
              | o_interval
              | o_shape
              | o_domain
              ;

prior_eq : prior_eq_opt EQUAL prior_eq_opt ';'
           {
             driver.copy_prior($1->at(0), $1->at(1), $1->at(2), $1->at(3),
                               $3->at(0), $3->at(1), $3->at(2), $3->at(3));
             delete $1;
             delete $3;
           }
         ;

prior_eq_opt : symbol '.' PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("par"));
                 $$->push_back($1);
                 $$->push_back(new string (""));
                 $$->push_back(new string (""));
               }
             | symbol '.' symbol '.' PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("par"));
                 $$->push_back($1);
                 $$->push_back(new string (""));
                 $$->push_back($3);
               }
             | STD '(' symbol ')' '.'  PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("std"));
                 $$->push_back($3);
                 $$->push_back(new string (""));
                 $$->push_back(new string (""));
               }
             | STD '(' symbol ')' '.' symbol '.' PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("std"));
                 $$->push_back($3);
                 $$->push_back(new string (""));
                 $$->push_back($6);
               }
             | CORR '(' symbol COMMA symbol ')' '.'  PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("corr"));
                 $$->push_back($3);
                 $$->push_back($5);
                 $$->push_back(new string (""));
               }
             | CORR '(' symbol COMMA symbol ')' '.' symbol '.' PRIOR
               {
                 $$ = new vector<string *>();
                 $$->push_back(new string ("corr"));
                 $$->push_back($3);
                 $$->push_back($5);
                 $$->push_back($8);
               }
             ;

options : symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_options($1, new string ("")); }
        | symbol '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_options($1, $3); }
        | STD '(' symbol ')' '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_std_options($3, new string ("")); }
        | STD '(' symbol ')' '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_std_options($3, $6); }
        | CORR '(' symbol COMMA symbol ')' '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_corr_options($3, $5, new string ("")); }
        | CORR '(' symbol COMMA symbol ')' '.' symbol '.' OPTIONS '(' options_options_list ')' ';'
          { driver.set_corr_options($3, $5, $8); }
        ;

options_options_list : options_options_list COMMA options_options
                     | options_options
                     ;

options_options : o_jscale
                | o_init
                | o_bounds
                ;

options_eq : options_eq_opt EQUAL options_eq_opt ';'
             {
               driver.copy_options($1->at(0), $1->at(1), $1->at(2), $1->at(3),
                                   $3->at(0), $3->at(1), $3->at(2), $3->at(3));
               delete $1;
               delete $3;
             }
           ;

options_eq_opt : symbol '.' OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("par"));
                   $$->push_back($1);
                   $$->push_back(new string (""));
                   $$->push_back(new string (""));
                 }
               | symbol '.' symbol '.' OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("par"));
                   $$->push_back($1);
                   $$->push_back(new string (""));
                   $$->push_back($3);
                 }
               | STD '(' symbol ')' '.'  OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("std"));
                   $$->push_back($3);
                   $$->push_back(new string (""));
                   $$->push_back(new string (""));
                 }
               | STD '(' symbol ')' '.' symbol '.' OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("std"));
                   $$->push_back($3);
                   $$->push_back(new string (""));
                   $$->push_back($6);
                 }
               | CORR '(' symbol COMMA symbol ')' '.'  OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("corr"));
                   $$->push_back($3);
                   $$->push_back($5);
                   $$->push_back(new string (""));
                 }
               | CORR '(' symbol COMMA symbol ')' '.' symbol '.' OPTIONS
                 {
                   $$ = new vector<string *>();
                   $$->push_back(new string ("corr"));
                   $$->push_back($3);
                   $$->push_back($5);
                   $$->push_back($8);
                 }
               ;

estimation : ESTIMATION ';'
             { driver.run_estimation(); }
           | ESTIMATION '(' estimation_options_list ')' ';'
             { driver.run_estimation(); }
           | ESTIMATION symbol_list ';'
             { driver.run_estimation(); }
           | ESTIMATION '(' estimation_options_list ')' symbol_list ';'
             { driver.run_estimation(); }
           ;

estimation_options_list : estimation_options_list COMMA estimation_options
                        | estimation_options
                        ;

estimation_options : o_datafile
                   | o_nobs
                   | o_first_obs
                   | o_prefilter
                   | o_presample
                   | o_lik_algo
                   | o_lik_init
                   | o_nograph
                   | o_nodisplay
                   | o_graph_format
                   | o_conf_sig
                   | o_mh_conf_sig
                   | o_mh_replic
                   | o_mh_drop
                   | o_mh_jscale
                   | o_optim
                   | o_mh_init_scale
                   | o_mode_file
                   | o_mode_compute
                   | o_mode_check
                   | o_mode_check_neighbourhood_size
                   | o_mode_check_symmetric_plots
                   | o_mode_check_number_of_points
                   | o_prior_trunc
                   | o_mh_mode
                   | o_mh_nblocks
                   | o_load_mh_file
                   | o_loglinear
		   | o_logdata
                   | o_nodiagnostic
                   | o_bayesian_irf
                   | o_dsge_var
                   | o_dsge_varlag
                   | o_irf
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_filtered_vars
                   | o_kalman_algo
                   | o_kalman_tol
                   | o_xls_sheet
                   | o_xls_range
                   | o_filter_step_ahead
                   | o_solve_algo
                   | o_constant
                   | o_noconstant
                   | o_mh_recover
                   | o_diffuse_filter
                   | o_plot_priors
                   | o_order
                   | o_aim_solver
                   | o_partial_information
                   | o_filter_covariance
                   | o_filter_decomposition
                   | o_selected_variables_only
                   | o_conditional_variance_decomposition
                   | o_cova_compute
                   | o_irf_shocks
                   | o_sub_draws
                   | o_sylvester
                   | o_sylvester_fixed_point_tol
                   | o_lyapunov
                   | o_lyapunov_fixed_point_tol
                   | o_lyapunov_doubling_tol
                   | o_dr
                   | o_dr_cycle_reduction_tol
                   | o_dr_logarithmic_reduction_tol
                   | o_dr_logarithmic_reduction_maxiter
                   | o_analytic_derivation
                   | o_ar
                   | o_endogenous_prior
                   | o_use_univariate_filters_if_singularity_is_detected
                   | o_qz_zero_threshold
                   | o_taper_steps
                   | o_geweke_interval
                   | o_mcmc_jumping_covariance
                   | o_irf_plot_threshold
                   | o_posterior_max_subsample_draws
                   | o_consider_all_endogenous
                   | o_consider_only_observed
                   ;

list_optim_option : QUOTED_STRING COMMA QUOTED_STRING
                    { driver.optim_options_string($1, $3); }
                  | QUOTED_STRING COMMA signed_number
                    { driver.optim_options_num($1, $3); }
                  ;

optim_options : list_optim_option
              | optim_options COMMA list_optim_option;
              ;

varobs : VAROBS { driver.check_varobs(); } varobs_list ';';

varobs_list : varobs_list symbol
              { driver.add_varobs($2); }
            | varobs_list COMMA symbol
              { driver.add_varobs($3); }
            | symbol
              { driver.add_varobs($1); }
            ;

observation_trends : OBSERVATION_TRENDS ';' trend_list END ';' { driver.set_trends(); };

trend_list : trend_list trend_element
           | trend_element
           ;

trend_element :  symbol '(' expression ')' ';' { driver.set_trend_element($1, $3); };

unit_root_vars : UNIT_ROOT_VARS symbol_list ';' { driver.set_unit_root_vars(); };

optim_weights : OPTIM_WEIGHTS ';' optim_weights_list END ';' { driver.optim_weights(); };

optim_weights_list : optim_weights_list symbol expression ';'
                     { driver.set_optim_weights($2, $3); }
                   | optim_weights_list symbol COMMA symbol expression ';'
                     { driver.set_optim_weights($2, $4, $5); }
                   | symbol expression ';'
                     { driver.set_optim_weights($1, $2); }
                   | symbol COMMA symbol expression ';'
                     { driver.set_optim_weights($1, $3, $4); }
                   ;

osr_params : OSR_PARAMS symbol_list ';' { driver.set_osr_params(); };


osr_options_list : osr_options_list COMMA osr_options
                 | osr_options
                 ;

osr_options : stoch_simul_primary_options
            | o_osr_maxit
            | o_osr_tolf
            ;

osr : OSR ';'
      { driver.run_osr(); }
    | OSR '(' osr_options_list ')' ';'
      { driver.run_osr(); }
    | OSR symbol_list ';'
      { driver.run_osr(); }
    | OSR '(' osr_options_list ')' symbol_list ';'
      {driver.run_osr(); }
    ;

dynatype : DYNATYPE '(' filename ')' ';'
           { driver.run_dynatype($3); }
         | DYNATYPE '(' filename ')' symbol_list ';'
           { driver.run_dynatype($3); }
         ;

dynasave : DYNASAVE '(' filename ')' ';'
           { driver.run_dynasave($3); }
         | DYNASAVE '(' filename ')' symbol_list ';'
           { driver.run_dynasave($3); }
         ;

load_params_and_steady_state : LOAD_PARAMS_AND_STEADY_STATE '(' filename ')' ';'
                               { driver.run_load_params_and_steady_state($3); }
                             ;

save_params_and_steady_state : SAVE_PARAMS_AND_STEADY_STATE '(' filename ')' ';'
                               { driver.run_save_params_and_steady_state($3); }
                             ;

identification : IDENTIFICATION ';'
                 { driver.run_identification(); }
               | IDENTIFICATION '(' identification_options_list ')' ';'
                 { driver.run_identification(); }
               ;

identification_options_list : identification_option COMMA identification_options_list
                            | identification_option
                            ;

identification_option : o_ar
                      | o_useautocorr
                      | o_load_ident_files
                      | o_prior_mc
                      | o_advanced
                      | o_max_dim_cova_group
                      | o_gsa_prior_range
                      | o_periods
                      | o_replic
                      | o_gsa_sample_file
                      | o_parameter_set
                      | o_lik_init
                      | o_kalman_algo
                      | o_nograph
                      | o_nodisplay
                      | o_graph_format
                      ;

model_comparison : MODEL_COMPARISON mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 | MODEL_COMPARISON '(' o_marginal_density ')' mc_filename_list ';'
                   { driver.run_model_comparison(); }
                 ;

filename : symbol
         | QUOTED_STRING
         ;

parallel_local_filename_list : filename
                               { driver.add_parallel_local_file($1); }
                             | parallel_local_filename_list COMMA filename
                               { driver.add_parallel_local_file($3); }
                             ;

mc_filename_list : filename
                   { driver.add_mc_filename($1); }
                 | filename '(' non_negative_number ')'
                   { driver.add_mc_filename($1, $3); }
                 | mc_filename_list filename
                   { driver.add_mc_filename($2); }
                 | mc_filename_list filename '(' non_negative_number ')'
                   { driver.add_mc_filename($2, $4); }
                 | mc_filename_list COMMA filename
                   { driver.add_mc_filename($3); }
                 | mc_filename_list COMMA filename '(' non_negative_number ')'
                   { driver.add_mc_filename($3, $5); }
                 ;

planner_objective : PLANNER_OBJECTIVE { driver.begin_planner_objective(); }
                    hand_side { driver.end_planner_objective($3); } ';';

ramsey_model : RAMSEY_MODEL ';'
                { driver.ramsey_model(); }
              | RAMSEY_MODEL '(' ramsey_model_options_list ')' ';'
                { driver.ramsey_model(); }
              | RAMSEY_MODEL symbol_list ';'
                { driver.ramsey_model(); }
              | RAMSEY_MODEL '(' ramsey_model_options_list ')' symbol_list ';'
                { driver.ramsey_model(); }
              ;

ramsey_policy : RAMSEY_POLICY ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY symbol_list ';'
                { driver.ramsey_policy(); }
              | RAMSEY_POLICY '(' ramsey_policy_options_list ')' symbol_list ';'
                { driver.ramsey_policy(); }
              ;

discretionary_policy : DISCRETIONARY_POLICY ';'
                       { driver.discretionary_policy(); }
                     | DISCRETIONARY_POLICY '(' discretionary_policy_options_list ')' ';'
                       { driver.discretionary_policy(); }
                     | DISCRETIONARY_POLICY symbol_list ';'
                       { driver.discretionary_policy(); }
                     | DISCRETIONARY_POLICY '(' discretionary_policy_options_list ')' symbol_list ';'
                       { driver.discretionary_policy(); }
                     ;

discretionary_policy_options_list : discretionary_policy_options_list COMMA discretionary_policy_options
                           | discretionary_policy_options
                           ;

discretionary_policy_options : ramsey_policy_options 
                             | o_discretionary_tol;
                             | o_dp_maxit;
                             ;

ramsey_model_options_list : ramsey_model_options_list COMMA ramsey_model_options
                           | ramsey_model_options
                           ;

ramsey_model_options :  o_planner_discount
                      | o_instruments
                      ;

ramsey_policy_options_list : ramsey_policy_options_list COMMA ramsey_policy_options
                           | ramsey_policy_options
                           ;

ramsey_policy_options : stoch_simul_primary_options
                      | o_planner_discount
                      | o_instruments
                      ;

write_latex_dynamic_model : WRITE_LATEX_DYNAMIC_MODEL ';'
                            { driver.write_latex_dynamic_model(); }
                          ;

write_latex_static_model : WRITE_LATEX_STATIC_MODEL ';'
                           { driver.write_latex_static_model(); }
                         ;

shock_decomposition : SHOCK_DECOMPOSITION ';'
                      {driver.shock_decomposition(); }
                    | SHOCK_DECOMPOSITION '(' shock_decomposition_options_list ')' ';'
                      { driver.shock_decomposition(); }
                    | SHOCK_DECOMPOSITION symbol_list ';'
                      { driver.shock_decomposition(); }
                    | SHOCK_DECOMPOSITION '(' shock_decomposition_options_list ')' symbol_list ';'
                      { driver.shock_decomposition(); }
                    ;

bvar_prior_option : o_bvar_prior_tau
                  | o_bvar_prior_decay
                  | o_bvar_prior_lambda
                  | o_bvar_prior_mu
                  | o_bvar_prior_omega
                  | o_bvar_prior_flat
                  | o_bvar_prior_train
                  ;

bvar_common_option : bvar_prior_option
                   | o_datafile
                   | o_xls_sheet
                   | o_xls_range
                   | o_first_obs
                   | o_presample
                   | o_nobs
                   | o_prefilter
                   | o_constant
                   | o_noconstant
                   ;

bvar_density_options_list : bvar_common_option COMMA bvar_density_options_list
                          | bvar_common_option
                          ;

bvar_density : BVAR_DENSITY INT_NUMBER ';'
               { driver.bvar_density($2); }
             | BVAR_DENSITY '(' bvar_density_options_list ')' INT_NUMBER ';'
               { driver.bvar_density($5); }
             ;

bvar_forecast_option : bvar_common_option
                     | o_forecast
                     | o_conf_sig
                     | o_bvar_replic
                     ;

bvar_forecast_options_list : bvar_forecast_option COMMA bvar_forecast_options_list
                           | bvar_forecast_option
                           ;

bvar_forecast : BVAR_FORECAST INT_NUMBER ';'
                { driver.bvar_forecast($2); }
              | BVAR_FORECAST '(' bvar_forecast_options_list ')' INT_NUMBER ';'
                { driver.bvar_forecast($5); }
              ;

sbvar_option : o_datafile
             | o_freq
             | o_initial_year
             | o_initial_subperiod
             | o_final_year
             | o_final_subperiod
             | o_data
             | o_vlist
             | o_vlistlog
             | o_vlistper
             | o_restriction_fname
             | o_nlags
             | o_cross_restrictions
             | o_contemp_reduced_form
             | o_real_pseudo_forecast
             | o_no_bayesian_prior
             | o_dummy_obs
             | o_nstates
             | o_indxscalesstates
             | o_alpha
             | o_beta
             | o_gsig2_lmdm
             | o_q_diag
             | o_flat_prior
             | o_ncsk
             | o_nstd
             | o_ninv
             | o_indxparr
             | o_indxovr
             | o_aband
             | o_indxap
             | o_apband
             | o_indximf
             | o_indxfore
             | o_foreband
             | o_indxgforhat
             | o_indxgimfhat
             | o_indxestima
             | o_indxgdls
             | o_eq_ms
             | o_cms
             | o_ncms
             | o_eq_cms
             | o_tlindx
             | o_tlnumber
             | o_cnum
             | o_forecast
             | o_coefficients_prior_hyperparameters;
             ;

sbvar_options_list : sbvar_option COMMA sbvar_options_list
                   | sbvar_option
                   ;

sbvar : SBVAR ';'
        { driver.sbvar(); }
      | SBVAR '(' sbvar_options_list ')' ';'
        { driver.sbvar(); }
      ;

ms_variance_decomposition_option : o_output_file_tag
                                 | o_file_tag
                                 | o_simulation_file_tag
                                 | o_filtered_probabilities
                                 | o_no_error_bands
                                 | o_error_band_percentiles
                                 | o_shock_draws
                                 | o_shocks_per_parameter
                                 | o_thinning_factor
                                 | o_free_parameters
                                 | o_regime
                                 | o_regimes
                                 | o_parameter_uncertainty
                                 | o_horizon
                                 ;

ms_variance_decomposition_options_list : ms_variance_decomposition_option COMMA ms_variance_decomposition_options_list
                                       | ms_variance_decomposition_option
                                       ;

ms_variance_decomposition : MS_VARIANCE_DECOMPOSITION ';'
                            { driver.ms_variance_decomposition(); }
                          | MS_VARIANCE_DECOMPOSITION '(' ms_variance_decomposition_options_list ')' ';'
                            { driver.ms_variance_decomposition(); }
                          ;

ms_forecast_option : o_output_file_tag
                   | o_file_tag
                   | o_simulation_file_tag
                   | o_data_obs_nbr
                   | o_error_band_percentiles
                   | o_shock_draws
                   | o_shocks_per_parameter
                   | o_thinning_factor
                   | o_free_parameters
                   | o_median
                   | o_regime
                   | o_regimes
                   | o_parameter_uncertainty
                   ;

ms_forecast_options_list : ms_forecast_option COMMA ms_forecast_options_list
                         | ms_forecast_option
                         ;

ms_forecast : MS_FORECAST ';'
              { driver.ms_forecast(); }
            | MS_FORECAST '(' ms_forecast_options_list ')' ';'
              { driver.ms_forecast(); }
            ;

ms_irf_option : o_output_file_tag
              | o_file_tag
              | o_simulation_file_tag
              | o_parameter_uncertainty
              | o_horizon
              | o_filtered_probabilities
              | o_error_band_percentiles
              | o_shock_draws
              | o_shocks_per_parameter
              | o_thinning_factor
              | o_free_parameters
              | o_median
              | o_regime
              | o_regimes
              ;

ms_irf_options_list : ms_irf_option COMMA ms_irf_options_list
                    | ms_irf_option
                    ;

ms_irf : MS_IRF ';'
         { driver.ms_irf(); }
       | MS_IRF '(' ms_irf_options_list ')' ';'
         { driver.ms_irf(); }
       | MS_IRF symbol_list ';'
         { driver.ms_irf(); }
       | MS_IRF '(' ms_irf_options_list ')' symbol_list ';'
         { driver.ms_irf(); }
       ;

ms_compute_probabilities_option : o_output_file_tag
                                | o_file_tag
                                | o_filtered_probabilities
                                | o_real_time_smoothed
                                ;

ms_compute_probabilities_options_list : ms_compute_probabilities_option COMMA ms_compute_probabilities_options_list
                                      | ms_compute_probabilities_option
                                      ;

ms_compute_probabilities : MS_COMPUTE_PROBABILITIES ';'
                           { driver.ms_compute_probabilities(); }
                         | MS_COMPUTE_PROBABILITIES '(' ms_compute_probabilities_options_list ')' ';'
                           { driver.ms_compute_probabilities(); }
                         ;

ms_compute_mdd_option : o_output_file_tag
                      | o_file_tag
                      | o_simulation_file_tag
                      | o_proposal_type
                      | o_proposal_lower_bound
                      | o_proposal_upper_bound
                      | o_proposal_draws
                      | o_use_mean_center
                      ;

ms_compute_mdd_options_list : ms_compute_mdd_option COMMA ms_compute_mdd_options_list
                            | ms_compute_mdd_option
                            ;

ms_compute_mdd : MS_COMPUTE_MDD ';'
                 { driver.ms_compute_mdd(); }
               | MS_COMPUTE_MDD '(' ms_compute_mdd_options_list ')' ';'
                 { driver.ms_compute_mdd(); }
               ;

ms_simulation_option : o_output_file_tag
                     | o_file_tag
                     | o_ms_mh_replic
                     | o_ms_drop
                     | o_thinning_factor
                     | o_adaptive_mh_draws
                     | o_save_draws
                     ;

ms_simulation_options_list : ms_simulation_option COMMA ms_simulation_options_list
                           | ms_simulation_option
                           ;

ms_simulation : MS_SIMULATION ';'
                { driver.ms_simulation(); }
              | MS_SIMULATION '(' ms_simulation_options_list ')' ';'
                { driver.ms_simulation(); }
              ;

ms_estimation_option : o_coefficients_prior_hyperparameters
                     | o_freq
                     | o_initial_year
                     | o_initial_subperiod
                     | o_final_year
                     | o_final_subperiod
                     | o_datafile
                     | o_xls_sheet
                     | o_xls_range
                     | o_nlags
                     | o_cross_restrictions
                     | o_contemp_reduced_form
                     | o_no_bayesian_prior
                     | o_alpha
                     | o_beta
                     | o_gsig2_lmdm
                     | o_specification
                     | o_output_file_tag
                     | o_file_tag
                     | o_no_create_init
                     | o_convergence_starting_value
                     | o_convergence_ending_value
                     | o_convergence_increment_value
                     | o_max_iterations_starting_value
                     | o_max_iterations_increment_value
                     | o_max_block_iterations
                     | o_max_repeated_optimization_runs
                     | o_function_convergence_criterion
                     | o_parameter_convergence_criterion
                     | o_number_of_large_perturbations
                     | o_number_of_small_perturbations
                     | o_number_of_posterior_draws_after_perturbation
                     | o_max_number_of_stages
                     | o_random_function_convergence_criterion
                     | o_random_parameter_convergence_criterion
                     ;

ms_estimation_options_list : ms_estimation_option COMMA ms_estimation_options_list
                           | ms_estimation_option
                           ;

ms_estimation : MS_ESTIMATION ';'
                { driver.ms_estimation(); }
              | MS_ESTIMATION '(' ms_estimation_options_list ')' ';'
                { driver.ms_estimation(); }
              ;

dynare_sensitivity : DYNARE_SENSITIVITY ';'
                     { driver.dynare_sensitivity(); }
                   | DYNARE_SENSITIVITY '(' dynare_sensitivity_options_list ')' ';'
                     { driver.dynare_sensitivity(); }
                   ;

dynare_sensitivity_options_list : dynare_sensitivity_option COMMA dynare_sensitivity_options_list
                                | dynare_sensitivity_option
                                ;

dynare_sensitivity_option : o_gsa_identification
                          | o_gsa_morris
                          | o_gsa_stab
                          | o_gsa_redform
                          | o_gsa_pprior
                          | o_gsa_prior_range
                          | o_gsa_ppost
                          | o_gsa_ilptau
                          | o_gsa_morris_nliv
                          | o_gsa_morris_ntra
                          | o_gsa_nsam
                          | o_gsa_load_redform
                          | o_gsa_load_rmse
                          | o_gsa_load_stab
                          | o_gsa_alpha2_stab
                          | o_gsa_ksstat
                          | o_gsa_logtrans_redform
                          | o_gsa_ksstat_redform
                          | o_gsa_alpha2_redform
                          | o_gsa_rmse
                          | o_gsa_lik_only
                          | o_gsa_pfilt_rmse
                          | o_gsa_istart_rmse
                          | o_gsa_alpha_rmse
                          | o_gsa_alpha2_rmse
                          | o_gsa_threshold_redform
                          | o_gsa_namendo
                          | o_gsa_namexo
                          | o_gsa_namlagendo
                          | o_gsa_var_rmse
                          | o_gsa_neighborhood_width
                          | o_gsa_pvalue_ks
                          | o_gsa_pvalue_corr
                          | o_datafile
                          | o_nobs
                          | o_first_obs
                          | o_prefilter
                          | o_presample
                          | o_nograph
                          | o_nodisplay
                          | o_graph_format
                          | o_conf_sig
                          | o_mh_conf_sig
                          | o_loglinear
                          | o_mode_file
                          | o_load_ident_files
                          | o_useautocorr
                          | o_ar
                          | o_kalman_algo
                          | o_lik_init
                          ;

shock_decomposition_options_list : shock_decomposition_option COMMA shock_decomposition_options_list
                                 | shock_decomposition_option
                                 ;

shock_decomposition_option : o_parameter_set
                           | o_datafile
                           ;

homotopy_setup: HOMOTOPY_SETUP ';' homotopy_list END ';'
               { driver.end_homotopy();};

homotopy_list : homotopy_item
              | homotopy_list homotopy_item
              ;

homotopy_item : symbol COMMA expression COMMA expression ';'
                { driver.homotopy_val($1, $3, $5);}
              | symbol COMMA expression ';'
                { driver.homotopy_val($1, NULL, $3);}
              ;

forecast: FORECAST ';' {driver.forecast();}
          | FORECAST '(' forecast_options ')' ';' {driver.forecast();}
          | FORECAST symbol_list ';' {driver.forecast();}
          | FORECAST '(' forecast_options ')' symbol_list ';' {driver.forecast();}
          ;

forecast_options: forecast_option
          | forecast_options COMMA forecast_option
          ;

forecast_option: o_periods
          | o_conf_sig
          | o_nograph
          | o_nodisplay
          | o_graph_format
          ;

conditional_forecast : CONDITIONAL_FORECAST '(' conditional_forecast_options ')' ';'
                       { driver.conditional_forecast(); }
                     ;

conditional_forecast_options : conditional_forecast_option
                             | conditional_forecast_options COMMA conditional_forecast_option
                             ;

conditional_forecast_option : o_periods
                            | o_replic
                            | o_conf_sig
                            | o_controlled_varexo
                            | o_parameter_set
                            ;

plot_conditional_forecast : PLOT_CONDITIONAL_FORECAST symbol_list ';'
                            { driver.plot_conditional_forecast(); }
                          | PLOT_CONDITIONAL_FORECAST '(' PERIODS EQUAL INT_NUMBER ')' symbol_list ';'
                            { driver.plot_conditional_forecast($5); }
                          ;

conditional_forecast_paths : CONDITIONAL_FORECAST_PATHS ';' conditional_forecast_paths_shock_list END ';'
                             { driver.conditional_forecast_paths(); }
                           ;

conditional_forecast_paths_shock_list : conditional_forecast_paths_shock_elem
                                      | conditional_forecast_paths_shock_list conditional_forecast_paths_shock_elem
                                      ;

conditional_forecast_paths_shock_elem : VAR symbol ';' PERIODS period_list ';' VALUES value_list ';'
                                        { driver.add_det_shock($2, true); }
                                      ;

steady_state_model : STEADY_STATE_MODEL ';' { driver.begin_steady_state_model(); }
                     steady_state_equation_list END ';' { driver.reset_data_tree(); }
                   ;

steady_state_equation_list : steady_state_equation_list steady_state_equation
                           | steady_state_equation
                           ;

steady_state_equation : symbol EQUAL expression ';'
                        { driver.add_steady_state_model_equal($1, $3); }
                      | '[' symbol_list ']' EQUAL expression ';'
                        { driver.add_steady_state_model_equal_multiple($5); }
                      ;

calib_smoother : CALIB_SMOOTHER ';'
                 { driver.calib_smoother(); }
               | CALIB_SMOOTHER '(' calib_smoother_options_list ')' ';'
                 { driver.calib_smoother(); }
               | CALIB_SMOOTHER symbol_list ';'
                 { driver.calib_smoother(); }
               | CALIB_SMOOTHER '(' calib_smoother_options_list ')' symbol_list ';'
                 { driver.calib_smoother(); }
               ;

calib_smoother_options_list : calib_smoother_option COMMA calib_smoother_options_list
                            | calib_smoother_option
                            ;

calib_smoother_option : o_filtered_vars
                      | o_filter_step_ahead
                      | o_datafile
                      ;

extended_path : EXTENDED_PATH ';'
                { driver.extended_path(); }
              | EXTENDED_PATH '(' extended_path_options_list ')' ';'            
                { driver.extended_path(); }
              ;

extended_path_options_list : extended_path_option COMMA extended_path_options_list
                           | extended_path_option
                           ;

extended_path_option : o_periods
                     | o_solver_periods
                     | o_extended_path_order
                     | o_hybrid
                     ;

model_diagnostics : MODEL_DIAGNOSTICS ';'
                    { driver.model_diagnostics(); }
                  ;

calibration_range : '[' signed_number_w_inf signed_number_w_inf ']'
                    {
                      $$ = new vector<string *>();
                      $$->push_back($2);
                      $$->push_back($3);
                    }
                  | '[' signed_number_w_inf COMMA signed_number_w_inf ']'
                    {
                      $$ = new vector<string *>();
                      $$->push_back($2);
                      $$->push_back($4);
                    }
                  | PLUS
                    {
                      $$ = new vector<string *>();
                      $$->push_back(new string("0"));
                      $$->push_back(new string("inf"));
                    }
                  | MINUS
                    {
                      $$ = new vector<string *>();
                      $$->push_back(new string("-inf"));
                      $$->push_back(new string("0"));
                    }
                  ;

moment_calibration : MOMENT_CALIBRATION ';' moment_calibration_list END ';'
                     { driver.end_moment_calibration(); }
                   ;

moment_calibration_list : moment_calibration_item
                        | moment_calibration_list moment_calibration_item
                        ;

moment_calibration_item : symbol COMMA symbol COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, new string("0"), $5); }
                        | symbol COMMA symbol '(' signed_integer ')' COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, $5, $8); }
                        | symbol COMMA symbol '(' signed_integer_range ')' COMMA calibration_range ';'
                          { driver.add_moment_calibration_item($1, $3, $5, $8); }
                        ;

irf_calibration : IRF_CALIBRATION ';' irf_calibration_list END ';'
                  { driver.end_irf_calibration(); }
                ;

irf_calibration_list : irf_calibration_item
                     | irf_calibration_list irf_calibration_item
                     ;

irf_calibration_item : symbol COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, new string("1"), $3, $5); }
                     | symbol '(' INT_NUMBER ')' COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, $3, $6, $8); }
                     | symbol '(' integer_range ')' COMMA symbol COMMA calibration_range ';'
                       { driver.add_irf_calibration_item($1, $3, $6, $8); }
                     ;

smoother2histval : SMOOTHER2HISTVAL ';'
                   { driver.smoother2histval(); }
                 | SMOOTHER2HISTVAL '(' smoother2histval_options_list ')' ';'
                   { driver.smoother2histval(); }
                 ;

smoother2histval_options_list : smoother2histval_option COMMA smoother2histval_options_list
                              | smoother2histval_option
                              ;

smoother2histval_option : o_infile
                        | o_invars
                        | o_period
                        | o_outfile
                        | o_outvars
                        ;

o_dr_algo : DR_ALGO EQUAL INT_NUMBER {
                                       if (*$3 == string("0"))
                                         driver.warning("dr_algo option is now deprecated, and may be removed in a future version of Dynare");
                                       else
                                         driver.error("dr_algo=1 option is no longer supported");
                                     };
o_solve_algo : SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("solve_algo", $3); };
o_simul_algo : SIMUL_ALGO EQUAL INT_NUMBER {
                                             if (*$3 == string("0"))
                                               driver.warning("simul_algo option is now deprecated, and may be removed in a future version of Dynare");
                                             else
                                               driver.error("simul_algo=1 option is no longer supported");
                                           };
o_stack_solve_algo : STACK_SOLVE_ALGO EQUAL INT_NUMBER { driver.option_num("stack_solve_algo", $3); };
o_endogenous_terminal_period : ENDOGENOUS_TERMINAL_PERIOD { driver.option_num("endogenous_terminal_period", "1"); };
o_linear : LINEAR { driver.linear(); };
o_order : ORDER EQUAL INT_NUMBER { driver.option_num("order", $3); };
o_replic : REPLIC EQUAL INT_NUMBER { driver.option_num("replic", $3); };
o_drop : DROP EQUAL INT_NUMBER { driver.option_num("drop", $3); };
o_ar : AR EQUAL INT_NUMBER { driver.option_num("ar", $3); };
o_nocorr : NOCORR { driver.option_num("nocorr", "1"); };
o_nofunctions : NOFUNCTIONS { driver.option_num("nofunctions", "1"); };
o_nomoments : NOMOMENTS { driver.option_num("nomoments", "1"); };
o_irf : IRF EQUAL INT_NUMBER { driver.option_num("irf", $3); };
o_irf_shocks : IRF_SHOCKS EQUAL '(' symbol_list ')' { driver.option_symbol_list("irf_shocks"); };
o_hp_filter : HP_FILTER EQUAL non_negative_number { driver.option_num("hp_filter", $3); };
o_hp_ngrid : HP_NGRID EQUAL INT_NUMBER { driver.option_num("hp_ngrid", $3); };
o_periods : PERIODS EQUAL INT_NUMBER { driver.option_num("periods", $3); };
o_solver_periods : SOLVER_PERIODS EQUAL INT_NUMBER { driver.option_num("ep.periods", $3); };
o_extended_path_order : ORDER EQUAL INT_NUMBER { driver.option_num("ep.stochastic.order", $3); };
o_hybrid : HYBRID { driver.option_num("ep.stochastic.hybrid_order", "2"); };
o_steady_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("steady.maxit", $3); };
o_simul_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("simul.maxit", $3); };
o_dp_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("dp.maxit", $3); };
o_osr_maxit : MAXIT EQUAL INT_NUMBER { driver.option_num("osr.maxit", $3); };
o_osr_tolf : TOLF EQUAL non_negative_number { driver.option_num("osr.tolf", $3); };
o_cutoff : CUTOFF EQUAL non_negative_number { driver.cutoff($3); };
o_markowitz : MARKOWITZ EQUAL non_negative_number { driver.option_num("markowitz", $3); };
o_minimal_solving_periods : MINIMAL_SOLVING_PERIODS EQUAL non_negative_number { driver.option_num("minimal_solving_periods", $3); };
o_mfs : MFS EQUAL INT_NUMBER { driver.mfs($3); };
o_simul : SIMUL; // Do nothing, only here for backward compatibility
o_simul_replic : SIMUL_REPLIC EQUAL INT_NUMBER { driver.option_num("simul_replic", $3); };
o_simul_seed : SIMUL_SEED EQUAL INT_NUMBER { driver.error("'simul_seed' option is no longer supported; use 'set_dynare_seed' command instead"); } ;
o_qz_criterium : QZ_CRITERIUM EQUAL non_negative_number { driver.option_num("qz_criterium", $3); };
o_qz_zero_threshold : QZ_ZERO_THRESHOLD EQUAL non_negative_number { driver.option_num("qz_zero_threshold", $3); };
o_file : FILE EQUAL filename { driver.option_str("file", $3); };
o_datafile : DATAFILE EQUAL filename { driver.option_str("datafile", $3); };
o_nobs : NOBS EQUAL vec_int
         { driver.option_vec_int("nobs", $3); }
       | NOBS EQUAL vec_int_number
         { driver.option_vec_int("nobs", $3); }
       ;
o_conditional_variance_decomposition : CONDITIONAL_VARIANCE_DECOMPOSITION EQUAL vec_int
                                       { driver.option_vec_int("conditional_variance_decomposition", $3); }
                                     | CONDITIONAL_VARIANCE_DECOMPOSITION EQUAL vec_int_number
                                       { driver.option_vec_int("conditional_variance_decomposition", $3); }
                                     ;
o_first_obs : FIRST_OBS EQUAL INT_NUMBER { driver.option_num("first_obs", $3); };
o_new_estimation_data_first_obs : FIRST_OBS EQUAL date_expr
                                  { driver.option_date("first_obs", $3); }
                                ;
o_last_obs : LAST_OBS EQUAL date_expr
             { driver.option_date("last_obs", $3); }
           ;
o_shift : SHIFT EQUAL signed_number { driver.option_num("shift", $3); };
o_shape : SHAPE EQUAL prior_distribution { driver.prior_shape = $3; };
o_mode : MODE EQUAL signed_number { driver.option_num("mode", $3); };
o_mean : MEAN EQUAL signed_number { driver.option_num("mean", $3); };
o_truncate : TRUNCATE EQUAL vec_value { driver.option_num("truncate", $3); };
o_stdev : STDEV EQUAL non_negative_number { driver.option_num("stdev", $3); };
o_jscale : JSCALE EQUAL non_negative_number { driver.option_num("jscale", $3); };
o_init : INIT EQUAL signed_number { driver.option_num("init", $3); };
o_bounds : BOUNDS EQUAL vec_value_w_inf { driver.option_num("bounds", $3); };
o_domain : DOMAINN EQUAL vec_value { driver.option_num("domain", $3); };
o_interval : INTERVAL EQUAL vec_value { driver.option_num("interval", $3); };
o_variance : VARIANCE EQUAL expression { driver.set_prior_variance($3); }
o_new_estimation_data_nobs : NOBS EQUAL INT_NUMBER { driver.option_num("nobs", $3); };
o_prefilter : PREFILTER EQUAL INT_NUMBER { driver.option_num("prefilter", $3); };
o_presample : PRESAMPLE EQUAL INT_NUMBER { driver.option_num("presample", $3); };
o_lik_algo : LIK_ALGO EQUAL INT_NUMBER { driver.option_num("lik_algo", $3); };
o_lik_init : LIK_INIT EQUAL INT_NUMBER { driver.option_num("lik_init", $3); };
o_nograph : NOGRAPH
            { driver.option_num("nograph","1"); }
          | GRAPH
            { driver.option_num("nograph", "0"); }
          ;
o_nodisplay : NODISPLAY { driver.option_num("nodisplay","1"); };
o_graph_format : GRAPH_FORMAT EQUAL allowed_graph_formats
                 { driver.process_graph_format_option(); }
               | GRAPH_FORMAT EQUAL '(' list_allowed_graph_formats ')'
                 { driver.process_graph_format_option(); }
               ;
allowed_graph_formats : EPS
                        { driver.add_graph_format("eps"); }
                      | FIG
                        { driver.add_graph_format("fig"); }
                      | PDF
                        { driver.add_graph_format("pdf"); }
                      | NONE
                        { driver.add_graph_format("none"); }
                      ;
list_allowed_graph_formats : allowed_graph_formats
                           | list_allowed_graph_formats COMMA allowed_graph_formats
                           ;

o_subsample_name : symbol EQUAL date_expr ':' date_expr
                   { driver.set_subsample_name_equal_to_date_range($1, $3, $5); }
                 ;
o_conf_sig : CONF_SIG EQUAL non_negative_number { driver.option_num("conf_sig", $3); };
o_mh_conf_sig : MH_CONF_SIG EQUAL non_negative_number { driver.option_num("mh_conf_sig", $3); };
o_mh_replic : MH_REPLIC EQUAL INT_NUMBER { driver.option_num("mh_replic", $3); };
o_posterior_max_subsample_draws : POSTERIOR_MAX_SUBSAMPLE_DRAWS EQUAL INT_NUMBER { driver.option_num("posterior_max_subsample_draws", $3); };
o_mh_drop : MH_DROP EQUAL non_negative_number { driver.option_num("mh_drop", $3); };
o_mh_jscale : MH_JSCALE EQUAL non_negative_number { driver.option_num("mh_jscale", $3); };
o_optim : OPTIM  EQUAL '(' optim_options ')';
o_mh_init_scale : MH_INIT_SCALE EQUAL non_negative_number { driver.option_num("mh_init_scale", $3); };
o_mode_file : MODE_FILE EQUAL filename { driver.option_str("mode_file", $3); };
o_mode_compute : MODE_COMPUTE EQUAL INT_NUMBER { driver.option_num("mode_compute", $3); };
               | MODE_COMPUTE EQUAL symbol { driver.option_str("mode_compute", $3); };
o_mode_check : MODE_CHECK { driver.option_num("mode_check.status", "1"); };
o_mode_check_neighbourhood_size : MODE_CHECK_NEIGHBOURHOOD_SIZE EQUAL signed_number_w_inf { driver.option_num("mode_check.neighbourhood_size", $3); };
o_mode_check_number_of_points : MODE_CHECK_NUMBER_OF_POINTS EQUAL INT_NUMBER { driver.option_num("mode_check.number_of_points", $3); };
o_mode_check_symmetric_plots : MODE_CHECK_SYMMETRIC_PLOTS EQUAL INT_NUMBER { driver.option_num("mode_check.symmetric_plots", $3); };
o_prior_trunc : PRIOR_TRUNC EQUAL non_negative_number { driver.option_num("prior_trunc", $3); };
o_mh_mode : MH_MODE EQUAL INT_NUMBER { driver.option_num("mh_mode", $3); };
o_mh_nblocks : MH_NBLOCKS EQUAL INT_NUMBER { driver.option_num("mh_nblck", $3); };
o_load_mh_file : LOAD_MH_FILE { driver.option_num("load_mh_file", "1"); };
o_loglinear : LOGLINEAR { driver.option_num("loglinear", "1"); };
o_logdata : LOGDATA { driver.option_num("logdata", "1"); };
o_nodiagnostic : NODIAGNOSTIC { driver.option_num("nodiagnostic", "1"); };
o_bayesian_irf : BAYESIAN_IRF { driver.option_num("bayesian_irf", "1"); };
o_dsge_var : DSGE_VAR EQUAL non_negative_number
             { driver.option_num("dsge_var", $3); }
           | DSGE_VAR EQUAL INF_CONSTANT
             { driver.option_num("dsge_var", "Inf"); }
           | DSGE_VAR
             { driver.option_str("dsge_var", "NaN"); }
           ;
o_dsge_varlag : DSGE_VARLAG EQUAL INT_NUMBER { driver.option_num("dsge_varlag", $3); };
o_tex : TEX { driver.option_num("TeX", "1"); };
o_forecast : FORECAST EQUAL INT_NUMBER { driver.option_num("forecast", $3); };
o_smoother : SMOOTHER { driver.option_num("smoother", "1"); };
o_moments_varendo : MOMENTS_VARENDO { driver.option_num("moments_varendo", "1"); };
o_filtered_vars : FILTERED_VARS { driver.option_num("filtered_vars", "1"); };
o_relative_irf : RELATIVE_IRF { driver.option_num("relative_irf", "1"); };
o_kalman_algo : KALMAN_ALGO EQUAL INT_NUMBER { driver.option_num("kalman_algo", $3); };
o_kalman_tol : KALMAN_TOL EQUAL non_negative_number { driver.option_num("kalman_tol", $3); };
o_marginal_density : MARGINAL_DENSITY EQUAL LAPLACE
                     { driver.option_str("mc_marginal_density", "laplace"); }
                   | MARGINAL_DENSITY EQUAL MODIFIEDHARMONICMEAN
                     { driver.option_str("mc_marginal_density", "modifiedharmonicmean"); }
                   ;
o_print : PRINT { driver.option_num("noprint", "0"); };
o_noprint : NOPRINT { driver.option_num("noprint", "1"); };
o_xls_sheet : XLS_SHEET EQUAL symbol { driver.option_str("xls_sheet", $3); };
o_xls_range : XLS_RANGE EQUAL range { driver.option_str("xls_range", $3); };
o_filter_step_ahead : FILTER_STEP_AHEAD EQUAL vec_int { driver.option_vec_int("filter_step_ahead", $3); };
o_taper_steps : TAPER_STEPS EQUAL vec_int { driver.option_vec_int("taper_steps", $3); };
o_geweke_interval : GEWEKE_INTERVAL EQUAL vec_value { driver.option_num("geweke_interval",$3); };
o_constant : CONSTANT { driver.option_num("noconstant", "0"); };
o_noconstant : NOCONSTANT { driver.option_num("noconstant", "1"); };
o_mh_recover : MH_RECOVER { driver.option_num("mh_recover", "1"); };
o_diffuse_filter: DIFFUSE_FILTER {driver.option_num("diffuse_filter", "1"); };
o_plot_priors: PLOT_PRIORS EQUAL INT_NUMBER {driver.option_num("plot_priors", $3); };
o_aim_solver: AIM_SOLVER {driver.option_num("aim_solver", "1"); };
o_partial_information : PARTIAL_INFORMATION {driver.option_num("partial_information", "1"); };
o_sub_draws: SUB_DRAWS EQUAL INT_NUMBER {driver.option_num("sub_draws",$3);};
o_planner_discount : PLANNER_DISCOUNT EQUAL expression { driver.declare_optimal_policy_discount_factor_parameter($3); };
o_sylvester : SYLVESTER EQUAL FIXED_POINT {driver.option_num("sylvester_fp", "1"); }
               | SYLVESTER EQUAL DEFAULT {driver.option_num("sylvester_fp", "0"); };
o_sylvester_fixed_point_tol : SYLVESTER_FIXED_POINT_TOL EQUAL non_negative_number {driver.option_num("sylvester_fixed_point_tol",$3);};
o_lyapunov : LYAPUNOV EQUAL FIXED_POINT {driver.option_num("lyapunov_fp", "1"); }
              | LYAPUNOV EQUAL DOUBLING {driver.option_num("lyapunov_db", "1"); }
              | LYAPUNOV EQUAL SQUARE_ROOT_SOLVER {driver.option_num("lyapunov_srs", "1"); }
              | LYAPUNOV EQUAL DEFAULT {driver.option_num("lyapunov_fp", "0");driver.option_num("lyapunov_db", "0"); driver.option_num("lyapunov_srs", "0");};
o_lyapunov_fixed_point_tol : LYAPUNOV_FIXED_POINT_TOL EQUAL non_negative_number {driver.option_num("lyapunov_fixed_point_tol",$3);};
o_lyapunov_doubling_tol : LYAPUNOV_DOUBLING_TOL EQUAL non_negative_number {driver.option_num("lyapunov_doubling_tol",$3);};
o_dr : DR EQUAL CYCLE_REDUCTION {driver.option_num("dr_cycle_reduction", "1"); }
       | DR EQUAL LOGARITHMIC_REDUCTION {driver.option_num("dr_logarithmic_reduction", "1"); }
       | DR EQUAL DEFAULT {driver.option_num("dr_cycle_reduction", "0"); driver.option_num("dr_logarithmic_reduction", "0");};
o_dr_cycle_reduction_tol : DR_CYCLE_REDUCTION_TOL EQUAL non_negative_number {driver.option_num("dr_cycle_reduction_tol",$3);};
o_dr_logarithmic_reduction_tol : DR_LOGARITHMIC_REDUCTION_TOL EQUAL non_negative_number {driver.option_num("dr_logarithmic_reduction_tol",$3);};
o_dr_logarithmic_reduction_maxiter : DR_LOGARITHMIC_REDUCTION_MAXITER EQUAL INT_NUMBER {driver.option_num("dr_logarithmic_reduction_maxiter",$3);};

o_bvar_prior_tau : BVAR_PRIOR_TAU EQUAL signed_number { driver.option_num("bvar_prior_tau", $3); };
o_bvar_prior_decay : BVAR_PRIOR_DECAY EQUAL non_negative_number { driver.option_num("bvar_prior_decay", $3); };
o_bvar_prior_lambda : BVAR_PRIOR_LAMBDA EQUAL signed_number { driver.option_num("bvar_prior_lambda", $3); };
o_bvar_prior_mu : BVAR_PRIOR_MU EQUAL non_negative_number { driver.option_num("bvar_prior_mu", $3); };
o_bvar_prior_omega : BVAR_PRIOR_OMEGA EQUAL INT_NUMBER { driver.option_num("bvar_prior_omega", $3); };
o_bvar_prior_flat : BVAR_PRIOR_FLAT { driver.option_num("bvar_prior_flat", "1"); };
o_bvar_prior_train : BVAR_PRIOR_TRAIN EQUAL INT_NUMBER { driver.option_num("bvar_prior_train", $3); };
o_bvar_replic : BVAR_REPLIC EQUAL INT_NUMBER { driver.option_num("bvar_replic", $3); };

o_gsa_identification : IDENTIFICATION EQUAL INT_NUMBER { driver.option_num("identification", $3); }; /*not in doc */
o_gsa_morris : MORRIS EQUAL INT_NUMBER { driver.option_num("morris", $3); };
o_gsa_stab : STAB  EQUAL INT_NUMBER { driver.option_num("stab", $3); };
o_gsa_redform : REDFORM  EQUAL INT_NUMBER { driver.option_num("redform", $3); };
o_gsa_pprior : PPRIOR  EQUAL INT_NUMBER { driver.option_num("pprior", $3); };
o_gsa_prior_range : PRIOR_RANGE  EQUAL INT_NUMBER { driver.option_num("prior_range", $3); };
o_gsa_ppost : PPOST  EQUAL INT_NUMBER { driver.option_num("ppost", $3); };
o_gsa_ilptau : ILPTAU EQUAL INT_NUMBER { driver.option_num("ilptau", $3); };
o_gsa_morris_nliv : MORRIS_NLIV EQUAL INT_NUMBER { driver.option_num("morris_nliv", $3); };
o_gsa_morris_ntra : MORRIS_NTRA EQUAL INT_NUMBER { driver.option_num("morris_ntra", $3); };
o_gsa_nsam : NSAM EQUAL INT_NUMBER { driver.option_num("Nsam", $3); }; /* not in doc ??*/
o_gsa_load_redform : LOAD_REDFORM EQUAL INT_NUMBER { driver.option_num("load_redform", $3); };
o_gsa_load_rmse : LOAD_RMSE EQUAL INT_NUMBER { driver.option_num("load_rmse", $3); };
o_gsa_load_stab : LOAD_STAB EQUAL INT_NUMBER { driver.option_num("load_stab", $3); };
o_gsa_alpha2_stab : ALPHA2_STAB EQUAL non_negative_number { driver.option_num("alpha2_stab", $3); };
o_gsa_ksstat : KSSTAT EQUAL non_negative_number { driver.option_num("ksstat", $3); };
o_gsa_logtrans_redform : LOGTRANS_REDFORM EQUAL INT_NUMBER { driver.option_num("logtrans_redform", $3); };
o_gsa_threshold_redform : THRESHOLD_REDFORM EQUAL vec_value_w_inf { driver.option_num("threshold_redform",$3); };
o_gsa_ksstat_redform : KSSTAT_REDFORM EQUAL non_negative_number { driver.option_num("ksstat_redform", $3); };
o_gsa_alpha2_redform : ALPHA2_REDFORM EQUAL non_negative_number { driver.option_num("alpha2_redform", $3); };
o_gsa_namendo : NAMENDO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namendo"); };
o_gsa_namlagendo : NAMLAGENDO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namlagendo"); };
o_gsa_namexo : NAMEXO EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("namexo"); };
o_gsa_rmse : RMSE EQUAL INT_NUMBER { driver.option_num("rmse", $3); };
o_gsa_lik_only : LIK_ONLY EQUAL INT_NUMBER { driver.option_num("lik_only", $3); };
o_gsa_var_rmse : VAR_RMSE EQUAL '(' symbol_list_ext ')' { driver.option_symbol_list("var_rmse"); };
o_gsa_pfilt_rmse : PFILT_RMSE EQUAL non_negative_number { driver.option_num("pfilt_rmse", $3); };
o_gsa_istart_rmse : ISTART_RMSE EQUAL INT_NUMBER { driver.option_num("istart_rmse", $3); };
o_gsa_alpha_rmse : ALPHA_RMSE EQUAL non_negative_number { driver.option_num("alpha_rmse", $3); };
o_gsa_alpha2_rmse : ALPHA2_RMSE EQUAL non_negative_number { driver.option_num("alpha2_rmse", $3); };
o_gsa_sample_file : GSA_SAMPLE_FILE EQUAL INT_NUMBER
                    { driver.option_num("gsa_sample_file", $3); }
                  | GSA_SAMPLE_FILE EQUAL filename
                    { driver.option_str("gsa_sample_file", $3); }
                  ;
o_gsa_neighborhood_width : NEIGHBORHOOD_WIDTH EQUAL non_negative_number { driver.option_num("neighborhood_width", $3); };
o_gsa_pvalue_ks : PVALUE_KS EQUAL  non_negative_number { driver.option_num("pvalue_ks", $3); };
o_gsa_pvalue_corr : PVALUE_CORR EQUAL  non_negative_number { driver.option_num("pvalue_corr", $3); };
o_load_ident_files : LOAD_IDENT_FILES EQUAL INT_NUMBER { driver.option_num("load_ident_files", $3); }
o_useautocorr : USEAUTOCORR EQUAL INT_NUMBER { driver.option_num("useautocorr", $3); }
o_prior_mc : PRIOR_MC EQUAL INT_NUMBER { driver.option_num("prior_mc", $3); }
o_advanced : ADVANCED EQUAL signed_integer { driver.option_num("advanced", $3); }
o_max_dim_cova_group : MAX_DIM_COVA_GROUP EQUAL INT_NUMBER { driver.option_num("max_dim_cova_group", $3); }

o_homotopy_mode : HOMOTOPY_MODE EQUAL INT_NUMBER {driver.option_num("homotopy_mode",$3); };
o_homotopy_steps : HOMOTOPY_STEPS EQUAL INT_NUMBER {driver.option_num("homotopy_steps",$3); };
o_homotopy_force_continue: HOMOTOPY_FORCE_CONTINUE EQUAL INT_NUMBER { driver.option_num("homotopy_force_continue",$3); };
o_nocheck : NOCHECK {driver.option_num("steadystate.nocheck","1"); };

o_controlled_varexo : CONTROLLED_VAREXO EQUAL '(' symbol_list ')' { driver.option_symbol_list("controlled_varexo"); };
o_parameter_set : PARAMETER_SET EQUAL PRIOR_MODE
                  { driver.option_str("parameter_set", "prior_mode"); }
                | PARAMETER_SET EQUAL PRIOR_MEAN
                  { driver.option_str("parameter_set", "prior_mean"); }
                | PARAMETER_SET EQUAL POSTERIOR_MEAN
                  { driver.option_str("parameter_set", "posterior_mean"); }
                | PARAMETER_SET EQUAL POSTERIOR_MODE
                  { driver.option_str("parameter_set", "posterior_mode"); }
                | PARAMETER_SET EQUAL POSTERIOR_MEDIAN
                  { driver.option_str("parameter_set", "posterior_median"); }
                | PARAMETER_SET EQUAL CALIBRATION
                  { driver.option_str("parameter_set", "calibration"); }
                ;
o_ms_drop : DROP EQUAL INT_NUMBER { driver.option_num("ms.drop", $3); };
o_ms_mh_replic : MH_REPLIC EQUAL INT_NUMBER { driver.option_num("ms.mh_replic", $3); };
o_freq : FREQ EQUAL INT_NUMBER
         { driver.option_num("ms.freq",$3); }
       | FREQ EQUAL MONTHLY
         { driver.option_num("ms.freq","12"); }
       | FREQ EQUAL QUARTERLY
         { driver.option_num("ms.freq","4"); }
       ;
o_initial_year : INITIAL_YEAR EQUAL INT_NUMBER {driver.option_num("ms.initial_year",$3); };
o_initial_subperiod : INITIAL_SUBPERIOD EQUAL INT_NUMBER {driver.option_num("ms.initial_subperiod",$3); };
o_final_year : FINAL_YEAR EQUAL INT_NUMBER {driver.option_num("ms.final_year",$3); };
o_final_subperiod : FINAL_SUBPERIOD EQUAL INT_NUMBER {driver.option_num("ms.final_subperiod",$3); };
o_data : DATA EQUAL filename { driver.option_str("ms.data", $3); };
o_vlist : VLIST EQUAL INT_NUMBER {driver.option_num("ms.vlist",$3); };
o_vlistlog : VLISTLOG EQUAL '(' symbol_list ')' {driver.option_symbol_list("ms.vlistlog"); };
o_vlistper : VLISTPER EQUAL INT_NUMBER {driver.option_num("ms.vlistper",$3); };
o_restriction_fname : RESTRICTION_FNAME EQUAL NAME
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname",$3);
                      }
                    | RESTRICTION_FNAME EQUAL UPPER_CHOLESKY
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname","upper_cholesky");
                      }
                    | RESTRICTION_FNAME EQUAL LOWER_CHOLESKY
                      {
                        driver.warning("restriction_fname is now deprecated, and may be removed in a future version of Dynare. Use svar_identification instead.");
                        driver.option_str("ms.restriction_fname","lower_cholesky");
                      }
                    ;
o_nlags : NLAGS EQUAL INT_NUMBER {driver.option_num("ms.nlags",$3); };
o_cross_restrictions : CROSS_RESTRICTIONS {driver.option_num("ms.cross_restrictions","1"); };
o_contemp_reduced_form : CONTEMP_REDUCED_FORM {driver.option_num("ms.contemp_reduced_form","1"); };
o_real_pseudo_forecast : REAL_PSEUDO_FORECAST EQUAL INT_NUMBER {driver.option_num("ms.real_pseudo_forecast",$3); };
o_no_bayesian_prior : NO_BAYESIAN_PRIOR {driver.option_num("ms.bayesian_prior","0"); };
o_dummy_obs : DUMMY_OBS EQUAL INT_NUMBER {driver.option_num("ms.dummy_obs",$3); };
o_nstates : NSTATES EQUAL INT_NUMBER {driver.option_num("ms.nstates",$3); };
o_indxscalesstates : INDXSCALESSTATES EQUAL INT_NUMBER {driver.option_num("ms.indxscalesstates",$3); };
o_alpha : ALPHA EQUAL non_negative_number {driver.option_num("ms.alpha",$3); };
o_beta : BETA EQUAL non_negative_number {driver.option_num("ms.beta",$3); };
o_gsig2_lmdm : GSIG2_LMDM EQUAL INT_NUMBER {driver.option_num("ms.gsig2_lmdm",$3); };
o_specification : SPECIFICATION EQUAL SIMS_ZHA
                  {driver.option_num("ms.specification","1"); }
                | SPECIFICATION EQUAL NONE
                  {driver.option_num("ms.specification","0"); }
                ;
o_q_diag : Q_DIAG EQUAL non_negative_number {driver.option_num("ms.q_diag",$3); };
o_flat_prior : FLAT_PRIOR EQUAL INT_NUMBER {driver.option_num("ms.flat_prior",$3); };
o_ncsk : NCSK EQUAL INT_NUMBER {driver.option_num("ms.ncsk",$3); };
o_nstd : NSTD EQUAL INT_NUMBER {driver.option_num("ms.nstd",$3); };
o_ninv : NINV EQUAL INT_NUMBER {driver.option_num("ms.ninv",$3); };
o_indxparr : INDXPARR EQUAL INT_NUMBER {driver.option_num("ms.indxparr",$3); };
o_indxovr : INDXOVR EQUAL INT_NUMBER {driver.option_num("ms.indxovr",$3); };
o_aband : ABAND EQUAL INT_NUMBER {driver.option_num("ms.aband",$3); };
o_indxap : INDXAP EQUAL INT_NUMBER {driver.option_num("ms.indxap",$3); };
o_apband : APBAND EQUAL INT_NUMBER {driver.option_num("ms.apband",$3); };
o_indximf : INDXIMF EQUAL INT_NUMBER {driver.option_num("ms.indximf",$3); };
o_indxfore : INDXFORE EQUAL INT_NUMBER {driver.option_num("ms.indxfore",$3); };
o_foreband : FOREBAND EQUAL INT_NUMBER {driver.option_num("ms.foreband",$3); };
o_indxgforhat : INDXGFOREHAT EQUAL INT_NUMBER {driver.option_num("ms.indxgforehat",$3); };
o_indxgimfhat : INDXGIMFHAT EQUAL INT_NUMBER {driver.option_num("ms.indxgimfhat",$3); };
o_indxestima : INDXESTIMA EQUAL INT_NUMBER {driver.option_num("ms.indxestima",$3); };
o_indxgdls : INDXGDLS EQUAL INT_NUMBER {driver.option_num("ms.indxgdls",$3); };
o_eq_ms : EQ_MS EQUAL INT_NUMBER {driver.option_num("ms.eq_ms",$3); };
o_cms : CMS EQUAL INT_NUMBER {driver.option_num("ms.cms",$3); };
o_ncms : NCMS EQUAL INT_NUMBER {driver.option_num("ms.ncms",$3); };
o_eq_cms : EQ_CMS EQUAL INT_NUMBER {driver.option_num("ms.eq_cms",$3); };
o_tlindx : TLINDX EQUAL INT_NUMBER {driver.option_num("ms.tlindx",$3); };
o_tlnumber : TLNUMBER EQUAL INT_NUMBER {driver.option_num("ms.tlnumber",$3); };
o_cnum : CNUM EQUAL INT_NUMBER {driver.option_num("ms.cnum",$3); };
o_k_order_solver : K_ORDER_SOLVER {driver.option_num("k_order_solver","1"); };
o_pruning : PRUNING { driver.option_num("pruning", "1"); };
o_chain : CHAIN EQUAL INT_NUMBER { driver.option_num("ms.chain",$3); };
o_restrictions : RESTRICTIONS EQUAL vec_of_vec_value
                 { driver.option_num("ms.restrictions",$3); }
               ;
o_duration : DURATION EQUAL non_negative_number
             { driver.option_num("ms.duration",$3); }
           | DURATION EQUAL vec_value_w_inf
             { driver.option_num("ms.duration",$3); }
           ;
o_number_of_regimes : NUMBER_OF_REGIMES EQUAL INT_NUMBER { driver.option_num("ms.number_of_regimes",$3); };
o_number_of_lags : NUMBER_OF_LAGS EQUAL INT_NUMBER { driver.option_num("ms.number_of_lags",$3); };
o_parameters : PARAMETERS EQUAL '[' symbol_list ']' { driver.option_symbol_list("ms.parameters"); };
o_coefficients : COEFFICIENTS { driver.option_str("ms.coefficients","svar_coefficients"); };
o_variances : VARIANCES { driver.option_str("ms.variances","svar_variances"); };
o_equations : EQUATIONS EQUAL vec_int
              { driver.option_vec_int("ms.equations",$3); }
            | EQUATIONS EQUAL vec_int_number
              { driver.option_vec_int("ms.equations",$3); }
            ;

o_instruments : INSTRUMENTS EQUAL '(' symbol_list ')' {driver.option_symbol_list("instruments"); };

o_ext_func_name : EXT_FUNC_NAME EQUAL filename { driver.external_function_option("name", $3); };
o_ext_func_nargs : EXT_FUNC_NARGS EQUAL INT_NUMBER { driver.external_function_option("nargs",$3); };
o_first_deriv_provided : FIRST_DERIV_PROVIDED EQUAL filename
                         { driver.external_function_option("first_deriv_provided", $3); }
                       | FIRST_DERIV_PROVIDED
                         { driver.external_function_option("first_deriv_provided", ""); }
                       ;
o_second_deriv_provided : SECOND_DERIV_PROVIDED EQUAL filename
                          { driver.external_function_option("second_deriv_provided", $3); }
                        | SECOND_DERIV_PROVIDED
                          { driver.external_function_option("second_deriv_provided", ""); }
                        ;
o_filter_covariance : FILTER_COVARIANCE 
                        { driver.option_num("filter_covariance","1");}
                      ;
o_filter_decomposition : FILTER_DECOMPOSITION
                           { driver.option_num("filter_decomposition","1");}
                         ;
o_selected_variables_only : SELECTED_VARIABLES_ONLY
                           { driver.option_num("selected_variables_only","1");}
                          ;
o_cova_compute : COVA_COMPUTE EQUAL INT_NUMBER
                 { driver.option_num("cova_compute",$3);}
               ;
o_output_file_tag : OUTPUT_FILE_TAG EQUAL filename {driver.option_str("ms.output_file_tag", $3); };
o_file_tag : FILE_TAG EQUAL filename { driver.option_str("ms.file_tag", $3); };
o_no_create_init : NO_CREATE_INIT { driver.option_num("ms.create_init", "0"); };
o_simulation_file_tag : SIMULATION_FILE_TAG EQUAL filename { driver.option_str("ms.simulation_file_tag", $3); };
o_coefficients_prior_hyperparameters : COEFFICIENTS_PRIOR_HYPERPARAMETERS EQUAL vec_value
                                       { driver.option_num("ms.coefficients_prior_hyperparameters",$3); };
o_convergence_starting_value : CONVERGENCE_STARTING_VALUE EQUAL non_negative_number
                               { driver.option_num("ms.convergence_starting_value",$3); };
o_convergence_ending_value : CONVERGENCE_ENDING_VALUE EQUAL non_negative_number
                             { driver.option_num("ms.convergence_ending_value",$3); };
o_convergence_increment_value : CONVERGENCE_INCREMENT_VALUE EQUAL non_negative_number
                                { driver.option_num("ms.convergence_increment_value",$3); };
o_max_iterations_starting_value : MAX_ITERATIONS_STARTING_VALUE EQUAL INT_NUMBER
                                  { driver.option_num("ms.max_iterations_starting_value",$3); };
o_max_iterations_increment_value : MAX_ITERATIONS_INCREMENT_VALUE EQUAL non_negative_number
                                   { driver.option_num("ms.max_iterations_increment_value",$3); };
o_max_block_iterations : MAX_BLOCK_ITERATIONS EQUAL INT_NUMBER
                         { driver.option_num("ms.max_block_iterations",$3); };
o_max_repeated_optimization_runs : MAX_REPEATED_OPTIMIZATION_RUNS EQUAL INT_NUMBER
                                   { driver.option_num("ms.max_repeated_optimization_runs",$3); };
o_function_convergence_criterion : FUNCTION_CONVERGENCE_CRITERION EQUAL non_negative_number
                                   { driver.option_num("ms.function_convergence_criterion",$3); };
o_parameter_convergence_criterion : PARAMETER_CONVERGENCE_CRITERION EQUAL non_negative_number
                                    { driver.option_num("ms.parameter_convergence_criterion",$3); };
o_number_of_large_perturbations : NUMBER_OF_LARGE_PERTURBATIONS EQUAL INT_NUMBER
                                  { driver.option_num("ms.number_of_large_perturbations",$3); };
o_number_of_small_perturbations : NUMBER_OF_SMALL_PERTURBATIONS EQUAL INT_NUMBER
                                  { driver.option_num("ms.number_of_small_perturbations",$3); };
o_number_of_posterior_draws_after_perturbation : NUMBER_OF_POSTERIOR_DRAWS_AFTER_PERTURBATION EQUAL INT_NUMBER
                                                 { driver.option_num("ms.number_of_posterior_draws_after_perturbation",$3); };
o_max_number_of_stages : MAX_NUMBER_OF_STAGES EQUAL INT_NUMBER
                         { driver.option_num("ms.max_number_of_stages",$3); };
o_random_function_convergence_criterion : RANDOM_FUNCTION_CONVERGENCE_CRITERION EQUAL non_negative_number
                                          { driver.option_num("ms.random_function_convergence_criterion",$3); };
o_random_parameter_convergence_criterion : RANDOM_PARAMETER_CONVERGENCE_CRITERION EQUAL non_negative_number
                                           { driver.option_num("ms.random_parameter_convergence_criterion",$3); };
o_thinning_factor : THINNING_FACTOR EQUAL INT_NUMBER { driver.option_num("ms.thinning_factor",$3); };
o_adaptive_mh_draws : ADAPTIVE_MH_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.adaptive_mh_draws",$3); };
o_save_draws : SAVE_DRAWS { driver.option_num("ms.save_draws","1"); };
o_proposal_draws : PROPOSAL_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.proposal_draws",$3); };
o_use_mean_center : USE_MEAN_CENTER { driver.option_num("ms.use_mean_center","1"); };
o_proposal_type : PROPOSAL_TYPE EQUAL INT_NUMBER { driver.option_num("ms.proposal_type",$3); }
o_proposal_lower_bound : PROPOSAL_LOWER_BOUND EQUAL signed_number { driver.option_num("ms.proposal_lower_bound",$3); }
o_proposal_upper_bound : PROPOSAL_UPPER_BOUND EQUAL signed_number { driver.option_num("ms.proposal_upper_bound",$3); }
o_parameter_uncertainty : PARAMETER_UNCERTAINTY { driver.option_num("ms.parameter_uncertainty","1"); };
o_horizon : HORIZON EQUAL INT_NUMBER { driver.option_num("ms.horizon",$3); };
o_filtered_probabilities : FILTERED_PROBABILITIES { driver.option_num("ms.filtered_probabilities","1"); };
o_real_time_smoothed : REAL_TIME_SMOOTHED { driver.option_num("ms.real_time_smoothed_probabilities","1"); };
o_no_error_bands : NO_ERROR_BANDS { driver.option_num("ms.error_bands","0"); };
o_error_band_percentiles : ERROR_BAND_PERCENTILES EQUAL vec_value { driver.option_num("ms.percentiles",$3); };
o_shock_draws : SHOCK_DRAWS EQUAL INT_NUMBER { driver.option_num("ms.shock_draws",$3); };
o_shocks_per_parameter : SHOCKS_PER_PARAMETER EQUAL INT_NUMBER { driver.option_num("ms.shocks_per_parameter",$3); };
o_free_parameters : FREE_PARAMETERS EQUAL vec_value { driver.option_num("ms.free_parameters",$3); };
o_median : MEDIAN { driver.option_num("ms.median","1"); }
         | MEDIAN EQUAL signed_number { driver.option_num("median", $3); };
o_regimes : REGIMES { driver.option_num("ms.regimes","1"); };
o_regime : REGIME EQUAL INT_NUMBER { driver.option_num("ms.regime",$3); };
o_data_obs_nbr : DATA_OBS_NBR EQUAL INT_NUMBER { driver.option_num("ms.forecast_data_obs",$3); };
o_discretionary_tol: DISCRETIONARY_TOL EQUAL non_negative_number { driver.option_num("discretionary_tol",$3); };
o_analytic_derivation : ANALYTIC_DERIVATION { driver.option_num("analytic_derivation", "1"); }
o_endogenous_prior : ENDOGENOUS_PRIOR { driver.option_num("endogenous_prior", "1"); }
o_use_univariate_filters_if_singularity_is_detected : USE_UNIVARIATE_FILTERS_IF_SINGULARITY_IS_DETECTED EQUAL INT_NUMBER { driver.option_num("use_univariate_filters_if_singularity_is_detected", $3); }
o_mcmc_jumping_covariance : MCMC_JUMPING_COVARIANCE EQUAL HESSIAN
                            { driver.option_str("MCMC_jumping_covariance", $3); }                                  | MCMC_JUMPING_COVARIANCE EQUAL PRIOR_VARIANCE
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          | MCMC_JUMPING_COVARIANCE EQUAL IDENTITY_MATRIX
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          | MCMC_JUMPING_COVARIANCE EQUAL filename
                            { driver.option_str("MCMC_jumping_covariance", $3); }
                          ;
o_irf_plot_threshold : IRF_PLOT_THRESHOLD EQUAL non_negative_number { driver.option_num("impulse_responses.plot_threshold", $3); };
o_consider_all_endogenous : CONSIDER_ALL_ENDOGENOUS { driver.option_str("endo_vars_for_moment_computations_in_estimation", "all_endogenous_variables"); };
o_consider_only_observed : CONSIDER_ONLY_OBSERVED { driver.option_str("endo_vars_for_moment_computations_in_estimation", "only_observed_variables"); };
o_no_homotopy : NO_HOMOTOPY { driver.option_num("no_homotopy", "1"); };

o_infile : INFILE EQUAL filename { driver.option_str("infile", $3); };
o_invars : INVARS EQUAL '(' symbol_list ')' { driver.option_symbol_list("invars"); };
o_period : PERIOD EQUAL INT_NUMBER { driver.option_num("period", $3); };
o_outfile : OUTFILE EQUAL filename { driver.option_str("outfile", $3); };
o_outvars : OUTVARS EQUAL '(' symbol_list ')' { driver.option_symbol_list("outvars"); };

range : symbol ':' symbol
        {
          $1->append(":");
          $1->append(*$3);
          delete $3;
          $$ = $1;
        };

integer_range : INT_NUMBER ':' INT_NUMBER
                {
                  $1->append(":");
                  $1->append(*$3);
                  delete $3;
                  $$ = $1;
                };
        
signed_integer_range : signed_integer ':' signed_integer
                       {
                         $1->append(":");
                         $1->append(*$3);
                         delete $3;
                         $$ = $1;
                       }
                     | MINUS '(' signed_integer ':' signed_integer ')'
                       {
                         $3->insert(0, "-(");
                         $3->append(":");
                         $3->append(*$5);
                         delete $5;
                         $3->append(")");
                         $$ = $3;
                       };

vec_int_number : INT_NUMBER { $$ = new vector<int>(); $$->push_back(atoi((*$1).c_str())); delete $1; };

vec_int_elem : vec_int_number
             | INT_NUMBER ':' INT_NUMBER
               {
                 $$ = new vector<int>();
                 for(int i=atoi((*$1).c_str()); i<=atoi((*$3).c_str()); i++)
                   $$->push_back(i);
                 delete $1;
                 delete $3;
               }
             ;

vec_int_1 : '[' vec_int_elem
            { $$ = $2;}
          | '[' COMMA vec_int_elem
            { $$ = $3;}
          | vec_int_1 vec_int_elem
            {
              $$ = $1;
              for (vector<int>::const_iterator it=$2->begin();
                   it!=$2->end(); it++)
                $1->push_back(*it);
              delete $2;
            }
          | vec_int_1 COMMA vec_int_elem
            {
              $$ = $1;
              for (vector<int>::const_iterator it=$3->begin();
                   it!=$3->end(); it++)
                $1->push_back(*it);
              delete $3;
            }
          ;

vec_int : vec_int_1 ']'
          { $$ = $1; }
        | vec_int_1 COMMA ']'
          { $$ = $1; }
        ;

vec_value_1 : '[' signed_number { $2->insert(0,"["); $$ = $2; }
            | '[' COMMA signed_number { $3->insert(0,"["); $$ = $3; }
            | vec_value_1 signed_number
              {
                $1->append(" ");
                $1->append(*$2);
                delete $2;
                $$ = $1;
              }
            | vec_value_1 COMMA signed_number
              {
                $1->append(" ");
                $1->append(*$3);
                delete $3;
                $$ = $1;
              }
            ;

vec_value : vec_value_1 ']' { $1->append("]"); $$ = $1; }
          | vec_value_1 COMMA ']' { $1->append("]"); $$ = $1; }
          ;

vec_value_list : vec_value_list COMMA vec_value
                 {
                   $1->append(",");
                   $1->append(*$3);
                   delete $3;
                   $$ = $1;
                 }
               | vec_value
                 { $$ = $1; }
               ;

vec_of_vec_value : '[' vec_value_list ']' { $$ = $2; }
                 | vec_value  { $$ = $1; };

vec_value_1_w_inf : '[' signed_number_w_inf
                    { $2->insert(0, "["); $$ = $2;}
                  | vec_value_1_w_inf signed_number_w_inf
                    {
                      $1->append(" ");
                      $1->append(*$2);
                      delete $2;
                      $$ = $1;
                    }
                  ;

vec_value_w_inf : vec_value_1_w_inf ']' { $1->append("]"); $$ = $1; };

symbol : NAME
       | ALPHA
       | BETA
       | NINV
       | ABAND
       | CMS
       | NCMS
       | CNUM
       | GAMMA
       | INV_GAMMA
       | INV_GAMMA1
       | INV_GAMMA2
       | NORMAL
       | UNIFORM
       | EPS
       | PDF
       | FIG
       | NONE
       | DR
       | PRIOR
       ;
%%

void
Dynare::parser::error(const Dynare::parser::location_type &l,
                      const string &m)
{
  driver.error(l, m);
}
