/*
 * Copyright (C) 2008-2017 Dynare Team
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

#include <iostream>

#include "ParsingDriver.hh"
#include "ModFile.hh"
#include "ConfigFile.hh"
#include "ExtendedPreprocessorTypes.hh"

void
main2(stringstream &in, string &basename, bool debug, bool clear_all, bool clear_global,
      bool no_tmp_terms, bool no_log, bool no_warn, bool warn_uninit, bool console,
      bool nograph, bool nointeractive, bool parallel, ConfigFile &config_file,
      WarningConsolidation &warnings, bool nostrict, bool check_model_changes,
      bool minimal_workspace, bool compute_xrefs, FileOutputType output_mode,
      LanguageOutputType language, int params_derivs_order
#if defined(_WIN32) || defined(__CYGWIN32__) || defined(__MINGW32__)
      , bool cygwin, bool msvc, bool mingw
#endif
      )
{
  ParsingDriver p(warnings, nostrict);

  // Do parsing and construct internal representation of mod file
  ModFile *mod_file = p.parse(in, debug);

  // Run checking pass
  mod_file->checkPass(nostrict);

  // Perform transformations on the model (creation of auxiliary vars and equations)
  mod_file->transformPass(nostrict);

  // Evaluate parameters initialization, initval, endval and pounds
  mod_file->evalAllExpressions(warn_uninit);

  // Do computations
  mod_file->computingPass(no_tmp_terms, output_mode, compute_xrefs, params_derivs_order);

  // Write outputs
  if (output_mode != none)
    mod_file->writeExternalFiles(basename, output_mode, language);
  else
    mod_file->writeOutputFiles(basename, clear_all, clear_global, no_log, no_warn, console, nograph,
                               nointeractive, config_file, check_model_changes, minimal_workspace, compute_xrefs
#if defined(_WIN32) || defined(__CYGWIN32__) || defined(__MINGW32__)
			       , cygwin, msvc, mingw
#endif
			       );

  delete mod_file;

  cout << "Preprocessing completed." << endl;
}
