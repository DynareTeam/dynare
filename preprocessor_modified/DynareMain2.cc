/*
 * Copyright (C) 2008-2013 Dynare Team
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

using namespace std;

#include <iostream>

#include "ParsingDriver.hh"
#include "ModFile.hh"
#include "ConfigFile.hh"

void
main2(stringstream &in, string &basename, bool debug, bool clear_all, bool no_tmp_terms, bool no_log, bool no_warn, bool warn_uninit, bool console,
      bool parallel, const string &parallel_config_file, const string &cluster_name, bool parallel_slave_open_mode,
      bool parallel_test, int dct
#if defined(_WIN32) || defined(__CYGWIN32__)
      , bool cygwin, bool msvc
#endif
      )
{
  WarningConsolidation warnings(no_warn);

  ParsingDriver p(warnings);

  // Do parsing and construct internal representation of mod file
  ModFile *mod_file = p.parse(in, debug);
  ConfigFile config_file(parallel, parallel_test, parallel_slave_open_mode, cluster_name, dct);
  config_file.getConfigFileInfo(parallel_config_file);

  // Run checking pass
  mod_file->checkPass();
  config_file.checkPass(warnings);

  // Perform transformations on the model (creation of auxiliary vars and equations)
  mod_file->transformPass();
  config_file.transformPass();

  // Evaluate parameters initialization, initval, endval and pounds
  mod_file->evalAllExpressions(warn_uninit);

  // Do computations
  mod_file->computingPass(no_tmp_terms);

  // Write outputs
  mod_file->writeOutputFiles(basename, clear_all, no_log, no_warn, console, config_file
#if defined(_WIN32) || defined(__CYGWIN32__)
                             , cygwin, msvc
#endif
                             );

  delete mod_file;

  cout << "Preprocessing completed." << endl
       << "Starting MATLAB/Octave computing." << endl;
}
