/*
 * Copyright (C) 2003-2013 Dynare Team
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
#include <sstream>
#include <fstream>

#include <cstdlib>
#include <cstring>
#ifndef PACKAGE_VERSION
# define PACKAGE_VERSION 4.
#endif
#include "macro/MacroDriver.hh"

#include <unistd.h>

/* Prototype for second part of main function
   Splitting main() in two parts was necessary because ParsingDriver.h and MacroDriver.h can't be
   included simultaneously (because of Bison limitations).
*/
void main2(stringstream &in, string &basename, bool debug, bool clear_all, bool no_tmp_terms, bool no_log, bool no_warn, bool warn_uninit, bool console,
           bool parallel, const string &parallel_config_file, const string &cluster_name, bool parallel_slave_open_mode,
           bool parallel_test, int dct
#if defined(_WIN32) || defined(__CYGWIN32__)
           , bool cygwin, bool msvc
#endif
           );

void
usage()
{
  cerr << "Dynare usage: dynare mod_file [debug] [noclearall] [savemacro[=macro_file]] [onlymacro] [nolinemacro] [notmpterms] [nolog] [warn_uninit]"
       << " [console] [parallel[=cluster_name]] [conffile=parallel_config_path_and_filename] [parallel_slave_open_mode] [parallel_test] "
       << " [-D<variable>[=<value>]]" << endl;
  cerr << " [parallel -1] for torque-queuing, [parallel] -0 for local parallelisation, using Distributed Computing Toolbox " 
#if defined(_WIN32) || defined(__CYGWIN32__)
       << " [cygwin] [msvc]"
#endif
       << endl;
  exit(EXIT_FAILURE);
}

int
main(int argc, char **argv)
{
  /*
    Redirect stderr to stdout.
    Made necessary because MATLAB/Octave can only capture stdout (but not
    stderr), in order to put it in the logfile (see issue #306)
  */
  dup2(STDOUT_FILENO, STDERR_FILENO);

  if (argc < 2)
    {
      cerr << "Missing model file!" << endl;
      usage();
    }

  bool clear_all = true;
  bool save_macro = false;
  string save_macro_file;
  bool debug = false;
  bool no_tmp_terms = false;
  bool only_macro = false;
  bool no_line_macro = false;
  bool no_log = false;
  bool no_warn = false;
  bool warn_uninit = false;
  bool console = false;
  int  dct = 0;
#if defined(_WIN32) || defined(__CYGWIN32__)
  bool cygwin = false;
  bool msvc = false;
#endif
  string parallel_config_file;
  bool parallel = false;
  string cluster_name;
  bool parallel_slave_open_mode = false;
  bool parallel_test = false;
  map<string, string> defines;

  // Parse options
  for (int arg = 2; arg < argc; arg++)
    {
      if (!strcmp(argv[arg], "debug"))
        debug = true;
      else if (!strcmp(argv[arg], "noclearall"))
        clear_all = false;
      else if (!strcmp(argv[arg], "onlymacro"))
        only_macro = true;
      else if (strlen(argv[arg]) >= 9 && !strncmp(argv[arg], "savemacro", 9))
        {
          save_macro = true;
          if (strlen(argv[arg]) > 9)
            {
              if (strlen(argv[arg]) == 10 || argv[arg][9] != '=')
                {
                  cerr << "Incorrect syntax for savemacro option" << endl;
                  usage();
                }
              save_macro_file = string(argv[arg] + 10);
            }
        }
      else if (!strcmp(argv[arg], "nolinemacro"))
        no_line_macro = true;
      else if (!strcmp(argv[arg], "notmpterms"))
        no_tmp_terms = true;
      else if (!strcmp(argv[arg], "nolog"))
        no_log = true;
      else if (!strcmp(argv[arg], "nowarn"))
        no_warn = true;
      else if (!strcmp(argv[arg], "warn_uninit"))
        warn_uninit = true;
      else if (!strcmp(argv[arg], "console"))
        console = true;
#if defined(_WIN32) || defined(__CYGWIN32__)
      else if (!strcmp(argv[arg], "cygwin"))
        cygwin = true;
      else if (!strcmp(argv[arg], "msvc"))
        msvc = true;
#endif
      else if (strlen(argv[arg]) >= 8 && !strncmp(argv[arg], "conffile", 8))
        {
          if (strlen(argv[arg]) <= 9 || argv[arg][8] != '=')
            {
              cerr << "Incorrect syntax for conffile option" << endl;
              usage();
            }
          parallel_config_file = string(argv[arg] + 9);
        }
      else if (!strcmp(argv[arg], "parallel_slave_open_mode"))
        parallel_slave_open_mode = true;
      else if (!strcmp(argv[arg], "parallel_test"))
        parallel_test = true;
      else if (strlen(argv[arg]) >= 8 && !strncmp(argv[arg], "parallel", 8))
        {
          parallel = true;
          if (strlen(argv[arg]) > 8)
            {
              if (strlen(argv[arg]) == 9 || argv[arg][8] != '=')
                {
                  cerr << "Incorrect syntax for parallel option" << endl;
                  usage();
                }
              cluster_name = string(argv[arg] + 9);
            }
        }
      else if (strlen(argv[arg]) >= 2 && !strncmp(argv[arg], "-D", 2))
        {
          if (strlen(argv[arg]) == 2)
            {
              cerr << "Incorrect syntax for command line define: the defined variable "
                   << "must not be separated from -D by whitespace." << endl;
              usage();
            }

          size_t equal_index = string(argv[arg]).find('=');
          if (equal_index != string::npos)
            {
              string key = string(argv[arg]).erase(equal_index).erase(0,2);
              defines[key] = string(argv[arg]).erase(0, equal_index+1);
            }
          else
            {
              string key = string(argv[arg]).erase(0,2);
              defines[key] = "1";
            }
        }
      else if (strlen(argv[arg]) == 2 && !strncmp(argv[arg], "-1", 2))
        {
           if (parallel)
               {
                dct=1;
               }
	   else
	       {
               usage();
	       }     
        }
      else if (strlen(argv[arg]) == 2 && !strncmp(argv[arg], "-2", 2))
        {
           if (parallel)
               {
                dct=2;
               }
           else
	       {
               usage();
	       }  
        }
      else
        {
          cerr << "Unknown option: " << argv[arg] << endl;
          usage();
        }
    }

  cout << "Starting Dynare (version " << PACKAGE_VERSION << ")." << endl
       << "Starting preprocessing of the model file ..." << endl;

  // Construct basename (i.e. remove file extension if there is one)
  string basename = argv[1];
  size_t pos = basename.find_last_of('.');
  if (pos != string::npos)
    basename.erase(pos);

  // Do macro processing
  MacroDriver m;

  stringstream macro_output;
  m.parse(argv[1], macro_output, debug, no_line_macro, defines);
  if (save_macro)
    {
      if (save_macro_file.empty())
        save_macro_file = basename + "-macroexp.mod";
      ofstream macro_output_file(save_macro_file.c_str());
      if (macro_output_file.fail())
        {
          cerr << "Cannot open " << save_macro_file << " for macro output" << endl;
          exit(EXIT_FAILURE);
        }
      macro_output_file << macro_output.str();
      macro_output_file.close();
    }

  if (only_macro)
    return EXIT_SUCCESS;

  // Do the rest
  main2(macro_output, basename, debug, clear_all, no_tmp_terms, no_log, no_warn, warn_uninit, console,
        parallel, parallel_config_file, cluster_name, parallel_slave_open_mode, parallel_test, dct
#if defined(_WIN32) || defined(__CYGWIN32__)
        , cygwin, msvc
#endif
        );

  return EXIT_SUCCESS;
}
