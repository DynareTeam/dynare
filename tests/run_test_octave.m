## Copyright (C) 2009-2016 Dynare Team
##
## This file is part of Dynare.
##
## Dynare is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## Dynare is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

## Implementation notes:
##
## Before every call to Dynare, the contents of the workspace is saved in
## 'wsOct', and reloaded after Dynare has finished (this is necessary since
## Dynare does a 'clear -all').

load_octave_packages

top_test_dir = getenv('TOP_TEST_DIR');
addpath(top_test_dir);
addpath('.'); % workaround for bug in Octave 4.2.0-rc2: https://savannah.gnu.org/bugs/?49334
addpath([top_test_dir filesep '..' filesep 'matlab']);

## Test Dynare Version
if !strcmp(dynare_version(), getenv("DYNARE_VERSION"))
    error("Incorrect version of Dynare is being tested")
endif

## Ask gnuplot to create graphics in text mode
graphics_toolkit gnuplot;
setenv("GNUTERM", "dumb");

## Test MOD files listed in Makefile.am
name = getenv("FILESTEM");
[directory, testfile, ext] = fileparts([top_test_dir '/' name]);
cd(directory);

printf("\n***  TESTING: %s ***\n", name);

cput = cputime;
save(['wsOct' testfile '.mat']);
try
  dynare([testfile ext])
  testFailed = false;
catch
  printMakeCheckOctaveErrMsg(getenv("FILESTEM"), lasterror);
  testFailed = true;
end_try_catch
top_test_dir = getenv('TOP_TEST_DIR');
name = getenv("FILESTEM");
[directory, testfile, ext] = fileparts([top_test_dir '/' name]);
load(['wsOct' testfile '.mat']);
ecput = cputime - cput;
delete(['wsOct' testfile '.mat']);

cd(top_test_dir);
fid = fopen([name '.o.trs'], 'w+');
if testFailed
  fprintf(fid,':test-result: FAIL\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 1\n');
  fprintf(fid,':list-of-failed-tests: %s\n', [name '.mod']);
else
  fprintf(fid,':test-result: PASS\n');
  fprintf(fid,':number-tests: 1\n');
  fprintf(fid,':number-failed-tests: 0\n');
  fprintf(fid,':list-of-passed-tests: %s\n', [name '.mod']);
end
fprintf(fid,':cputime: %f\n', ecput);
fclose(fid);

## Local variables:
## mode: Octave
## End:
