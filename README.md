# Dynare

Described on the homepage: <http://www.dynare.org/>

Most users should use the precompiled package available for your OS, also
available via the Dynare homepage: <http://www.dynare.org/download/dynare-stable>.

# License

Most of the source files are covered by the GNU General Public Licence version
3 or later (there are some exceptions to this, see [license.txt](license.txt) in
Dynare distribution for specifics).

# Building Dynare From Source

Here, we explain how to build from source:
- Dynare, including preprocessor and MEX files for MATLAB and Octave
- Dynare++
- all the associated documentation (PDF and HTML)

This source can be retrieved in three forms:
- via git, at <https://github.com/DynareTeam/dynare.git>
- using the stable source archive of the latest Dynare version (currently 4.4) from <http://www.dynare.org/download/dynare-stable/>
- using a source snapshot of the unstable version, from <http://www.dynare.org/download/dynare-unstable/source-snapshot>

Note that if you obtain the source code via git, you will need to install more tools (see below).

The first section of this page gives general instructions, which apply to all platforms. Then some specific platforms are discussed.

**NB**: Here, when we refer to 32-bit or 64-bit, we refer to the type of MATLAB installation, not the type of Windows installation. It is perfectly possible to run a 32-bit MATLAB on a 64-bit Windows: in that case, instructions for Windows 32-bit should be followed. To determine the type of your MATLAB installation, type:
```matlab
>> computer
```
at the MATLAB prompt: if it returns `PCWIN`, then you have a 32-bit MATLAB; if it returns `PCWIN64`, then you have a 64-bit MATLAB.

**Contents**

1. [**General Instructions**](#general-instructions)
1. [**Debian or Ubuntu**](#debian-or-ubuntu)
1. [**Fedora**](#fedora)
1. [**Windows**](#windows)
1. [**Mac OS X**](#mac-os-x)

## General Instructions

### Prerequisites

A number of tools and libraries are needed in order to recompile everything. You don't necessarily need to install everything, depending on what you want to compile.

- A POSIX compliant shell and an implementation of Make (mandatory)
- The [GNU Compiler Collection](http://gcc.gnu.org/), with gcc, g++ and gfortran (mandatory)
- MATLAB (if you want to compile MEX for MATLAB)
- [GNU Octave](http://www.octave.org), with the development headers (if you want to compile MEX for Octave)
- [Boost libraries](http://www.boost.org), version 1.36 or later
- [Bison](http://www.gnu.org/software/bison/), version 2.5 or later (only if you get the source through Git)
- [Flex](http://flex.sourceforge.net/), version 2.5.4 or later (only if you get the source through Git)
- [Autoconf](http://www.gnu.org/software/autoconf/), version 2.62 or later (only if you get the source through Git) (see [Installing an updated version of Autoconf in your own directory, in GNU/Linux](http://www.dynare.org/DynareWiki/AutoMake))
- [Automake](http://www.gnu.org/software/automake/), version 1.11.2 or later (only if you get the source through Git) (see [Installing an updated version of AutoMake in your own directory, in GNU/Linux](http://www.dynare.org/DynareWiki/AutoMake))
- [CWEB](http://www-cs-faculty.stanford.edu/%7Eknuth/cweb.html), with its tools `ctangle` and `cweave` (only if you want to build Dynare++ and get the source through Git)
- An implementation of BLAS and LAPACK: either [ATLAS](http://math-atlas.sourceforge.net/), [OpenBLAS](http://xianyi.github.com/OpenBLAS/), Netlib ([BLAS](http://www.netlib.org/blas/), [LAPACK](http://www.netlib.org/lapack/)) or [MKL](http://software.intel.com/en-us/intel-mkl/) (only if you want to build Dynare++)
- An implementation of [POSIX Threads](http://en.wikipedia.org/wiki/POSIX_Threads) (optional, for taking advantage of multi-core)
- [MAT File I/O library](http://sourceforge.net/projects/matio/) (if you want to compile Markov-Switching code, the estimation DLL, k-order DLL and Dynare++)
- [SLICOT](http://www.slicot.org) (if you want to compile the Kalman steady state DLL)
- [GSL library](http://www.gnu.org/software/gsl/) (if you want to compile Markov-Switching code)
- A decent LaTeX distribution (if you want to compile PDF documentation). The following extra components may be needed:
  - [Eplain](http://www.tug.org/eplain/) TeX macros (only if you want to build Dynare++ source documentation)
  - [Beamer](http://latex-beamer.sourceforge.net/) (for some PDF presentations)
- For building the reference manual:
  - [GNU Texinfo](http://www.gnu.org/software/texinfo/)
  - [Texi2HTML](http://www.nongnu.org/texi2html) and [Latex2HTML](http://www.latex2html.org), if you want nice mathematical formulas in HTML output
  - [Doxygen](http://www.stack.nl/%7Edimitri/doxygen/) (if you want to build Dynare preprocessor source documentation)
- For Octave, the development libraries corresponding to the UMFPACK packaged with Octave

### Preparing the sources

If you have downloaded the sources from an official source archive or the source snapshot, just unpack it.

If you want to use Git, do the following from a terminal:

    git clone --recursive http://github.com/DynareTeam/dynare.git
    cd dynare
    autoreconf -si

The last line runs Autoconf and Automake in order to prepare the build environment (this is not necessary if you got the sources from an official source archive or the source snapshot).

### Configuring the build tree

Simply launch the configure script from a terminal:
```
./configure
```
If you have MATLAB, you need to indicate both the MATLAB location and version. For example, on GNU/Linux:
```
./configure --with-matlab=/usr/local/MATLAB/R2013a MATLAB_VERSION=8.1
```
Note that the MATLAB version can also be specified via the MATLAB family product release (R2009a, R2008b, ...).

**NB**: For MATLAB versions strictly older than 7.1, you need to explicitly give the MEX extension, via `MEXEXT` variable of the configure script (for example, `MEXEXT=dll` for Windows with MATLAB \< 7.1).

Alternatively, you can disable the compilation of MEX files for MATLAB with the `--disable-matlab` flag, and MEX files for Octave with `--disable-octave`.

You may need to specify additional options to the configure script, see the platform specific instructions below.

Note that if you don't want to compile with debugging information, you can specify the `CFLAGS` and `CXXFLAGS` variables to configure, such as:
```
./configure CFLAGS="-O3" CXXFLAGS="-O3"
```
If you want to give a try to the parallelized versions of some mex files (`A_times_B_kronecker_C` and `sparse_hessian_times_B_kronecker_C` used to get the reduced form of the second order approximation of the model) you can add the `--enable-openmp` flag, for instance:
```
./configure --with-matlab=/usr/local/matlab78 MATLAB_VERSION=7.8 --enable-openmp
```
If the configuration goes well, the script will tell you which components are correctly configured and will be built.

### Bulding

Binaries and Info documentation are built with:
```
make
```
PDF and HTML documentation are respectively built with:
```
make pdf
make html
```
The testsuites can be run with:
```
make check
```
## Debian or Ubuntu

All the prerequisites are packaged.

The easiest way to install the pre-requisites in Debian is to use Debian's dynare package and do:
```
apt-get build-dep dynare
```
Alternatively, if you want to build everything, manually install the following packages:

- `build-essential` (for gcc, g++ and make)
- `liboctave-dev` or `octave3.2-headers` (will install ATLAS)
- `libboost-graph-dev`
- `libgsl0-dev`
- `libmatio-dev`
- `libslicot-dev` and `libslicot-pic`
- `libsuitesparse-dev`
- `flex`
- `bison`
- `autoconf`
- `automake`
- `texlive`
- `texlive-publishers` (for Econometrica bibliographic style)
- `texlive-extra-utils` (for CWEB)
- `texlive-formats-extra` (for Eplain)
- `texlive-latex-extra` (for fullpage.sty)
- `latex-beamer`
- `texinfo`
- `texi2html`, `latex2html`
- `doxygen`

## Fedora

**NB**: Documentation still in progress…
- `octave-devel`
- `boost-devel`
- `gsl-devel`
- `matio-devel`
- `flex`
- `bison`
- `autoconf`
- `automake`
- `texlive`
- `texinfo`
- `texi2html`, `latex2html`
- `doxygen`

## Windows

The following instructions are compatible with MATLAB or with Octave/MinGW (as downloadable [here](http://www.dynare.org/download/octave)).

### Setting up the Compilation Environment

- First, you need to setup a Cygwin environment, following the instructions at <http://www.cygwin.com>. You can install either the 32-bit or the 64-bit version. If you opt for the latter, you need to replace `c:\cygwin` by `c:\cygwin64` in the following.
- Install the following packages:
    - `make`
    - `bison`
    - `flex`
    - `git`, `gettext`
    - `autoconf`
    - `automake`
    - `texi2html`
    - `texlive`, `texlive-collection-latexextra`, `texlive-collection-formatsextra`, `texlive-collection-publishers`, `texlive-collection-fontsrecommended`, `texlive-collection-fontsextra`, `texlive-collection-bibtexextra`, `texlive-collection-genericrecommended`, `texlive-collection-mathextra`, `texlive-collection-binextra`
    - `texinfo`
    - `doxygen`
    - `mingw64-i686-gcc-core`, `mingw64-i686-gcc-g++`, `mingw64-i686-gcc-fortran` (if you have Octave/MinGW or if you have MATLAB 32-bit)
    - `mingw64-x86_64-gcc-core`, `mingw64-x86_64-gcc-g++`, `mingw64-x86_64-gcc-fortran` (if you have MATLAB 64-bit)
- Second, install precompiled librairies for BLAS, LAPACK, Boost and GSL:
    - If you have Octave or MATLAB 32-bit, download [dynare-mingw32-libs.zip](http://www.dynare.org/build/dynare-mingw32-libs.zip), and uncompress it in `c:\cygwin\usr\local\lib\mingw32`
    - If you have MATLAB 64-bit, download [dynare-mingw64-libs.zip](http://www.dynare.org/build/dynare-mingw64-libs.zip), and uncompress it in `c:\cygwin\usr\local\lib\mingw64`

*Remark*: You need to make sure that Cygwin’s git is used and not a potentially installed msysgit. The latter typically happens when one installs msysgit and allows it to set a system path. This will result in wrong line endings and cryptic error messages à la "syntax error near unexpected token `fi'". In that case it might be necessary to uninstall msysgit and reinstall it without setting a system path.

### Compiling the preprocessor, Dynare++, the MEX for MATLAB and the documentation

Install the Dynare source tree, let’s say in `c:\cygwin\home\user\dynare` (by either uncompressing a source archive, or by cloning the git repository with `git clone --recursive http://github.com/DynareTeam/dynare.git`).

Launch a Cygwin shell, and enter the Dynare source tree:
```
cd dynare
```
If you retrieved the source from Git, don't forget to do:
```
autoreconf -i -s
```
Then, configure the package.

- If your MATLAB is 32-bit, let's say version R2008b installed in `c:\Program Files\MATLAB\R2008b`
```
./configure --host=i686-w64-mingw32 --with-boost=/usr/local/lib/mingw32/boost --with-blas=/usr/local/lib/mingw32/blas/libopenblas.a --with-lapack=/usr/local/lib/mingw32/lapack/liblapack.a --with-gsl=/usr/local/lib/mingw32/gsl --with-matio=/usr/local/lib/mingw32/matio --with-slicot=/usr/local/lib/mingw32/slicot --with-matlab=/cygdrive/c/Progra~1/MATLAB/R2008b MATLAB_VERSION=R2008b --disable-octave
```
- If your MATLAB is 64-bit:
```
./configure --host=x86_64-w64-mingw32 --with-boost=/usr/local/lib/mingw64/boost --with-blas=/usr/local/lib/mingw64/blas/libopenblas.a --with-lapack=/usr/local/lib/mingw64/lapack/liblapack.a --with-gsl=/usr/local/lib/mingw64/gsl --with-matio=/usr/local/lib/mingw64/matio --with-slicot=/usr/local/lib/mingw64/slicot --with-matlab=/cygdrive/c/Progra~1/MATLAB/R2008b MATLAB_VERSION=R2008b --disable-octave
```
A few remarks:

- Note that here we use `Progra~1` (the 8.3 filename) instead of `Program Files`. This is because spaces in filenames confuse the configuration scripts.
- If you don’t have MATLAB, then drop the `--with-matlab` and `MATLAB_VERSION` options
- If your MATLAB is 32-bit and your Windows is 64-bit, you need to explicitly give the MEX extension, with `MEXEXT=mexw32`

Then compile everything with:
```
make all pdf html
```
This should build:

- Dynare preprocessor
- Dynare MEX files for MATLAB (provided you gave the MATLAB path to configure)
- Dynare++
- Part of the documentation

### Compiling the MEX for Octave (MinGW package)

Launch a Cygwin shell, and enter the Dynare source tree for Octave MEX:

    cd dynare/mex/build/octave

Configure and make:

    ./configure MKOCTFILE=/usr/local/lib/mingw32/mkoctfile-win --with-boost=/usr/local/lib/mingw32/boost --with-gsl=/usr/local/lib/mingw32/gsl --with-matio=/usr/local/lib/mingw32/matio --with-slicot=/usr/local/lib/mingw32/slicot-underscore
    make

## Mac OS X

- Install the Xcode Command Line Tools:
    - Download "Command Line Tools (OS X 10.X) for Xcode," where 10.X corresponds to your OS X version, from https://developer.apple.com/downloads/index.action
- Install [Homebrew](http://mxcl.github.io/homebrew/) by following the instructions on the website
- Tap [Homebrew Science](https://github.com/Homebrew/homebrew-science) by opening Terminal and typing:
    - ```brew tap homebrew/science```
- Install the following brews:
    - ```brew install gcc```
    - ```brew install automake```
    - ```brew install gsl```
    - ```brew install bison```
    - ```brew install boost```
    - ```brew install libmatio --with-hdf5```
    - ```brew install slicot --with-default-integer-8```
- Force link ```flex``` and ```bison```
    - ```brew link --force flex```
    - ```brew link --force bison```
- **(Optional)** To compile Dynare mex files for use on Octave:
    - ```brew install octave```
- **(Optional)** To compile Dynare's documentation, first install the latest version of [MacTeX](http://www.tug.org/mactex/). Then install `doxygen`, `latex2html` and `texi2html` via Homebrew with the following commands:
    - ```brew install doxygen```
    - ```brew install texinfo```
    - ```brew install latex2html```
    - ```brew install texi2html```
- **(On OS X 10.7 Only)** Copy [FlexLexer.h](http://www.dynare.org/build/FlexLexer.h) into the `preprocessor` directory (there was an error in the `FlexLexer.h` file distributed with 10.7)
- Ensure `/usr/local/bin` is at the beginning of your path:
    - ```PATH=/usr/local/bin:$PATH```
- Finally, switch to the root dynare directory and compile dynare
    - `autoreconf -si`
    - `./configure --with-matlab=/Applications/MATLAB_R2014b.app MATLAB_VERSION=8.4`
    - `make`
    - `make pdf TEXI2DVI=/usr/local/Cellar/texinfo/5.2/bin/texi2dvi`, where you replace everything after the equal sign with the path to the `texi2dvi` installed by homebrew when you installed `texinfo`.
