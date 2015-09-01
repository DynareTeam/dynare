<a name="logo"/>
<div align="center">
<a href="http://www.dynare.org/" target="_blank">
<img src="http://www.dynare.org/img/dynare.png" alt="Dynare Logo"></img>
</a>
</div>

# Dynare

[![Join the chat at https://gitter.im/DynareTeam/dynare](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/DynareTeam/dynare?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

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
- `gfortran`
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

We no longer support compilation on Windows. To use the unstable version of Dynare on a Windows system, please download it from the [Dynare website](http://www.dynare.org/download/dynare-unstable).

## Mac OS X

- Install the Xcode Command Line Tools:
    - Download "Command Line Tools (OS X 10.X) for Xcode," where 10.X corresponds to your OS X version, from https://developer.apple.com/downloads/index.action
- Install the latest version of [MacTeX](http://www.tug.org/mactex/)
- Install [Homebrew](http://mxcl.github.io/homebrew/) by following the instructions on the website
- Tap [Homebrew Science](https://github.com/Homebrew/homebrew-science) by opening Terminal and typing:
- ```brew tap homebrew/science```
- **(Optional)** To compile Dynare mex files for use on Octave:
    - ```brew install octave```
- To see the available options for compiling Dynare, type:
    - ```brew info dynare```
- Install Dynare via a command of the form:
    - (basic) ```brew install dynare --HEAD --without-check```
    - (with Matlab mex) ```brew install dynare --HEAD --without-check --with-matlab=/Applications/MATLAB_R2015a.app --with-matlab-version=8.5```
- **NB**: If compiling Dynare documentation, add ```--with-doc``` to the installation command
- **NB**: If not compiling Dynare mex files for Octave, add ```--without-octave``` to the installation command
- **NB**: To compile the latest stable version of dynare, follow the same instructions as above, omitting the ```--HEAD``` argument
- **NB**: To update a ```--HEAD``` install of dynare you need to uninstall it then install it again: ```brew uninstall dynare; brew install dynare --HEAD```.
- **NB**: If you want to maintain a separate git directory of dynare, you can do a ```--HEAD``` install of dynare, then uninstall it. This will have the effect of bringing in all the dependencies you will need to then compile dynare from your git directory. Then, change to the git directory and type:
    - ```autoreconf -si; ./configure --with-matlab=/Applications/MATLAB_R2015a.app MATLAB_VERSION=R2015a```, adjusting the Matlab path and version to accord with your version
- Once compilation is done, open Matlab and type the last line shown when you type ```brew info dynare``` in the Terminal window. With the typical Homebrew setup, this is:
    - ```addpath /usr/local/opt/dynare/lib/dynare/matlab```
