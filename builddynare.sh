PATH=/usr/local/bin:/usr/bin:/usr/lib/lapack
rm apt-cyg
lynx -source rawgit.com/transcode-open/apt-cyg/master/apt-cyg > apt-cyg
install apt-cyg /usr/local/bin
apt-cyg install nano zip unzip wget make libtool bison flex git gettext autoconf automake texi2html texlive texlive-collection-latexextra texlive-collection-formatsextra texlive-collection-publishers texlive-collection-fontsrecommended texlive-collection-fontsextra texlive-collection-bibtexextra texlive-collection-genericrecommended texlive-collection-mathextra texlive-collection-binextra texinfo doxygen mingw64-i686-gcc-core mingw64-i686-gcc-g++ mingw64-i686-gcc-fortran mingw64-x86_64-gcc-core mingw64-x86_64-gcc-g++ mingw64-x86_64-gcc-fortran
cd /cygdrive/c/cygwin64/
wget -N http://www.cygwin.com/setup-x86_64.exe
chmod +x ./setup-x86_64.exe
./setup-x86_64.exe --no-desktop --no-shortcuts --no-startmenu --quiet-mode
mkdir -p /usr/local/lib/mingw64
cd /usr/local/lib/mingw64
wget http://www.dynare.org/build/dynare-mingw64-libs.zip
unzip dynare-mingw64-libs.zip
rm dynare-mingw64-libs.zip
cd ~
git clone http://git.savannah.gnu.org/r/gsl.git
cd ~/gsl
git pull
autoreconf -si
./configure CFLAGS="-O3" CXXFLAGS="-O3" --host=x86_64-w64-mingw32 --prefix=/usr/local/lib/mingw64/gsl
make
make install
cd /cygdrive/c/dynare/dynare/
git fetch upstream --recurse-submodules
git checkout master
git merge upstream/master
git commit -m "merge"
git push
autoreconf -si
./configure CFLAGS="-O3 -DCUDA=1" CXXFLAGS="-O3 -DCUDA=1" --host=x86_64-w64-mingw32 --enable-openmp --with-boost=/usr/local/lib/mingw64/boost --with-blas=/usr/local/lib/mingw64/blas/libopenblas.a --with-lapack=/usr/local/lib/mingw64/lapack/liblapack.a --with-gsl=/usr/local/lib/mingw64/gsl --with-matio=/usr/local/lib/mingw64/matio --with-slicot=/usr/local/lib/mingw64/slicot --with-matlab=/cygdrive/c/Progra~1/MATLAB/R2015a MATLAB_VERSION=R2015a --disable-octave
make all html
