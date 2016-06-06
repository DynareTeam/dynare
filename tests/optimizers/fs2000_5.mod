@#include "fs2000.common.inc"

estimation(mode_compute=5,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0);
estimation(mode_compute=5,order=1, datafile='../fs2000/fsdat_simul', nobs=192, mh_replic=0, optim=('Hessian',2));
