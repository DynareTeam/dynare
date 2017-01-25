//Test mh_recover function for RW-MH when load_mh_file was also used
//Previous chain did fill last mh_file so mh_recover does not need to fill that file

@#include "fs2000.common.inc"

options_.MaxNumberOfBytes=2000*11*8/4;
estimation(order=1, datafile='../fsdat_simul',nobs=192, loglinear, mh_replic=1000, mh_nblocks=2, mh_jscale=0.8);
estimation(order=1,mode_compute=0,mode_file=fs2000_recover_3_mode, datafile='../fsdat_simul',nobs=192, loglinear, load_mh_file,mh_replic=1000, mh_nblocks=2, mh_jscale=0.8);
copyfile([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat'],[M_.dname '_mh1_blck1.mat'])
copyfile([M_.dname filesep 'metropolis' filesep M_.dname '_mh3_blck2.mat'],[M_.dname '_mh3_blck2.mat'])
copyfile([M_.dname filesep 'metropolis' filesep M_.dname '_mh4_blck2.mat'],[M_.dname '_mh4_blck2.mat'])
delete([M_.dname filesep 'metropolis' filesep M_.dname '_mh3_blck2.mat'])
delete([M_.dname filesep 'metropolis' filesep M_.dname '_mh4_blck2.mat'])

estimation(order=1, datafile='../fsdat_simul',mode_compute=0,mode_file=fs2000_recover_3_mode, nobs=192, loglinear, mh_replic=2000, mh_nblocks=2, mh_jscale=0.8,mh_recover);

%check first unaffected chain
temp1=load([M_.dname '_mh1_blck1.mat']);
temp2=load([M_.dname filesep 'metropolis' filesep M_.dname '_mh1_blck1.mat']);

if max(max(abs(temp1.x2-temp2.x2)))>1e-10
    error('Draws of unaffected chain are not the same')
end

%check second, affected chain with last unaffected file
temp1=load([M_.dname '_mh3_blck2.mat']);
temp2=load([M_.dname filesep 'metropolis' filesep M_.dname '_mh3_blck2.mat']);

if max(max(abs(temp1.x2-temp2.x2)))>1e-10
    error('Draws of affected chain''s unaffected files are not the same')
end

if ~isoctave
    %check second, affected chain with affected file
    temp1=load([M_.dname '_mh4_blck2.mat']);
    temp2=load([M_.dname filesep 'metropolis' filesep M_.dname '_mh4_blck2.mat']);

    if max(max(abs(temp1.x2-temp2.x2)))>1e-10
        error('Draws of affected chain''s affected files are not the same')
    end
end