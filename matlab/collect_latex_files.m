function collect_latex_files
% function collect_LaTeX_Files;
% Creates TeX-File embedding all eps-loaders created for current mod-file
% 
% Inputs: none
% 
% Notes: 
%   - The packages loaded enable pdflatex to run
%   - The _dynamic and _static TeX-model files are not included as they are standalone TeX-files

% Copyright (C) 2015-16 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

global M_
%% Write header
f_name_binder=[M_.fname,'_TeX_binder.tex'];
fid=fopen(f_name_binder,'w+');
fprintf(fid,'%s \n','\documentclass[12pt]{article}');
fprintf(fid,'%s \n','\usepackage[margin=2cm]{geometry}');
fprintf(fid,'%s \n','\usepackage{psfrag}');
fprintf(fid,'%s \n','\usepackage{graphicx}');
fprintf(fid,'%s \n','\usepackage{epstopdf}');
fprintf(fid,'%s \n','\usepackage{longtable,booktabs}');
fprintf(fid,'%s \n','\usepackage{amsmath,amsfonts}');
fprintf(fid,'%s \n','\usepackage{breqn}');
fprintf(fid,'%s \n','\usepackage{float,morefloats,caption}');
fprintf(fid,'%s \n','\begin{document}');

%% Root directory
TeX_Files=dir([M_.fname,'*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_dynamic.tex']) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_static.tex']) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_original.tex'])
        fprintf(fid,'%s \n',['\include{',f_name,'}']);    
    end
end

%% Output directory
TeX_Files=dir([M_.dname filesep 'Output' filesep  M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/Output' '/',f_name,'}']);    
    end
end

%% graphs directory
TeX_Files=dir([M_.dname filesep 'graphs' filesep  M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/graphs' '/',f_name,'}']);    
    end
end

%% Identification directory
TeX_Files=dir([M_.dname filesep 'identification' filesep  M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/identification' '/',f_name,'}']);    
    end
end


%% Identification/Output directory
TeX_Files=dir([M_.dname filesep 'identification' filesep 'Output' filesep M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/identification/Output' '/',f_name,'}']);    
    end
end

%% GSA directory
TeX_Files=dir([M_.dname filesep 'gsa' filesep  M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/gsa' '/',f_name,'}']);    
    end
end

%% GSA/Output directory
TeX_Files=dir([M_.dname filesep 'gsa' filesep 'Output' filesep  M_.fname '*.tex']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/gsa/Output' '/',f_name,'}']);    
    end
end


dirinfo_parent = dir([M_.dname filesep 'gsa' filesep 'redform*']);
dirinfo_parent(~[dirinfo_parent.isdir]) = [];  %remove non-directories
tf = ismember( {dirinfo_parent.name}, {'.', '..'});
dirinfo_parent(tf) = [];  %remove current and parent directory.
numsubdir_level1 = length(dirinfo_parent);
for level1_iter = 1:numsubdir_level1
    dirinfo_subfolder = dir([M_.dname filesep 'gsa' filesep dirinfo_parent(level1_iter).name]);
    dirinfo_subfolder(~[dirinfo_subfolder.isdir]) = [];  %remove non-directories
    tf = ismember( {dirinfo_subfolder.name}, {'.', '..'});
    dirinfo_subfolder(tf) = [];  %remove current and parent directory.
    numsubdir_level2 = length(dirinfo_subfolder);
    for level2_iter = 1:numsubdir_level2
        TeX_Files=dir([M_.dname filesep 'gsa' filesep dirinfo_parent(level1_iter).name filesep  dirinfo_subfolder(level2_iter).name filesep M_.fname '*.tex']);
        for ii=1:length(TeX_Files)
            [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
            if ~strcmp(TeX_Files(ii).name,f_name_binder)
                fprintf(fid,'%s \n',['\include{', M_.dname '/gsa/',dirinfo_parent(level1_iter).name '/'  dirinfo_subfolder(level2_iter).name ,'/',f_name,'}']);
            end
        end
        TeX_Files=dir([M_.dname filesep 'gsa' filesep dirinfo_parent(level1_iter).name filesep  dirinfo_subfolder(level2_iter).name filesep 'Output' filesep  M_.fname '*.tex']);
        for ii=1:length(TeX_Files)
            [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
            if ~strcmp(TeX_Files(ii).name,f_name_binder)
                fprintf(fid,'%s \n',['\include{', M_.dname '/gsa/', dirinfo_parent(level1_iter).name '/'  dirinfo_subfolder(level2_iter).name, '/Output' '/',f_name,'}']);
            end
        end
    end
end




%% Write footer
fprintf(fid,'%s \n','\end{document}');

fclose(fid);
