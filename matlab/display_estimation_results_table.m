function oo_=display_estimation_results_table(xparam1,stdh,M_,options_,estim_params_,bayestopt_,oo_,pnames,table_title,field_name)
%function oo_=display_results_table(xparam1,stdh,M_,estim_params_,bayestopt_,oo_,pnames,table_title,field_name)
% Display estimation results on screen and write them to TeX-file
% 
% INPUTS 
%   o xparam1       [double]   (p*1) vector of estimate parameters.
%   o stdh          [double]   (p*1) vector of estimate parameters.
%   o M_                        Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).          
%   o estim_params_             Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   o options_                  Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   o bayestopt_                Matlab's structure describing the priors (initialized by dynare, see @ref{bayesopt_}).
%   o oo_                       Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%   o pnames        [string]    Character Array storing the names for prior distributions     
%   o table_title   [string]    Title of the Table     
%   o field_name    [string]    String storing the name of the fields for oo_ where the parameters are stored
%  
% OUTPUTS 
%   o oo_                       Matlab's structure gathering the results
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2013 Dynare Team
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

nvx = estim_params_.nvx;  % Variance of the structural innovations (number of parameters).
nvn = estim_params_.nvn;  % Variance of the measurement innovations (number of parameters).
ncx = estim_params_.ncx;  % Covariance of the structural innovations (number of parameters).
ncn = estim_params_.ncn;  % Covariance of the measurement innovations (number of parameters).
np  = estim_params_.np ;  % Number of deep parameters.
nx  = nvx+nvn+ncx+ncn+np; % Total number of parameters to be estimated.

disp(' ')
disp(['RESULTS FROM ' upper(table_title) ' ESTIMATION'])
LaTeXtitle=strrep(table_title,' ','_');
tstath = abs(xparam1)./stdh;

header_width = row_header_width(M_,estim_params_,bayestopt_);
if strcmp(field_name,'posterior')
    tit1 = sprintf('%-*s %7s %8s %7s %4s %6s\n',header_width-2,' ','prior mean', ...
                   'mode','s.d.','prior','pstdev');
else
    tit1 = sprintf('%-*s %10s %7s %6s\n',header_width-2,' ','Estimate','s.d.','t-stat');
end
if np
    ip = nvx+nvn+ncx+ncn+1;
    disp('parameters')
    disp(tit1)
    for i=1:np
        name = bayestopt_.name{ip};
        if strcmp(field_name,'posterior')
            fprintf('%-*s %7.3f %8.4f %7.4f %4s %6.4f \n', ...
                     header_width,name, ...
                     bayestopt_.p1(ip),xparam1(ip),stdh(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.p2(ip));
        else
            fprintf('%-*s %8.4f %7.4f %7.4f \n', ...
                 header_width,name,xparam1(ip),stdh(ip),tstath(ip));
        end
        eval(['oo_.' field_name '_mode.parameters.' name ' = xparam1(ip);']);
        eval(['oo_.' field_name '_std.parameters.' name ' = stdh(ip);']);
        ip = ip+1;
    end
    disp(' ')
end
if nvx
    ip = 1;
    disp('standard deviation of shocks')
    disp(tit1)
    for i=1:nvx
        k = estim_params_.var_exo(i,1);
        name = deblank(M_.exo_names(k,:));
        if strcmp(field_name,'posterior')
            fprintf('%-*s %7.3f %8.4f %7.4f %4s %6.4f \n', ...
                     header_width,name,bayestopt_.p1(ip),xparam1(ip), ...
                     stdh(ip),pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.p2(ip));
        else
            fprintf('%-*s %8.4f %7.4f %7.4f \n',header_width,name,xparam1(ip),stdh(ip),tstath(ip));
        end
        M_.Sigma_e(k,k) = xparam1(ip)*xparam1(ip);
        eval(['oo_.' field_name '_mode.shocks_std.' name ' = xparam1(ip);']);
        eval(['oo_.' field_name '_std.shocks_std.' name ' = stdh(ip);']);
        ip = ip+1;
    end
    disp(' ')
 end
 if nvn
    disp('standard deviation of measurement errors')
    disp(tit1)
    ip = nvx+1;
    for i=1:nvn
        name = deblank(options_.varobs(estim_params_.nvn_observable_correspondence(i,1),:));
        if strcmp(field_name,'posterior')           
            fprintf('%-*s %7.3f %8.4f %7.4f %4s %6.4f \n', ...
                     header_width,name,bayestopt_.p1(ip), ...
                     xparam1(ip),stdh(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), ...
                     bayestopt_.p2(ip));
        else
            fprintf('%-*s %8.4f %7.4f %7.4f \n',header_width,name,xparam1(ip),stdh(ip),tstath(ip))            
        end
        eval(['oo_.' field_name '_mode.measurement_errors_std.' name ' = xparam1(ip);']);
        eval(['oo_.' field_name '_std.measurement_errors_std.' name ' = stdh(ip);']);
        ip = ip+1;
    end
    disp(' ')
 end

if ncx
    disp('correlation of shocks')
    disp(tit1)
    ip = nvx+nvn+1;
    for i=1:ncx
        k1 = estim_params_.corrx(i,1);
        k2 = estim_params_.corrx(i,2);
        name = [deblank(M_.exo_names(k1,:)) ',' deblank(M_.exo_names(k2,:))];
        NAME = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];
        if strcmp(field_name,'posterior')           
            fprintf('%-*s %7.3f %8.4f %7.4f %4s %6.4f \n', ...
                     header_width,name,bayestopt_.p1(ip),xparam1(ip),stdh(ip),  ...
                     pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.p2(ip));
        else
            fprintf('%-*s %8.4f %7.4f %7.4f \n', header_width,name,xparam1(ip),stdh(ip),tstath(ip));            
        end
        M_.Sigma_e(k1,k2) = xparam1(ip)*sqrt(M_.Sigma_e(k1,k1)*M_.Sigma_e(k2,k2));
        M_.Sigma_e(k2,k1) = M_.Sigma_e(k1,k2);
        eval(['oo_.' field_name '_mode.shocks_corr.' NAME ' = xparam1(ip);']);
        eval(['oo_.' field_name '_std.shocks_corr.' NAME ' = stdh(ip);']);
        ip = ip+1;
    end
    disp(' ')
end

if ncn
    disp('correlation of measurement errors')
    disp(tit1)
    ip = nvx+nvn+ncx+1;
    for i=1:ncn
        k1 = estim_params_.corrn(i,1);
        k2 = estim_params_.corrn(i,2);
        name = [deblank(M_.endo_names(k1,:)) ',' deblank(M_.endo_names(k2,:))];
        NAME = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
        if strcmp(field_name,'posterior')                 
            fprintf('%-*s %7.3f %8.4f %7.4f %4s %6.4f \n', ...
                     header_width,name,bayestopt_.p1(ip),xparam1(ip),stdh(ip), ...
                     pnames(bayestopt_.pshape(ip)+1,:), bayestopt_.p2(ip));
        else
            fprintf('%-*s %8.4f %7.4f %7.4f \n',header_width,name,xparam1(ip),stdh(ip),tstath(ip));            
        end
        eval(['oo_.' field_name '_mode.measurement_errors_corr.' NAME ' = xparam1(ip);']);
        eval(['oo_.' field_name '_std.measurement_errors_corr.' NAME ' = stdh(ip);']);
        ip = ip+1;
    end
    disp(' ')
end

OutputDirectoryName = CheckPath('Output',M_.dname);

if any(bayestopt_.pshape > 0) && options_.TeX %% Bayesian estimation (posterior mode) Latex output
    if np
        filename = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_1.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (parameters)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcccc} \n');
        fprintf(fidTeX,'\\caption{Results from posterior maximization (parameters)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:1}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:1}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{6}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+ncx+ncn+1;
        for i=1:np
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    M_.param_names_tex(estim_params_.param_vals(i,1),:),...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                    bayestopt_.p1(ip),...
                    bayestopt_.p2(ip),...
                    xparam1(ip),...
                    stdh(ip));
            ip = ip + 1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvx
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_2.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of structural shocks)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcccc} \n');
        fprintf(fidTeX,'\\caption{Results from posterior maximization (standard deviation of structural shocks)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:2}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:2}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. & Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{6}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = 1;
        for i=1:nvx
            k = estim_params_.var_exo(i,1);
            fprintf(fidTeX,[ '$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                    deblank(M_.exo_names_tex(k,:)),...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)),...
                    bayestopt_.p1(ip),...
                    bayestopt_.p2(ip),...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvn
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_3.TeX'];
        fidTeX  = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (standard deviation of measurement errors)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcccc} \n');
        fprintf(fidTeX,'\\caption{Results from posterior maximization (standard deviation of measurement errors)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:3}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:3}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{6}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+1;
        for i=1:nvn
            idx = strmatch(options_.varobs(estim_params_.nvn_observable_correspondence(i,1),:),M_.endo_names);
            fprintf(fidTeX,'$%s$ & %4s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    deblank(M_.endo_names_tex(idx,:)), ...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip),...
                    xparam1(ip),...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncx
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_4.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of structural shocks)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcccc} \n');
        fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of structural shocks)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:4}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:4}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{6}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+1;
        for i=1:ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            fprintf(fidTeX,[ '$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n'],...
                    [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))], ...
                    deblank(pnames(bayestopt_.pshape(ip)+1,:)), ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncn
        TeXfile = [OutputDirectoryName '/' M_.fname '_Posterior_Mode_5.TeX'];
        fidTeX = fopen(TeXfile,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,'%% RESULTS FROM POSTERIOR MAXIMIZATION (correlation of measurement errors)\n');
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcccc} \n');
        fprintf(fidTeX,'\\caption{Results from posterior parameters (correlation of measurement errors)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:5}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,'\\label{Table:Posterior:5}\\\\\n');
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Prior distribution & Prior mean  & Prior s.d. &  Posterior mode & s.d. \\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{6}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+ncx+1;
        for i=1:ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            fprintf(fidTeX,'$%s$ & %s & %7.3f & %6.4f & %8.4f & %7.4f \\\\ \n',...
                    [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))], ...
                    pnames(bayestopt_.pshape(ip)+1,:), ...
                    bayestopt_.p1(ip), ...
                    bayestopt_.p2(ip), ...
                    xparam1(ip), ...
                    stdh(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
elseif all(bayestopt_.pshape == 0) && options_.TeX %% MLE and GMM Latex output
    if np        
        filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_1.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,['%% RESULTS FROM ' table_title ' MAXIMIZATION (parameters)\n']);
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcc} \n');
        fprintf(fidTeX,['\\caption{Results from ' table_title ' maximization (parameters)}\n ']);
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':1}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':1}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{4}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+ncx+ncn+1;
        for i=1:np
            fprintf(fidTeX,'$%s$ & %8.4f & %7.4f & %7.4f\\\\ \n',...
                    M_.param_names_tex(estim_params_.param_vals(i,1),:),...
                    xparam1(ip),...
                    stdh(ip),...
                    tstath(ip));                    
            ip = ip + 1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvx
        filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_2.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,['%% RESULTS FROM ' table_title ' MAXIMIZATION (standard deviation of structural shocks)\n']);
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcc} \n');
        fprintf(fidTeX,['\\caption{Results from ' table_title ' maximization (standard deviation of structural shocks)}\n ']);
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':2}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':2}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{4}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = 1;
        for i=1:nvx
            k = estim_params_.var_exo(i,1);
            fprintf(fidTeX,[ '$%s$ & %8.4f & %7.4f & %7.4f\\\\ \n'],...
                    deblank(M_.exo_names_tex(k,:)),...
                    xparam1(ip), ...
                    stdh(ip),...
                    tstath(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if nvn
        filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_3.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,['%% RESULTS FROM ' table_title ' MAXIMIZATION (standard deviation of measurement errors)\n']);
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcc} \n');
        fprintf(fidTeX,['\\caption{Results from ' table_title ' maximization (standard deviation of measurement errors)}\n ']);
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':3}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':3}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{4}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+1;
        for i=1:nvn
           idx = strmatch(options_.varobs(estim_params_.nvn_observable_correspondence(i,1),:),M_.endo_names);
           fprintf(fidTeX,'$%s$ & %8.4f & %7.4f & %7.4f \\\\ \n',...
                    deblank(M_.endo_names_tex(idx,:)), ...
                    xparam1(ip),...
                    stdh(ip),...
                    tstath(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncx
        filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_4.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,['%% RESULTS FROM ' table_title ' MAXIMIZATION (correlation of structural shocks)\n']);
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcc} \n');
        fprintf(fidTeX,['\\caption{Results from ' table_title ' maximization (correlation of structural shocks)}\n ']);
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':4}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':4}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{4}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+1;
        for i=1:ncx
            k1 = estim_params_.corrx(i,1);
            k2 = estim_params_.corrx(i,2);
            fprintf(fidTeX,[ '$%s$  & %8.4f & %7.4f & %7.4f \\\\ \n'],...
                    [deblank(M_.exo_names_tex(k1,:)) ',' deblank(M_.exo_names_tex(k2,:))], ...
                    xparam1(ip), ...
                    stdh(ip),...
                    tstath(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end
    if ncn
        filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '_Mode_5.TeX'];
        fidTeX = fopen(filename,'w');
        fprintf(fidTeX,'%% TeX-table generated by dynare_estimation (Dynare).\n');
        fprintf(fidTeX,['%% RESULTS FROM ' table_title ' MAXIMIZATION (correlation of measurement errors)\n']);
        fprintf(fidTeX,['%% ' datestr(now,0)]);
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,' \n');
        fprintf(fidTeX,'\\begin{center}\n');
        fprintf(fidTeX,'\\begin{longtable}{l|lcc} \n');
        fprintf(fidTeX,['\\caption{Results from ' table_title ' maximization (correlation of measurement errors)}\n ']);
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':5}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endfirsthead \n');
        fprintf(fidTeX,'\\caption{(continued)}\n ');
        fprintf(fidTeX,['\\label{Table:' LaTeXtitle ':5}\\\\\n']);
        fprintf(fidTeX,'\\hline\\hline \\\\ \n');
        fprintf(fidTeX,'  & Mode & s.d. & t-stat\\\\ \n');
        fprintf(fidTeX,'\\hline \\endhead \n');
        fprintf(fidTeX,'\\hline \\multicolumn{4}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n');
        fprintf(fidTeX,'\\hline \\hline \\endlastfoot \n');
        ip = nvx+nvn+ncx+1;
        for i=1:ncn
            k1 = estim_params_.corrn(i,1);
            k2 = estim_params_.corrn(i,2);
            fprintf(fidTeX,'$%s$  & %8.4f & %7.4f & %7.4f \\\\ \n',...
                    [deblank(M_.endo_names_tex(k1,:)) ',' deblank(M_.endo_names_tex(k2,:))], ...
                    xparam1(ip), ...
                    stdh(ip),...
                    tstath(ip));
            ip = ip+1;
        end
        fprintf(fidTeX,'\\end{longtable}\n ');    
        fprintf(fidTeX,'\\end{center}\n');
        fprintf(fidTeX,'%% End of TeX file.\n');
        fclose(fidTeX);
    end

end