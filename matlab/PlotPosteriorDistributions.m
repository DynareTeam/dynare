function oo_ = PlotPosteriorDistributions(estim_params_, M_, options_, bayestopt_, oo_)

% function PlotPosteriorDistributions()
% plots posterior distributions
%
% INPUTS
%    estim_params_   [structure] 
%    M_              [structure]
%    options_        [structure] 
%    bayestopt_      [structure]
%    oo_             [structure]
%    
% OUTPUTS
%    oo_             [structure]  
%    
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2005-2016 Dynare Team
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

OutputDirectoryName = CheckPath('Output',M_.dname);

TeX     = options_.TeX;
nblck   = options_.mh_nblck;
nvx     = estim_params_.nvx;
nvn     = estim_params_.nvn;
ncx     = estim_params_.ncx;
ncn     = estim_params_.ncn;
np      = estim_params_.np ;
npar    = nvx+nvn+ncx+ncn+np;

MaxNumberOfPlotPerFigure = 9;% The square root must be an integer!
nn = sqrt(MaxNumberOfPlotPerFigure);

figurename = 'Priors and posteriors';

if TeX && any(strcmp('eps',cellstr(options_.graph_format)))    
    fidTeX = fopen([OutputDirectoryName '/' M_.fname '_PriorsAndPosteriors.tex'],'w');
    fprintf(fidTeX,'%% TeX eps-loader file generated by PlotPosteriorDistributions.m (Dynare).\n');
    fprintf(fidTeX,['%% ' datestr(now,0) '\n']);
    fprintf(fidTeX,' \n');
end

figunumber = 0;
subplotnum = 0;

for i=1:npar
    subplotnum = subplotnum+1;
    if subplotnum == 1
        figunumber = figunumber+1;
        hfig=dyn_figure(options_,'Name',figurename);
    end
    [nam,texnam] = get_the_name(i,TeX,M_,estim_params_,options_);
    if subplotnum == 1
        NAMES = nam;
        if TeX
            TeXNAMES = texnam;
        end
    else
        NAMES = char(NAMES,nam);
        if TeX
            TeXNAMES = char(TeXNAMES,texnam);
        end
    end
    [x2,f2,abscissa,dens,binf2,bsup2] = draw_prior_density(i,bayestopt_);
    top2 = max(f2); 
    if i <= nvx
        name = deblank(M_.exo_names(estim_params_.var_exo(i,1),:));  
        x1 = oo_.posterior_density.shocks_std.(name)(:,1);
        f1 = oo_.posterior_density.shocks_std.(name)(:,2);
        oo_.prior_density.shocks_std.(name)(:,1) = x2;
        oo_.prior_density.shocks_std.(name)(:,2) = f2;
        if ~options_.mh_posterior_mode_estimation
            pmod = oo_.posterior_mode.shocks_std.(name);
        end
    elseif i <= nvx+nvn
        name = options_.varobs{estim_params_.nvn_observable_correspondence(i-nvx,1)};
        x1 = oo_.posterior_density.measurement_errors_std.(name)(:,1);
        f1 = oo_.posterior_density.measurement_errors_std.(name)(:,2);
        oo_.prior_density.measurement_errors_std.(name)(:,1) = x2;
        oo_.prior_density.measurement_errors_std.(name)(:,2) = f2;
        if ~options_.mh_posterior_mode_estimation
            pmod = oo_.posterior_mode.measurement_errors_std.(name);
        end     
    elseif i <= nvx+nvn+ncx
        j = i - (nvx+nvn);
        k1 = estim_params_.corrx(j,1);
        k2 = estim_params_.corrx(j,2);
        name = [deblank(M_.exo_names(k1,:)) '_' deblank(M_.exo_names(k2,:))];  
        x1 = oo_.posterior_density.shocks_corr.(name)(:,1);
        f1 = oo_.posterior_density.shocks_corr.(name)(:,2);
        oo_.prior_density.shocks_corr.(name)(:,1) = x2;
        oo_.prior_density.shocks_corr.(name)(:,2) = f2;
        if ~options_.mh_posterior_mode_estimation
            pmod = oo_.posterior_mode.shocks_corr.(name);  
        end
    elseif i <= nvx+nvn+ncx+ncn
        j = i - (nvx+nvn+ncx);
        k1 = estim_params_.corrn(j,1);
        k2 = estim_params_.corrn(j,2);
        name = [deblank(M_.endo_names(k1,:)) '_' deblank(M_.endo_names(k2,:))];
        x1 = oo_.posterior_density.measurement_errors_corr.(name)(:,1);
        f1 = oo_.posterior_density.measurement_errors_corr.(name)(:,2);
        oo_.prior_density.measurement_errors_corr.(name)(:,1) = x2;
        oo_.prior_density.measurement_errors_corr.(name)(:,2) = f2;
        if ~options_.mh_posterior_mode_estimation
            pmod = oo_.posterior_mode.measurement_errors_corr.(name);
        end
    else
        j = i - (nvx+nvn+ncx+ncn);
        name = deblank(M_.param_names(estim_params_.param_vals(j,1),:));
        x1 = oo_.posterior_density.parameters.(name)(:,1);
        f1 = oo_.posterior_density.parameters.(name)(:,2);
        oo_.prior_density.parameters.(name)(:,1) = x2;
        oo_.prior_density.parameters.(name)(:,2) = f2;
        if ~options_.mh_posterior_mode_estimation
            pmod = oo_.posterior_mode.parameters.(name);
        end
    end
    top1 = max(f1);
    top0 = max([top1;top2]);
    binf1 = x1(1);
    bsup1 = x1(end);
    borneinf = min(binf1,binf2);
    bornesup = max(bsup1,bsup2);
    subplot(nn,nn,subplotnum)
    hh = plot(x2,f2,'-k','linewidth',2);
    set(hh,'color',[0.7 0.7 0.7]);
    hold on;
    plot(x1,f1,'-k','linewidth',2);
    if ~options_.mh_posterior_mode_estimation
        plot( [pmod pmod], [0.0 1.1*top0], '--g', 'linewidth', 2);
    end
    box on;
    axis([borneinf bornesup 0 1.1*top0]);
    title(nam,'Interpreter','none');
    hold off;
    drawnow
    if subplotnum == MaxNumberOfPlotPerFigure || i == npar;
        dyn_saveas(hfig,[OutputDirectoryName '/' M_.fname '_PriorsAndPosteriors' int2str(figunumber)],options_);
        if TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            for j = 1:size(NAMES,1)
                fprintf(fidTeX,'\\psfrag{%s}[1][][0.5][0]{%s}\n',deblank(NAMES(j,:)),deblank(TeXNAMES(j,:)));
            end    
            fprintf(fidTeX,'\\centering\n');
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s/%s_PriorsAndPosteriors%s}\n',options_.figures.textwidth*min(subplotnum/nn,1),OutputDirectoryName,M_.fname,int2str(figunumber));
            fprintf(fidTeX,'\\caption{Priors and posteriors.}');
            fprintf(fidTeX,'\\label{Fig:PriorsAndPosteriors:%s}\n',int2str(figunumber));
            fprintf(fidTeX,'\\end{figure}\n');
            fprintf(fidTeX,' \n');
            if i == npar
                fprintf(fidTeX,'%% End of TeX file.\n');
                fclose(fidTeX);
            end
        end
        subplotnum = 0;
    end
end