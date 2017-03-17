function indcorr = stab_map_2(x,alpha2, pvalue_crit, fnam, dirname,xparam1,figtitle)
% function stab_map_2(x, alpha2, pvalue, fnam, dirname,xparam1)
%
% Written by Marco Ratto
% Joint Research Centre, The European Commission,
% marco.ratto@ec.europa.eu 

% Copyright (C) 2011-2016 European Commission
% Copyright (C) 2011-2016 Dynare Team
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

%global bayestopt_ estim_params_ dr_ options_ ys_ fname_
global bayestopt_ estim_params_ options_ oo_ M_

npar=size(x,2);
nsam=size(x,1);
ishock= npar>estim_params_.np;
nograph = options_.nograph;
if nargin<4,
  fnam='';
end
if nargin<5,
  dirname='';
  nograph=1;
end
if nargin<6,
  xparam1=[];
end
if nargin<7,
  figtitle=fnam;
end

ys_ = oo_.dr.ys;
dr_ = oo_.dr;
fname_ = M_.fname;
nshock = estim_params_.nvx;
nshock = nshock + estim_params_.nvn;
nshock = nshock + estim_params_.ncx;
nshock = nshock + estim_params_.ncn;

[c0, pvalue] = corrcoef(x);
c00=tril(c0,-1);
fig_nam_=[fname_,'_',fnam,'_corr_'];
fig_nam_tex_table=strrep([fnam,'_corr'],' ','_');
fig_nam_=strrep(fig_nam_,' ','_');

ifig=0;
j2=0;
if ishock==0
  npar=estim_params_.np;
  if ~isempty(xparam1),
      xparam1=xparam1(nshock+1:end);
  end
else
  npar=estim_params_.np+nshock;
end
skipline();
title_string=['Correlation analysis for ',fnam];
title_string_tex=['Correlation analysis for ',strrep(fnam,'_','\\_')];

indcorr = [];
entry_iter=1;
for j=1:npar,
    i2=find(abs(c00(:,j))>alpha2);
    if length(i2)>0,
        for jx=1:length(i2),
            if pvalue(j,i2(jx))<pvalue_crit,
                indcorr = [indcorr; [j i2(jx)]];
                j2=j2+1;
                if ishock,
                    if options_.TeX
                        [param_name_temp1, param_name_tex_temp1]= get_the_name(j,options_.TeX,M_,estim_params_,options_);
                        param_name_tex_temp1 = strrep(param_name_tex_temp1,'$','');
                        [param_name_temp2, param_name_tex_temp2]= get_the_name(i2(jx),options_.TeX,M_,estim_params_,options_);
                        param_name_tex_temp2 = strrep(param_name_tex_temp2,'$','');
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        tmp_name_tex=(['[',param_name_tex_temp1,',',param_name_tex_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                        name_tex{entry_iter,1}=tmp_name_tex;
                    else
                        [param_name_temp1]= get_the_name(j,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2]= get_the_name(i2(jx),options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                    end
                else
                    if options_.TeX
                        [param_name_temp1, param_name_tex_temp1]= get_the_name(j+nshock,options_.TeX,M_,estim_params_,options_);
                        param_name_tex_temp1 = strrep(param_name_tex_temp1,'$','');
                        [param_name_temp2, param_name_tex_temp2]= get_the_name(i2(jx)+nshock,options_.TeX,M_,estim_params_,options_);
                        param_name_tex_temp2 = strrep(param_name_tex_temp2,'$','');
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        tmp_name_tex=(['[',param_name_tex_temp1,',',param_name_tex_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                        name_tex{entry_iter,1}=tmp_name_tex;
                    else
                        [param_name_temp1]= get_the_name(j+nshock,options_.TeX,M_,estim_params_,options_);
                        [param_name_temp2]= get_the_name(i2(jx)+nshock,options_.TeX,M_,estim_params_,options_);
                        tmp_name=(['[',param_name_temp1,',',param_name_temp2,']']);
                        name{entry_iter,1}=tmp_name;
                    end
                end                
                data_mat(entry_iter,1)=c0(i2(jx),j);
                entry_iter=entry_iter+1;
                
                if ~nograph,
                    if mod(j2,12)==1,
                    ifig=ifig+1;
                    hh=dyn_figure(options_,'name',[figtitle,' sample bivariate projection ', num2str(ifig)]);
                end
                subplot(3,4,j2-(ifig-1)*12)
                %             bar(c0(i2,j)),
                %             set(gca,'xticklabel',bayestopt_.name(i2)),
                %             set(gca,'xtick',[1:length(i2)])
                %plot(stock_par(ixx(nfilt+1:end,i),j),stock_par(ixx(nfilt+1:end,i),i2(jx)),'.k')
                %hold on,
                plot(x(:,j),x(:,i2(jx)),'.')
                if ~isempty(xparam1)
                    hold on, plot(xparam1(j),xparam1(i2(jx)),'ro')
                end
                %             xlabel(deblank(estim_params_.param_names(j,:)),'interpreter','none'),
                %             ylabel(deblank(estim_params_.param_names(i2(jx),:)),'interpreter','none'),
                if ishock,
                    xlabel(bayestopt_.name{j},'interpreter','none'),
                    ylabel(bayestopt_.name{i2(jx)},'interpreter','none'),
                else
                    xlabel(bayestopt_.name{j+nshock},'interpreter','none'),
                    ylabel(bayestopt_.name{i2(jx)+nshock},'interpreter','none'),
                end
                title(['cc = ',num2str(c0(i2(jx),j))])
                if (mod(j2,12)==0) && j2>0,
                    dyn_saveas(hh,[dirname,filesep,fig_nam_,int2str(ifig)],options_);
                    if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
                        fidTeX = fopen([dirname,filesep,fig_nam_,int2str(ifig),'.tex'],'w');
                        fprintf(fidTeX,'%% TeX eps-loader file generated by stab_map_2.m (Dynare).\n');
                        fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
                        fprintf(fidTeX,'\\begin{figure}[H]\n');
                        fprintf(fidTeX,'\\centering \n');
                        fprintf(fidTeX,'\\includegraphics[width=0.8\\textwidth]{%s}\n',strrep([dirname,'/',fig_nam_,int2str(ifig)],'\','/'));
                        fprintf(fidTeX,'\\caption{%s.}',[figtitle,' sample bivariate projection ', num2str(ifig)]);
                        fprintf(fidTeX,'\\label{Fig:%s:%u}\n',fig_nam_,ifig);
                        fprintf(fidTeX,'\\end{figure}\n\n');
                        fprintf(fidTeX,'%% End Of TeX file. \n');
                        fclose(fidTeX);
                    end
                end
                end
            end
            
        end
    end
    if ~nograph && (j==(npar)) && j2>0 && (mod(j2,12)~=0),
        dyn_saveas(hh,[dirname,filesep,fig_nam_,int2str(ifig)],options_);
        if options_.TeX && any(strcmp('eps',cellstr(options_.graph_format)))
            fidTeX = fopen([dirname,filesep,fig_nam_,int2str(ifig),'.tex'],'w');
            fprintf(fidTeX,'%% TeX eps-loader file generated by stab_map_2.m (Dynare).\n');
            fprintf(fidTeX,['%% ' datestr(now,0) '\n\n']);
            fprintf(fidTeX,'\\begin{figure}[H]\n');
            fprintf(fidTeX,'\\centering \n');
            fprintf(fidTeX,'\\includegraphics[width=%2.2f\\textwidth]{%s}\n',options_.figures.textwidth*min((j2-(ifig-1)*12)/3,1),strrep([dirname,'/',fig_nam_,int2str(ifig)],'\','/'));
            fprintf(fidTeX,'\\caption{%s.}',[figtitle,' sample bivariate projection ', num2str(ifig)]);
            fprintf(fidTeX,'\\label{Fig:%s:%u}\n',fig_nam_,ifig);
            fprintf(fidTeX,'\\end{figure}\n\n');
            fprintf(fidTeX,'%% End Of TeX file. \n');
            fclose(fidTeX);
        end
    end    
end

if j2==0,
    disp(['No correlation term with pvalue <', num2str(pvalue_crit),' and |corr. coef.| >',num2str(alpha2),' found for ',fnam])
else
    headers=strvcat('Parameters','corrcoef');
    dyntable(options_,title_string,headers,char(name),data_mat, 0, 7, 3);
    if options_.TeX
        dyn_latex_table(M_,options_,title_string_tex,fig_nam_tex_table,headers,char(name_tex),data_mat,0,7,3);
    end
end
%close all
