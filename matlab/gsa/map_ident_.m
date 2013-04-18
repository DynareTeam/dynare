function map_ident_(OutputDirectoryName,opt_gsa)

% Copyright (C) 2012 Dynare Team
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

global bayestopt_ M_ options_ estim_params_ oo_

% opt_gsa = options_.opt_gsa;
fname_ = M_.fname;
nliv   = opt_gsa.morris_nliv;
ntra   = opt_gsa.morris_ntra;
itrans = opt_gsa.trans_ident;

np = estim_params_.np;
if opt_gsa.load_ident_files,
  gsa_flag=0;
else
  gsa_flag=-2;
end

pnames = M_.param_names(estim_params_.param_vals(:,1),:);

filetoload=[OutputDirectoryName '/' fname_ '_prior'];
load(filetoload,'lpmat','lpmat0','istable','T','yys','nspred','nboth','nfwrd')
if ~isempty(lpmat0),
  lpmatx=lpmat0(istable,:);
else
  lpmatx=[];
end
Nsam = size(lpmat,1);
nshock = size(lpmat0,2);
npT = np+nshock;

fname_ = M_.fname;

if opt_gsa.load_ident_files==0,
  % th moments
%     options_.ar = min(3,options_.ar);

  mss = yys(bayestopt_.mfys,:);
  mss = teff(mss(:,istable),Nsam,istable);
  yys = teff(yys(oo_.dr.order_var,istable),Nsam,istable);
  if exist('T'),
      [vdec, cc, ac] = mc_moments(T, lpmatx, oo_.dr);
  else
      return,
  end


  if opt_gsa.morris==2,
   pdraws = dynare_identification(options_.options_ident,[lpmatx lpmat(istable,:)]);
%    [pdraws, TAU, GAM] = dynare_identification(options_.options_ident,[lpmatx lpmat(istable,:)]);
    if max(max(abs(pdraws-[lpmatx lpmat(istable,:)])))==0,
      disp(['Sample check OK ', num2str(max(max(abs(pdraws-[lpmatx lpmat(istable,:)]))))]),
      clear pdraws;
    end
%     for j=1:length(istable), gas(:,j)=[vech(cc(:,:,j)); vec(ac(:,:,j))];  end
%     if ~isempty(mss),
%     gas = [mss(istable,:)'; gas];
%     end
%     if max(max(abs(GAM-gas)))<=1.e-8,
%       disp(['Moments check OK ',num2str(max(max(abs(GAM-gas))))]),
      clear GAM gas
%     end
  end
  if opt_gsa.morris~=1 & M_.exo_nbr>1,
    ifig=0;
    for j=1:M_.exo_nbr,
      if mod(j,6)==1
        hh=dyn_figure(options_,'name',['Variance decomposition shocks']);
        ifig=ifig+1;
        iplo=0;
      end
      iplo=iplo+1;
      subplot(2,3,iplo)
      myboxplot(squeeze(vdec(:,j,:))',[],'.',[],10)
      %     boxplot(squeeze(vdec(:,j,:))','whis',10,'symbol','.r')
      set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:size(options_.varobs,1)])
      set(gca,'xlim',[0.5 size(options_.varobs,1)+0.5])
      set(gca,'ylim',[-2 102])
      for ip=1:size(options_.varobs,1),
        text(ip,-4,deblank(options_.varobs(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
      end
      xlabel(' ')
      ylabel(' ')
      title(M_.exo_names(j,:),'interpreter','none')
      if mod(j,6)==0 | j==M_.exo_nbr,
        dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_vdec_exo_',int2str(ifig)],options_);  
      end
    end
  end
  for j=1:size(cc,1)
    cc(j,j,:)=stand_(squeeze(log(cc(j,j,:))))./2;
  end
  [vdec, j0, ir_vdec, ic_vdec] = teff(vdec,Nsam,istable);
  [cc, j0, ir_cc, ic_cc] = teff(cc,Nsam,istable);
  [ac, j0, ir_ac, ic_ac] = teff(ac,Nsam,istable);

  [nr1, nc1, nnn] = size(T);
    endo_nbr = M_.endo_nbr;
    nstatic = M_.nstatic;
    nspred = M_.nspred;
    iv = (1:endo_nbr)';
    ic = [ nstatic+(1:nspred) endo_nbr+(1:size(oo_.dr.ghx,2)-nspred) ]';
  
    dr.ghx = T(:, [1:(nc1-M_.exo_nbr)],1);
    dr.ghu = T(:, [(nc1-M_.exo_nbr+1):end], 1);
    [Aa,Bb] = kalman_transition_matrix(dr,iv,ic,M_.exo_nbr);
%     bayestopt_.restrict_var_list, ...
%     bayestopt_.restrict_columns, ...
%     bayestopt_.restrict_aux, M_.exo_nbr);
  A = zeros(size(Aa,1),size(Aa,2)+size(Aa,1),length(istable));
  % Sig(estim_params_.var_exo(:,1))=lpmatx(1,:).^2;
  if ~isempty(lpmatx),
      set_shocks_param(lpmatx(1,:));
  end
  A(:,:,1)=[Aa, triu(Bb*M_.Sigma_e*Bb')];
  for j=2:length(istable),
    dr.ghx = T(:, [1:(nc1-M_.exo_nbr)],j);
    dr.ghu = T(:, [(nc1-M_.exo_nbr+1):end], j);
    [Aa,Bb] = kalman_transition_matrix(dr, iv, ic, M_.exo_nbr);
%       bayestopt_.restrict_var_list, ...
%       bayestopt_.restrict_columns, ...
%       bayestopt_.restrict_aux, M_.exo_nbr);
    if ~isempty(lpmatx),
        set_shocks_param(lpmatx(j,:));
    end
    A(:,:,j)=[Aa, triu(Bb*M_.Sigma_e*Bb')];
  end
  clear T;
  clear lpmatx;

  [nr,nc,nn]=size(A);
  io=bayestopt_.mf2;
  % T1=A(io,1:nr,:);
  % ino=find(~ismember([1:nr],io));
  % T2=A(ino,1:nr,:);
  R=A(:,nr+1:nc,:);
%   [tadj, iff] = gsa_speed(A(1:nr,1:nr,:),R,io,0.5);
%   [tadj, j0, ir_tadj, ic_tadj] = teff(tadj,Nsam,istable);
%   [iff, j0, ir_if, ic_if] = teff(iff,Nsam,istable);


  [yt, j0]=teff(A,Nsam,istable);
  yt = [yys yt];
  if opt_gsa.morris==2,
%     iii=find(std(yt(istable,:))>1.e-8);
%     if max(max(abs(TAU-yt(istable,iii)')))<= 1.e-8,
%       err = max(max(abs(TAU-yt(istable,iii)')));
%       disp(['Model check OK ',num2str(err)]),
      clear TAU A
%     end
  else
    clear A,
  end
  % [yt1, j01]=teff(T1,Nsam,istable);
  % [yt2, j02]=teff(T2,Nsam,istable);
  % [ytr, j0r]=teff(R,Nsam,istable);
  %
  % yt=[yt1 yt2 ytr];
  save([OutputDirectoryName,'/',fname_,'_main_eff'],'ac','cc','vdec','yt','mss')
else
  if opt_gsa.morris==2,
%    [pdraws, TAU, GAM] = dynare_identification([1:npT]); %,[lpmatx lpmat(istable,:)]);
%    [pdraws, TAU, GAM] = dynare_identification(options_.options_ident);
   pdraws = dynare_identification(options_.options_ident);
  end
  load([OutputDirectoryName,'/',fname_,'_main_eff'],'ac','cc','vdec','yt','mss')
end

%   for j=1:nr,
%     for i=1:nc,
%       y0=squeeze(A(j,i,:));
%       if max(y0)-min(y0)>1.e-10,
%         j0=j0+1;
%         y1=ones(size(lpmat,1),1)*NaN;
%         y1(istable,1)=y0;
%         yt(:,j0)=y1;
%       end
%     end
%   end
%   yt = yt(:,j0);

if opt_gsa.morris==1,
  %OutputDir = CheckPath('gsa/screen');
  if ~isempty(vdec),
  if opt_gsa.load_ident_files==0,
  SAMorris = [];
  for i=1:size(vdec,2),
    [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], vdec(:,i),nliv);
  end
  SAvdec = squeeze(SAMorris(:,1,:))';
  save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAvdec','vdec','ir_vdec','ic_vdec')
  else
    load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAvdec','vdec','ir_vdec','ic_vdec')
  end
  
  hh = dyn_figure(options_);
%   boxplot(SAvdec,'whis',10,'symbol','r.')
  myboxplot(SAvdec,[],'.',[],10)
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('All variance decomposition')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_morris_vdec'],options_);
  else
  save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'vdec')
    
  end

%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['EET variance decomposition observed variables']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_vdec==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAvdec(iv,:),[],'.',[],3)
%       else
%         plot(SAvdec(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_vdec_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_vdec_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_vdec_varobs_',int2str(ifig)]);
%       close(gcf)
%     end
%   end
% 
%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['EET variance decomposition shocks']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_vdec==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAvdec(iv,:),[],'.',[],3)
%       else
%         plot(SAvdec(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-2,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(M_.exo_names(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==M_.exo_nbr,
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_vdec_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_vdec_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_vdec_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end


  if opt_gsa.load_ident_files==0,
  SAMorris = [];
  ccac = [mss cc ac];
  for i=1:size(ccac,2),
    [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], [ccac(:,i)],nliv);
  end
  SAcc = squeeze(SAMorris(:,1,:))';
  SAcc = SAcc./(max(SAcc')'*ones(1,npT));
  save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAcc','cc','ir_cc','ic_cc','-append')
  save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'ac','ir_ac','ic_ac','-append')
  else
    load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAcc','cc','ir_cc','ic_cc')
    load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'ac','ir_ac','ic_ac')
  end
  
  hh=dyn_figure(options_);
%   boxplot(SAcc,'whis',10,'symbol','r.')
  myboxplot(SAcc,[],'.',[],10)
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('EET All moments')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_morris_moments'],options_);
%   close(gcf),

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MORRIS FOR DERIVATIVES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if opt_gsa.load_ident_files==0,
%     for j=1:npT,
%   SAMorris = [];
%   ddd=NaN(size(lpmat,1),size(JJ,1));
%   ddd(istable,:) = squeeze(JJ(:,j,:))';
%   for i=1:size(ddd,2),
%     [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], [ddd(:,i)],nliv);
%   end
%   SAddd(:,:,j) = squeeze(SAMorris(:,1,:))';
%   SAddd(:,:,j) = SAddd(:,:,j)./(max(SAddd(:,:,j)')'*ones(1,npT));
%   sad(:,j) = median(SAddd(find(~isnan(squeeze(SAddd(:,1,j)))),:,j))'; 
%     end
%   save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAddd','sad','-append')
%   else
%     load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAddd','sad')
%   end
%   figure,
%   contourf(sad,10), colorbar
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%   set(gca,'yticklabel',' ','fontsize',10,'ytick',[1:npT])
%   for ip=1:npT,
%     text(ip,0.9,['D(',bayestopt_.name{ip},')'],'rotation',90,'HorizontalAlignment','right','interpreter','none')
%     text(0.9,ip,[bayestopt_.name{ip}],'rotation',0,'HorizontalAlignment','right','interpreter','none')
%   end
%   [m,im]=max(sad);
%   iii = find((im-[1:npT])==0);
%   disp('Most identified params')
%   disp(bayestopt_.name(iii))
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF MORRIS FOR DERIVATIVES  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['EET cross-correlations']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_cc==j);
%     iv = [iv; find(ic_cc==j)];
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAcc(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAcc(iv,:),[],'.',[],3)
%       else
%         plot(SAcc(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_cc_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_cc_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_cc_',int2str(ifig)]);
%       close(gcf),
%     end
%   end


%   if opt_gsa.load_ident_files==0,
%   SAMorris = [];
%   for i=1:size(ac,2),
%     [SAmeas, SAMorris(:,:,i)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], ac(:,i),nliv);
%   end
%   %end
%   SAac = squeeze(SAMorris(:,1,:))';
%   save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAac','ac','ir_ac','ic_ac','-append')
%   else
%     load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAac','ac','ir_ac','ic_ac')
%   end
%   figure,
% %   boxplot(SAac,'whis',10,'symbol','r.')
%   myboxplot(SAac,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%   set(gca,'xlim',[0.5 npT+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:npT,
%     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title('EET All auto-correlations')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_ac'])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_ac']);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_ac']);
%   close(gcf),

%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['EET auto-correlations']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_ac==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAac(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAac(iv,:),[],'.',[],3)
%       else
%         plot(SAac(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_ac_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_ac_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_ac_',int2str(ifig)]);
%       close(gcf),
%     end
%   end

%   if opt_gsa.load_ident_files==0,
%   js=0;
%   %for j=1:size(tadj,1),
%   SAMorris = [];
%   for i=1:size(tadj,2),
%     js=js+1;
%     [SAmeas, SAMorris(:,:,js)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], tadj(:,i),nliv);
%   end
%   %end
%   SAM = squeeze(SAMorris(nshock+1:end,1,:));
%   for j=1:js,
%     SAtadj(:,j)=SAM(:,j)./(max(SAM(:,j))+eps);
%   end
%   SAtadj = SAtadj';
%   save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAtadj','tadj','ir_tadj','ic_tadj','-append')
%   else
%     load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAtadj','tadj','ir_tadj','ic_tadj')
%   end
%   if opt_gsa.load_ident_files==0,
%   js=0;
%   SAMorris = [];
%   for i=1:size(iff,2),
%     js=js+1;
%     [SAmeas, SAMorriss(:,:,js)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], iff(:,i),nliv);
%   end
%   SAM = squeeze(SAMorriss(nshock+1:end,1,:));
%   for j=1:js,
%     SAIF(:,j)=SAM(:,j)./(max(SAM(:,j))+eps);
%   end
%   SAIF = SAIF';
%   save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAIF','iff','ir_if','ic_if','-append')
%   else
%     load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAIF','iff','ir_if','ic_if')
%   end
%   figure,
%   %bar(SAtadj),
% %   boxplot(SAtadj,'whis',10,'symbol','r.')
%   myboxplot(SAtadj,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   set(gca,'ylim',[0 1])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title('All half-life')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_tadj'])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_tadj']);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_tadj']);
%   close(gcf),

%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['EET speed of adjustment observed variables']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_tadj==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAtadj(iv,:),[],'.',[],3)
%       else
%         plot(SAtadj(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_tadj_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_tadj_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_tadj_varobs_',int2str(ifig)]);
%       close(gcf),
%     end
%   end

%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['EET speed of adjustment shocks']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_tadj==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAtadj(iv,:),[],'.',[],3)
%       else
%         plot(SAtadj(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(M_.exo_names(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==M_.exo_nbr,
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_tadj_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_tadj_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_tadj_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end

%   figure,
%   %bar(SAIF),
% %   boxplot(SAIF,'whis',10,'symbol','r.')
%   myboxplot(SAIF,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   set(gca,'ylim',[0 1])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   ylabel('Elementary Effects')
%   title('Steady state gains (impact factors)')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_gain'])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_gain']);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_gain']);
%   close(gcf),
  %figure, bar(SAIF'), title('All Gain Relationships')
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['EET steady state gain observed series']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_if==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAIF(iv,:),'whis',10,'symbol','r.');
%         myboxplot(SAIF(iv,:),[],'.',[],10)
%       else
%         plot(SAIF(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_gain_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_gain_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_gain_varobs_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['EET steady state gain shocks']);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_if==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAIF(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAIF(iv,:),[],'.',[],3)
%       else
%         plot(SAIF(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(M_.exo_names(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==M_.exo_nbr,
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_gain_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_gain_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_gain_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end


  if opt_gsa.load_ident_files==0,
  SAMorris = [];
  for j=1:j0,
    [SAmeas, SAMorris(:,:,j)] = Morris_Measure_Groups(npT, [lpmat0 lpmat], yt(:,j),nliv);
  end

%   SAM = squeeze(SAMorris(nshock+1:end,1,:));
  SAM = squeeze(SAMorris(1:end,1,:));
  for j=1:j0
    SAnorm(:,j)=SAM(:,j)./max(SAM(:,j));
    irex(j)=length(find(SAnorm(:,j)>0.01));
  end
  [dum, irel]=sort(irex);

%   SAMmu = squeeze(SAMorris(nshock+1:end,2,:));
  SAMmu = squeeze(SAMorris(1:end,2,:));
  for j=1:j0
    SAmunorm(:,j)=SAMmu(:,j)./max(SAM(:,j));  % normalised w.r.t. mu*
  end
%   SAMsig = squeeze(SAMorris(nshock+1:end,3,:));
  SAMsig = squeeze(SAMorris(1:end,3,:));
  for j=1:j0
    SAsignorm(:,j)=SAMsig(:,j)./max(SAMsig(:,j));
  end
    save([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAnorm','SAmunorm','SAsignorm','-append')
  else
    load([OutputDirectoryName,'/',fname_,'_morris_IDE'],'SAnorm','SAmunorm','SAsignorm')
  end
  hh=dyn_figure(options_); %bar(SAnorm(:,irel))
%   boxplot(SAnorm','whis',10,'symbol','r.')
  myboxplot(SAnorm',[],'.',[],10)
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:npT,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('Elementary effects in the model')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_morris_par'],options_);

  hh=dyn_figure(options_); %bar(SAmunorm(:,irel))
%   boxplot(SAmunorm','whis',10,'symbol','r.')
  myboxplot(SAmunorm',[],'.',[],10)
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  set(gca,'ylim',[-1 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:npT,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('\mu in the model')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_morrismu_par'],options_);

  hh=dyn_figure(options_); %bar(SAsignorm(:,irel))
%   boxplot(SAsignorm','whis',10,'symbol','r.')
  myboxplot(SAsignorm',[],'.',[],10)
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  set(gca,'ylim',[0 1])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  xlabel(' ')
  for ip=1:npT,
    text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title('\sigma in the model')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_morrissig_par'],options_);

  %     figure, bar(SAnorm(:,irel)')
  %     set(gca,'xtick',[1:j0])
  %     set(gca,'xlim',[0.5 j0+0.5])
  %     title('Elementary effects relationships')
  %     saveas(gcf,[OutputDirectoryName,'/',fname_,'_morris_redform'])
  %     eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_morris_redform']);
  %     eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_morris_redform']);

elseif opt_gsa.morris==3,
    return
    
  np=estim_params_.np;
  na=(4*np+1)*opt_gsa.Nsam;
  for j=1:j0,
    [idex(j,:), yd(j,:)] = spop_ide(lpmat, yt(:,j), opt_gsa.Nsam, 5-1);
  end
  iok=find(~isnan(yt(1:opt_gsa.Nsam,1)));
  yr=NaN*ones(size(lpmat,1),j0);
  for j=1:j0,
    ys(j,:)=yd(j,:)./max(yd(j,:));
    [dum, is]=sort(yt(iok,j));
    yr(iok(is),j)=[1:length(iok)]'./length(iok);
    yr(istable(length(iok)+1:end),j) = interp1(yt(iok,j),yr(iok,j),yt(istable(length(iok)+1:end),j),'','extrap');
    ineg=find(yr(:,j)<0);
    if any(ineg),
      [dum, is]=sort(yr(ineg,j));
      yr(ineg(is),j)=-[length(ineg):-1:1]./length(iok);

    end
    [idex_r(j,:), yd_r(j,:)] = spop_ide(lpmat, yr(:,j), opt_gsa.Nsam, 5-1);
    ys_r(j,:)=yd_r(j,:)./max(yd_r(j,:));

  end,
  figure, bar((idex.*ys)./opt_gsa.Nsam), title('Relationships')
  figure, bar((idex.*ys)'./opt_gsa.Nsam), title('Parameters')
  figure, bar((idex_r.*ys_r)./opt_gsa.Nsam), title('Relationships rank')
  figure, bar((idex_r.*ys_r)'./opt_gsa.Nsam), title('Parameters rank')
  [v0,d0]=eig(corrcoef(yt(iok,:)));
  ee=diag(d0);
  ee=ee([end:-1:1])./j0;
  i0=length(find(ee>0.01));
  v0=v0(:,[end:-1:1]);
  for j=1:i0,
    [idex_pc(j,:), yd_pc(j,:)] = spop_ide(lpmat, yt*v0(:,j), opt_gsa.Nsam, 5-1);
  end
  for j=1:i0,
    ys_pc(j,:)=yd_pc(j,:)./max(yd_pc(j,:));
  end,
  figure, bar((idex_pc.*ys_pc)./opt_gsa.Nsam), title('Relationships PCA')
  figure, bar((idex_pc.*ys_pc)'./opt_gsa.Nsam), title('Parameters PCA')

  [vr,dr]=eig(corrcoef(yr(iok,:)));
  er=diag(dr);
  er=er([end:-1:1])./j0;
  ir0=length(find(er>0.01));
  vr=vr(:,[end:-1:1]);
  for j=1:ir0,
    [idex_pcr(j,:), yd_pcr(j,:)] = spop_ide(lpmat, yr*vr(:,j), opt_gsa.Nsam, 5-1);
  end
  for j=1:ir0,
    ys_pcr(j,:)=yd_pcr(j,:)./max(yd_pcr(j,:));
  end,
  figure, bar((idex_pcr.*ys_pcr)./opt_gsa.Nsam), title('Relationships rank PCA')
  figure, bar((idex_pcr.*ys_pcr)'./opt_gsa.Nsam), title('Parameters rank PCA')
  
elseif opt_gsa.morris==2,   % ISKREV staff
  return,

  
else,  % main effects analysis
  
  if itrans==0,
    fsuffix = '';
  elseif itrans==1,
    fsuffix = '_log';
  else
    fsuffix = '_rank';
  end
  
  imap=[1:npT];

  if isempty(lpmat0),
      x0=lpmat(istable,:);
  else
      
      x0=[lpmat0(istable,:), lpmat(istable,:)];
  end
  nrun=length(istable);
  nest=min(250,nrun);
  nfit=min(1000,nrun);
  
%   opt_gsa.load_ident_files=0;
  
%   if opt_gsa.load_ident_files==0,
%   try 
%     EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAvdec','vdec','ir_vdec','ic_vdec');
%   catch
%     EET=[];
%   end
%   SAvdec=zeros(size(vdec,2),npT);
% 
%   for j=1:size(vdec,2),
%     if itrans==0,
%       y0 = vdec(istable,j);
%     elseif itrans==1,
%       y0 = log_trans_(vdec(istable,j));
%     else
%       y0 = trank(vdec(istable,j));
%     end
%     if ~isempty(EET),
% %       imap=find(EET.SAvdec(j,:));
% %       [dum, isort]=sort(-EET.SAvdec(j,:));
%       imap=find(EET.SAvdec(j,:) >= (0.1.*max(EET.SAvdec(j,:))) );
%     end
%   gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
%       2, [],[],[],0,[OutputDirectoryName,'/map_vdec',fsuffix,int2str(j)], pnames);
%   if nfit>nest,
%     gsa_(j) = gsa_sdp(y0(1:nfit), x0(1:nfit,imap), ...
%         -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_vdec',fsuffix,int2str(j)], pnames);
%   end
%     
%     SAvdec(j,imap)=gsa_(j).si;
%     imap_vdec{j}=imap;
%   end
%   save([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_vdec','SAvdec','vdec','ir_vdec','ic_vdec','-append')
%   else
%   load([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_vdec','SAvdec','vdec','ir_vdec','ic_vdec')
%   end
%   figure,
% %   boxplot(SAvdec,'whis',10,'symbol','r.')
%   myboxplot(SAvdec,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 npT+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:npT,
%     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects variance decomposition ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_vdec',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix]);
%   close(gcf),
% 
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['Main effects observed variance decomposition ',fsuffix]);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_vdec==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAvdec(iv,:),'whis',10,'symbol','r.');
%         myboxplot(SAvdec(iv,:),[],'.',[],10)
%       else
%         plot(SAvdec(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_vdec',fsuffix,'_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix,'_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix,'_varobs_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['Main effects shocks variance decomposition ',fsuffix]);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_vdec==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAvdec(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAvdec(iv,:),[],'.',[],10)
%       else
%         plot(SAvdec(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',3,'xtick',[1:np])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%       set(gca,'fontsize',10)
%     end
%     title(M_.exo_names(j,:),'interpreter','none','fontsize',10)
%     if mod(j,6)==0 | j==M_.exo_nbr
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_vdec',fsuffix,'_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix,'_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_vdec',fsuffix,'_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end

  if opt_gsa.load_ident_files==0,
  try 
    EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAcc','ir_cc','ic_cc');
  catch
    EET=[];
  end
  ccac = stand_([mss cc ac]);
  [pcc, dd] = eig(cov(ccac(istable,:)));
  [latent, isort] = sort(-diag(dd));
  latent = -latent;
  figure, bar(latent)
  title('Eigenvalues in PCA')
  pcc=pcc(:,isort);
  ccac = ccac*pcc;
%   npca = min(40, max(find(cumsum(latent)./length(latent)<0.99))+1);
  npca = max(find(cumsum(latent)./length(latent)<0.99))+1;
  siPCA = (EET.SAcc'*abs(pcc'))';
%   siPCA = siPCA./(max(siPCA')'*ones(1,npT)).*(latent*ones(1,npT));
  siPCA = siPCA./(max(siPCA')'*ones(1,npT));
%   siPCA = sum(siPCA,1);
%   siPCA = siPCA./max(siPCA);
  SAcc=zeros(size(ccac,2),npT);
  for j=1:npca, %size(ccac,2),
    if itrans==0,
      y0 = ccac(istable,j);
    elseif itrans==1,
      y0 = log_trans_(ccac(istable,j));
    else
      y0 = trank(ccac(istable,j));
    end
    if ~isempty(EET),
%       imap=find(EET.SAvdec(j,:));
%       [dum, isort]=sort(-EET.SAvdec(j,:));
      imap=find(siPCA(j,:) >= (0.1.*max(siPCA(j,:))) );
%       imap=find(EET.SAcc(j,:) >= (0.1.*max(EET.SAcc(j,:))) );
    end
  gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
      2, [],[],[],0,[OutputDirectoryName,'/map_cc',fsuffix,int2str(j)], pnames);
%   if nfit>nest,
%     gsa_(j) = gsa_sdp(y0(1:nfit), x0(1:nfit,imap), ...
%         -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_cc',fsuffix,int2str(j)], pnames);
%   end
    SAcc(j,imap)=gsa_(j).si;
    imap_cc{j}=imap;

  end
  save([OutputDirectoryName,'/map_cc',fsuffix],'gsa_')
  save([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_cc','SAcc','ccac','-append')
  else
    load([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_cc','SAcc','ccac')
    
  end
%   figure,
% %   boxplot(SAcc,'whis',10,'symbol','r.')
%   myboxplot(SAcc,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 npT+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:npT,
%     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   ylabel(' ')
%   title(['Main effects moments''s PCA ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_cc',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_moments',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_moments',fsuffix]);
%   close(gcf),

%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['Main effects cross-covariances ',fsuffix]);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_cc==j);
%     iv = [iv; find(ic_cc==j)];
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAcc(iv,:),'whis',10,'symbol','r.');
%         myboxplot(SAcc(iv,:),[],'.',[],10)
%       else
%         plot(SAcc(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%       set(gca,'fontsize',10)
%     end
%     title(options_.varobs(j,:),'interpreter','none','fontsize',10)
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_cc',fsuffix,'_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_cc',fsuffix,'_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_cc',fsuffix,'_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
%   if opt_gsa.load_ident_files==0,
%   try 
%     EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAac','ir_ac','ic_ac');
%   catch
%     EET=[];
%   end
%   SAac=zeros(size(ac,2),npT);
%   for j=1:size(ac,2),
%     if itrans==0,
%       y0 = ac(istable,j);
%     elseif itrans==1,
%       y0 = log_trans_(ac(istable,j));
%     else
%       y0 = trank(ac(istable,j));
%     end
%     if ~isempty(EET),
%       imap=find(EET.SAac(j,:) >= (0.1.*max(EET.SAac(j,:))) );
%     end
% %     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
% %       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
%   gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
%       2, [],[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
%   if nfit>nest,
%     gsa_(j) = gsa_sdp(y0(1:nfit), x0(1:nfit,imap), ...
%         -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_ac',fsuffix,int2str(j)], pnames);
%   end
%     SAac(j,imap)=gsa_(j).si;
%     imap_ac{j}=imap;
% 
%   end
%   save([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_ac','SAac','ac','ir_ac','ic_ac','-append')
%   else
%   load([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_ac','SAac','ac','ir_ac','ic_ac')
%   end
% 
%   figure,
% %   boxplot(SAac,'whis',10,'symbol','r.')
%   myboxplot(SAac,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 npT+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects 1 lag auto-covariances ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_ac',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_ac',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_ac',fsuffix]);
%   close(gcf),
% 
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['Main effects auto-covariances ',fsuffix]);
%       ifig=ifig+1;
%       iplo = 0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_ac==j);
%     %iv = [iv; find(ic_ac==j)];
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAac(iv,:),'whis',10,'symbol','r.');
%         myboxplot(SAac(iv,:),[],'.',[],10)
%       else
%         plot(SAac(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 npT+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:npT,
%         text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%       set(gca,'fontsize',10)
%     end
%     title(options_.varobs(j,:),'interpreter','none','fontsize',10)
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_ac',fsuffix,'_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_ac',fsuffix,'_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_ac',fsuffix,'_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
  
%   x0=x0(:,nshock+1:end);
  imap=[1:npT];

%   if opt_gsa.load_ident_files==0,
%   try 
%     EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAtadj','ir_tadj','ic_tadj');
%     ny=size(EET.SAtadj,1);
%   catch
%     EET=[];
%   end
%   SAtadj=zeros(size(tadj,2),np);
%   for j=1:size(tadj,2),
%     if itrans==0,
%       y0 = tadj(istable,j);
%     elseif itrans==1,
%       y0 = log_trans_(tadj(istable,j));
%     else
%       y0 = trank(tadj(istable,j));
%     end
%     if ~isempty(EET),
%       if size(tadj,2)~=ny,
%         jj=find(EET.ir_tadj==ir_tadj(j));
%         jj=jj(find(EET.ic_tadj(jj)==ic_tadj(j)));
%         if ~isempty(jj),
%           imap=find(EET.SAtadj(jj,:) >= (0.1.*max(EET.SAtadj(jj,:))) );
%         else
%           imap=[1:np];
%         end
%       else
%         imap=find(EET.SAtadj(j,:) >= (0.1.*max(EET.SAtadj(j,:))) );
%       end
%     end
% %     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
% %       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
%   gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
%       2, [],[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
%   if nfit>nest,
%     gsa_(j) = gsa_sdp(y0(1:nfit), x0(1:nfit,imap), ...
%         -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_tadj',fsuffix,int2str(j)], pnames);
%   end
%     SAtadj(j,imap)=gsa_(j).si;
%     imap_tadj{j}=imap;
% 
%   end
%   save([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_tadj','SAtadj','tadj','ir_tadj','ic_tadj','-append')
%   else
%   load([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_tadj','SAtadj','tadj','ir_tadj','ic_tadj')
%   end
% 
%   figure,
% %   boxplot(SAtadj,'whis',10,'symbol','r.')
%   myboxplot(SAtadj,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects speed of adjustment ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_tadj',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix]);
%   close(gcf),
% 
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['Main effects observed speed adjustment ',fsuffix]);
%       ifig=ifig+1;
%       iplo = 0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_tadj==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAtadj(iv,:),[],'.',[],10)
%       else
%         plot(SAtadj(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_tadj',fsuffix,'_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix,'_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix,'_varobs_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['Main effects shocks speed of adjustment ',fsuffix]);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_tadj==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAtadj(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAtadj(iv,:),[],'.',[],10)
%       else
%         plot(SAtadj(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(M_.exo_names(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==M_.exo_nbr,
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_tadj',fsuffix,'_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix,'_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_tadj',fsuffix,'_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
% 
%   if opt_gsa.load_ident_files==0,
%   try 
%     EET=load([OutputDirectoryName,'/SCREEN/',fname_,'_morris_IDE'],'SAIF','ir_if','ic_if');
%   catch
%     EET=[];
%   end
%   SAif=zeros(size(iff,2),np);
%   for j=1:size(iff,2),
%     if itrans==0,
%       y0 = iff(istable,j);
%     elseif itrans==1,
%       y0 = log_trans_(iff(istable,j));
%     else
%       y0 = trank(iff(istable,j));
%     end
%     if ~isempty(EET),
%       imap=find(EET.SAIF(j,:) >= (0.1.*max(EET.SAIF(j,:))) );
%     end
% %     gsa_(j) = gsa_sdp_dyn( y0, lpmat(istable,:), ...
% %       gsa_flag, [],[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
%   gsa_(j) = gsa_sdp(y0(1:nest), x0(1:nest,imap), ...
%       2, [],[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
%   if nfit>nest,
%     gsa_(j) = gsa_sdp(y0(1:nfit), x0(1:nfit,imap), ...
%         -2, gsa_(j).nvr*nest^3/nfit^3,[],[],0,[OutputDirectoryName,'/map_if',fsuffix,int2str(j)], pnames);
%   end
%     SAif(j,imap)=gsa_(j).si;
%     imap_if{j}=imap;
%     
%   end
%   save([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_if','SAif','iff','ir_if','ic_if','-append')
%   else
%   load([OutputDirectoryName,'/',fname_,'_main_eff'],'imap_if','SAif','iff','ir_if','ic_if')
%   end
% 
%   figure,
% %   boxplot(SAif,'whis',10,'symbol','r.')
%   myboxplot(SAif,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects impact factors ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_if',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_if',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_if',fsuffix]);
%   close(gcf),
% 
%   ifig = 0;
%   for j=1:size(options_.varobs,1)
%     if mod(j,6)==1
%       figure('name',['Main effects observed impact factors ',fsuffix]);
%       ifig=ifig+1;
%       iplo = 0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ir_if==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAif(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAif(iv,:),[],'.',[],10)
%       else
%         plot(SAif(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(options_.varobs(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==size(options_.varobs,1)
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_if',fsuffix,'_varobs_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_if',fsuffix,'_varobs_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_if',fsuffix,'_varobs_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
% 
%   ifig = 0;
%   for j=1:M_.exo_nbr,
%     if mod(j,6)==1
%       figure('name',['Main effects shocks impact factors ',fsuffix]);
%       ifig=ifig+1;
%       iplo=0;
%     end
%     iplo=iplo+1;
%     subplot(3,2,iplo)
%     iv = find(ic_if==j);
%     if ~isempty(iv)
%       if length(iv)>1
% %         boxplot(SAif(iv,:),'whis',3,'symbol','r.');
%         myboxplot(SAif(iv,:),[],'.',[],10)
%       else
%         plot(SAif(iv,:),'r.');
%       end
%       set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%       set(gca,'xlim',[0.5 np+0.5])
%       ydum = get(gca,'ylim');
%       set(gca,'ylim',[0 ydum(2)])
%       for ip=1:np,
%         text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%       end
%       xlabel(' ')
%     end
%     title(M_.exo_names(j,:),'interpreter','none')
%     if mod(j,6)==0 | j==M_.exo_nbr
%       saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_if',fsuffix,'_exo_',int2str(ifig)])
%       eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_if',fsuffix,'_exo_',int2str(ifig)]);
%       eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_if',fsuffix,'_exo_',int2str(ifig)]);
%       close(gcf),
%     end
%   end
%   SAmom = [SAvdec' SAcc' SAac']';
%   SAdyn = [SAtadj' SAif']';
%   SAall = [SAmom(:,nshock+1:end)' SAdyn']';
%   
%   figure,
%   %   boxplot(SAtadj,'whis',10,'symbol','r.')
%   myboxplot(SAmom,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
%   set(gca,'xlim',[0.5 npT+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:npT,
%     %     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects theoretical moments ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_moments',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_moments',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_moments',fsuffix]);
% %   close(gcf),
%   
%   figure,
%   %   boxplot(SAtadj,'whis',10,'symbol','r.')
%   myboxplot(SAdyn,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects short-long term dynamics ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_dynamics',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_dynamics',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_dynamics',fsuffix]);
% %   close(gcf),
% 
%   figure,
%   %   boxplot(SAtadj,'whis',10,'symbol','r.')
%   myboxplot(SAall,[],'.',[],10)
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:np])
%   set(gca,'xlim',[0.5 np+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:np,
%     text(ip,-0.02,deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
% %     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Main effects all ',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_map_ALL',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_map_ALL',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_map_ALL',fsuffix]);
% %   close(gcf),

%   for j=1:size(SAall,1),
%     SAallN(j,:)=SAall(j,:)./max(SAall(j,:));
%   end
%   SAmean=mean(SAallN);
%   for j=1:size(SAmom,1),
%     SAmomN(j,:)=SAmom(j,1:nshock)./max(SAmom(j,1:nshock));
%   end
%   SAmomN(find(isnan(SAmomN)))=0;
%   SAmeanexo=mean(SAmomN(:,1:nshock));

%   figure, bar(latent'*SAcc),
  hh=dyn_figure(options_);
  bar(sum(SAcc)),
  set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:npT])
  set(gca,'xlim',[0.5 npT+0.5])
  ydum = get(gca,'ylim');
  set(gca,'ylim',[0 ydum(2)])
  set(gca,'position',[0.13 0.2 0.775 0.7])
  for ip=1:npT,
    text(ip,-0.02*(ydum(2)),bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
    %     text(ip,-0.02,bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
  end
  xlabel(' ')
  title(['Identifiability indices in the ',fsuffix,' moments.'],'interpreter','none')
  dyn_saveas(hh,[OutputDirectoryName,'/',fname_,'_ident_ALL',fsuffix],options_);
  
%   figure, bar(SAmeanexo),
%   set(gca,'xticklabel',' ','fontsize',10,'xtick',[1:nshock])
%   set(gca,'xlim',[0.5 nshock+0.5])
%   ydum = get(gca,'ylim');
%   set(gca,'ylim',[0 ydum(2)])
%   set(gca,'position',[0.13 0.2 0.775 0.7])
%   for ip=1:nshock,
%     %     text(ip,-0.02*(ydum(2)),deblank(pnames(ip,:)),'rotation',90,'HorizontalAlignment','right','interpreter','none')
%     text(ip,-0.02*(ydum(2)),bayestopt_.name{ip},'rotation',90,'HorizontalAlignment','right','interpreter','none')
%   end
%   xlabel(' ')
%   title(['Identifiability indices for shocks',fsuffix],'interpreter','none')
%   saveas(gcf,[OutputDirectoryName,'/',fname_,'_ident_SHOCKS',fsuffix])
%   eval(['print -depsc2 ' OutputDirectoryName '/' fname_ '_ident_SHOCKS',fsuffix]);
%   eval(['print -dpdf ' OutputDirectoryName '/' fname_ '_ident_SHOCKS',fsuffix]);
end

return
