var pie y R g z YGR INFL INT;
varexo e_R e_g e_z;
parameters tau kap psi1 psi2 rhoR rhog rhoz rr_steady pi_steady gam_steady std_R std_g std_z;

tau        = 2;
kap        = 0.15;

psi1       = 1.5;
psi2       = 1;

rhoR       = 0.96;
rhog       = 0.95;
rhoz       = 0.65;

rr_steady  = 0.4;
pi_steady  = 4.0;
gam_steady = 0.5;
std_R = 0.2;
std_g = 0.8;
std_z = 0.45;

model;
y = y(+1) + g - g(+1) - 1/tau*(R-pie(+1)-z(+1));
pie = 1/(1+rr_steady/400)*pie(+1)+kap*(y-g);
R = rhoR*R(-1)+(1-rhoR)*psi1*pie+(1-rhoR)*psi2*(y-y(-1)+z)+e_R*std_R;
g = rhog*g(-1)+e_g*std_g;
z = rhoz*z(-1)+e_z*std_z;
YGR = gam_steady+100*(y-y(-1)+z);
INFL = pi_steady+400*pie;
INT = pi_steady+rr_steady+4*gam_steady+400*R;
end;

estimated_params;
tau, gamma_pdf, 2, 0.5; 
kap, beta_pdf, 0.2, 0.1;
psi1, gamma_pdf, 1.5, 0.25;
psi2, gamma_pdf, 0.5, 0.25;
rhoR, beta_pdf, 0.5, 0.2;
rhog, beta_pdf, 0.8, 0.1;
rhoz, beta_pdf, 0.66, 0.15;
rr_steady, gamma_pdf, 0.8, 0.5;
pi_steady, gamma_pdf, 0.8, 0.5;
gam_steady, gamma_pdf, 0.8, 0.5;
std_R, inv_gamma_pdf, 0.05, inf;
std_g, inv_gamma_pdf, 0.05, inf;
std_z, inv_gamma_pdf, 0.05, inf;
end;

steady;
check;

varobs YGR INFL INT;

shocks;
var e_g;
stderr 1;
var e_R;
stderr 1;
var e_z;
stderr 1;
end;
options_.TeX=1;
identification;

identification(advanced=1,max_dim_cova_group=3,prior_mc=250);

collect_latex_files;
if system(['pdflatex -halt-on-error -interaction=batchmode ' M_.fname '_TeX_binder.tex'])
    error('TeX-File did not compile.')
end
