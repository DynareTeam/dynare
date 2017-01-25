var z dw dx dy dc1 dc2;
varexo e_w e_x e_y e_z;

parameters rho_w rho_x rho_y rho_z a1 a2 a3 b c;

model(linear);
dw = rho_w*dw(-1)+a1*(dc1(-1))+e_w;
dx = rho_x*dx(-1)+a2*(dc1(-1))+e_x;
dy = rho_y*dy(-1)+a3*(dc2(-1))+e_y;
z = rho_z*z(-1)+dw-dx+e_z;
dc1 = dc1(-1)+dx-b*dy-c*dw;
dc2 = dc2(-1)+dx-b*dy;
end;

estimated_params;
rho_w, normal_pdf, 0.5,0.2;
rho_x, normal_pdf, 0.5,0.2;
rho_y, normal_pdf, 0.5,0.2;
rho_z, normal_pdf, 0.8,0.2;

a1, normal_pdf, 0.1,0.2;
a2, normal_pdf,  -0.1,0.2;
a3, normal_pdf,  0.1,0.2;
b , normal_pdf,  1,0.2;
c , normal_pdf,  1,0.2;

stderr e_w, uniform_pdf,,, 0.01, 0.1;
stderr e_x, uniform_pdf,,, 0.01, 0.1;
stderr e_y, uniform_pdf,,, 0.01, 0.1;
stderr e_z, uniform_pdf,,, 0.01, 0.1;

stderr z, normal_pdf, 0, 0.1;
end;

varobs dw dx dy z;
       
estimation(datafile=data,first_obs=1000,nobs=200,mh_replic=0,filtered_vars,smoothed_state_uncertainty);

//checking smoother consistency
X = oo_.SmoothedVariables;
S = [X.z X.dw X.dx X.dy X.dc1 X.dc2];
X = oo_.SmoothedShocks;
E = [X.e_w X.e_x X.e_y X.e_z];
A = oo_.dr.ghx;
B = oo_.dr.ghu;
err = zeros(6,200);
for t=2:200;
    err(:,t) = S(t,:)'-A*S(t-1,:)'-B*E(t,:)';
end;
if max(max(abs(err))) > 1e-10;
   error('Test fails');
end;

d=load('data');
dat = [d.dw d.dx d.dy d.z];
X = oo_.SmoothedMeasurementErrors;
ME = [X.dw X.dx X.dy X.z];
if max(max(abs(dat(1000:1199,:)-S(:,[2:4 1])-ME))) > 1e-10;
   error('Test fails');
end;
