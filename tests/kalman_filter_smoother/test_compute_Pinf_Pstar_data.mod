// 3 independent local trend models
// test that the three local trends processes are still independent in the smoother

var	
	beta1, beta2, beta3,
	mu1, mu2, mu3,
	psi1, psi2, psi3,
	tren1, tren2, tren3,
	slop1, slop2, slop3,
	cyc1, cyc2, cyc3,
	v1, v2, v3;

varexo
	e_irr1, e_irr2, e_irr3,
	e_lev1, e_lev2, e_lev3,
	e_slp1, e_slp2, e_slp3,
	e_cyc1, e_cyc2, e_cyc3;

parameters
	rho,
	sdirr1, sdirr2, sdirr3,
	sdlev1, sdlev2, sdlev3,
	sdslp1, sdslp2, sdslp3,
	sdcyc1, sdcyc2, sdcyc3;
	
	rho    = 0.75;   
	sdirr1 = 0.005;  
	sdirr2 = 0.005;  
	sdirr3 = 0.005;  
	sdlev1 = 0.0001; 
	sdlev2 = 0.0001; 
	sdlev3 = 0.0001; 
	sdslp1 = 0.0001; 
	sdslp2 = 0.0001; 
	sdslp3 = 0.0001;
	sdcyc1 = 0.005; 
	sdcyc2 = 0.005; 
	sdcyc3 = 0.005; 

model(linear);

	mu1	= mu1(-1) + beta1(-1) + sdlev1*e_lev1;
	beta1 =	beta1(-1) + sdslp1*e_slp1;
	psi1 = rho*psi1(-1) + sdcyc1*e_cyc1;
	
	mu2	= mu2(-1) + beta2(-1) + sdlev2*e_lev2;
	beta2 =	beta2(-1) + sdslp2*e_slp2;
	psi2 = rho*psi2(-1) + sdcyc2*e_cyc2;
	
	mu3	= mu3(-1) + beta3(-1) + sdlev3*e_lev3;
	beta3 =	beta3(-1) + sdslp3*e_slp3;
	psi3 = rho*psi3(-1) + sdcyc3*e_cyc3;
	
	tren1 =	mu1(-1);
	tren2 =	mu2(-1);
	tren3 =	mu3(-1);
	
	slop1 =	beta1(-2);
	slop2 =	beta2(-2);
	slop3 =	beta3(-2);
	
	cyc1 = psi1(-1);
	cyc2 = psi2(-1);
	cyc3 = psi3(-1);
	
	v1 = tren1 + cyc1 + sdirr1*e_irr1;
	v2 = tren2 + cyc2 + sdirr2*e_irr2;
	v3 = tren3 + cyc3 + sdirr3*e_irr3;
	
end;

shocks;
	var e_irr1; stderr 1;
	var e_irr2; stderr 1;
	var e_irr3; stderr 1;
	var e_lev1; stderr 1;
	var e_lev2; stderr 1;
	var e_lev3; stderr 1;
	var e_slp1; stderr 1;
	var e_slp2; stderr 1;
	var e_slp3; stderr 1;
	var e_cyc1; stderr 1;
	var e_cyc2; stderr 1;
	var e_cyc3; stderr 1;
end;

stoch_simul(order=1,irf=20,periods=500);

save Data.mat v1 v2 v3;
