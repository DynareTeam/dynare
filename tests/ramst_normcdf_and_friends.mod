// Tests the normcdf(), normpdf() and erf() functions, in the static M-file, and in a dynamic C-file

var c k t u v sin_x cos_x tan_x asin_y acos_y atan_y;
varexo x y;

parameters alph gam delt bet aa;
alph=0.5;
gam=0.5;
delt=0.02;
bet=0.05;
aa=0.5;

model(use_dll);
c + k - aa*x*k(-1)^alph - (1-delt)*k(-1);
c^(-gam) - (1+bet)^(-1)*(aa*alph*x(+1)*k^(alph-1) + 1 - delt)*c(+1)^(-gam);
t = normcdf(x, 2, 3);
u = normpdf(x, 1, 0.5);
v = erf(x);
sin_x = sin(x);
cos_x = cos(x);
tan_x = tan(x);
asin_y = asin(y);
acos_y = acos(y);
atan_y = atan(y);
end;

initval;
x = 1;
y = 0;
k = ((delt+bet)/(1.0*aa*alph))^(1/(alph-1));
c = aa*k^alph-delt*k;
t = 0;
u = 0;
v = 0;
end;

steady;

check;

shocks;
var x;
periods 1;
values 1.2;
var y;
periods 1,2;
values 0.5, 1;
end;

simul(periods=20);

if(abs(oo_.steady_state(5) - erf(1)) > 1e-10)
   error('Test failed in static M-file for erf')
end

if (abs(oo_.endo_simul(5, 2) - erf(1.2)) > 1e-10)
   error('Test failed in dynamic M-file for erf')
end

if(abs(oo_.steady_state(4) - normpdf(1, 1, 0.5)) > 1e-10)
   error('Test failed in static M-file for normpdf')
end

if (abs(oo_.endo_simul(4, 2) - normpdf(1.2, 1, 0.5)) > 1e-10)
   error('Test failed in dynamic M-file for normpdf')
end

if (abs(oo_.steady_state(3) - normcdf(1, 2, 3)) > 1e-10)
   error('Test failed in static M-file for normcdf')
end

if (abs(oo_.endo_simul(3, 2) - normcdf(1.2, 2, 3)) > 1e-10)
   error('Test failed in dynamic M-file for normcdf')
end


sin_x_pos=strmatch('sin_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(sin_x_pos) - sin(1)) > 1e-10)
   error('Test failed in static M-file for sin')
end

if (abs(oo_.endo_simul(sin_x_pos, 2) - sin(1.2)) > 1e-10)
   error('Test failed in dynamic M-file for sin')
end

cos_x_pos=strmatch('cos_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(cos_x_pos) - cos(1)) > 1e-10)
   error('Test failed in static M-file for cos')
end

if (abs(oo_.endo_simul(cos_x_pos, 2) - cos(1.2)) > 1e-10)
   error('Test failed in dynamic M-file for cos')
end

tan_x_pos=strmatch('tan_x',M_.endo_names,'exact');
if(abs(oo_.steady_state(tan_x_pos) - tan(1)) > 1e-10)
   error('Test failed in static M-file for tan')
end

if (abs(oo_.endo_simul(tan_x_pos, 2) - tan(1.2)) > 1e-10)
   error('Test failed in dynamic M-file for tan')
end

asin_y_pos=strmatch('asin_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(asin_y_pos) - asin(0)) > 1e-10)
   error('Test failed in static M-file for asin')
end

if (abs(oo_.endo_simul(asin_y_pos,1) - asin(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(asin_y_pos, 2) - asin(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(asin_y_pos, 3) - asin(1)) > 1e-10)
   error('Test failed in dynamic M-file for asin')
end

acos_y_pos=strmatch('acos_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(acos_y_pos) - acos(0)) > 1e-10)
   error('Test failed in static M-file for acos')
end

if (abs(oo_.endo_simul(acos_y_pos,1) - acos(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(acos_y_pos, 2) - acos(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(acos_y_pos, 3) - acos(1)) > 1e-10)
   error('Test failed in dynamic M-file for acos')
end

atan_y_pos=strmatch('atan_y',M_.endo_names,'exact');
if(abs(oo_.steady_state(atan_y_pos) - atan(0)) > 1e-10)
   error('Test failed in static M-file for atan')
end

if (abs(oo_.endo_simul(atan_y_pos,1) - atan(0)) > 1e-10) || ...
        (abs(oo_.endo_simul(atan_y_pos, 2) - atan(0.5)) > 1e-10) || ...
        (abs(oo_.endo_simul(atan_y_pos, 3) - atan(1)) > 1e-10)
   error('Test failed in dynamic M-file for atan')
end