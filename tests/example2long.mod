var y, c, k, a, h, b, w;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;

phi   = 0.1;

model;
c*theta*h^(1+psi)=(1-alpha)*y;
k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(h^(1-alpha));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + e;
b = tau*a(-1)+rho*b(-1) + u;
w = y
@#for i in 1:50
  + beta*(y(@{i})
@#endfor
@#for i in 1:50
  )
@#endfor
;
end;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
h = 0.29175631001732;
k = 11.08360443260358;
w = y*(1-beta^51)/(1-beta);
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e; stderr 0.009;
var u; stderr 0.009;
var e, u = phi*0.009*0.009;
end;

stoch_simul;

if ~exist('example1long_results.mat','file');
   error('example1long must be run first');
end;

oo1 = load('example1long_results','oo_');

dr0 = oo1.oo_.dr;
dr = oo_.dr;

d_ghx = max(abs(dr0.ghx(:) - dr.ghx(:)));
d_ghu = max(abs(dr0.ghu(:) - dr.ghu(:)));
d_ghxx = max(abs(dr0.ghxx(:) - dr.ghxx(:)));
d_ghxu = max(abs(dr0.ghxu(:) - dr.ghxu(:)));
d_ghuu = max(abs(dr0.ghuu(:) - dr.ghuu(:)));
d_ghs2 = max(abs(dr0.ghs2 - dr.ghs2));

skipline()
disp(sprintf('ghx max. abs.  diff. is %s.', num2str(d_ghx)))
disp(sprintf('ghu max. abs.  diff. is %s.', num2str(d_ghu)))
disp(sprintf('ghxx max. abs.  diff. is %s.', num2str(d_ghxx)))
disp(sprintf('ghxu max. abs.  diff. is %s.', num2str(d_ghxu)))
disp(sprintf('ghuu max. abs.  diff. is %s.', num2str(d_ghuu)))
disp(sprintf('ghs2 max. abs.  diff. is %s.', num2str(d_ghs2)))
skipline()

epsilon = 1e-12;

if d_ghx>epsilon
   error('error in ghx')
end

if d_ghu>epsilon
   error('error in ghu')
end

if d_ghxx>epsilon
   error('error in ghxx')
end

if d_ghxu>epsilon
   error('error in ghxu')
end

if d_ghuu>epsilon
   error('error in ghuu')
end

if d_ghs2>epsilon
   error('error in ghs2')
end