var a b;
varexo epsil ;
parameters betta;
betta = 0.97;

model;
a = betta*a(+1)+epsil;
b = a;
end;

initval;
a=0;
end;
steady;

shocks;
var epsil; stderr .2;
end;

steady;
check;

stoch_simul(periods=0, irf=30, order=1);
stoch_simul(periods=2000, irf=30, order=1);
stoch_simul(periods=0, irf=30, order=2);
stoch_simul(periods=2000, irf=30, order=2);
stoch_simul(periods=0, irf=30, order=3);
stoch_simul(periods=2000, irf=30, order=3);


stoch_simul(periods=0, irf=30, order=1,hp_filter=1600);
stoch_simul(periods=2000, irf=30, order=1,hp_filter=1600);
stoch_simul(periods=0, irf=30, order=2,hp_filter=1600);
stoch_simul(periods=2000, irf=30, order=2,hp_filter=1600);
stoch_simul(periods=0, irf=30, order=3,hp_filter=1600);
stoch_simul(periods=2000, irf=30, order=3,hp_filter=1600);
