verbatim;
if options_.loglinear
Y_obs=log(Y_obs(options_.first_obs:end));
P_obs=log(P_obs(options_.first_obs:end));
junk2=log(junk2(options_.first_obs:end));
else
Y_obs=Y_obs(options_.first_obs:end);
P_obs=P_obs(options_.first_obs:end);
junk2=junk2(options_.first_obs:end);    
end
end;
