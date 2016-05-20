//This mod-file tests the value of the objective function in the presence of non-stationary variables

Var


p_U_Hat            

pi_U_Hat               
r_Hat              
r_Gap
r_Eff
realrate_U_Hat
realrate_U_Eff
realrate_U_Gap     


c_U_Gap            
c_U_Eff
c_U_Hat
g_U_Gap            
g_U_Eff
g_U_Hat

y_U_Gap            
y_U_Eff
y_U_Hat

l_U_Gap             
l_U_Eff
l_U_Hat

a_U_Hat            



p_A_Hat             
pc_A_Hat           

tt_A_Hat           
tt_A_Gap
tt_A_Eff

pi_A_Hat           
pic_A_Hat         
c_A_Gap             
c_A_Eff
c_A_Hat

g_A_Gap             
g_A_Eff
g_A_Hat

realrate_A_Hat      

l_A_Gap             
l_A_Eff
l_A_Hat

muW_A_Hat          
y_A_Gap             
y_A_Eff
y_A_Hat

a_A_Hat             
mc_A_Hat           

tau_A_Hat           
tau_A_Gap
tau_A_Eff

Kdebt_A        
                    





p_B_Hat             
pc_B_Hat            

tt_B_Hat            
tt_B_Gap
tt_B_Eff

pi_B_Hat           
pic_B_Hat          
c_B_Gap             
c_B_Eff
c_B_Hat

g_B_Gap             
g_B_Eff
g_B_Hat

realrate_B_Hat      

l_B_Gap             
l_B_Eff
l_B_Hat

muW_B_Hat          
y_B_Gap             
y_B_Eff
y_B_Hat

a_B_Hat             
mc_B_Hat          

tau_B_Hat          
tau_B_Gap
tau_B_Eff

Kdebt_B        


;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF EXOGENOUS VARIABLES%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Varexo

eps_muW_A       
eps_a_A         

eps_muW_B       
eps_a_B         

; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%DECLARATION OF PARAMETERS%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters
n               
alpha           
beta            
rho             
lambdaA         
lambdaB         
sigma          
chi             
psi             

epsilon         
muP             
gamma           
ell             
phi_GY          
phi_caps       

theta_A         
               
theta_B         
                
phi_A           
phi_B          


%%%%%%
//Rule Parameters

ruleR_U_y      
ruleR_U_pi     
ruleR_U_debt

ruleG_A_y
ruleG_A_tt
ruleG_A_debt


ruleG_B_y
ruleG_B_tt
ruleG_B_debt



ruleT_A_y
ruleT_A_tt
ruleT_A_debt



ruleT_B_y
ruleT_B_tt
ruleT_B_debt

%%%%%%



rho_a_A         
rho_a_B	
rho_muW_A      
rho_muW_B
rho_mp          



r_SS

muW_SS                  
tau_A_SS                
DebtOutput_A_SS         


tau_B_SS
DebtOutput_B_SS    


;


%%%%%%%%%%%%%%%%%%%%%%%%
%%SET PARAMETER VALUES%%
%%%%%%CALIBRATION%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
n = 0.5; 
alpha = 0.4;    
beta = 0.99 ;   
rho = -log(beta);
lambdaA = 1-(1-n)*alpha;
lambdaB = 1-n*alpha;
sigma = 0.4;   
chi = 3;    
psi = 0.4;    

epsilon = 11;  
muP = 1/(epsilon - 1);
gamma	= 4.5;
ell = 11;      

phi_GY = 0.25; 
phi_caps = alpha*(gamma-(1-alpha)*(-gamma+sigma));

theta_A = 0.75;     
theta_B = 0.75;     
phi_A = (1 - theta_A*beta)*(1 - theta_A)/(theta_A*(1 + epsilon*chi)); 
phi_B = (1 - theta_B*beta)*(1 - theta_B)/(theta_B*(1 + epsilon*chi));  




//Rule Parameters

ruleR_U_y = 0.5;     
ruleR_U_pi = 1.5;     
ruleR_U_debt = 0;


ruleG_A_y = 0;
ruleG_A_tt = 0;
ruleG_A_debt = -0.05;


ruleG_B_y = ruleG_A_y;
ruleG_B_tt = 0;
ruleG_B_debt = ruleG_A_debt;


ruleT_A_y = 0;
ruleT_A_tt = 0;
ruleT_A_debt = 0.05;


ruleT_B_y = ruleT_A_y;
ruleT_B_tt = 0;
ruleT_B_debt = ruleT_A_debt;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rho_a_A = 0.85;    
rho_a_B = 0.85;
rho_muW_A = 0;
rho_muW_B = 0;
rho_mp = 0;      



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% STEADY STATE VALUES %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
muW_SS = 1/(ell - 1);

r_SS = log(1/beta);

DebtOutput_A_SS = 4*0.6;
                      
DebtOutput_B_SS = 4*0.6; 


tau_A_SS = phi_GY + r_SS*DebtOutput_A_SS; 
tau_B_SS = phi_GY + r_SS*DebtOutput_B_SS; 







model(linear);


p_U_Hat = n*p_A_Hat + (1-n)*p_B_Hat; 
pi_U_Hat = n*pi_A_Hat + (1-n)*pi_B_Hat; 


c_U_Gap = c_U_Gap(+1)- sigma*(r_Gap - pi_U_Hat(+1));


realrate_U_Gap = realrate_U_Hat - realrate_U_Eff; 
realrate_U_Eff = ((1 + chi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*(a_U_Hat(+1) - a_U_Hat); 

realrate_U_Hat = r_Hat - pi_U_Hat(+1);


r_Eff = ((1 + chi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*(a_U_Hat(+1) - a_U_Hat);
r_Hat = r_Gap + r_Eff;



// Interest rate rule
r_Gap = ruleR_U_y*y_U_Gap(-1) + ruleR_U_pi*pi_U_Hat(-1) + ruleR_U_debt*(n*Kdebt_A(-1) + (1-n)*Kdebt_B(-1));




c_U_Gap = n*c_A_Gap + (1-n)*c_B_Gap; 
c_U_Eff = n*c_A_Eff + (1-n)*c_B_Eff;
c_U_Hat = n*c_A_Hat + (1-n)*c_B_Hat;

g_U_Gap = n*g_A_Gap + (1-n)*g_B_Gap; 
g_U_Eff = n*g_A_Eff + (1-n)*g_B_Eff;
g_U_Hat = n*g_A_Hat + (1-n)*g_B_Hat;


y_U_Gap = n*y_A_Gap + (1-n)*y_B_Gap; 
y_U_Eff = n*y_A_Eff + (1-n)*y_B_Eff;
y_U_Hat = n*y_A_Hat + (1-n)*y_B_Hat;

l_U_Gap = n*l_A_Gap + (1-n)*l_B_Gap; 
l_U_Eff = n*l_A_Eff + (1-n)*l_B_Eff;
l_U_Hat = n*l_A_Hat + (1-n)*l_B_Hat;

a_U_Hat = n*a_A_Hat + (1-n)*a_B_Hat; 


tt_A_Hat = (p_U_Hat - p_A_Hat)/(1-n);                   

pi_A_Hat = p_A_Hat - p_A_Hat(-1); 			
pic_A_Hat = pc_A_Hat - pc_A_Hat(-1); 		


pic_A_Hat = pi_A_Hat + (1-lambdaA)*(tt_A_Hat - tt_A_Hat(-1)); 
			
		

c_A_Eff = (((1 + chi)*sigma)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 - alpha)*(1 + chi)*sigma)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_A_Hat - a_U_Hat);

g_A_Eff = (((1 + chi)*psi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 + chi)*psi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_A_Hat - a_U_Hat);

y_A_Eff = (((1 + chi)*(phi_GY*psi + (1 - phi_GY)*sigma))/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 + chi)*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha))))/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_A_Hat - a_U_Hat);

l_A_Eff = y_A_Eff - a_A_Hat;


realrate_A_Hat = r_Hat - pic_A_Hat(+1);


tt_A_Gap = tt_A_Hat - tt_A_Eff;
tt_A_Eff = ((1 + chi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*((a_A_Hat - a_U_Hat)/(1-n)); 


y_A_Gap = l_A_Gap; 


pi_A_Hat = beta*pi_A_Hat(+1)+ phi_A*mc_A_Hat;

mc_A_Hat = tau_A_Gap/(1 - tau_A_SS) +(chi + 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))* y_A_Gap - (phi_GY /((1- phi_GY)*(phi_caps+sigma*(1-alpha))))*g_A_Gap + (1/(sigma*(1- phi_GY))- 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))* y_U_Gap - phi_GY *(1/(sigma*(1 - phi_GY)) - 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))*g_U_Gap;


tau_A_Hat = tau_A_Gap + tau_A_Eff;
tau_A_Eff = - (1 - tau_A_SS)*muW_A_Hat;


y_A_Gap = (1- phi_GY)*c_U_Gap + (1 - phi_GY)*(sigma*(1 - alpha) + phi_caps)*(1-n)*tt_A_Gap + phi_GY*g_A_Gap;



Kdebt_A - ((1 + r_SS)/r_SS)*(tau_A_SS - phi_GY)*r_Hat= (1/beta)*(Kdebt_A(-1) - ((1 + r_SS)/r_SS)*(tau_A_SS - phi_GY)*pi_A_Hat + (phi_GY*g_A_Gap - tau_A_SS*y_A_Gap - tau_A_Gap)) + (1/beta)*(phi_GY*g_A_Eff - tau_A_SS*y_A_Eff + (1 - tau_A_SS)*muW_A_Hat) + (1/beta)*((1 + r_SS)/r_SS)*(tau_A_SS - phi_GY)*alpha*(1-n)*(tt_A_Gap(-1) -(1/(1 + r_SS))*tt_A_Gap) + (1/beta)*((1 + r_SS)/r_SS)*(tau_A_SS - phi_GY)*alpha*(1-n)*(tt_A_Eff(-1) -(1/(1 + r_SS))*tt_A_Eff); 



//FISCAL RULES

g_A_Gap = ruleG_A_y*(y_A_Gap(-1) - y_U_Gap(-1)) + ruleG_A_tt*tt_A_Gap(-1) + ruleG_A_debt*Kdebt_A(-1); 

tau_A_Gap = ruleT_A_y*(y_A_Gap(-1) - y_U_Gap(-1)) + ruleT_A_tt*tt_A_Gap(-1)+ ruleT_A_debt*Kdebt_A(-1); 


c_A_Hat = c_A_Gap + c_A_Eff;
g_A_Hat = g_A_Gap + g_A_Eff;
y_A_Hat = y_A_Gap + y_A_Eff;
l_A_Hat = l_A_Gap + l_A_Eff;







tt_B_Hat = (p_U_Hat - p_B_Hat)/n;                  

pi_B_Hat = p_B_Hat - p_B_Hat(-1); 			
pic_B_Hat = pc_B_Hat - pc_B_Hat(-1); 		

pic_B_Hat = pi_B_Hat + (1-lambdaB)*(tt_B_Hat - tt_B_Hat(-1)); 
			

c_B_Eff = (((1 + chi)*sigma)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 - alpha)*(1 + chi)*sigma)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_B_Hat - a_U_Hat);

g_B_Eff = (((1 + chi)*psi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 + chi)*psi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_B_Hat - a_U_Hat);

y_B_Eff = (((1 + chi)*(phi_GY*psi + (1 - phi_GY)*sigma))/(1 + chi*(phi_GY*psi + (1 - phi_GY)*sigma)))*a_U_Hat + (((1 + chi)*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha))))/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*(a_B_Hat - a_U_Hat);

l_B_Eff = y_B_Eff - a_B_Hat;


realrate_B_Hat = r_Hat - pic_B_Hat(+1);


c_A_Gap = c_B_Gap + sigma*(1 - alpha)*((1-n)*tt_A_Gap-n*(tt_B_Gap));


tt_B_Gap = tt_B_Hat - tt_B_Eff;
tt_B_Eff = ((1 + chi)/(1 + chi*(phi_GY*psi + (1 - phi_GY)*(gamma + (sigma - gamma)*(1 - 2*alpha + alpha*alpha)))))*((a_B_Hat - a_U_Hat)/n); 


y_B_Gap = l_B_Gap; //l_B_Gap = y_B_Gap + z_B , up to a first order approximation z_B=0


pi_B_Hat = beta*pi_B_Hat(+1)+ phi_B*mc_B_Hat;

mc_B_Hat = tau_B_Gap/(1 - tau_B_SS) +(chi + 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))* y_B_Gap - (phi_GY /((1- phi_GY)*(phi_caps+sigma*(1-alpha))))*g_B_Gap + (1/(sigma*(1- phi_GY))- 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))* y_U_Gap - phi_GY *(1/(sigma*(1 - phi_GY)) - 1/((1- phi_GY)*(phi_caps+sigma*(1-alpha))))*g_U_Gap;


tau_B_Hat = tau_B_Gap + tau_B_Eff;
tau_B_Eff = - (1 - tau_B_SS)*muW_B_Hat;


y_B_Gap = (1- phi_GY)*c_U_Gap + (1 - phi_GY)*(sigma*(1 - alpha) + phi_caps)*n*tt_B_Gap + phi_GY*g_B_Gap;


Kdebt_B - ((1 + r_SS)/r_SS)*(tau_B_SS - phi_GY)*r_Hat= (1/beta)*(Kdebt_B(-1) - ((1 + r_SS)/r_SS)*(tau_B_SS - phi_GY)*pi_B_Hat + (phi_GY*g_B_Gap - tau_B_SS*y_B_Gap - tau_B_Gap)) + (1/beta)*(phi_GY*g_B_Eff - tau_B_SS*y_B_Eff + (1 - tau_B_SS)*muW_B_Hat) + (1/beta)*((1 + r_SS)/r_SS)*(tau_B_SS - phi_GY)*alpha*n*(tt_B_Gap(-1) -(1/(1 + r_SS))*tt_B_Gap) + (1/beta)*((1 + r_SS)/r_SS)*(tau_B_SS - phi_GY)*alpha*n*(tt_B_Eff(-1) -(1/(1 + r_SS))*tt_B_Eff); 




//FISCAL RULES

g_B_Gap = ruleG_B_y*(y_B_Gap(-1) - y_U_Gap(-1)) + ruleG_B_tt*tt_B_Gap(-1) + ruleG_B_debt*Kdebt_B(-1); // function of domestic inflation since government spending is oriented only to domestic production  

tau_B_Gap = ruleT_B_y*(y_B_Gap(-1) - y_U_Gap(-1)) + ruleT_B_tt*tt_B_Gap(-1) + ruleT_B_debt*Kdebt_B(-1); // function of domenstic inflation(!?)



c_B_Hat = c_B_Gap + c_B_Eff;
g_B_Hat = g_B_Gap + g_B_Eff;
y_B_Hat = y_B_Gap + y_B_Eff;
l_B_Hat = l_B_Gap + l_B_Eff;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Exogenous processes %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

a_A_Hat 	= rho_a_A* a_A_Hat(-1) + eps_a_A;
a_B_Hat 	= rho_a_B* a_B_Hat(-1) + eps_a_B;

muW_A_Hat 	= rho_muW_A*muW_A_Hat(-1) + eps_muW_A;
muW_B_Hat	= rho_muW_B*muW_B_Hat(-1) + eps_muW_B;


end;



shocks;

var eps_a_B; stderr 1;
var eps_a_A; stderr 1;

var eps_muW_B;  stderr 1;
var eps_muW_A;  stderr 1;

end;


steady; 
//check; //do not use this command under optimal policy!!!

//stoch_simul(order=1, irf=100);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimized Simple Rule (osr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   


optim_weights;
                                                    //for country A

pi_A_Hat (1/2)*n*(epsilon/phi_A);
c_A_Gap (1/2)*n*(1 - phi_GY)*(1/sigma + (1 - phi_GY)*chi);
g_A_Gap (1/2)*n*phi_GY*(1/psi + phi_GY*chi);
tt_A_Gap (1/2)*n*(1 - phi_GY)*(gamma*alpha*(alpha - 2) + 2*phi_caps + (1 - phi_GY)*(phi_caps^2)*chi)*((1-n)^2); 

c_A_Gap, g_A_Gap (1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*chi;
c_A_Gap, tt_A_Gap (1/2)*(1/2)*n*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*(1-n);
g_A_Gap, tt_A_Gap (1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*phi_caps*chi*(1-n);

g_A_Gap, c_A_Gap (1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*chi;
tt_A_Gap, c_A_Gap (1/2)*(1/2)*n*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*(1-n);
tt_A_Gap, g_A_Gap (1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*phi_caps*chi*(1-n);



                                                    //now, the same for country B
pi_B_Hat (1/2)*(1-n)*(epsilon/phi_B);
c_B_Gap (1/2)*(1-n)*(1 - phi_GY)*(1/sigma + (1 - phi_GY)*chi);
g_B_Gap (1/2)*(1-n)*phi_GY*(1/psi + phi_GY*chi);
tt_B_Gap (1/2)*(1-n)*(1 - phi_GY)*(gamma*alpha*(alpha - 2) + 2*phi_caps + (1 - phi_GY)*(phi_caps^2)*chi)*(n^2); 

c_B_Gap, g_B_Gap (1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*chi;
c_B_Gap, tt_B_Gap (1/2)*(1/2)*(1-n)*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*n;
g_B_Gap, tt_B_Gap (1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*phi_caps*chi*n;

g_B_Gap, c_B_Gap (1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*chi;
tt_B_Gap, c_B_Gap (1/2)*(1/2)*(1-n)*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*n;
tt_B_Gap, g_B_Gap (1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*phi_caps*chi*n;

end;



osr_params 


ruleR_U_y      
ruleR_U_pi    
ruleR_U_debt




ruleG_A_y
ruleG_A_tt
ruleG_A_debt



ruleG_B_y
ruleG_B_tt
ruleG_B_debt



ruleT_A_y
ruleT_A_tt
ruleT_A_debt



ruleT_B_y
ruleT_B_tt
ruleT_B_debt

;

osr(irf=5,maxit=10000, nograph);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

asymptotic_loss=((1/2)*n*(epsilon/phi_A))...
        *oo_.var(strmatch('pi_A_Hat',M_.endo_names,'exact'),strmatch('pi_A_Hat',M_.endo_names,'exact'))+...
            ((1/2)*n*(1 - phi_GY)*(1/sigma + (1 - phi_GY)*chi))...
        *oo_.var(strmatch('c_A_Gap',M_.endo_names,'exact'),strmatch('c_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*n*phi_GY*(1/psi + phi_GY*chi))...
        *oo_.var(strmatch('g_A_Gap',M_.endo_names,'exact'),strmatch('g_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*n*(1 - phi_GY)*(gamma*alpha*(alpha - 2) + 2*phi_caps + (1 - phi_GY)*(phi_caps^2)*chi)*((1-n)^2))...
        *oo_.var(strmatch('tt_A_Gap',M_.endo_names,'exact'),strmatch('tt_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*chi)...
        *oo_.var(strmatch('c_A_Gap',M_.endo_names,'exact'),strmatch('g_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*(1-n))...
        *oo_.var(strmatch('c_A_Gap',M_.endo_names,'exact'),strmatch('tt_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*phi_caps*chi*(1-n))...
        *oo_.var(strmatch('g_A_Gap',M_.endo_names,'exact'),strmatch('tt_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*chi)...
        *oo_.var(strmatch('g_A_Gap',M_.endo_names,'exact'),strmatch('c_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*(1-n))...
        *oo_.var(strmatch('tt_A_Gap',M_.endo_names,'exact'),strmatch('c_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*n*2*phi_GY*(1 - phi_GY)*phi_caps*chi*(1-n))...
        *oo_.var(strmatch('tt_A_Gap',M_.endo_names,'exact'),strmatch('g_A_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1-n)*(epsilon/phi_B))...
        *oo_.var(strmatch('pi_B_Hat',M_.endo_names,'exact'),strmatch('pi_B_Hat',M_.endo_names,'exact'))+...
            ((1/2)*(1-n)*(1 - phi_GY)*(1/sigma + (1 - phi_GY)*chi))...
        *oo_.var(strmatch('c_B_Gap',M_.endo_names,'exact'),strmatch('c_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1-n)*phi_GY*(1/psi + phi_GY*chi))...
        *oo_.var(strmatch('g_B_Gap',M_.endo_names,'exact'),strmatch('g_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1-n)*(1 - phi_GY)*(gamma*alpha*(alpha - 2) + 2*phi_caps + (1 - phi_GY)*(phi_caps^2)*chi)*(n^2))...
        *oo_.var(strmatch('tt_B_Gap',M_.endo_names,'exact'),strmatch('tt_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*chi)...
        *oo_.var(strmatch('c_B_Gap',M_.endo_names,'exact'),strmatch('g_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*n)...
        *oo_.var(strmatch('c_B_Gap',M_.endo_names,'exact'),strmatch('tt_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*phi_caps*chi*n)...
        *oo_.var(strmatch('g_B_Gap',M_.endo_names,'exact'),strmatch('tt_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*chi)...
        *oo_.var(strmatch('g_B_Gap',M_.endo_names,'exact'),strmatch('c_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*(1 - phi_GY)*(2*alpha + 2*(1 - phi_GY)*phi_caps*chi)*n)...
        *oo_.var(strmatch('tt_B_Gap',M_.endo_names,'exact'),strmatch('c_B_Gap',M_.endo_names,'exact'))+...
            ((1/2)*(1/2)*(1-n)*2*phi_GY*(1 - phi_GY)*phi_caps*chi*n)...
        *oo_.var(strmatch('tt_B_Gap',M_.endo_names,'exact'),strmatch('g_B_Gap',M_.endo_names,'exact'));

if abs(asymptotic_loss-oo_.osr.objective_function)>1e-10
    error('Objective Function with Non-stationary Variables is wrong')
end
        
        