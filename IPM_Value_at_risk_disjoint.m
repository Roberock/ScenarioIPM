function [c,ceq]=IPM_Value_at_risk_disjoint(x,Y,X,Degree_up,Degree_low,BasisType,Alpa_level,Std)
% this time we assume error in y (gaussian error)
% Y=expected value of the measurements
% X= ordinata (fixed)
% x= coefficient of the polynomial basys
% Degree_up,Degree_low degrees of the upper and lower poly bounds
% Alpa_level  = alpha level we want to compute CVAR for e.g. 0.5 , 0.95 
% Std = error in the y measurement
ceq=[]; 
fl = compute_Basis(X,x(Degree_up+1:Degree_up+Degree_low),Degree_low,BasisType) ;
fu = compute_Basis(X,x(1:Degree_up),Degree_up,BasisType); 
ValatRisk=@(mu,sig,alpha)  mu+norminv(alpha)*sig;

% generate samples 
c_lb= fl-ValatRisk(Y,Std(X,Y),1-Alpa_level) ; % compute lower bound:  VAR^-_alpha(Y)>f_lower  
c_ub=ValatRisk(Y,Std(X,Y),Alpa_level)-fu;% compute upper bound: VAR^+_alpha(Y)<f_upper 
c_upperdominance=fl-fu; % fu>=fl 
c=[c_lb';c_ub';c_upperdominance']; 
end
