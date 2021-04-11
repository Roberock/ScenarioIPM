function [c,ceq,g_val]=IPM_SoftConstraint(theta,Y,X,Degree_up,Degree_low,BasisType,Len)
% theta=[u0,u1....udu,l0,l1,...,ldl,s1,s2,....,sN] coefficents of the bounds
% and slack varibales si
% si={0,1} binary...if si=0 then ith constraint disappear form the program
% Nout= number of outliers to be removed
% Y dependent variable samples (N samples)
% X explanatory variable samples (N samples)
% Degree_up=du number of coefficents for the upper bound
% Degree_low=dl number of coefficents for the lower bound
% BasisType=type of basis functions (e.g. 'poly', 'fourier' etc)

Deg=Degree_up+Degree_low;
N=length(Y); % number of samples 
 
fu=compute_Basis(X,theta(1:Degree_up),Degree_up,BasisType,Len);
fl=compute_Basis(X,theta(Degree_up+1:Deg),Degree_low,BasisType,Len);
% g_val=max(fl-Y,Y-fu);
g_val=[];

Zeta=theta(Deg+1:Deg+N); 
c_lb = (fl-Y-Zeta)' ; % yi<=fu(xi)  for all i
c_ub = (Y-fu-Zeta)'; % yi>=fl(xi)  for all i
 
X_support=linspace(min(X),max(X),5000);
fu= compute_Basis(X_support,theta(1:Degree_up),Degree_up,BasisType,Len);
fl= compute_Basis(X_support,theta(Degree_up+1:Deg),Degree_low,BasisType,Len);
c_upperdominance=max(fl-fu); % fu>=fl


 c=[c_lb;c_ub;c_upperdominance];
%c=[max(c_lb,c_ub); c_upperdominance]; 
ceq=[]; % no equality constraints

end

