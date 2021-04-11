function [c,ceq]=IPM_constraints_no_exceptions(theta,Y,X,Degree_up,Degree_low,BasisType,Len)
% theta=[u0,u1....udu,l0,l1,...,ldl] coefficents of the IPM lower and upper bounds
% Y dependent variable samples (N samples)
% X explanatory variable samples (N samples)
% Degree_up=du number of coefficents for the upper bound 
% Degree_low=dl number of coefficents for the lower bound 
% BasisType=type of basis functions (e.g. 'poly', 'fourier' etc)
ceq=[];  
fl = compute_Basis(X,theta(Degree_up+1:Degree_up+Degree_low),Degree_low,BasisType,Len) ;
fu = compute_Basis(X,theta(1:Degree_up),Degree_up,BasisType,Len); 
c_lb = (fl-Y);  % f_lower(xi)=<y_i 
c_ub = (Y-fu);  % f_upper(xi)>=y_i 
C=[c_lb';c_ub']; %C=max(max(c_lb,c_ub)); 

X_linspace=linspace(min(X),max(X),5000);
fl = compute_Basis(X_linspace,theta(Degree_up+1:Degree_up+Degree_low),Degree_low,BasisType,Len) ;
fu = compute_Basis(X_linspace,theta(1:Degree_up),Degree_up,BasisType,Len); 
c_upperdominance=max((fl-fu)); % fu>=fl 
c=[C;c_upperdominance'];
end
