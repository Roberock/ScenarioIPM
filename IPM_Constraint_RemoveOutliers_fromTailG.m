function [c,ceq]=IPM_Constraint_RemoveOutliers_fromTailG(theta,Y,X,Nout,Degree_up,Degree_low,BasisType,Len)
% theta=[u0,u1....udu,l0,l1,...,ldl,s1,s2,....,sN] coefficents of the bounds
% and slack varibales si
% si={0,1} binary...if si=0 then ith constraint disappear form the program
% Nout= number of outliers to be removed
% Y dependent variable samples (N samples)
% X explanatory variable samples (N samples)
% Degree_up=du number of coefficents for the upper bound
% Degree_low=dl number of coefficents for the lower bound
% BasisType=type of basis functions (e.g. 'poly', 'fourier' etc)
%%
N=length(Y); % number of samples
ceq=[];
fl = compute_Basis(X,theta(Degree_up+1:Degree_up+Degree_low),Degree_low,BasisType,Len) ;
fu = compute_Basis(X,theta(1:Degree_up),Degree_up,BasisType,Len);
c_lb = (fl-Y);  % f_lower(xi)=<y_i 
c_ub = (Y-fu);  % f_upper(xi)>=y_i 
%C=max(c_ub,c_lb);
%Csorted=sort(C,'ascend');
%C=Csorted(1:(N-Nout))';
c_lbsorted=sort(c_lb,'ascend');
c_ubsorted=sort(c_ub,'ascend');
C=[c_lbsorted(1:N-round(Nout/2))';c_ubsorted(1:N-round(Nout/2))'];
 
X_linspace=linspace(min(X),max(X),5000);
fl = compute_Basis(X_linspace,theta(Degree_up+1:Degree_up+Degree_low),Degree_low,BasisType,Len) ;
fu = compute_Basis(X_linspace,theta(1:Degree_up),Degree_up,BasisType,Len);
c_upperdominance=max((fl-fu));
c=[C;c_upperdominance'];
 
%ceq_Number_of_outliers=sum(SVal)-Nout;  % sum(s) = Nout 
% Equality constraints
%  ceq_one_or_zero=theta(Degtot+1:Degtot+N).*(1-theta(Degtot+1:Degtot+N)); % s(1-s)=0 ....s=1 or s=0
% ceq= [ceq_one_or_zero ceq_Number_of_outliers]; % we must keep only N-Nout constraints


end
