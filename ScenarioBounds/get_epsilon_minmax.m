function   [Ub,Lb,Confidence]=get_epsilon_minmax(X,nd,beta_targ)
% X = empirical cost 
% beta_targ = target confidence level
% nd = number of optimization variables
Xs=sort(X,'ascend'); %(ordered statisitc of the empirical costs)
N=length(Xs); % number of samples
Lb=zeros(length(Xs),1); % Lower bound
Ub=ones(length(Xs),1);  % Upper bound
for x=1:N 
    [epsL, epsU] = epsLU_fast(N-x,N,beta_targ);
    %epsU = getWaitandJudgeEpsilon_fast(N-x+1,N,beta_targ);
    %epsL = getWaitandJudgeEpsilon_fast(N-x,N,beta_targ);
    Lb(x)= 1-epsU;  
     if x>N-nd+1
    Ub(x)= 1;  
     else
    Ub(x)= 1-epsL;  
     end 
end
% Lb=[0 Lb Lb(end)];
% Ub=[Ub(1) Ub(1:N-nd) ones(1,nd+1)];
Confidence=1-2*beta_targ;

end