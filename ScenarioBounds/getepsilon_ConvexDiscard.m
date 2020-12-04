function Epsilon_target=getepsilon_ConvexDiscard(N,beta_targ,k,Nd)
%% This function gives a-priori sample-and-discard upper bound on the risk
% it is applicable prior solving a convex scenario optimization program

% reference
% [] [] Campi, M.C., Garatti, S. A Sampling-and-Discarding Approach to
% Chance-Constrained Optimization: Feasibility and Optimality.
% J Optim Theory Appl 148, 257–280 (2011). https://doi.org/10.1007/s10957-010-9754-6

% beta_target: target confidence level
% k: Number of discarded samples
% N: number of samples
% Nd optimization variables

%% Example: how to use?
% N=50000; % 50k samples
% k=30; % we discard 30 befoe running the optimize  e.g. the samples making the data non-linearly separable
% beta_target=10^-6;
% Nd= 15; % number of desing variables (number of hyperparameters)
% for a linear SVM the Nd is equal to the number of featrues (+1 due to the bias ?)
% Epsilon_targt=getepsilon_ConvexDiscard(N,beta_target,k,Nd)
%% NOTE N musst be > Nd
if Nd>=N || k>=N
    warning('getConfidence_ConvexDiscard(N,beta_targ,k,Nd), we need Nd<N and k<N')
    Epsilon_target=NaN;
    return
end
%% STAT
epsilon=linspace(10^-5,1-10^-5,10^3);
Beta=0;  % confidence parameter
for j=0:(k+Nd-1)
    a=N-j+1;
    b=j+1;
    Beta=Beta+betapdf(1-epsilon,a,b);
end
% for j=0:(k+Nd-1)
%     beta=beta+nchoosek(N,j).*epsilon.^j.*(1-epsilon).^(N-j);
% end
 
Beta=Beta*1/((N+1)*beta(N-k+1,k+1)); % equivalent Beta= Beta*nchoosek(k+Nd-1,k) but more stable 
Epsilon_target=epsilon(find(Beta<=beta_targ,1,'first'));
end

 
% BETAPdf=0;
% for j=0:(Nd-1)
%     a=N-j+1
%     b=j+1
%     BETAPdf=BETAPdf+betapdf(1-epsilon,a,b);
% end