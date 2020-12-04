clear variables;
clc;
close all;
%% Comparison of the bounds
addpath([pwd '\ScenarioBounds'])
Nd=100;
SN=95;
N=5000;
N_discard=5;
bet=10^-6;

N_linspace=[200:10:5000];
for i=1:length(N_linspace)
    N=N_linspace(i);
    Epsilon_apriori=getepsilon_ConvexDiscard(N,bet,0,Nd);
    Epsilon_apriori_discard=getepsilon_ConvexDiscard(N,bet,N_discard,Nd);
    Epsilon_aposteriori = getWaitandJudgeEpsilon_fast(SN,N,bet);
    Epsilon_aposteriori_discard = getWaitandJudgeEpsilon_fast(SN+N_discard,N,bet);
    [epsL, epsU] = epsLU_fast(SN,N,bet);
    epsilon_nonconvex=getConfidence_nonconvex(SN,N,bet);
    %[Ub,Lb,Confidence]=get_epsilon_minmax(X,nd,bet);
    Results(i,:)=[Epsilon_apriori , Epsilon_apriori_discard , Epsilon_aposteriori , Epsilon_aposteriori_discard , epsL, epsU ,epsilon_nonconvex];
end

%% Compare results of the IPM programs
N_train=2000; % available scenarios (for training the IPM)
N_test= 5*10^4; % scenarios for testing
% Generate samples for training
DGMtype = 2 ;% the data geneating mechanism
[X,Y,X_test,Y_test,Xsupport,Ysupport]= DGMs(DGMtype,N_train,N_test);
beta_target=10^-6;

options=optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-6,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e4,...
    'Display','off','Algorithm','sqp');

% build model
IPM=IPM_Model('degreeup', 80,'degreelow', 80, 'basis','rbf', 'Length',45,...
    'Xsup',Xsupport, 'Ysup',Ysupport,'options',options);

for TesProg=1:4
    
    
    if TesProg==1 % Hard-constrained no exceptions
        Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);% Desing Method with no uncertainty 
    elseif TesProg==2  % Hard-constrained with 10 discarded samples
        Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,10);% Desing Method with no uncertainty
    elseif TesProg==3  % Hard-constrained with 100 discarded samples
        Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,100);% Desing Method with no uncertainty
    elseif TesProg==4 % Min-max layer
        Design=IPM.DesingIPM_MinMaxLayer(X,Y);% Desing Method with no uncertainty
    elseif TesProg==5 % Soft-constrained with rho=10^5
        Design=IPM.DesingIPM_SoftConstraied(X,Y,10^5);
    elseif TesProg==6 % Soft-constrained with rho=10^2
        Design=IPM.DesingIPM_SoftConstraied(X,Y,10^2);
    end
 
    
    
end