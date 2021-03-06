clear variables
clc
close all
%% add folders to the path
addpath([ pwd '\ScenarioBounds'])
addpath([ pwd '\datasets']) 
addpath([ pwd '\IPM_functions']) 
%% START:
%% 1) Generate samples from an unknown process y=f(x)+uncertainty
Nsamples=100;
X=unifrnd(-3,10,[1,Nsamples]);
Y= sin(X.*rand([1,Nsamples])).^2+X.*exprnd(2,[1,Nsamples]);  
%% 2) construct IPM object
Xsupport=[min(X),max(X)]; % define data support X
Ysupport=[min(Y),max(Y)]; % define data support Y

options=optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-4,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e4,...
    'Display','iter','Algorithm','sqp'); 

IPM=IPM_Model('degreeup', 6,'degreelow',6 ,...
    'basis','poly', 'Length',[],...
    'Xsup',Xsupport, 'Ysup',Ysupport,'options',options);
    
%% 3) run different training programs
%% %% %%      1) hard-constrained solver  
%% %% %%      1.1)   Reliability bounds ( A-priori, data-independent, generalization error bounds )
beta_target=10^-5;
Ndiscarded=0; % number of discarded samples
Nd=12; %  number of optimization variables 
epsilon_apriori=getepsilon_ConvexDiscard(Nsamples,beta_target,Ndiscarded,Nd);
RelBound_apriori=[1-epsilon_apriori,1]; % minimum reliability for the IPM

%% %% %%      1.2) optimize the model
Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);
Area(1)=Design.Area; % area between the IPM bounds (to be minimized)
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);  % points to plot
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y)); % de-normalize 

%% %% %%      1.3) plot results
subplot(2,2,1)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on; grid on;
scatter(X,Y,'.b')
title('hard-constrained')
%% %% %%      1.4) Reliability bounds ( A-posteriori, data-deopendent, generalization error bounds )
SN=Design.Generalization.NSupportConstraints;
epsilon_aposteriori=getWaitandJudgeEpsilon_fast(SN,Nsamples,beta_target);
RelBound=[1-epsilon_aposteriori,1];

%% %% %%     2) hard-constrained with discarded samples
Ndiscarded=round(Nsamples/10); % remove 10% of the saples
Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,Ndiscarded);
Area(2)=Design.Area;

X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)'; % optimize the model
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
%% %% %%     2.1) plot
subplot(2,2,2)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on;  grid on;
scatter(X,Y,'.b')
title('hard-constrained discarded')

%% %% %%     2.2) Reliability bounds ( A-posteriori) 
SN=Design.Generalization.NSupportConstraints+Ndiscarded;
epsilon_aposteriori=getWaitandJudgeEpsilon_fast(SN,Nsamples,beta_target);
RelBound=[1-epsilon_aposteriori,1];

%% 3) minimax layer
Design=IPM.DesingIPM_MinMaxLayer(X,Y);
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)'; % optimize the model
I_predict1=[I_predict(:,1)-Design.EmpiricalCosts(1) I_predict(:,2)+Design.EmpiricalCosts(1)]; % maximum empirical layer
I_predict_denorm1= IPM.De_Normalize(I_predict1,min(Y),max(Y));
I_predict50=[I_predict(:,1)-Design.EmpiricalCosts(50) I_predict(:,2)+Design.EmpiricalCosts(50)]; % maximum empirical layer
I_predict_denorm50= IPM.De_Normalize(I_predict50,min(Y),max(Y)); 

%% 3.1) plot
subplot(2,2,3)
plot(X_linpace,I_predict_denorm1,'k','LineWidth',1.5);
hold on;  grid on;
plot(X_linpace,I_predict_denorm50,'--k','LineWidth',1.5);
scatter(X,Y,'.b')
title('hard-constrained minmax')

%% 4) soft constrained layer
CostofViolations=0.1;
Design=IPM.DesingIPM_SoftConstraied(X,Y,CostofViolations); % optimize the model
 
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));

%% 4.1) plot
subplot(2,2,4)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on; grid on;
scatter(X,Y,'.b')
title('soft-constrained')
