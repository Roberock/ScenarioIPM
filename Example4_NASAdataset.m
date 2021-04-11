clear variables
clc
close all
%% add folders to the path
addpath([ pwd '\ScenarioBounds'])
addpath([ pwd '\datasets']) 
addpath([ pwd '\IPM_functions']) 

%% NASA Example
load('NASA_data_example.mat') 
% Data from  https://uqtools.larc.nasa.gov/nasa-uq-challenge-problem-2020/ 
%  
% [] The NASA langley challenge on optimization under uncertainty, LG
% Crespo, SP Kenny, Mechanical Systems and Signal Processing 152, 107405

X=T(:)'; Y=ReferenceSignal(:)'; 

Nsamples=2000; % number of samples
Perm=randperm(length(X), Nsamples); 
X=X(Perm); Y=Y(Perm);  % take a subset of samples for training

%% construct IPM object
Xsupport=[min(X),max(X)];
Ysupport=[min(Y),max(Y)];
options=optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-4,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e4,...
    'Display','iter','Algorithm','sqp');

IPM=IPM_Model('degreeup',14,'degreelow',14 ,...
    'basis','rbf', 'Length',20,...
    'Xsup',Xsupport, 'Ysup',Ysupport,'options',options);
%% optimize IPM
Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));


%% plot 
figure(1)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on; grid on;
scatter(X,Y,'*r')
plot(T(1,:),ReferenceSignal,':b','LineWidth',0.1)
xlabel('Time x(t) [sec]')
ylabel('System controlled response z_2(t)')
legend('f_u(x;\theta)','f_l(x;\theta)','Training samples','Validation Samples') 
 
%% failure bounds for discrete-time (point) realization
beta_target=10^-6;
SN=Design.Generalization.NSupportConstraints;
epsilon=getWaitandJudgeEpsilon_fast(SN,Nsamples,beta_target);
RelBound_hard=[1-epsilon,1];
[Reliability_hard,Y_target_up,Y_target_low] = IPM.Validation_Samples(T(:)',ReferenceSignal(:)',X_linpace,I_predict_denorm);
 
%% analyze failure per time-series (process)
IDX_failedprocess=reshape(Reliability_hard.IdxFailed,100,5001);
PF_perProcess_hard=mean(any(IDX_failedprocess(:,50:end),2));
 