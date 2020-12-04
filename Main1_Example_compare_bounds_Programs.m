clear variables;
clc; close all;
%% Comparison of the bounds
addpath([pwd '\ScenarioBounds']) % add path with the .m for the scenario-based bound computations
Nd=100; SN=95;
N_discard=5;
bet=10^-6;
N_linspace=[400:50:5000];
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
figure(1)
plot(N_linspace,Results,'DisplayName','Results')
xlabel('number of samples N'); grid on;
ylabel('Bound \epsilon on the probability of contraint violation')
legend('a-priori no exceptions','a-priori with N_odiscarded samples','a-posteriori (Wait-and-judge)','a-posteriori with N_o discarded',...
    'lower \epsilon (soft-constrained)','upper \epsilon (soft-constrained)', 'non-convex')
%% Compare results of the IPM programs
N_train=1000; % available scenarios (for training the IPM)
N_test= 1*10^4; % scenarios for testing
% Generate samples  
DGMtype = 1 ;% the data geneating mechanism
[X,Y,X_test,Y_test,Xsupport,Ysupport]= DGMs(DGMtype,N_train,N_test);
 
if DGMtype==1 || DGMtype==3
Deg_u=10;Deg_l=Deg_u;  
ntheta=Deg_u+Deg_l;   
Length=2; BasisType='poly';
elseif DGMtype==2
Deg_u=80;Deg_l=Deg_u;  
ntheta=Deg_u+Deg_l;   
Length=45; BasisType='rbf';    
end
     
options=optimoptions('fmincon','MaxIterations',1e4,'ConstraintTolerance',1e-6,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',2e4,...
    'Display','off','Algorithm','sqp');

% build model
IPM=IPM_Model('degreeup', Deg_u,'degreelow', Deg_l, 'basis',BasisType, 'Length',Length,...
    'Xsup',Xsupport, 'Ysup',Ysupport,'options',options);

beta_target=10^-6;
X_linpace=linspace(min(X),max(X),IPM.IntegrationKnots); 
Noutliers=200;
EmpiricalLevel=20; 
PLOT=1; % set = 1 to show the plots
%% optimize IPM bounds using different scenario programs
for TesProg=1:6  
    if TesProg==1 % Hard-constrained no-exceptions
        Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);% Desing Method with no uncertainty
    elseif TesProg==2  % Hard-constrained with 10 discarded samples
        Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,Noutliers,ThetaOptimal{1});% Desing Method with no uncertainty
    elseif TesProg==3  % Hard-constrained with 100 discarded samples
        Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,Noutliers*2,ThetaOptimal{1});% Desing Method with no uncertainty
    elseif TesProg==4 % Min-max layer
        Design=IPM.DesingIPM_MinMaxLayer(X,Y,ThetaOptimal{1}(1:Deg_u));% Desing Method with no uncertainty
    elseif TesProg==5 % Soft-constrained with rho=0.02
        Design=IPM.DesingIPM_SoftConstraied(X,Y,0.02,ThetaOptimal{1});
    elseif TesProg==6 % Soft-constrained with rho=0.01
        Design=IPM.DesingIPM_SoftConstraied(X,Y,0.01,ThetaOptimal{1});
    end
 
    % Assess generalization and reliability bounds
    if TesProg==1 % Hard-constrained no exceptions
        SN=Design.Generalization.NSupportConstraints;
        epsilon=getWaitandJudgeEpsilon_fast(SN,N_train,beta_target);
        RelBound=[1-epsilon,1];
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
        I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
    elseif TesProg==2  % Hard-constrained with Noutliers discarded samples
        SN=Design.Generalization.NSupportConstraints+Design.Generalization.Nout;
        epsilon=getWaitandJudgeEpsilon_fast(SN,N_train,beta_target);
        RelBound=[1-epsilon,1];
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
        I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
    elseif TesProg==3  % Hard-constrained with Noutliers*10 discarded samples
        SN=Design.Generalization.NSupportConstraints+Design.Generalization.Nout;
        epsilon=getWaitandJudgeEpsilon_fast(SN,N_train,beta_target);
        RelBound=[1-epsilon,1];
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
        I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
    elseif TesProg==4 % Min-max layer
        [Ub_eps,Lb_eps]=get_epsilon_minmax(Design.EmpiricalCosts,ntheta,beta_target);
        RelBound=[Lb_eps(end-EmpiricalLevel+1),Ub_eps(end-EmpiricalLevel+1)]; 
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup);% predict method 
        I_predict_denorm=IPM.De_Normalize([I_predict+ repmat([Design.EmpiricalCosts(EmpiricalLevel);-Design.EmpiricalCosts(EmpiricalLevel)],1,size(I_predict,2))],min(Y),max(Y))';
    elseif TesProg==5 % Soft-constrained with rho=10^6
        SN=Design.Generalization.NSupportConstraints;
        [epsL, epsU] = epsLU_fast(SN,N_train,beta_target);
        epsilon=[epsU, epsL];
        RelBound=1-epsilon;
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
        I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
    elseif TesProg==6 % Soft-constrained with rho=1
        SN=Design.Generalization.NSupportConstraints;
        [epsL, epsU] = epsLU_fast(SN,N_train,beta_target);
        epsilon=[epsU, epsL];
        RelBound=1-epsilon;
        I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
        I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
    end
    
    % True Reliability 
    [Reliability,Y_target_up,Y_target_low] = IPM.Validation_Samples(X_test,Y_test,X_linpace,I_predict_denorm);
    
    % Save results
    Reli_estimate=1-Reliability.PF_all;
    Area_normalized_trapz=Design.Area;
%     Area_denorm_average=mean(I_predict_denorm(:,1)-I_predict_denorm(:,2));
%     Area_denorm_trapz=trapz(X_linpace,I_predict_denorm(:,1)-I_predict_denorm(:,2));
%     Areas=[Area_normalized_trapz,Area_denorm_average,Area_denorm_trapz];
    CPUtime=Design.ComputationalTime/60;  % in minutes
    ResultMatrix(TesProg,:)=[CPUtime,Area_normalized_trapz,Reli_estimate,RelBound]; 
    SupporNumber(TesProg)=SN; 
    ThetaOptimal{TesProg}=Design.OptTheta;
    IPM_objects{TesProg}=Design;
    
    %  plot 
    if PLOT 
       subplot(3,2,TesProg)
       scatter(X,Y,'o','MarkerFaceColor','b','MarkerFaceAlpha',0.3); box on; grid on; hold on; 
        plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
        idx_Failed=Reliability.IdxFailed;
        scatter(X_test(idx_Failed),Y_test(idx_Failed),'xr','MarkerEdgeAlpha',0.05')
        
    end
end