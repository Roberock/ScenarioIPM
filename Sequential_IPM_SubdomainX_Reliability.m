%% Try to Build an IPM with different level of reliability on the X axis

%
 clear variables; close all;
PlotL=1;
%% Generate samples for training
 
DGMtype = 2 ;% the data geneating mechanism
N_train=5000; % available scenarios
N_test=2e5; 
[X_all,Y_all,X_test,Y_test,Xsupport,Ysupport]= DGMs(DGMtype,N_train,N_test);

%scatter(X_test,Y_test,'xr','MarkerEdgeAlpha',0.05);
hold on; grid on; box on;
scatter(X_all,Y_all,'.b');
ylabel('Y'); xlabel('X')

 
DEG=5;
Nintervals=150;
EPS=1e-6;
Xdomain=[linspace(min(X_all)-EPS,max(X_all)+EPS,Nintervals+1)];

for i=1:Nintervals
     % data for training
    Somainidx=(X_all>=Xdomain(i)  & X_all<Xdomain(i+1));
    N=sum(Somainidx);
    X=X_all(Somainidx);
    Y=Y_all(Somainidx);
    % data for testing
    Somainidx=(X_test>=Xdomain(i)  & X_test<Xdomain(i+1));   N_test_domain=sum(Somainidx);
    X_test_domain=X_test(Somainidx);
    Y_test_domain=Y_test(Somainidx); 
    
    %% define IPM and train parameters
    % Scenario_Pf_constraint compute the constraints defined as w(d,delta_i)>\zeta_i
    SelectAlgorithm='sqp'; % 'interior-point' ,'sqp-legacy',  'trust-region-reflective' 'sqp' 'active-set'
    options = optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-6,...
        'StepTolerance',1e-6,'MaxFunctionEvaluations',1e3,'Display','off','Algorithm',SelectAlgorithm);
    Deg_u=DEG; % selct degree of the upper bound
    Deg_l=Deg_u;% selct degree of the lower bound
    BasisType='rbf'; 
    IPM=IPM_Model('degreeup', Deg_u,'degreelow',Deg_l ,...
        'basis',BasisType, 'Length',1,'Xsup',Xsupport, 'Ysup',Ysupport, 'options',options);
  
    Design_IPM=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);% Desing Method with no uncertainty
    %  Design_IPM=IPM.DesingIPM_Outliers_noMeasureUnc(X,Y,1);% Desing Method with soft constraints

    %% post process 
    Time(i)=Design_IPM.ComputationalTime;
    X_opt=Design_IPM.OptTheta;
    c=Design_IPM.Generalization.ConstraintsVal;
   % [c_test,~]=IPM_Constraint_No_measurmentuncertainty(X_opt,Y_test_domain,X_test_domain,Deg_u,Deg_l,BasisType);
%     Support_set_low=find(c(1:N)>-1e-3);
%     Support_set_up=find(c(1+N:end)>-1e-3);
%     Support_Set_p1=union(Support_set_low,Support_set_up);
   SupportCardinality(i)=Design_IPM.Generalization.NSupportConstraints;
 
    % compute scenario reliability-confidence
    epsilon=0.001:0.001:0.99;
    k=0; %number of discarded samples
    Nd=Deg_u+Deg_l; % number of desing variables
    beta_target=1e-6;
    epsilon_domain(i)=getepsilon_ConvexDiscard(N,beta_target,k,Nd); 
    OUT= getWaitandJudgeEpsilon_fast(Design_IPM.Generalization.NSupportConstraints,N,beta_target);
    epsilon_domain_posteriori(i) =OUT(end);
  
    I_predict_test=IPM.Predict(X_test_domain,Design_IPM.OptThetalow,Design_IPM.OptThetaup);% predict method
    I_predict_test= IPM.De_Normalize(I_predict_test,min(Y),max(Y))';
    idx_Failed= Y_test_domain>I_predict_test(:,1)' |  Y_test_domain<I_predict_test(:,2)';
    Pf(i)=mean(idx_Failed);
    if PlotL==1
        %plot
        figure(2)
        x_predict=min(X):0.01:max(X); 
        I_predict=IPM.Predict(x_predict,Design_IPM.OptThetalow,Design_IPM.OptThetaup);% predict method
        scatter(X,Y,'b.'); 
        box on; grid on; hold on;
        scatter(X_test_domain(idx_Failed),Y_test_domain(idx_Failed),'r*'); box on; grid on; hold on;
        plot(x_predict,IPM.De_Normalize(I_predict,min(Y),max(Y))','k','LineWidth',1.5); 
    end
end
%plot(MFRFX(Data_SuspensionArm.Cracktype==crackIdx,50:5:end-50)')
OUT=getWaitandJudgeEpsilon_fast(sum(SupportCardinality),N_train,beta_target);
epsilon_domain_posteriori_all=OUT(end); 
epsilon_domain_priori_all=getepsilon_ConvexDiscard(N_train,beta_target,0,Nd*Nintervals); 
TrueViolationProbability=mean(Pf);
display(['is the a posteriori bound  ' num2str(epsilon_domain_posteriori_all) '>=' num2str(TrueViolationProbability) ' ?'])
display([' is the a priori bound ' num2str(epsilon_domain_priori_all) '>=' num2str(TrueViolationProbability) ' ?'])