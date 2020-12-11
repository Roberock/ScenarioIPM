#  Classdef Interval Predictor Model

Interval Predictive Model (IPM) gives an interval-valued characterization of the uncertainty affecting the stochastic process
This repository presents a matlab class to optimized the bounding functions defining an IPM.
The reliability of the optmized predictor (probability that future samples will fall outside from the predictive bounds) is formally bounded thanks to scenario theory

<p align="center">
  <img src="./figs/IPM_Example.png" alt="Size Limit CLI" width="650">
</p>
 
 
```
   an IPM is a rule I(x;theta) which assign to an input x an interval for a dependent quantity y.
    two bounding functions deinfe
    
    I(x;theta)=[fl(x;theta);fu(x;theta)]
    
    where \theta are the fitting coefficients to be optimized and,
    fl(x;theta) is a linear combination of \theta and basis functions. 
```    
```
   IPM are built given-data from a data set of sampels of the process
   D_{N}={(x1,y1),(x2,y2),..., (xN,yN)} assumed to be iid at each i
```
##  Different methods are available in this class to optimize I(x;theta)
### 1) full-data enclosure (no exeption) and minimization of the area between the bounding functions

    %    min_{\theta} Area(theta) s.t.
    %    f_u(xi;theta)>=y_i and  f_l(xi;theta)<=y_i for all i=1,..,D_{N}
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
### 2) discard No samples from the data base D_{N} and minimization of the area for the remaining scenarios
    
    %    min_{\theta} Area(theta) s.t.
    %    f_u(xi;theta)>=y_i and  f_l(xi;theta)<=y_i for all i in D_{N-No}
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
### 3) Minimax layer IPM: minimizes the maximum distance between a regression function and the samples y_i

    min_{\theta} max_{i=1,..,N} |f_reg(x_i;theta)-y_i|
     f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
### 4) Soft-constrained IPM: minizes a combination of area metric and cost of violations given by a parameter rho>0
     
     min_{\theta} Area(theta)+rho \sum_{i=1}^N \zeta_i s.t.
     f_u(xi;theta)>=y_i-\zeta_i and  f_l(xi;theta)<=y_i-\zeta_i for all x_i,y_i \in D_{N}
     \zeta_i>=0 for i=1,...,N  (slack variables)
     f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
    
### 5) CVaR method (To be updated)

     min_{\theta} Area(theta)
     s.t.    f_u(xi;theta)>=CVAR(y_i,alpha) for all i=1,..,D_{N}
            f_l(xi;theta)<=CVAR(y_i,1-alpha) for all i=1,..,D_{N}
    
   
 <p align="center">
  <img src="./figs/IPM_Example_hard_vs_soft_Constraints.png" alt="Size Limit CLI" width="650">
</p>
   
   
### A simple Example :
  ```
  %inputs 
  
    Xdn vector of explanatory variables (1-Dimensional [1xNsamples])
    Ydn vector of dependent variables (1-Dimensional [1xNsamples])
   % OUTPUTS: a structure named Design containing the following fields
    
     Design.Area: The area betweeen the optimized (accuracy)
     Design.OptTheta: optimized fitting coefficients defining the bounds
     Design.Generalization: Compelxity of the solution and other
     propreties needed to evaluate scenario-based reliabiity/error bounds
```

### Example Scenaro-based Reliability Bounds:

% Generate scenarios 
Nsamples=400;
X=unifrnd(-3,10,[1,Nsamples]);
Y= (X.*rand([1,Nsamples])).^2+X.*exprnd(2,[1,Nsamples]);


% prepare IPM object
Xsupport=[min(X),max(X)];
Ysupport=[min(Y),max(Y)];
options=optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-4,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e4,...
    'Display','iter','Algorithm','sqp');

IPM=IPM_Model('degreeup', 6,'degreelow',6 ,...
    'basis','poly', 'Length',[],...
    'Xsup',Xsupport, 'Ysup',Ysupport,'options',options);
    
    
    %% hard-constrained solver
Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);

X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
%plot
subplot(2,2,1)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on; grid on;
scatter(X,Y,'.b')
title('hard-constrained')
%% hard-constrained with discarded samples
Ndiscarded=round(Nsamples/10); % remove 10% of the saples
Design=IPM.DesingIPM_hard_constrained_discarded(X,Y,Ndiscarded);

X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)'; % optimize the model
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));
%plot
subplot(2,2,2)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on;  grid on;
scatter(X,Y,'.b')
title('hard-constrained discarded')
%% minimax layer
Design=IPM.DesingIPM_MinMaxLayer(X,Y);
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)'; % optimize the model
I_predict1=[I_predict(:,1)-Design.EmpiricalCosts(1) I_predict(:,2)+Design.EmpiricalCosts(1)]; % maximum empirical layer
I_predict_denorm1= IPM.De_Normalize(I_predict1,min(Y),max(Y));
I_predict50=[I_predict(:,1)-Design.EmpiricalCosts(50) I_predict(:,2)+Design.EmpiricalCosts(50)]; % maximum empirical layer
I_predict_denorm50= IPM.De_Normalize(I_predict50,min(Y),max(Y)); 

%plot
subplot(2,2,3)
plot(X_linpace,I_predict_denorm1,'k','LineWidth',1.5);
hold on;  grid on;
plot(X_linpace,I_predict_denorm50,'--k','LineWidth',1.5);
scatter(X,Y,'.b')
title('hard-constrained minmax')
%% soft constrained layer
CostofViolations=0.1;
Design=IPM.DesingIPM_SoftConstraied(X,Y,CostofViolations); % optimize the model
 
X_linpace=linspace(Xsupport(1),Xsupport(2),10^4);
I_predict=IPM.Predict(X_linpace,Design.OptThetalow,Design.OptThetaup)';% predict method
I_predict_denorm= IPM.De_Normalize(I_predict,min(Y),max(Y));

%plot
subplot(2,2,4)
plot(X_linpace,I_predict_denorm,'k','LineWidth',1.5);
hold on; grid on;
scatter(X,Y,'.b')
title('soft-constrained')


