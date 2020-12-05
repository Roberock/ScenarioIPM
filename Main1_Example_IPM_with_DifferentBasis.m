
clc
close all
clear variables
%%  Data Gen. Mechanism
Xa=40; Xb=25; % X support to be examined
Xsupport=[-Xa,Xb-Xa];
my_DGM_fun=@(X)(X.^2).*cos(X)-sin(3*X).*exp(-X.^2)-X-cos(X.^2)+X.*(3*randn(1,length(X)));
% training data set
N_train=100;
X=[-Xa+Xb*rand(1,N_train)]; % generate x values uniformly in [Xa,Xb]
Y=my_DGM_fun(X);  % random dependent variable y=f(x)
% testing data set
N_test=10000;
X_test=[-Xa+Xb*rand(1,N_test)]; % generate x values uniformly in [Xa,Xb]
Y_test=my_DGM_fun(X_test);  % random dependent variable y=f(x)
Ysupport=[min(Y_test),max(Y_test)];
% display data
figure(1)
scatter(X_test,Y_test,'dr','filled','MarkerFaceAlpha',0.1,'MarkerEdgeAlpha',0.1)
hold on;grid on;xlabel('X'); ylabel('Y')
scatter(X,Y,'b','filled') 
legend('Samples for Testing','Samples for Training')

%% Plot 3 example of sum of the basis functios terms in x
Degrees=[2 3 5];
Lengths=[1 5 20]; % length parameters for the radail basis functions
IntegrationKnots=5000; % number of points in x to compute the basis
X_support=linspace(-1,1,IntegrationKnots); % points where the basis is evaluated
BasisType={'rbf','rbf','rbf'};
Len =1;
Theta=ones(1,Degrees(1)); % the fitting coefficients for the function f(x_j)=sum_i theta_i*Basis_i(x_j)
f1= compute_Basis(X_support,Theta,Degrees(1),BasisType{1},Lengths(1));% upper boud  
Theta=ones(1,Degrees(2)); % the fitting coefficients for the function f(x_j)=sum_i theta_i*Basis_i(x_j)
f2= compute_Basis(X_support,Theta,Degrees(2),BasisType{2},Lengths(2));% upper boud
Theta=ones(1,Degrees(3)); % the fitting coefficients for the function f(x_j)=sum_i theta_i*Basis_i(x_j)
f3= compute_Basis(X_support,Theta,Degrees(3),BasisType{3},Lengths(3));% upper boud
 
 figure(2)
 plot(X_support,f1); hold on;
 plot(X_support,f2)
 plot(X_support,f3)
 xlabel('x'); ylabel('BasisFunction')
 legend(['Basisi Type' BasisType{1} ' Degree=' num2str(Degrees(1)) ' Length parameter=' num2str(Lengths(1))],...
     ['Basisi Type' BasisType{2} ' Degree=' num2str(Degrees(2)) ' Length parameter=' num2str(Lengths(2))],...
     ['Basisi Type' BasisType{3} ' Degree=' num2str(Degrees(3)) ' Length parameter=' num2str(Lengths(3))])
 %% Example with 6 basis type and differen rbf critical length
for i=1:6
     %% Define IPM model and options 
    display(['Evaluate Kernel Type number ' num2str(i)])
    if i==1; BASISTYPE='poly'; Dup=5; Dlow=Dup; Color='k'; Length=[];
    elseif i==2; BASISTYPE='poly'; Dup=10; Dlow=Dup; Color='k'; Length=[];
    elseif i==3; BASISTYPE='rbf'; Dup=5;Dlow=Dup; Color='r'; Length=2;
    elseif i==4; BASISTYPE='rbf'; Dup=10; Dlow=Dup;Color='r'; Length=2;
    elseif i==5; BASISTYPE='rbf'; Dup=10;Dlow=Dup; Color='r'; Length=10;
    elseif i==6; BASISTYPE='fourier'; Dup=10; Color='b'; Length=[];
    end 

    options=optimoptions('fmincon','MaxIterations',1e3,'ConstraintTolerance',1e-6,...
    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e3,...
    'Display','off','Algorithm','sqp'); 

    IPM=IPM_Model('degreeup', Dup,'degreelow',Dlow ,...
        'basis',BASISTYPE, 'Length',Length,...
     'Xsup',Xsupport, 'Ysup',Ysupport,'options',options); 
  
    %% Train IPM:
    % Design method: Full data enclosure with no exceptions  
    Design=IPM.DesingIPM_hard_constrained_noExeptions(X,Y);% Desing Method with no uncertainty 
    %% display results
    figure(3)
    subplot(2,3,i)
    X_support=linspace(min(X),max(X),IPM.IntegrationKnots);
    %X_pred_nor=Normalize(X_predict);
    I_predict=IPM.Predict(X_support,Design.OptThetalow,Design.OptThetaup);% predict method
    scatter(X,Y,'.');
    box on; grid on; hold on;
    plot(X_support,IPM.De_Normalize(I_predict,min(Y),max(Y))',Color,'LineWidth',1.5);
    title([BASISTYPE ' degree ' num2str(Dup) ' len=' num2str(Length)])
end
