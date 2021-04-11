function [X,Y,X_test,Y_test,Xsupport,Ysupport]= DGMs(DGMtype,N_train,N_test)
   %rng('default'); rng(1); % for reproducibility
if DGMtype==1  % trigonometic function with random noise
    Xa=10; Xb=15; % X support to be examined
    Xsupport=[-Xa,Xb-Xa];
    my_fun=@(X)(X.^2).*cos(X)-sin(3*X).*exp(-X.^2)-X-cos(X.^2)+X.*(3*randn(1,length(X))); 
    X=[-Xa+Xb*rand(1,N_train)]; % generate x values uniformly in [Xa,Xb]
    Y=my_fun(X);  % random dependent variable y=f(x)
    % for testing
    X_test=[-Xa+Xb*rand(1,N_test)]; % generate x values uniformly in [Xa,Xb]
    Y_test=my_fun(X_test);  % random dependent variable y=f(x) 
    Ysupport=[min(Y_test),max(Y_test)];
elseif DGMtype==2 % crkced suspension arm data generating mechanism
    [X,Y]=SuspensionArm_DGM(N_train,1);  % for training
    [X_test,Y_test]=SuspensionArm_DGM(N_test,1);  % for testing
    Xsupport=[1,180];
    Ysupport=[min(Y_test),max(Y_test)];
elseif DGMtype==3 % crkced suspension arm data generating mechanism
    Xa=20; Xb=45; % X support to be examined
    Xsupport=[-Xa,Xb-Xa];
    X=[-Xa+Xb*rand(1,N_train)]; %-(1+3*rand(1,40))]; %generate x values
    my_DGM=@(X)(X.^2).*cos(X)-sin(3*X).*exp(-X.^2)-X-cos(X.^2)+X+lognrnd(0,1.3,[1,length(X)]); %generate y values
    Y=my_DGM(X);
    Noise=  wblrnd(1,2,[1,N_train]); % random measurement noise
    Y=[Y.*Noise]; 
    % testing data
    X_test =[-Xa+Xb*rand(1,N_test)];
    Y_test =my_DGM(X_test);
    Noise= wblrnd(1,2,[1,N_test]); % random measurement noise
    Y_test=[Y_test.*Noise]; 
    Ysupport=[min(Y_test),max(Y_test)];
end

end