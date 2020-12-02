classdef IPM_Model
    %IPM_Model: Construct and traing an Interval Predictor Model (IPM)
    % INPUTS
    % Xdn vector of explanatory variables (1-Dimensional [1xNsamples])
    % Ydn vector of dependent variables (1-Dimensional [1xNsamples])
    
    properties
        Deg_u % degree of the upper bound (e.g. edgree of the polynomial basis)
        Deg_l % degree of the lower bound (e.g. edgree of the polynomial basis)
        BasisType % the type of basis function 'poly', 'fourier', 'rbf', ....(consider adding more, e.g, 'B-spline', 'Wavelet')
        Len % this is a length parameter defining the spread of 'radial basisi functions'
        options % the options to solve the MATLAB fmincon Program
        Xsupport % the support of the input variables
        Ysupport % the support of the dependent variables
    end
    
    properties (Constant=true, Hidden=false)
        IntegrationKnots=5000  % Number of integration knots in X we used to approx integrals (needed for the area metric)
    end
    
    methods
        %% CONSTRUCTOR
        function IPMobject = IPM_Model(varargin)  % method to define an object of type Scenario_Based_Reliability
            if nargin==0;  return; end  % Create an empty object
            for k=1:2:length(varargin)
                switch lower(varargin{k})
                    case {'deg_u','degreeup','dup'}; IPMobject.Deg_u=varargin{k+1};
                    case {'deg_l','degreelow','dlow'}; IPMobject.Deg_l=varargin{k+1};
                    case {'basistype','basis'};  IPMobject.BasisType=varargin{k+1};
                    case {'len','length'};  IPMobject.Len=varargin{k+1};
                    case {'xsup','xsupport'};  IPMobject.Xsupport=varargin{k+1};
                    case {'ysup','ysupport'};  IPMobject.Ysupport=varargin{k+1};
                    case {'options','optimoptions'};  IPMobject.options=varargin{k+1};
                end
            end
            if isempty(IPMobject.Deg_u); IPMobject.Deg_u=3; end % if is empty, set default degree = 3
            if isempty(IPMobject.Deg_l);IPMobject.Deg_l=3;end % if is empty, set default degree = 3
            if isempty(IPMobject.BasisType);IPMobject.BasisType='poly';end % if is empty, set to the default 'poly'
            if isempty(IPMobject.Len);IPMobject.Len=2; end% if is empty, set to the default value of 2
            if isempty(IPMobject.options) % default optimization option if not provided
                IPMobject.options=optimoptions('fmincon','MaxIterations',1e4,'ConstraintTolerance',1e-8,...
                    'StepTolerance',1e-6,'MaxFunctionEvaluations',1e4,'Display','off','Algorithm','sqp');  end
        end     %of constructor
        
        %% IPM Training Methods
        function [IPM] = DesingIPM_hard_constrained_noExeptions(obj,Xdn,Ydn)
            
            X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2));  %normalize in [-1 +1]
            Y=obj.Normalize(Ydn,min(Ydn),max(Ydn));    %normalize in [-1 +1]
            
            % AREA Function: approximate integral of the difference between two curves
            [~,~,Area_fun]=Def_fun_basis_objective(obj,X);
            
            % Define constraints and initial guess
            X0=[0.001*ones(1,obj.Deg_u),-0.001*ones(1,obj.Deg_l)]; % first guess
            LB=[-inf*ones(1,obj.Deg_u),-inf*ones(1,obj.Deg_l)];
            UB=[+inf*ones(1,obj.Deg_u),+inf*ones(1,obj.Deg_l)];
            ConstrainFun= @(x) IPM_constraints_no_exceptions(x,Y,X,obj.Deg_u,obj.Deg_l,obj.BasisType,obj.Len);
            
            % Tune IPM parameters (TO DO: replace with linprog to speed up calculations)
            tic
            [Theta_opt,Area,exitflag,output,lambda]= fmincon(@(x) Area_fun(x),... % objective function
                X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                @(x) ConstrainFun(x),obj.options); % non linear constraints and options
            
            G_low=lambda.ineqnonlin(1:length(Y)); % lagrangian multipliers for the lower bound constraints
            G_up=lambda.ineqnonlin(length(Y)+1:length(Y)*2); % lagrangian multipliers for the upper bound constraints
            
            % Save results
            IPM.ComputationalTime=toc; % computatioanl time to desing the IPM
            IPM.Area=Area; % estimator of the area in between the bounds
            IPM.OptTheta= Theta_opt; % optimzed design varibales
            IPM.OptThetaup=Theta_opt(1:obj.Deg_u); % upper bound  optimal degree
            IPM.OptThetalow=Theta_opt(obj.Deg_u+1:obj.Deg_u+obj.Deg_l); % lower bound optimal degree
            IPM.fminconOutput_Structure.output=output;% structure output of the fmincon
            IPM.fminconOutput_Structure.exitflag=exitflag; % exit flat
            
            IPM.Data.Y=Ydn;  IPM.Data.X=Xdn;
            IPM.Data.Xsupport=[min(X),max(X)]; % support of the explanatory variable
            
            % generealization info
            IPM.Generalization.IdxSupportConstraints_upperbound= find(G_up>0);
            IPM.Generalization.IdxSupportConstraints_lowerbound= find(G_low>0);
            IPM.Generalization.IdxSupportall_union=union( find(G_up>0),  find(G_low>0));
            IPM.Generalization.NSupportConstraints=sum(lambda.ineqnonlin>0);
            IPM.Generalization.ConstraintsVal= ConstrainFun(Theta_opt);
            IPM.Generalization.Nout=[];
        end
        %% DesingIPM by removing Nout outlier samples
        function [IPM] = DesingIPM_hard_constrained_discarded(obj,Xdn,Ydn,Nout,Baseline_design)
            % Nout  number of samples discarded from the set of scenario constraints
            X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2));
            Y=obj.Normalize(Ydn,min(Ydn),max(Ydn));
            
            % AREA Function: approximate integral of the difference between two curves
            [~,~,Area_fun]=Def_fun_basis_objective(obj,X);
            
            if nargin>4 % if we provide a baseline solution for the IPM
                X0=[Baseline_design];
            else
                X0=[0.001*ones(1,obj.Deg_u),-0.001*ones(1,obj.Deg_l)];
            end
            
            LB=[-inf*ones(1,obj.Deg_u),-inf*ones(1,obj.Deg_l)];
            UB=[+inf*ones(1,obj.Deg_u),+inf*ones(1,obj.Deg_l)];
            ConstrainFun= @(x) IPM_Constraint_RemoveOutliers_fromTailG(x,Y,X,Nout,obj.Deg_u,obj.Deg_l,obj.BasisType,obj.Len);
            % Tune IPM parameters (TO DO: replace with linprog to speed up calculations)
            tic
            [Theta_opt,Area,exitflag,output,lambda]= fmincon(@(x) Area_fun(x),... % objective function
                X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                @(x) ConstrainFun(x),obj.options);
            
            IPM.ComputationalTime=toc; % computatioanl time to desing the IPM
            IPM.Area=Area; % estimator of the area in between the bounds
            IPM.OptTheta= Theta_opt; % optimzed design varibales
            IPM.OptThetaup=Theta_opt(1:obj.Deg_u); % upper bound  optimal degree
            IPM.OptThetalow=Theta_opt(obj.Deg_u+1:obj.Deg_u+obj.Deg_l); % lower bound optimal degreedegree
            IPM.Outlier= find(round(Theta_opt(obj.Deg_u+obj.Deg_l+1:end))==0); % the outliers
            IPM.Inlier= find(round(Theta_opt(obj.Deg_u+obj.Deg_l+1:end))==1); % the inliers
            IPM.fminconOutput_Structure.output=output; % structure output of the fmincon
            IPM.fminconOutput_Structure.exitflag=exitflag;
            
            % save data info
            IPM.Data.Y=Ydn;  IPM.Data.X=Xdn;
            IPM.Data.Xsupport=[min(X),max(X)]; % support of the explanatory variable
            
            % Structure for the analysis of Generealization bounds
            c_low=compute_Basis(X,Theta_opt(obj.Deg_u+1:obj.Deg_u+obj.Deg_l),obj.Deg_l,obj.BasisType,obj.Len)-Y;
            c_up=Y-compute_Basis(X,Theta_opt(1:obj.Deg_u),obj.Deg_u,obj.BasisType,obj.Len);
            IPM.Generalization.IdxSupportConstraints_upperbound= find(c_up>=0);
            IPM.Generalization.IdxSupportConstraints_lowerbound= find(c_low>=0);
            IPM.Generalization.IdxSupportall_union=find(lambda.ineqnonlin>0);
            IPM.Generalization.NSupportConstraints=sum(c_up>=0 | c_low>=0);
            %IPM.Generalization.NSupportConstraints=sum(lambda.ineqnonlin>0);
            IPM.Generalization.Nout=Nout;
            IPM.Generalization.ConstraintsVal= ConstrainFun(Theta_opt);
        end
        
        %% DesingIPM: minmax layer (L-infinity norm minimization)
        function [IPM] = DesingIPM_MinMaxLayer(obj,Xdn,Ydn,Baseline_design)
            % Xdn vector of explanatory variables (1-Dimensional [1xNsamples])
            % Ydn vector of dependent variables (1-Dimensional [1xNsamples])
            X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2)); % normalize in [-1 1]
            Y=obj.Normalize(Ydn,min(Ydn),max(Ydn)); % normalize in [-1 1]
            
            if nargin>3 % if we provide a baseline solution for the IPM
                X0=[Baseline_design];
            else;  X0=ones(1,obj.Deg_u);% initial guess
            end
            LB= -inf*ones(1,obj.Deg_u); % lower bounds
            UB= +inf*ones(1,obj.Deg_u); % upper bounds
            w=@(x)  abs(Y-compute_Basis(X,x(1:obj.Deg_u),obj.Deg_u,obj.BasisType,obj.Len)) ; % Errors
            
            % Tune IPM parameters (TO DO: replace with linprog to speed up calculations)
            tic
            [Theta_opt,Minmax_val,exitflag,output]= fmincon(@(x) max(w(x)),... % objective function
                X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                [],obj.options); % non linear constraints and options
            IPM.ComputationalTime=toc; % computatioanl time to desing the
            
            % Save results
            X_linspace=linspace(min(X),max(X),obj.IntegrationKnots);
            f_regression=@(x) compute_Basis(X_linspace,x(1:obj.Deg_u),obj.Deg_u,obj.BasisType,obj.Len);% upper boud
            Area_fun= @(x)  trapz(X_linspace, (f_regression(x)+Minmax_val-(f_regression(x)-Minmax_val)));
            %[~,~,Area_fun]=Def_fun_basis_objective(obj,X);
            IPM.Area=Area_fun(Theta_opt); % estimator of the area in between the bounds
            IPM.Minmax_val=Minmax_val;
            IPM.OptTheta= Theta_opt; % optimzed design varibales
            IPM.OptThetaup=Theta_opt;
            IPM.OptThetalow=Theta_opt;
            IPM.fminconOutput_Structure.output=output;% structure output of the fmincon
            IPM.fminconOutput_Structure.exitflag=exitflag;
            % save data info
            IPM.Data.Y=Ydn; IPM.Data.X=Xdn;
            IPM.Data.Xsupport=[min(X),max(X)]; % support of the explanatory variable
            % Structure for the analysis of Generealization bounds
            IPM.EmpiricalCosts=sort(w(Theta_opt),'descend');
            IPM.Generalization.NSupportConstraints=[];
            IPM.Generalization.Nout=[];
            IPM.Generalization.ConstraintsVal= [];
        end 
        %% DesingIPM: soft-constrained optimization (slack variables relaxin sample constraints)
        function [IPM] = DesingIPM_SoftConstraied(obj,Xdn,Ydn,Rho,Baseline_design)
            % Rho: Regularization parameter
            % X vector of explanatory variables (1-Dimensional [1xNsamples])
            % Y vector of dependent variables (1-Dimensional [1xNsamples])
            %             X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2));
            %             Y=obj.Normalize(Ydn,obj.Ysupport(1),obj.Ysupport(2));
            X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2)); % normalize in [-1 1]
            Y=obj.Normalize(Ydn,min(Ydn),max(Ydn)); % normalize in [-1 1]
            
            N=length(Y);
            [~,~,Area_fun]=Def_fun_basis_objective(obj,X);
            Slack_fun= @(x) sum(x(obj.Deg_u+obj.Deg_l+1:end));
            if nargin>4 % if we provide a baseline solution for the IPM
                X0=[Baseline_design, zeros(1,N)];
            else
                X0=[0.001*ones(1,obj.Deg_u), -0.001*ones(1,obj.Deg_l), zeros(1,N)];
            end
            LB=[-inf*ones(1,obj.Deg_u), -inf*ones(1,obj.Deg_l), zeros(1,N)];
            UB=[+inf*ones(1,obj.Deg_u), +inf*ones(1,obj.Deg_l), +inf*ones(1,N)];
            
            ConstrainFun= @(x)  IPM_SoftConstraint(x,Y,X,obj.Deg_u,obj.Deg_l,obj.BasisType,obj.Len);
            objfun=@(x) Area_fun(x)+Rho*Slack_fun(x);
            % Tune IPM parameters (TO DO: replace with linprog to speed up calculations)
            tic
            [Theta_opt,~,exitflag,output,lambda]= fmincon(@(x) objfun(x),... % objective function
                X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                @(x) ConstrainFun(x),obj.options);
            
            IPM.ComputationalTime=toc; % computatioanl time to desing the IPM
            IPM.Area=Area_fun(Theta_opt); % estimator of the area in between the bounds
            IPM.OptTheta= Theta_opt; % optimzer output
            IPM.OptThetaup=Theta_opt(1:obj.Deg_u); % upper bound  optimal degree
            IPM.OptThetalow=Theta_opt(obj.Deg_u+1:obj.Deg_u+obj.Deg_l); % lower bound optimal degree
            IPM.Zeta= Theta_opt(obj.Deg_u+obj.Deg_l+1:end); % slack variable value
            IPM.fminconOutput_Structure.output=output; % structure output of the fmincon
            IPM.fminconOutput_Structure.exitflag=exitflag;
            % save data info
            IPM.Data.Y=Ydn;  IPM.Data.X=Xdn;
            IPM.Data.Xsupport=[min(X),max(X)]; % support of the explanatory variable
            % Generealization info
            IPM.Generalization.NSupportConstraints=max(sum(IPM.Zeta>0),sum(lambda.ineqnonlin>0));
            IPM.Generalization.NSupportConstraints=sum(lambda.ineqnonlin>0);%sum(g_val>=0);
            IPM.Generalization.Nout=[];
            IPM.Generalization.ConstraintsVal= ConstrainFun(Theta_opt);
        end 
        %% DesingIPM CVAR METHOD (To be updated)
        function [IPM] = DesingIPM_CVAR_VAR(obj,Xdn,Ydn,Alpa_level,Error_fun,VarorCvar)
            % Alpa_level alpha level of the error distribution we whant to include in the bounds
            % Error_fun a function of the error (assumed gaussian),
            % e.g. Error_fun= @(x,y) 2; %(constant error)%or Error_fun= @(x,y) abs(y.^2)*(0.01);
            X=obj.Normalize(Xdn,obj.Xsupport(1),obj.Xsupport(2));
            Y=obj.Normalize(Ydn,min(Ydn),max(Ydn));
            
            [~,~,Area_fun]=Def_fun_basis_objective(obj,X);
            X0=[0.001*ones(1,obj.Deg_u),-0.001*ones(1,obj.Deg_l)]; % initial guess
            LB=[-inf*ones(1,obj.Deg_u),-inf*ones(1,obj.Deg_l)];
            UB=[+inf*ones(1,obj.Deg_u),+inf*ones(1,obj.Deg_l)];
            % optimize IPM desing
            tic
            if strcmp(VarorCvar,'var')
                [Theta_opt,Area,exitflag,output]= fmincon(@(x) Area_fun(x),... % objective function
                    X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                    @(x) IPM_Value_at_risk_disjoint(x,Y,X,obj.Deg_u,obj.Deg_l,obj.BasisType,Alpa_level,Error_fun,obj.Len),obj.options);
            elseif strcmp(VarorCvar,'cvar')
                [Theta_opt,Area,exitflag,output]= fmincon(@(x) mean(f_up(x)-f_low(x)),... % objective function
                    X0,[],[],[],[], LB, UB,...% initial guess and design bounds
                    @(x) IPM_CVaRConstraint_y_unc(x,Y,X,obj.Deg_u,obj.Deg_l,obj.BasisType,Alpa_level,Error_fun,obj.Len),obj.options);
            end
            IPM.ComputationalTime=toc; % computatioanl time to desing the IPM
            IPM.Area=Area; % estimator of the area in between the bounds
            IPM.OptTheta= Theta_opt; % optimzed design varibales
            IPM.OptThetaup=Theta_opt(1:obj.Deg_u); % upper bound  optimal degree
            IPM.OptThetalow=Theta_opt(obj.Deg_u+1:obj.Deg_u+obj.Deg_l); % lower bound optimal degreedegree
            IPM.fminconOutput_Structure.output=output;% structure output of the fmincon
            IPM.fminconOutput_Structure.exitflag=exitflag;
            % save data info
            IPM.Data.Y=Ydn;  IPM.Data.X=Xdn;
            IPM.Data.Xsupport=[min(X),max(X)]; % support of the explanatory variable
        end
        %% Function definition
        function [f_up,f_low,objfun]=Def_fun_basis_objective(obj,X)
            X_linspace=linspace(min(X),max(X),obj.IntegrationKnots);
            f_up=@(x) compute_Basis(X_linspace,x(1:obj.Deg_u),obj.Deg_u,obj.BasisType,obj.Len);% upper boud
            f_low=@(x) compute_Basis(X_linspace,x(obj.Deg_u+1:obj.Deg_u+obj.Deg_l),obj.Deg_l,obj.BasisType,obj.Len);  % lower boud
            %objfun=@(x) mean(f_up(x)-f_low(x)); % area metric as the average of intervals
            objfun= @(x) trapz(X_linspace,f_up(x)-f_low(x)); % trapezoidal integration method
        end
        %% Map: {xdn_i,ydn_i} -->  {x_i,y_i} where x_i,y_i \in [-1,1]
        function Xn=Normalize(obj,X,MIN,MAX)
            Xn=obj.mapfun(X,MIN,MAX,-1,+1);
        end
        function Xdn=De_Normalize(obj,X,MIN,MAX)
            Xdn=obj.mapfun(X,-1,+1,MIN,MAX);
        end
        function output = mapfun(~,value,fromLow,fromHigh,toLow,toHigh)
            output = (value - fromLow) .* (toHigh - toLow) ./ (fromHigh - fromLow) + toLow;
        end
        
        %% IPM bounds methods
        function [Bounds] = Predict(obj,x_predict,OptThetalow,OptThetaup)
            % x_predict a vector or a point where to evaluate the bounds
            % OptDeglow,OptDegup degree of the basis to be evaluated
            %x_predict=min(X):0.01:max(X);
            x_predict=obj.Normalize(x_predict,obj.Xsupport(1),obj.Xsupport(2)); %% FIX
            F_up_predict=compute_Basis(x_predict,OptThetaup,obj.Deg_u,obj.BasisType,obj.Len);
            F_low_predict=compute_Basis(x_predict,OptThetalow,obj.Deg_l,obj.BasisType,obj.Len);
            Bounds=[F_up_predict;F_low_predict];
        end
        %% IPM bounds methods
        function [Reliability,Y_target_up,Y_target_low] = Validation_Samples(~,X_test,Y_test,X_bnd_dn,Y_bnd_dn)
            YUb=Y_bnd_dn(:,1);
            YLb=Y_bnd_dn(:,2);
            [Y_target_up,Y_target_low]=deal(zeros(length(X_test),1));
            for i=1:length(X_test)
                idx_1=find(X_bnd_dn<X_test(i),1,'last');
                idx_2=find(X_bnd_dn>X_test(1),1,'first');
                if isempty(idx_1)
                    idx_1=find(X_bnd_dn==min(X_bnd_dn));
                end
                if  isempty(idx_2)
                    idx_2=find(X_bnd_dn==max(X_bnd_dn));
                end
                X_1=  X_bnd_dn(idx_1); X_2=  X_bnd_dn(idx_2);
                Y_1u=  YUb(idx_1);   Y_1l=  YLb(idx_1);
                Y_2u=  YUb(idx_2);   Y_2l=  YLb(idx_2);
                DeltaYu=Y_1u-Y_2u;
                DeltaYl=Y_1l-Y_2l;
                DeltaX=X_1-X_2;
                Y_target_up(i)=Y_1u-(X_1-X_test(i))*(DeltaYu)/ (DeltaX);
                Y_target_low(i) =Y_1l-(X_1-X_test(i))*(DeltaYl)/ (DeltaX);
            end
            Reliability.IdxFailed_up= Y_test'>Y_target_up;
            Reliability.PF_up  = mean(Reliability.IdxFailed_up); % failure probability upper bound
            Reliability.IdxFailed_low= Y_test'<Y_target_low;
            Reliability.PF_low = mean(Reliability.IdxFailed_low);  % failure probability lower bound
            Reliability.IdxFailed=(Reliability.IdxFailed_up | Reliability.IdxFailed_low);
            Reliability.PF_all =mean(Reliability.IdxFailed); % failure probability interval
        end
    end
end

