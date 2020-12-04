# Interval Predictor Model

Interval Predictive Model gives an interval-valued characterization of the uncertainty affecting the stochastic process
This repository presents a matlab class to optimized the bounding functions defining an IPM 


    % an IPM is a rule I(x;theta) which assign to an input x an interval
    % for a dependent quantity y.
    % two bounding functions deinfe I(x;theta)=[fl(x;theta);fu(x;theta)]
    % where \theta are the fitting coefficients to be optimized and,
    % fl(x;theta) is a linear combination of \theta and basis functions.
    
    % IPM are built given-data from a set of uncertain scenarios
    % D_{N}={(x1,y1),(x2,y2),..., (xN,yN)} assumed to be iid
    
    % Different methods are available in this class to optimize I(x;theta)
    % 1) full-data enclosure (no exeption) and minimization of the area
    %    between the bounding functions
    %    min_{\theta} Area(theta) s.t.
    %    f_u(xi;theta)>=y_i and  f_l(xi;theta)<=y_i for all i=1,..,D_{N}
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
    % 2) discard No samples from the data base D_{N} and minimization of the
    %    area for the remaining scenarios
    %    min_{\theta} Area(theta) s.t.
    %    f_u(xi;theta)>=y_i and  f_l(xi;theta)<=y_i for all i in D_{N-No}
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
    % 3) Minimax layer IPM: minimizes the maximum distance between a
    %    regression function and the samples y_i
    %     min_{\theta} max_{i=1,..,N} |f_reg(x_i;theta)-y_i|
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
    % 4) Soft-constrained IPM: minizes a combination of area metric and
    % cost of violations given by a parameter rho>0
    %    min_{\theta} Area(theta)+rho \sum_{i=1}^N \zeta_i s.t.
    %    f_u(xi;theta)>=y_i-\zeta_i and  f_l(xi;theta)<=y_i-\zeta_i for all x_i,y_i \in D_{N}
    %     \zeta_i>=0 for i=1,...,N  (slack variables)
    %    f_u(x;theta)>=f_l(x;theta) for all x (up bound dominance)
    
    
    % 5) CVaR method (To be updated)
    % min_{\theta} Area(theta)
    % s.t.    f_u(xi;theta)>=CVAR(y_i,alpha) for all i=1,..,D_{N}
    %         f_l(xi;theta)<=CVAR(y_i,1-alpha) for all i=1,..,D_{N}
    % INPUTS:
    % Xdn vector of explanatory variables (1-Dimensional [1xNsamples])
    % Ydn vector of dependent variables (1-Dimensional [1xNsamples])
    % OUTPUTS:Design: the output structure
    
    % Design.Area: The area betweeen the optimized (accuracy)
    % Design.OptTheta: optimized fitting coefficients defining the bounds
    % Design.Generalization: Compelxity of the solution and other
    % propreties needed to evaluate scenario-based reliabiity/error bounds
