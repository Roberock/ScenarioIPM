function   [Bound,Basis,Basisx]=Compute_Basis_and_Bound(x,Theta,Degree,BasesType,Leng,varargin)
% INPUTS
% x: samples of the dependent variable points {x_i| i=1,...,N}
% Degree: degree of the basis function (number of fitting coefficents theta)
% Theta: fitting coefficents of the function fu(x;theta) or fl(x;theta)
% BasesType: type of basis function {'poly','rbf','fourier',...}
% Leng: critical length for the rbf basis

% OUTPUTS
% Bound: the value fu(x;theta) of the bound
% Basisx: value of the basis in each x_i
% summedX: % Expected value of basis
%% default values
if nargin<=3; BasesType='poly';Leng=2; end
Basisx=zeros(Degree,length(x));  %Preallocate Basis for speed
%% Compute the basis
for i=1:Degree 
    
    if strcmp(BasesType,'poly') %(polynomial basis- default option)
        if i==1;  Basisx(i,:)=  ones(size(Basisx(i,:)));
        else; Basisx(i,:)=  x.^(i-1); % Compute polynomial Basis
        end
        
    elseif  strcmp(BasesType,'fourier')%(fourier basis)
        if i==1; Basisx(i,:)=  ones(size(Basisx(i,:)));
        else
            if mod(i,2) % if i is odd
                Basisx(i,:)= cos(pi.*tanh(x).*(i-1)); % Compute Fourier Basis
            else % if i is even
                Basisx(i,:)=  sin(pi.*tanh(x).*(i-1)); % Compute Fourier Basis
            end
        end
        
    elseif  strcmp(BasesType,'rbf') %(radial basis)
        % C(i) radius for the basis
        if i==1; Basisx(i,:)=  ones(size(Basisx(i,:)));
        else
            Centers=linspace(min(x),max(x),Degree); % konts for the basis
            Basisx(i,:)=  exp(-Leng*(x-Centers(i)).^2/(Centers(2)-Centers(1))); % Compute Gausiann Basis
        end
        
        %     elseif  strcmp(BasesType,'SE') % Square Expnential basis
        %         D=size(x,1);sigmaL10 = 0.1*ones(D,1);sigmaL20 = 0.1; sigmaF10 = 1; sigmaF20 = 1;
        %         theta0   = [log(sigmaL10);log(sigmaL20);log(sigmaF10);log(sigmaF20)];
        %         Cov=mykernel(x',x',theta0);
        %         Basisx(i,:)= C(i).*sum(Cov);   
        %     elseif  strcmp(BasesType,'SE') % Spline basis
     
    end
    
    
end
Basis=sum((Basisx),2);
if isempty(Theta)
    Bound=[];
else
    Bound= sum((Basisx).*repmat(Theta',1,size(Basisx,2)),2);
end
end
