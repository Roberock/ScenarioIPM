function   [SummedBasis]=compute_Basis(x,Theta,Degree,BasesType,Leng,varargin)
% INPUTS
% x dependent variable points (x_i points with i=1,...,N)
% NPoly= degree of the basis
% Theta=coefficents of the function
% OUTPUTS
% Basisx= value of the basis in each x_i
% summedX= % Expected value of basis
%%
%NPoly=4;
if nargin<=3; BasesType='poly';Leng=2; end
Basisx=zeros(Degree,length(x));  %Preallocate Basis for speed
% summedX=zeros(NTerms,1);
for i=1:Degree
    
    if strcmp(BasesType,'poly') %(default option)
        if i==1 
            Basisx(i,:)= Theta(i)*ones(size(Basisx(i,:)));
        else
            Basisx(i,:)= Theta(i)*x.^(i-1); % Compute polynomial Basis
        end
        
    elseif  strcmp(BasesType,'fourier')
        if i==1
            Basisx(i,:)= Theta(i)*ones(size(Basisx(i,:)));
        else
        if mod(i,2) % if i is odd
            Basisx(i,:)= Theta(i)*cos(pi.*tanh(x).*(i-1)); % Compute Fourier Basis
        else % if i is even
            Basisx(i,:)= Theta(i)*sin(pi.*tanh(x).*(i-1)); % Compute Fourier Basis
        end
        end 
        
    elseif  strcmp(BasesType,'rbf')
        % C(i) radius for the basis
        if i==1; Basisx(i,:)= Theta(i)*ones(size(Basisx(i,:)));
        else
            Centers=[linspace(min(x),max(x),(Degree+1))]; %konts for the basis
            Basisx(i,:)= Theta(i).*exp(-Leng*(x-Centers(i)).^2/(Centers(2)-Centers(1))); % Compute Gausiann Basis
           %  Basisx(i,:)= Theta(i).*exp(-Leng*(x-Centers(i)).^2/(Centers(2)-Centers(1))); % Compute Gausiann Basis
        end
        
        %     elseif  strcmp(BasesType,'SE')
        %         D=size(x,1);sigmaL10 = 0.1*ones(D,1);sigmaL20 = 0.1; sigmaF10 = 1; sigmaF20 = 1;
        %         theta0   = [log(sigmaL10);log(sigmaL20);log(sigmaF10);log(sigmaF20)];
        %         Cov=mykernel(x',x',theta0);
        %         Basisx(i,:)= C(i).*sum(Cov);  % Compute Squared Exp?
        %
    end
    
    
end
SummedBasis=sum(Basisx);
end
