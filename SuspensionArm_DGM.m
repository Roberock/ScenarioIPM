function [X,Y]=SuspensionArm_DGM(N,crackIdx)
% Data Generating Mechaism for the FRF response of the suspension are in Y
% direction
  addpath([pwd '\datasets\FRF of a damaged suspension arm\'])
  Temp=load('MFRFX654321.mat');
  load('Data_SuspensionArm.mat')
  MFRFX=Temp.MFRFX;
  MFRFX=MFRFX(Data_SuspensionArm.Cracktype==crackIdx,:);
  GoodSimulationIndx=find(~isnan(MFRFX(:,1)));
  MFRFX=MFRFX(GoodSimulationIndx,50:5:end-50);
  ReferenceFRFY=mean(MFRFX);
  Std_vals=std(MFRFX);
  
   Std_vals=std(MFRFX)*4;
%   if crackIdx==1
%   ReferenceFRFY=MFRFX(1,50:5:end-50);
%   elseif crackIdx==2
%   ReferenceFRFY=MFRFX(600,50:5:end-50); 
%   elseif crackIdx==3
%    ReferenceFRFY=MFRFX(1200,50:5:end-50);
%    elseif crackIdx==4    
%    ReferenceFRFY=MFRFX(1800,50:5:end-50);   
%    elseif crackIdx==5  
%   ReferenceFRFY=MFRFX(2050,50:5:end-50);   
%      elseif crackIdx==6  
%   ReferenceFRFY=MFRFX(2700,50:5:end-50);
%   end 
  
 Xref=1:1:length(ReferenceFRFY);
[X,Y]=deal(zeros(1,N));
for i=1:N
    Idx=Xref(randi(length(Xref)-1))+1;
    rndX=unifrnd(0,1);
    X(i)=Idx-rndX;
   % Y(i)=normrnd(ReferenceFRFY(Idx),ReferenceFRFY(Idx)/10);
    InterpolMean=ReferenceFRFY(Idx)-rndX*((ReferenceFRFY(Idx)-ReferenceFRFY(Idx-1)));
    InterpolSTD=Std_vals(Idx)-rndX*((Std_vals(Idx)-Std_vals(Idx-1)));
    Y(i)=normrnd(InterpolMean,InterpolSTD);
    if Y(i)<0
        Y(i)=0;
    end
end

end

 