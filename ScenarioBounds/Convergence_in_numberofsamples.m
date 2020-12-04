clc
clear all
m=10:1:1000;
BetaTarget=10^-3;
ScenarioBound=zeros(1,length(m));
for n=1:length(m)
    Ndiscarded=0;
    Ndim=10;
    epsilon=linspace(0.0001,0.9999,10^4);
    beta=getConfidence_ConvexDiscard(m(n),epsilon,Ndiscarded,Ndim); 
    %epsilon_nnconvex(n)=getConfidence_nonconvex(Ndim,m(n),BetaTarget);
     ScenarioBound(n)=epsilon(find(beta<=BetaTarget,1,'first')); 
end

subplot(2,2,1)
plot(m,ScenarioBound./(log(m)./m),':k','LineWidth',3)
xlabel('m')
legend('b(m)/(log(m)/m)')
grid on
set(gca,'FontSize',20)
 
subplot(2,2,2)
plot(m,ScenarioBound./(1./m),'r')
xlabel('m')
legend( 'b(m)/(1/m)' )
grid on
set(gca,'FontSize',20)
 
subplot(2,2,3)
plot(m,(ScenarioBound./(log(log(m)))),'b')
xlabel('m')
legend( 'b(m)/(log(log(m)))' )
grid on
set(gca,'FontSize',20)
 
subplot(2,2,4)
plot(m,(ScenarioBound./(log(m).^(1/2))),'g','LineWidth',2)
xlabel('m')
legend( 'b(m)/(log(m)^0^.^5) ')
grid on
set(gca,'FontSize',20)