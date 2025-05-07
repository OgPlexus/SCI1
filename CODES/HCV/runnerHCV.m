clc;clear all;close all;
% beta11=0.072/2;beta12=0.072/2;beta21=0.072/2;beta22=0.072/2;
% alpha11=4/2;alpha12=4/2;alpha21=4/2;alpha22=4/2;
% beta11=0.072/2;beta12=0;beta21=0;beta22=0.072/2;
% alpha11=2;alpha12=0;alpha21=0;alpha22=2;
% beta11=0.072/2;beta12=0.072/2*0.1;beta21=0.072/2*0.1;beta22=0.072/2;
% alpha11=4/2;alpha12=4/2*0.1;alpha21=4/2*0.1;alpha22=4/2;
beta11=0.072/2;beta22=0.072/2;
alpha11=4/2;alpha22=4/2;


ku1=0.143; ku=[ku1;ku1];
piS1=47/2; piS=[piS1;piS1];
muH1=2.7*10^(-4);muH=[muH1;muH1];
ki1=0.143;ki=[ki1;ki1];
epsN1=1.17;epsN=[epsN1;epsN1];
piN1=3.14/2*10^4;piN=[piN1;piN1];
tauH1=0.011;tauH=[tauH1;tauH1];
omgH1=0.006; omgH=[omgH1;omgH1];
phiH1=4.7*10^(-3);phiH=[phiH1;phiH1];
a=[0 0.01 0.05 1];%a=linspace(0,4,1000);
tauH=[0.04;0.0001];
%phiH=[0.05;5*10^(-8)];
%epsN=[2;0.33];
for i=1:length(a)
    %%tauH=[tauH1;a(i)*tauH1];
beta12=0.072/2*a(i);beta21=0.072/2*a(i);
alpha12=4/2*a(i);alpha21=4/2*a(i);
betaTran=[beta11, beta12;beta21,beta22];
alphaTran=[alpha11,alpha12;alpha21,alpha22];
%epsN=[epsN1;a(i)*epsN1];
%phiH=[phiH1;a(i)*phiH1];
%tauH=[tauH1;a(i)*tauH1];
para=[ku piS muH ki epsN piN tauH omgH phiH];
[H W]=SCIhcv(betaTran,alphaTran,para);
Ro(i)=sqrt(eigs(H*W,1));
Ro1(i)=sqrt(H(1,1)*W(1,1));
Ro2(i)=sqrt(H(2,2)*W(2,2));
valC1(i)=(sqrt(eigs(H*W,1))-sqrt(H(2,2)*W(2,2)))/Ro(i);
valC2(i)=(sqrt(eigs(H*W,1))-sqrt(H(1,1)*W(1,1)))/Ro(i);
end
%%



%%
plot(a,valC1,'LineWidth',2)
hold on
plot(a,valC2,'LineWidth',2)
plot(a,Ro,'k--','LineWidth',2)
xlabel('$\epsilon_2/\epsilon_1$ with $\epsilon_1=1.17$','Interpreter','latex');
ylabel('SCI:$\mathcal{C}_j$ or $\mathcal{R}_0^{HCV(SD)}$','Interpreter','latex');
legend('$\mathcal{C}_1$','$\mathcal{C}_2$','$\mathcal{R}_0^{HCV(SD)}$','location','best')