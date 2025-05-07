clc;clear all; close all;
w1=0.5; 
w2=0.5;
w0=0;
syms a a1 N1;
N=400;
N2=N-N1;
a1=[(N1-N2)/N,-0.99,-0.5,0,0.5,0.99];
%% %Fig 2
close all
markers = {'none', 'x', '*', '|', 'o', '>'};
 for i=1:length(a1)
     marker = markers{i};
 a=a1(i);
A1= (1+w0).*(1+a)./N1;
A2= (1-w0).*(1-a)./N2;
sqPartC=sqrt((A1.*w1+A2.*w2).^2+4*A1.*A2.*(1-w1-w2));
A=(N/4)*((A1*w1+A2*w2)+sqPartC);
fplot(N1,A,[50,350],'Marker',marker,'MarkerSize',8,'MarkerEdgeColor','black',LineWidth=2);hold on
 end
xlabel('$N_1$')
ylabel('$\mathcal{A}$')
ylim([0.5,4])
legend('$a=\frac{N_1-N_2}{N}$','$a=-0.99$','$a=-0.5$','$a=0$','$a=0.5$','$a=0.99$','location','best');