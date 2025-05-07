clc;clear all;close all;
N1=300;% group 1 population size
w1=0.5; 
w2=0.5;
w0=0;
syms a;
N=400;% Total population size
N2=N-N1; %% group 2 population size
A1= (1+w0).*(1+a)./N1;
A2= (1-w0).*(1-a)./N2;
sqPartC=sqrt((A1.*w1+A2.*w2).^2+4*A1.*A2.*(1-w1-w2));
A=(N/4)*((A1*w1+A2*w2)+sqPartC);
beta=1;gam=1;
C1=A1./(A1+A2);
C2=A2./(A1+A2);
%% %fig 3
close all; clc;
astar=(N1-N2)/N;
figure;fplot(a,C1,[-1,1],'Marker','none',LineWidth=2);hold on
fplot(a,C2,[-1,1],'Marker','x','MarkerEdgeColor','black',LineWidth=2);hold on
fplot(a,A,[-1,1],'Marker','o','MarkerEdgeColor','black',LineWidth=2);hold on
plot([astar, astar], [0, 2], 'r--', 'LineWidth', 2); 
%plot([astar, astar], [0, 1.2], 'r--', 'LineWidth', 2);% case: N1=N2=200
xlabel('$a$')
ylabel('$\mathcal{R}_0$, $\mathcal{C}_i$')
%legend('$\mathcal{C}_1$','$\mathcal{C}_2$','$\mathcal{R}_0^{SIR(SD)}$','a = $\frac{N_1-N_2}{N}$','location','bestoutside','Orientation','horizontal');
%text(astar, 0.8, 'a = $\frac{N_1-N_2}{N}$', 'Rotation', 90,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'black', 'Interpreter','latex');% case: N1=N2=200
text(astar, 1.4, 'a = $\frac{N_1-N_2}{N}$', 'Rotation', 90,'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Color', 'black', 'Interpreter','latex');
