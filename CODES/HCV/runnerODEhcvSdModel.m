% Define time span
close all; clear all; clc
tspan = [0 9000]; 

% parameters
pi_S1 = 47/2;
pi_S2 = 47/2;
phi1 = 4.7*10^(-3);
phi2 = 4.7*10^(-3);
% phi1 = 5*10^(-2);%
% phi2 = 5*10^(-8);%
beta_11 = 0.072/2;
beta_12 = 0;%beta_12 = 0.072/2;%
beta_21 = 0;%beta_21 = 0.072/2;%
beta_22 = 0.072/2;
omega = 0.006;
tau1 = 0.011;%
tau2 = 0.011;%
% tau1 = 0.04;
% tau2 = 0.0001;
mu = 2.7*10^(-4);
pi_N1 = 3.14/2*10^4;
pi_N2 = 3.14/2*10^4;
alpha_11 = 2;
alpha_12 = 0;%alpha_12 = 2;%
alpha_21 = 0;%alpha_21 = 2;%
alpha_22 = 2;
k_u = 0.143;
k_i = 0.143;
% epsilon1 = 1.17;
% epsilon2 = 1.17;
epsilon1 = 2;
epsilon2 = 0.33;%
% Extracting variables
S1_0 = 170000/2;
S2_0 = 170000/2;
I_E1_0 = 1;
I_E2_0 = 1;
I_L1_0 = 0;
I_L2_0 = 0;
N_u1_0 = 220000/2;
N_u2_0 = 220000/2;
N_i1_0 = 1;
N_i2_0 = 1;

% Define initial conditions
y0 = [S1_0, S2_0, I_E1_0, I_E2_0, I_L1_0, I_L2_0, N_u1_0, N_u2_0, N_i1_0, N_i2_0];
%%
a=[0 0.01 0.05 1]
figure; hold on;
for i=1:length(a)
    beta_12=a(i)*beta_11; beta_21=a(i)*beta_22;
    alpha_12=a(i)*alpha_22; alpha_21=a(i)*alpha_11;% Define parameters  
params = [pi_S1, pi_S2, phi1, beta_11, beta_12, beta_21, beta_22, omega, tau1, mu, pi_N1, pi_N2, alpha_11, alpha_12, alpha_21, alpha_22, k_u, k_i, epsilon1, epsilon2,phi2,tau2];

% Solve the ODEs
[t, y] = ode45(@(t, y) HCV_SD_ODEs(t, y, params), tspan, y0);
% resultstau{i}=y;
% timetau{i}=t;
semilogy(t, y(:, 3),'LineWidth',2)
drawnow
end
% save('resultstau.mat', 'resultstau', 'timetau');
ylabel('Community 2: $I_{E_2}$ Population');
set(gca, 'YScale', 'log');
years = 5:5:ceil(max(t)/365); % Calculate number of years
xticks(years*365); % Set ticks at the start of each year
xticklabels(years); % Set tick labels as years
%xtickangle(45);
% Label the x-axis
xlabel('Time (years)');
 legend('Isolated','$1\%$','$5\%$','Homogeneous mixing',Location='southoutside',Orientation='horizontal');

%%
% Plot the results (example)
figure;
semilogy(t, y(:, 3), t, y(:, 4),'LineWidth',2);
xlabel('Time(days)');
ylabel('Population');
legend('$I_{E_1}$', '$I_{E_2}$',Location='best');
title(sprintf('$I_{E_1}$ Population dynamics over time:  $\\epsilon_1=2,\\ \\epsilon_1=0.33$'),Interpreter="latex");
