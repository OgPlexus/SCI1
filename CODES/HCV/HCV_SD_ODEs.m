function dydt = HCV_SD_ODEs(t, y, params)
% HCV_SD_ODEs: Function representing the system of ordinary differential equations for the HCV(SD) model

% Extracting parameters
pi_S1 = params(1);
pi_S2 = params(2);
phi1 = params(3);
phi2 = params(21);
beta_11 = params(4);
beta_12 = params(5);
beta_21 = params(6);
beta_22 = params(7);
omega = params(8);
tau1 = params(9);
tau2 = params(22);
mu = params(10);
pi_N1 = params(11);
pi_N2 = params(12);
alpha_11 = params(13);
alpha_12 = params(14);
alpha_21 = params(15);
alpha_22 = params(16);
k_u = params(17);
k_i = params(18);
epsilon1 = params(19);
epsilon2 = params(20);

% Extracting variables
S1 = y(1);
S2 = y(2);
I_E1 = y(3);
I_E2 = y(4);
I_L1 = y(5);
I_L2 = y(6);
N_u1 = y(7);
N_u2 = y(8);
N_i1 = y(9);
N_i2 = y(10);

% Calculating derivatives
dydt = zeros(10, 1);
dydt(1) = pi_S1 + phi1*I_E1 - (beta_11*S1*N_i1)/(N_i1 + N_u1) - (beta_12*S1*N_i2)/(N_i2 + N_u2) - mu*S1;
dydt(2) = pi_S2 + phi2*I_E2 - (beta_21*S2*N_i1)/(N_i1 + N_u1) - (beta_22*S2*N_i2)/(N_i2 + N_u2) - mu*S2;
dydt(3) = (beta_11*S1*N_i1)/(N_i1 + N_u1) + (beta_12*S1*N_i2)/(N_i2 + N_u2) - (omega + tau1 + mu + phi1)*I_E1;
dydt(4) = (beta_21*S2*N_i1)/(N_i1 + N_u1) + (beta_22*S2*N_i2)/(N_i2 + N_u2) - (omega + tau2 + mu + phi2)*I_E2;
dydt(5) = omega*I_E1 - (tau1 + mu)*I_L1;
dydt(6) = omega*I_E2 - (tau2 + mu)*I_L2;
dydt(7) = pi_N1 - alpha_11*(I_E1 + I_L1)*N_u1/(N_i1 + N_u1) - alpha_12*(I_E2 + I_L2)*N_u1/(N_i1 + N_u1) - k_u*N_u1 + epsilon1*N_i1;
dydt(8) = pi_N2 - alpha_21*(I_E1 + I_L1)*N_u2/(N_i2 + N_u2) - alpha_22*(I_E2 + I_L2)*N_u2/(N_i2 + N_u2) - k_u*N_u2 + epsilon2*N_i2;
dydt(9) = alpha_11*(I_E1 + I_L1)*N_u1/(N_i1 + N_u1) + alpha_12*(I_E2 + I_L2)*N_u1/(N_i1 + N_u1) - k_i*N_i1 - epsilon1*N_i1;
dydt(10) = alpha_21*(I_E1 + I_L1)*N_u2/(N_i2 + N_u2) + alpha_22*(I_E2 + I_L2)*N_u2/(N_i2 + N_u2) - k_i*N_i2 - epsilon2*N_i2;

end
