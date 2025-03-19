% Task 2
clear; close all; clc;

%  Discretize the model using the Euler method:

a = 0;
b = 50;
step=0.01;
n = (b-a)/step;
t = linspace(0,1,n);

% SIR model state variables:
S = zeros(size(t));
I = zeros(size(t));
R = zeros(size(t));

% Initial values:
S(1) =  0.9999;
I(1) = 0.0001;
R(1) = 0;

% Parameters:
beta  = 1.1;
gamma = 0.4;

% Discretization:
for i = 2:size(t,2)
    % x(k+1):
    S(i) = S(i-1) + step * -beta*S(i-1)*I(i-1);
    
    I(i) = I(i-1) + step * (beta*S(i-1)*I(i-1) - gamma* I(i-1));
    
    R(i) = R(i-1) + step * gamma * I(i-1); 
    
end

figure(1);
grid on;
hold on;
plot(t,S,'Color','red');
plot(t,I,'Color','blue');
plot(t,R,'Color','green');
hold off;

legend('S(t)','I(t)','R(t)')

sir = readmatrix('sir.csv');

S = sir(:,1);
I = sir(:,2);
R = sir(:,3);


%%

clc;
sir = readmatrix('sir.csv');

S = sir(:,1);
I = sir(:,2);
R = sir(:,3);


s = -step*S.*I;
s_x = [S,s];
s_x = s_x(1:end-1,1:2);

r_x = I(1:end-1,1);
r_y = R(2:end,1);
p_hat_r = pinv(r_x' * r_x) * (r_x' * r_y);
p_hat_s = pinv(s_x' * s_x) * (s_x' * S(2:end));

% Simulate the system with the estimated parameters:

% Values calculated for estimated values:
S_est = zeros(size(S));
I_est = zeros(size(I));
R_est = zeros(size(R));

S_est(1) = 0.99;
I_est(1) = 0.01;
R_est(1) = 0;

beta_est  = p_hat_s(2);
gamma_est = p_hat_r(1);

for i = 2:size(S,1)
    % x(k+1):
    S_est(i) = S_est(i-1) + step * -beta_est*S_est(i-1)*I_est(i-1);
    
    I_est(i) = I_est(i-1) + step * (beta_est*S_est(i-1)*I_est(i-1) - gamma_est* I_est(i-1));
    
    R_est(i) = R_est(i-1) + step * gamma_est * I_est(i-1); 
    
end




% Original values obtained from the file:
figure(2);
grid on;
hold on;
plot(linspace(0,1,size(S,1)),S,'Color','red','LineWidth',1);
plot(linspace(0,1,size(I,1)),I,'Color','blue','LineWidth',1);
plot(linspace(0,1,size(R,1)),R,'Color','green','LineWidth',1);

plot(linspace(0,1,size(S_est,1)),S_est,'Color','red','LineStyle',':','LineWidth',2);
plot(linspace(0,1,size(I_est,1)),I_est,'Color','blue','LineStyle',':','LineWidth',2);
plot(linspace(0,1,size(R_est,1)),R_est,'Color','green','LineStyle',':','LineWidth',2);
hold off;
legend('S(t)','I(t)','R(t)','S\_est(t)','I\_est(t)','R\_est(t)')




