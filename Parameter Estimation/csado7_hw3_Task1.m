% Task 1 Homework 1 part:
clear; close all; clc;

% 1. Load the downloaded file into the workspace:
u = load("paramest_hw1_task1.mat");
u = [0; 0; u.u];
N = 100;

% 2. Calculate the output y(k) of the system:
% epsilon is expected value 0 and σ = 1.2.
sigma = 1.2;
epsilon = sigma * randn(N+2, 1);
epsilon(1:2) = 0;

y = zeros(N+2,1)+2;

X = zeros(N,4);
% p = [1.5, -0.7, -0.25, 1.64];

for k = 3:N+2
    y(k) = 1.5*y(k-1)-0.7*y(k-2)-0.25*u(k-1)+1.64*u(k-2)-1.6*epsilon(k-1)+0.75*epsilon(k-2);
    x = [y(k-1),y(k-2),u(k-1),u(k-2)];
    X(k,1:4) = x;
end


% LSQ Estimator:
p_hat = pinv(X' * X) * (X' * y);
disp(p_hat)
% Some result I have got for p_hat = [1.2860, -0.4961, -0.2532, 1.6504]

% Task 1, second exercise: 
% Simulate the system using the estimated parameters:
y_hat = zeros(N+2,1)+2;
for k = 3:N+2
    y_hat(k) = p_hat(1)*y_hat(k-1)+p_hat(2)*y_hat(k-2)+p_hat(3)*u(k-1)+p_hat(4)*u(k-2); 
end


% Plot the output alongside the results from Homework 1:
figure(1);
grid on;
hold on;
plot((1:N+2),y,'Color','red')
plot((1:N+2),y_hat,'Color','blue')
legend('y(k) with real parameters','y(k) with estimated parameters')



