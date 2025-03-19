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

for k = 3:N+2
    y(k) = 1.5*y(k-1)-0.7*y(k-2)-0.25*u(k-1)+1.64*u(k-2)-1.6*epsilon(k-1)+0.75*epsilon(k-2);
end


% 3. Create a figure showing the output of the system:
figure;
grid on;
hold on;
plot((1:N+2),y)

% 4. Apply a moving average transformation:
M = movmean(y,5);

% 5. Plot the transformed data on the same figure along y:

hold on;
plot((1:size(M,1)),M,'Color','red')

% 6. Approximate the transformed data with a cubic spline:
xx = 0:0.05:N;
yy = spline((1:size(M,1)),M,xx);
plot(xx,yy,'Color','green')
legend('y(k)','Moving Average','Spline')

% 7. Calculate the first, second and third derivative of the spline:

figure;
hold on;
grid on;

first_diff = diff(yy);
second_diff= diff(first_diff);
third_diff = diff(second_diff);

xx_fd = xx(1:end-1);
plot(xx_fd,first_diff,'Color','blue');

xx_sd = xx(1:end-2);
plot(xx_sd,second_diff,'Color','red');

xx_td = xx(1:end-3);
plot(xx_td,third_diff,'Color','green');

legend('First Derivative','Second Derivative','Third Derivative')

