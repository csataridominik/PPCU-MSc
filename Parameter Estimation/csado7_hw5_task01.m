% Homework 05 Task 01:
clear; close all; clc;

% 1. Load the downloaded file into the workspace:
vessel = load("vessel.mat");

c123 = [1,2,0;vessel.c123]; % Give the t = 0 conectrations
t = vessel.t;
V = 50;
c01 = 1;
c02 = 2;

% 1. Discretize the dynamical model with the Euler method:


y = zeros(size(t,1),3);
X = zeros(size(c123));
X2 = zeros(size(c123));
X3 = zeros(size(c123));
%phi = [F1, F2,kr];

for i = 2:size(t,1)
    
    % This part is for discretization:
    y(i,1) = (c123(i,1)-c123(i-1,1))/0.1;
    y(i,2) = (c123(i,2)-c123(i-1,2))/0.1;
    y(i,3) = (c123(i,3)-c123(i-1,3))/0.1;
    
end

for i = 1:size(t,1)
    
    X(i,:) = [0.1/V-c123(i,1)/V,-c123(i,1)/V,-c123(i,1)*c123(i,2)];
    X2(i,:) = [0.2/V-c123(i,2)/V,-c123(i,2)/V,-c123(i,1)*c123(i,2)];
    X3(i,:) = [-c123(i,3)/V,-c123(i,3)/V,c123(i,1)*c123(i,2)];
    
     
end



x_sum = [X(3:end,:);X2(3:end,:);X3(3:end,:)];
y_sum = [y(2:end,1);y(2:end,2);y(2:end,3)];
% 2. Give the model in predictive form for least squares estimation:


% This is the leas square estimation result:
theta_hat_sum = inv(x_sum'*x_sum)*x_sum'*y_sum;




%%
% 3. Estimate the unknown parameter values:
% P_prev = inv(X' * X); as mentioned in  the lecture. But with more data
% incoming I decided to initialize it like this:
P_prev = eye(3);

theta_pred_prev = P_prev*(transpose(X(2:end,:))*y(:,1));
theta_sum = zeros(size(t,1)-1,3);
%%

alpha = 1;
lambda = 0.95;

for j = 2:size(t,1)-1

    for k = 1:3
        if k == 1
            phi = [c01/V-c123(j,1)/V,-c123(j,1)/V,-c123(j,1)*c123(j,2)]';
        elseif k == 2
            phi = [c02/V-c123(j,2)/V,-c123(j,2)/V,-c123(j,1)*c123(j,2)]';
        elseif k == 3
            phi = [-c123(j,1)/V,-c123(j,2)/V,c123(j,1)*c123(j,2)]';
        end
        act = j+(k-1) * size(t,1);
        
        P = 1/lambda*(P_prev-(P_prev*phi*phi'*P_prev)/(lambda/alpha + phi' * P_prev * phi));
        P_prev = P;
        theta_pred = theta_pred_prev + alpha * P * phi * (y(j,k)-phi'*theta_pred_prev);
        theta_pred_prev = theta_pred;
        theta_sum(j,:) = theta_pred;
    end

end


%%
% Plot the time-dependent estimates of the parameters:

% Setting lambda this high means it is steadier, however it needs mroe
% data, taht is why I chose to show the plot in range [100:500]:
figure(1)
hold on;
plot(100:499,theta_sum(100:end,1))
plot(100:499,theta_sum(100:end,2))
plot(100:499,theta_sum(100:end,3))


dc1 = zeros(size(t));
dc2 = zeros(size(t));
dc3 = zeros(size(t));


for k = 1:size(t,1)-1
    curr_theta = theta_sum(k,:);
    dc1(k) = curr_theta(1)*c01/V - curr_theta(3)*c123(k,1)*c123(k,2) -(curr_theta(1)+curr_theta(2))*c123(k,1)/V;
    dc2(k) = curr_theta(1)*c02/V - curr_theta(3)*c123(k,1)*c123(k,2) -(curr_theta(1)+curr_theta(2))*c123(k,2)/V;
    dc3(k) = curr_theta(3)*c123(k,1)*c123(k,2) - (curr_theta(1)+curr_theta(2))*c123(k,3)/V;
end

sim_c1 = zeros(size(t,1)+1,1);
sim_c2 = zeros(size(t,1)+1,1);
sim_c3 = zeros(size(t,1)+1,1);

for s = 1:size(t,1)
    sim_c1(s+1) = sim_c1(s) + 0.1 * dc1(s);
    sim_c2(s+1) = sim_c2(s) + 0.1 * dc2(s);
    sim_c3(s+1) = sim_c3(s) + 0.1 * dc3(s);
end


figure(2)
hold on;
plot(t,c123(2:end,1))
plot(t,c123(2:end,2))
plot(t,c123(2:end,3))
plot(t,sim_c1(2:end),'--')
plot(t,sim_c2(2:end),'--')
plot(t,sim_c3(2:end),'--')
legend('c1','c2','c3','c1\_pred','c2\_pred','c3\_pred')
