% Task 1:
close all; clc; clear;

% 1. Let the simulation time horizon contain N = 1000 samples:
N = 1000;
sigma1 = sqrt(1.2);
sigma2 = sqrt(3);

epsilon1 = sigma1 * randn(N+1, 1);
epsilon1(1) = 0;
epsilon2 = sigma2 * randn(N+1, 1);
epsilon2(1) = 0;


% With this line we can check, that variance is 1.2 and 3 respectively:
% stats = [mean(epsilon1) std(epsilon1) var(epsilon1)];
% disp(stats)

theta1 = 0.3;
theta2 = 0.7;

u = 1 * randn(N+1, 1);

y = zeros(N+1,1);

for k = 2:N+1
    if k > 500
        y(k) = theta1 * y(k-1) + theta2 * u(k-1) + epsilon1(k);
    else
        y(k) = theta1 * y(k-1) + theta2 * u(k-1) + epsilon2(k);
    end
end

% 3. Plot the log-likelihood surface:
n = 200; 
theta1_range = linspace(0,1,n);
theta2_range = linspace(0,1,n);

% Creating a surface matrix values for each theta:
Z = zeros(size(theta1_range,1),size(theta1_range,1));

for i = 1:n
    a = theta1_range(i);
    for j = 1:n
        b = theta2_range(j);
        Z(i,j) = log_likelihood(u,y,a,b);
    end
end


figure(1);
surf(theta1_range,theta1_range,Z)
xlabel('Theta1')
ylabel('Theta2')
zlabel('Log Likelihood')

% 4. Estimate the parameters:
% I have estimated the parameters with the method mentioned in the seminar.
% I have saved the values for log_likelihood, than I searched for minimum
% of the surface:
% minimum of loglikelihood matrix 'Z':
[M,I] = min(Z,[],"all");
[row,col] = ind2sub(size(Z),I);
disp(theta1_range(row))
disp(theta2_range(col))

% 5. Run 100 simulations of the model:
%% Here I run 100 simulations for the task above:
% Note, that for Z, I have used a smaller n, which would result in worse
% estimation, but faster run time.
est_log = zeros(100,2);
est_lse = zeros(100,2);

for sim = 1:100
    N = 1000;
    sigma1 = sqrt(1.2);
    sigma2 = sqrt(3);
    
    epsilon1 = sigma1 * randn(N+1, 1);
    epsilon1(1) = 0;
    epsilon2 = sigma2 * randn(N+1, 1);
    epsilon2(1) = 0;
    
    
    % With this line we can check, that variance is 1.2 and 3 respectively:
    % stats = [mean(epsilon1) std(epsilon1) var(epsilon1)];
    % disp(stats)
    
    theta1 = 0.3;
    theta2 = 0.7;
    
    y = zeros(N+1,1);
    
    for k = 2:N+1
        if k > 500
            y(k) = theta1 * y(k-1) + theta2 * u(k-1) + epsilon1(k);
        else
            y(k) = theta1 * y(k-1) + theta2 * u(k-1) + epsilon2(k);
        end
    end
    
    n = 100; 
    theta1_range = linspace(0,1,n);
    theta2_range = linspace(0,1,n);
    
    % Creating a surface matrix values for each theta:
    Z = zeros(size(theta1_range,1),size(theta1_range,1));
    
    for i = 1:n
        a = theta1_range(i);
        for j = 1:n
            b = theta2_range(j);
            Z(i,j) = log_likelihood(u,y,a,b);
        end
    end
    
    % 4. Estimate the parameters:
    [M,I] = min(Z,[],"all");
    [row,col] = ind2sub(size(Z),I);


    % lse estimations:
    X = zeros(size(u,1)-1,2);
    
    X(:,1) = y(1:end-1);
    X(:,2) = u(1:end-1);
    
    % Because of some dimensionality error, I did
    % the multiplication one-by-one.
    temp = [X(:,1)' * y(2:end) ; X(:,2)' * y(2:end)];
    % Lse estimation is stored in Lse:
    p_hat = pinv(X' * X) * temp;

    % Storing the reuslts for both estimations:
    est_log(sim,:) = [theta1_range(row), theta2_range(col)];
    est_lse(sim,:) = p_hat;

end

% With this we can see the results are simular, both are close to the real
% values of the parameters. This part was not mandatory. The figure can be 
% seen in the lab summary.
figure(2)
plot(est_lse)
grid on;
hold on;
plot(est_log)
plot(0.7*ones(100,1),'Color','black',LineWidth=2)
plot(0.3*ones(100,1),'Color','black',LineWidth=2)
legend('est\_log theta1','est\_log theta2','est\_lse theta1','est\_lse theta2')


% Compute and compare the covariance matrices:
disp("This is the covariance matrix of est_log:")
disp(cov(est_log))
disp("This is the covariance matrix of est_lse:")
disp(cov(est_lse))

% It is apparent by Figure 2, that both method change the same direction
% whe the e(k) changing. this means, that by checking the following:
% disp(cov(est_lse(:,1),est_log(:,1)))
% dispcov(est_lse(:,2),est_log(:,2)))
% We would still get low variance.

