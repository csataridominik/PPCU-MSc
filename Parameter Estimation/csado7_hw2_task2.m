clear; close all; clc;

% Task 2:
% 1. Generate n data points using the equation:

% Number of samples:
N = 200;

% Noise:
sigma = 1.0;
epsilon = sigma * randn(N, 1);
% vector of arbitary chosen parameters for the vector product:
p = [1 2 4];

% result vector
y = zeros(N,1);

X = zeros(N,3);

for i = 1:N
    x = [1, i*randn(), i*randn()];
    X(i,1:3) = x;
    y(i) = x * transpose(p) + epsilon(i);
end

% 2. Plot the generated data points:
figure;
scatter3(X(:,1),X(:,2),y)
xlabel('X')
ylabel('Y')
zlabel('Z')
grid on;
hold on;

% 3. Use least square estimation to estimate the parameters:
% LSQ Estimator task 1:
p_hat = pinv(X' * X) * (X' * y);
% we need to use pinv() which is pseudo inverse in case of a not invertible
% matrix, or if we have too many data points, it is a long calculation to
% take.

% 3. You should get something similar to your values chosen in 1(to p):
disp([0.9776, 2.0046, 3.9999])
% These are some values I've got.


% 4. Calculate the loss function:

L = (y- X * p_hat)' * (y - X* p_hat);


%% 5.  Show how the loss function depends on the amplitude of the noise:
clc; clear;
% Fixed number of samples:
N = 200;

% Noise:
sigma = [1,2, 5,8, 10,15,20, 40,50,80,100];


% vector of arbitary chosen parameters for the vector product:
p = [1 2 4];

% result vector
y = zeros(N,1);

X = zeros(N,3);

% This is where I am going to store the losses:
Losses = zeros(size(sigma,1));

for j = 1:size(sigma,2)
    epsilon = sigma(j) * randn(N, 1);

    for i = 1:N
        x = [1, i*randn(), i*randn()];
        X(i,1:3) = x;
        y(i) = x * transpose(p) + epsilon(i);
    end

    p_hat = pinv(X' * X) * (X' * y);
    L = (y- X * p_hat)' * (y - X* p_hat);

    Losses(j) = L;

end

figure;
plot(sigma,Losses)
hold on;

scatter(sigma,Losses)
xlabel('sigma value')
ylabel('Loss')
grid on;
hold on;



%% 6. Show how the loss function depends on the sample size:
clc; clear;

% Vector of possible number of samples:
N = [10, 50, 100, 200, 500, 800 ,1200,2000,5000,7500, 10000];

% Fixed noise:
sigma = 2.0;

% vector of arbitary chosen parameters for the vector product:
p = [1 2 4];


% This is where I am going to store the losses:
Losses = zeros(size(N,1));

for j = 1:size(N,2)
    epsilon = sigma * randn(N(j), 1);

    y = zeros(N(j),1);
    
    X = zeros(N(j),3);


    for i = 1:N(j)
        x = [1, i*randn(), i*randn()];
        X(i,1:3) = x;
        y(i) = x * transpose(p) + epsilon(i);
    end

    p_hat = pinv(X' * X) * (X' * y);
    L = (y- X * p_hat)' * (y - X* p_hat);

    Losses(j) = L;

end

figure;
plot(N,Losses)
hold on;
scatter(N,Losses)
xlabel('Number of samples')
ylabel('Loss')
grid on;
hold on;
%% showing the loss in function of the noise amplitude and sample size:
% changing both sigma and N, neither is fixed.
clc; clear;


% Vector of possible number of samples and sigmas:
N = [10, 50, 100, 200, 500, 800 ,1200,2000,5000,7500, 10000];
sigma = [1,2, 5,8, 10,15,20, 40,50,80,100];


% vector of arbitary chosen parameters for the vector product:
p = [1 2 4];


% This is where I am going to store the losses:
Losses = zeros(size(N,1));

for j = 1:size(N,2)
    epsilon = sigma(j) * randn(N(j), 1);
    
    y = zeros(N(j),1);
    
    X = zeros(N(j),3);


    for i = 1:N(j)
        x = [1, i*randn(), i*randn()];
        X(i,1:3) = x;
        y(i) = x * transpose(p) + epsilon(i);
    end

    p_hat = pinv(X' * X) * (X' * y);
    L = (y- X * p_hat)' * (y - X* p_hat);

    Losses(j) = L;

end

figure;
plot3(N,sigma,Losses)
hold on;
scatter3(N,sigma,Losses)
xlabel('Number of samples')
ylabel('Simga')
zlabel('Loss')
grid on;
hold on;

% What is the expected value of θ_hatLS? How does the covariance matrix of θ_hatLS depend on Σϵ?
% - The expected value for the estimated parameters is to be close
% to the parameters(so p ~= p_hat in my case).
