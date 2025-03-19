% 1. Load the downloaded file into the workspace:
clear; close all; clc;

X = load("paramest_hw1_task2.mat").X;
N=1000;

% 2. Calculate the mean value and covariance matrix of the data:
M = mean(X);
C = cov(X);

% 3. Plot the multivariate probability density function having the previously calculated mean and covariance:
x_range = min(X(:,1)):0.1:max(X(:,1));
y_range = min(X(:,2)):0.1:max(X(:,2));

[Xx, Yy] = meshgrid(x_range, y_range);

pdf = mvnpdf([Xx(:), Yy(:)], M, C);
pdf = reshape(pdf, size(Xx));
pdf = pdf / sum(pdf(:));

figure;
hold on;
mesh(Xx,Yy,pdf)
title('Multivariate Probability Density Function');

% 4. Check the positive definiteness of the covariance matrix:
disp('Covariance Matrix:')
disp(C)

try chol(A)
    disp('Matrix is symmetric positive definite.')
catch ME
    disp('Matrix is not symmetric positive definite')
end

