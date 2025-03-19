% 1. Task:
close all; clear; clc;


% Values for initializtaion:
N = 1000; % With higher values or more simulations and its average it shows better results.
y = zeros(N,1);
sigma = 0.1;
epsilon = sigma * randn(N, 1);

theta1 = -0.12;
theta2 = -0.42;

theta_est = [0, -1];
sigma_est = 0.1;

x = -1:0.01:1;

% 2. Create grid:
[Theta1, Theta2] = meshgrid(x, x);

% 3. Prior and Normalization step:
prior = mvnpdf([Theta1(:), Theta2(:)], [0 -1], [0.5 0.5]);
prior = reshape(prior, size(Theta1));
prior = prior / sum(prior(:));


% 4. Iterate over the data series:
for k = 3:N

    % 1. Create AR:
    y(k) = theta1 * y(k-1) + theta2 * y(k-2) + epsilon(k);


    % Likelihood from seminar example:
    p_y_Dk_prev = (y(k) - Theta1 * y(k-1) - Theta2 * y(k-2)).^2/(2*sigma_est^2);
    p_y_Dk_prev = exp(-p_y_Dk_prev);

    posterior = prior .* p_y_Dk_prev;
    prior = posterior / sum(posterior(:));
end


% 5. Plot 2D posterior probability values:
figure(1);
contourf(Theta1, Theta2, posterior, 50, 'LineStyle', 'none');
colorbar;
xlabel('\theta_1');
ylabel('\theta_2');
title('Posterior Probability Density');

[max_prob, idx] = max(posterior(:));
[max_theta1_idx, max_theta2_idx] = ind2sub(size(posterior), idx);
max_theta1 = Theta1(max_theta1_idx, max_theta2_idx);
max_theta2 = Theta2(max_theta1_idx, max_theta2_idx);

% 6. The msot probable values:
disp('Candidate theta with highest probability:');
disp(max_theta1)
disp(max_theta2)



