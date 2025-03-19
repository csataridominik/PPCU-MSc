% Task 1
clear; clc; close all;
% 1. Load data:
T  = readtable("data.xlsx");
y  = T.y;
y1 = T.y1;
y2 = T.y2;
u  = T.u;

% 2. Estimate the system with LSQ:

x_g = zeros(size(y,1),2);
x_n1 = zeros(size(y,1),2);
x_n2 = zeros(size(y,1),3);

for k = 2 : size(y,1)
    
    x_g(k,:) = [y(k-1),u(k-1)];
    
    x_n1(k,:) = [y1(k-1),y(k-1)];
    
    x_n2(k,:) = [y1(k-1),y(k-1),y2(k-1)];
    
end

params12 = pinv(x_g' * x_g) * (x_g' * y);
params45 = pinv(x_n1' * x_n1) * (x_n1' * y1);
params789 = pinv(x_n2' * x_n2) * (x_n2' * y2);

disp('Parameters Estimated with LSQ:')
disp(params12)
disp(params45)
disp(params789)


% 3. Simulate the system with the estimated parameters:
y_est = zeros(size(y));
y1_est = zeros(size(y));
y2_est = zeros(size(y));

for k = 2 : size(y,1) 

    y_est(k) = params12(1)*y_est(k-1)+params12(2)*u(k-1);
    y1_est(k) = params45(1)*y1_est(k-1)+params45(2)*y_est(k-1);
    y2_est(k) = params789(1)*y1_est(k-1)+params789(2)*y_est(k-1)+params789(3)*y2_est(k-1);

end


% Plotting the results:

figure(1);
grid on;

subplot(3,3,1:3)
plot(1:size(y,1),y,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y_est,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y(k)','y\_est(k)');
title("Plot of y(k) and the estimation of it.")

subplot(3,3,4:6)
plot(1:size(y,1),y1,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y1_est,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y1(k)','y1\_est(k)');
title("Plot of y1(k) and the estimation of it.")

subplot(3,3,7:9)
plot(1:size(y,1),y2,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y2_est,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y2(k)','y2\_est(k)');

title("Plot of y2(k) and the estimation of it.")
hold off;


%% 4.
close all;

% This part is for estimating params12:
Pz = y_est *pinv(y_est'*y_est)*y_est';
Beta_2sls_12 = pinv(x_g'*Pz*x_g) * x_g' * Pz * y;

% This part is for estimating params45:
Pz = y1_est *pinv(y1_est'*y1_est)*y1_est';
Beta_2sls_45 = pinv(x_n1'*Pz*x_n1) * x_n1' * Pz * y1;

% This part is for estimating params789:
Pz = y2_est *pinv(y2_est'*y2_est)*y2_est';
Beta_2sls_789 = pinv(x_n2'*Pz*x_n2) * x_n2' * Pz * y2;

y_est_IV = zeros(size(y));
y1_est_IV = zeros(size(y));
y2_est_IV = zeros(size(y));

for k = 2 : size(y,1) 
    y_est_IV(k)  = Beta_2sls_12(1)*y_est_IV(k-1)+Beta_2sls_12(2)*u(k);
    y1_est_IV(k) = Beta_2sls_45(1)*y1_est_IV(k-1)+Beta_2sls_45(2)*y_est_IV(k-1);
    y2_est_IV(k) = Beta_2sls_789(1)*y1_est_IV(k-1)+Beta_2sls_789(2)*y_est_IV(k-1)*Beta_2sls_789(3)*y2_est_IV(k-1);
    

end

figure(2)
subplot(3,3,1:3)
plot(1:size(y,1),y,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y_est_IV,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y(k)','y\_est\_IV(k)');
title("Plot of y(k) and the estimation of it.")

subplot(3,3,4:6)
plot(1:size(y,1),y1,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y1_est_IV,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y1(k)','y1\_est\_IV(k)');
title("Plot of y1(k) and the estimation of it.")

subplot(3,3,7:9)
plot(1:size(y,1),y2,'Color','green','LineWidth',1);
hold on;
plot(1:size(y,1),y2_est_IV,'Color','[1, 0.5, 0]','LineStyle',':','LineWidth',2);
legend('y2(k)','y2\_est\_IV(k)');

title("Plot of y2(k) and the estimation of it.")
hold off;



