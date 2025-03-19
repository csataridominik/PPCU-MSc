% Project 2024.06.07. | Dominik Csatári | FV1TW4 | Parameter Estimation
clear;clc; close all;

% The input values are defined in the documentation:
% https://www.mathworks.com/help/simbio/ug/defining-reaction-rates-with-enzyme-kinetics.html
% The data simulation process was ChatGPT assisted.

% Create a SimBiology model:
model = sbiomodel('EnzymeKinetics');

% Add compartments:
compartment = addcompartment(model, 'Compartment');

% Add species:
E_original = addspecies(compartment, 'E', 'InitialAmount', 4);
S_original = addspecies(compartment, 'S', 'InitialAmount', 8);
C_original = addspecies(compartment, 'C', 'InitialAmount', 0);
P_original = addspecies(compartment, 'P', 'InitialAmount', 0);

% Add parameters:
k1 = addparameter(model, 'k1', 2);
k_r = addparameter(model, 'k_r', 1);
k2 = addparameter(model, 'k2', 1.5);

% Reactions
reaction1 = addreaction(model, 'E + S <-> C');
reaction2 = addreaction(model, 'C -> P + E');

% Set kinetics adn parameters:
kinetics1 = addkineticlaw(reaction1, 'MassAction');
kinetics2 = addkineticlaw(reaction2, 'MassAction');

addparameter(kinetics1, 'k1', 'Value', k1.Value, 'ConstantValue', true);
addparameter(kinetics1, 'k_r', 'Value', k_r.Value, 'ConstantValue', true);
set(kinetics1, 'ParameterVariableNames', {'k1', 'k_r'});

addparameter(kinetics2, 'k2', 'Value', k2.Value, 'ConstantValue', true);
set(kinetics2, 'ParameterVariableNames', 'k2');

% I chose to set simulation time to 5, as after that the simulation comes
% to a steady-state, it does not change.
configset = getconfigset(model);
set(configset, 'StopTime', 5);

% Setting the timestep:
set(configset.SolverOptions, 'AbsoluteTolerance', 1e-8);
set(configset.SolverOptions, 'RelativeTolerance', 1e-5);

set(configset.SolverOptions, 'MaxStep', 0.1); 

% Simulate:
simdata = sbiosimulate(model);

% Plot the results
sbioplot(simdata);



%% Save and export model:
%sbmlexport(model, 'EnzymeKinetics.sbproj');
sbiosaveproject allSimBioModels
%save('EnzymeKinetics.mat', 'model');

%% 1. Estimation LSQ Estimator:

data = simdata.Data;
E_original_data = simdata.Data(:,1);
S_original_data = simdata.Data(:,2);
C_original_data = simdata.Data(:,3);
P_original_data = simdata.Data(:,4);
t = simdata.Time;

n = size(simdata.Data(:,1),1);


y = zeros(n,4);
% Discretize the dynamical model with Euler method:
% (Note: Need to calculate ts at each timestap as Simulink uses adaptive stepping
% mechanism.)
for i=2:n
    ts = t(i)-t(i-1);
    y(i,1) = (E_original_data(i)-E_original_data(i-1))/ts ;
    y(i,2) = (S_original_data(i)-S_original_data(i-1))/ts;
    y(i,3) = (C_original_data(i)-C_original_data(i-1))/ts ;
    y(i,4) = (P_original_data(i)-P_original_data(i-1))/ts ;

end

% Add noise:

simulations = 100;
noises = 5;

EstimatedParameters = zeros(simulations,noises,3);

for noise=1:noises
    sigma = 0.2*noise;
    
    for sim=1:(simulations)
        E = E_original_data + sigma * randn(n,1);
        S = S_original_data + sigma * randn(n,1);
        C = C_original_data + sigma * randn(n,1);
        P = P_original_data + sigma * randn(n,1);

        % 1. predicito method is LSQ Estimator:
        % Predicting k2 here:
        k2_est = pinv(C' * C) * (C' * y(:,4));
        
        % Predicting k1 here:
        y_p2 = y(:,2)+y(:,1)-y(:,4);
        X = [-2*E.*S, 2*C];
        
        k1_and_k1_r_est = pinv(X' * X) * (X' * y_p2);
        k1_est = k1_and_k1_r_est(1);
        
        % Predicting k1_r here:
        y_p3 = y(:,3)+y(:,4)-k1_est*E.*S;
        X = -C;
        k1_r_est = pinv(X' * X) * (X' * y_p3);
        EstimatedParameters(sim,noise,:) = [k1_est,k1_and_k1_r_est(2),k2_est];
    end
end

for noise = 1:noises
    disp(['Average estimated parameters at noise level ', num2str(noise), ' are:'])
    disp(['k1: ', num2str(mean(EstimatedParameters(:,noise,1)))])
    disp(['k1_r: ', num2str(mean(EstimatedParameters(:,noise,2)))])
    disp(['k2: ', num2str(mean(EstimatedParameters(:,noise,3)))])
    disp('')
end


%% 1. Estimation MLE Estimator:

noises = 3;

EstimatedParameters_mle = zeros(noises,3);

sigmas = [0.1,0.2,0.3];
for h=1:3
    sigma = sigmas(h);
    disp(['Estimated parameters with MLE at noise level ', num2str(h), ' are:'])
    
    E = E_original_data + sigma * randn(n,1);
    S = S_original_data + sigma * randn(n,1);
    C = C_original_data + sigma * randn(n,1);
    P = P_original_data + sigma * randn(n,1);
    
    range = 500;
    theta1_range = linspace(0,5,range);
    theta2_range = linspace(0,5,range);
    
    Z = zeros(size(theta1_range,1),size(theta1_range,1));
    Z_k2 = zeros(size(theta1_range,1),1);
    
    for i = 1:range
        a = theta1_range(i);
        for j = 1:range
            b = theta2_range(j);
            Z(i,j) = log_likelihood(E,S,C,y,a,b,sigma );
        end
        Z_k2(i) = log_likelihood_k2(C,y,a,sigma );
    end
    
    [M,I] = max(Z,[],"all");
    [row,col] = ind2sub(size(Z),I);
    
    [M_k2,I_k2] = max(Z_k2,[],"all");
    [row_k2,col_k2] = ind2sub(size(Z_k2),I_k2);
    
    disp(['k1: ', num2str(theta1_range(row))])
    disp(['k1_r: ', num2str(theta2_range(col))])
    disp(['k2: ', num2str(theta1_range(col_k2))])
    EstimatedParameters_mle(h,:) = [theta1_range(row),theta2_range(col),theta1_range(col_k2)];
end

%% Generating data with the obtained parameters:
% For LSQ first, for all the noise levels:

% Noise level 1-5:
% estiamte() is a function at the alst sectio of this code,
% estimating data based on the parameters:
for it=1:5
    Estimated_data = estimate(n,t,mean(EstimatedParameters(:,it,1)),mean(EstimatedParameters(:,it,2)),mean(EstimatedParameters(:,it,3)));
    plotting([S_original_data,E_original_data,C_original_data,P_original_data], Estimated_data,n,it+1, ['LSQ noise level ',num2str(it)]);
end

% For MLE, for all the noise levels:
for it=1:3
    Estimated_data = estimate(n,t,EstimatedParameters_mle(it,1),EstimatedParameters_mle(it,2),EstimatedParameters_mle(it,3));
    plotting([S_original_data,E_original_data,C_original_data,P_original_data], Estimated_data,n,it+6, ['MLE noise level ',num2str(it)]);
end

%% Closing the generated figures:
%close all;
%% Functions:

function loglikelihood = log_likelihood(E,S,C,y,a,b,sigma)
    loglikelihood = 0;
    sigma1 = sqrt(sigma);
    
    % y_p2 = y(:,2)+y(:,1)-y(:,4);
    % X = [-2*E.*S, 2*C];
    for k = 2:size(y,1)
            loglikelihood = loglikelihood + (1/(sigma1*sqrt(2*pi))) +...
                log(exp(-((y(k,2)+y(k,1)-y(k,4)+a*2*E(k)*S(k)-2*b*C(k))^2/(2*sigma1^2))));
    end

end

function loglikelihood = log_likelihood_k2(C,y,a,sigma)
    loglikelihood = 0;
    sigma1 = sqrt(sigma);
    
    for k = 2:size(y,1)
            loglikelihood = loglikelihood + (1/(sigma1*sqrt(2*pi))) +...
                log(exp(-((y(k,4)-a*C(k))^2/(2*sigma1^2))));
    end

end

% This function is for plotting the original data (with no noise) and
% plotting the data generated with the calculated parameters. It also 
% calculates the squarred losses between them for each instances:
function plotting(Original, Estimated,n,fig_num,Current_process)
    
    figure(fig_num);
    grid on;
    hold on;
    plot(linspace(0,5,n),Original(:,1),'Color','red','LineWidth',3);
    plot(linspace(0,5,n),Original(:,2),'Color','green','LineWidth',3);
    plot(linspace(0,5,n),Original(:,3),'Color','blue','LineWidth',3);
    plot(linspace(0,5,n),Original(:,4),'Color','magenta','LineWidth',3);
    
    plot(linspace(0,5,n),Estimated(:,1),'Color','red','LineStyle',':','LineWidth',3);
    plot(linspace(0,5,n),Estimated(:,2),'Color','green','LineStyle',':','LineWidth',3);
    plot(linspace(0,5,n),Estimated(:,3),'Color','blue','LineStyle',':','LineWidth',3);
    plot(linspace(0,5,n),Estimated(:,4),'Color','magenta','LineStyle',':','LineWidth',3);
    title(Current_process)
    
    hold off;
    legend('S(t)','E(t)','C(t)','P(t)','S\_est(t)','E\_est(t)','C\_est(t)','P\_est(t)')
    % At this part it calculates the losses:
    sq_loss = zeros(4,1); 
    for j = 1:n
        
        sq_loss(1) = sq_loss(1) + (Original(j,1) - Estimated(j,1))^2;
        sq_loss(2) = sq_loss(2) + (Original(j,2) - Estimated(j,2))^2;
        sq_loss(3) = sq_loss(3) + (Original(j,3) - Estimated(j,3))^2;
        sq_loss(4) = sq_loss(4) + (Original(j,4) - Estimated(j,4))^2;
        
    end

    disp(['Squared losses for process ', Current_process,':'])
    disp(sq_loss)

end

% This function generates data with the estimated parameters
function Estimated_data = estimate(n,t,k1_est,k1_r_est,k2_est)

    Estimated_data = zeros(n,4);
    Estimated_data(1,1) = 8;
    Estimated_data(1,2) = 4;
    
    for h=1:n-1
        
        dS = -k1_est * Estimated_data(h,2)*Estimated_data(h,1) + k1_r_est*Estimated_data(h,3);
        dE = -k1_est * Estimated_data(h,2)*Estimated_data(h,1) + (k1_r_est+k2_est)*Estimated_data(h,3);
        dC = k1_est * Estimated_data(h,2)*Estimated_data(h,1) - (k1_r_est+k2_est)*Estimated_data(h,3);
        dP = k2_est*Estimated_data(h,3);
  
        % Times stepping:
        ts = t(h+1)-t(h);

        Estimated_data(h+1,1) = Estimated_data(h,1) + ts*dS;
        Estimated_data(h+1,2) = Estimated_data(h,2) + ts*dE;
        Estimated_data(h+1,3) = Estimated_data(h,3) + ts*dC;
        Estimated_data(h+1,4) = Estimated_data(h,4) + ts*dP;
    end
end





