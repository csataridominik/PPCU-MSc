close all; clear;clc;
% 1. Load the data:
load('keydyn_db.mat');


% 2. Normalize the data:
for i = 1:length(recordings)
    data = recordings{i, 2};
    
    data(2, :) = data(2, :) / 1000;
    % Normalize key codes
    data(3, :) = data(3, :) / 255;
    recordings{i, 2} = data;
end

% 3. Feature extraction:
features_mean = zeros(24,4);
features_euclidean_distance = zeros(24,4,40);

for i = 1:length(recordings)
    data = recordings{i, 2};
    hl = diff( data(2, :));
    il = diff(data(2, data(1, :) == 68)) - diff(data(2, data(1, :) == 85));
    pl = diff(data(2, data(1, :) == 68));
    rl = diff(data(2, data(1, :) == 85));
    
    features_mean(i,:) = [mean(il), mean(pl), mean(rl),mean(hl)];
    features_euclidean_distance(i,:,1:40) = [il(1:40); pl(1:40); rl(1:40); hl(1:40)];
    
end
% Step 4: Classification
k = 8;
[idx, centroids] = kmeans(features_mean, k);

figure(1);
scatter3(features_mean(:,1),features_mean(:,2),features_mean(:,3), 'filled');
xlabel('mean il');
ylabel('mean pl');
zlabel('mean rl');
title('3D Plot of Features');

figure(2);
hold on;
for i = 1:8
    scatter3(features_mean(idx==i,1), features_mean(idx==i,2), features_mean(idx==i,3), '.');
end

%% This part is just for indices matching and displaying the rates of the k search matching:
cell = {};
rates = zeros(8,1);

for i = 1:8
    indices = find(idx == i);

    temp = [0, recordings{indices(1), 1}];
    for j = 1:length(indices)
        
         content = find(temp(:,2) == recordings{indices(j), 1});
         if isempty(content)
             temp = [temp; 1, recordings{indices(j), 1}];
         else
             temp(content,1) = temp(content,1) +1;
         end
    end
    
    % Greedy choice:
    [max_value, max_index] = max(temp(:,1));
    disp("Success rate for each cluster:")
    disp(max_value/length(temp))

end

%% Euclidean Distance from references:


references = zeros(8,1);


rates = zeros(8,2);
it = 1;
for i = 1:8
    reference = squeeze(features_euclidean_distance(it,:,1:40));
    it = it+3;
    
    distances = zeros(23, 1);
    for k = 1:23
        
        distances(k) = norm(squeeze(features_euclidean_distance(k,:,:)) - reference);
    end
    % Get the 2 closest ones:
    [sorted_distances, sorted_indices] = sort(distances);
    
    original = recordings{it-3,1};
    closest_vector_1 = recordings{sorted_indices(1),1};
    closest_vector_2 = recordings{sorted_indices(2),1};
    if closest_vector_1 == original && closest_vector_2 == original 
        rates(i,1) = 1;
    elseif closest_vector_1 == original || closest_vector_2 == original 
        rates(i,1) = 2/3;
    else
         rates(i,1) = 1/3;
    end
    rates(i,2) = original;

end

disp("With Euclidean distance the error rates for the indecies are:")
disp("Success Rate: "+rates(:,1) +"  |   Index: "+num2str(rates(:,2)))




