% Exercise 1:
close all; clear; clc;

N = 40 * 10;
W = 92;
H = 112;


% Create Gamma:
Gamma = zeros(W*H, N-1);

% Store the images with dir command:
folder = dir("Cambridge_FaceDB_jpg\Cambridge_FaceDB_jpg\*.jpg");

% Select a random test image and delete from Gamma:
% I did it so I choose the image beforehand and I do not ad it to Gamma.
p = randperm(N,1);



for i = 1:size(folder, 1)
    
    temp = imread("Cambridge_FaceDB_jpg\Cambridge_FaceDB_jpg\"+folder(i).name);

    if size(temp,3) == 3
        temp = rgb2gray(temp);
    end
    
    temp = temp(:);

    if i == p
        test_img = temp;
    else    
        Gamma(:,i) = temp;
    end
end


% Calculate the average image:
Psi = mean( Gamma , 2 );

% Create matrix A as the difference between 
% the individual images and the mean image:

A = Gamma - Psi;

% PCA-base creation:

L = A' * A;
[V, D] = eig(L, 'vector');
[D, idx] = sort(D, 'descend');
V = V(:, idx);
Uraw = A*V;

n = vecnorm(Uraw);

for i = 1: size(n,2)
    Uraw(:,i) = Uraw(:,i)/n(i);
end


% Here I can set different values for K:
K = 320;
% K = 50;

U = Uraw(:,1:K);

% Calculating the PCA-projection:
Y = U' * A;

% For the test image as well:
Phi = double(test_img) - Psi;
y_test = U' * Phi;

% Select the closest value:
distance = abs(Y - y_test);
difference = vecnorm(distance);

% Find the smallest element of this vector:
[closest,closest_ind] = min(difference);


% Reconstruct the found element:
Phi_rec = U*Y(:,closest_ind);
recognised_test_img = Phi_rec + Psi;

%
recognised_test_img = reshape(recognised_test_img,[H,W]);
test_img = reshape(test_img,[H,W]);
Psi = reshape(Psi,[H,W]);

figure(1);
subplot(1,3,1);
imagesc(recognised_test_img);
colormap('gray')
title('Test Image')


subplot(1,3,2);
imagesc(Psi);
colormap('gray')
title('Mean Face')

subplot(1,3,3 );
imagesc(test_img);
colormap('gray')
title('Reconstructed from the closest PCA projection, K='+string(K))

% 12. Plot the first 10 eigenfaces:
%%
figure(2);
for j = 1:10
    curr = reshape(U(:,j),[H,W]);

    subplot(2,5,j);
    imagesc(curr);
    colormap('gray')
    
end







