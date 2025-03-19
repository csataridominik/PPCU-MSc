% Task 1:
close all; clc; clear;


T = readtable("pharmacokinetics.xlsx");
x1 = [0.27; T.Var1];
x2 = [0.2; T.Var2];
x3 = [0.2; T.Var3];
x4 = [0.2; T.Var4];

Y = zeros((size(x1,1)-1)*2,1);
X = zeros(size(Y,1),10);

it = 1;
it2 = 2;
for k=2:size(T.Var1,1)+1

    Y(it) = x1(k)+(x1(k) + x1(k-1))*x1(k);
    

    X(it,1:5) = [-x1(k) + x1(k-1), x1(k)*x2(k)*0.1 - x1(k)*x1(k)*0.1, ...
                (-x1(k) + x1(k-1))*x3(k),x2(k)*0.1 - x1(k)*0.1, ...
                x3(k)*x2(k)*0.1 - x3(k)*x1(k)*0.1];
    it = it+2;

    Y(it2) = x3(k)+(x3(k) + x3(k-1))*x3(k);
    
    
    X(it2,6:10) = [-x3(k) + x3(k-1), x4(k)*x3(k)*0.1 - x3(k)*x3(k)*0.1, ...
                (-x3(k) + x3(k-1))*x1(k),x4(k)*0.1 - x3(k)*0.1, ...
                x1(k)*x4(k)*0.1 - x3(k)*x1(k)*0.1];
    it2 = it2+2;

end

D = 4.5;
C = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0];


temp = inv([2*X'*X, C'; C, 0]);
Theta_cons = temp*[2*X'*Y ;D];

kc = Theta_cons(1);
ka = Theta_cons(6);


% In this part I calculate alpha2 and beta2.
% These are jsut placehodler vectors for the results of eq. 2 and 3.
X_ = zeros(size(x2));
X__ = zeros(size(x4));

% I calculate in other X's alpha1 and beta1 again hence _a and _b:
X_a = zeros(size(x2));
X_b = zeros(size(x4));

Y_a = zeros(size(x2));
Y_b = zeros(size(x2));

for i = 1:size(x2)
    X_(i) = x1(i)-x2(i); 
    X__(i) = x3(i)-x4(i); 
    
    X_a(i) = x2(i) - x1(i);
    X_b(i) = x4(i) - x3(i);

    Y_a(i) = x1(i) + (ka*10*x1(i))/(kc*ka+kc*x3(i)+ka*x1(i)); 
    Y_b(i) = x3(i) + (kc*10*x3(i))/(kc*ka+kc*x3(i)+ka*x1(i)); 
end

alpha1 = pinv(X_a' * X_a) * (X_a' * Y_a);
beta1 = pinv(X_b' * X_b) * (X_b' * Y_b);

alpha2 = pinv(X_' * X_) * (X_' * x2);
beta2 = pinv(X__' * X__) * (X__' * x4);




disp("The value of kc is: ");
disp(kc);
disp("The value of ka is: ");
disp(ka);
disp("The value of alpha1 is: ");
disp(alpha1);
disp("The value of beta1 is: ");
disp(beta1);
disp("The value of alpha2 is: ");
disp(alpha2);
disp("The value of beta2 is: ");
disp(beta2);

% In this part I simulate x1-x4 and calculate squared loss function.

Sim_X = zeros(size(x1,1),4);

sq_loss = zeros(4,1); 


for j = 1:size(Sim_X,1)
    
    Sim_X(j,1) = alpha1*(x2(j)-x1(j))-(ka*10*x1(j))/(kc*ka+kc*x3(j)+ka*x1(j));
    Sim_X(j,2) = alpha2*(x1(j)-x2(j));

    Sim_X(j,3) = beta1*(x4(j)-x3(j))-(kc*10*x3(j))/(kc*ka+kc*x3(j)+ka*x1(j));
    Sim_X(j,4) = beta2*(x3(j)-x4(j));


    % Here I calculate the losses:
    sq_loss(1) = sq_loss(1) + (Sim_X(j,1) - x1(j))^2;
    sq_loss(2) = sq_loss(2) + (Sim_X(j,2) - x1(j))^2;
    sq_loss(3) = sq_loss(3) + (Sim_X(j,3) - x1(j))^2;
    sq_loss(4) = sq_loss(4) + (Sim_X(j,4) - x1(j))^2;

    
end



figure(1);
grid on;
hold on;
plot(linspace(0,51,size(x1,1)),x1,'Color','red','LineWidth',1);
plot(linspace(0,51,size(x1,1)),x2,'Color','blue','LineWidth',1);
plot(linspace(0,51,size(x1,1)),x3,'Color','green','LineWidth',1);
plot(linspace(0,51,size(x1,1)),x4,'Color','magenta','LineWidth',1);

plot(linspace(0,51,size(x1,1)),Sim_X(:,1),'Color','red','LineStyle',':','LineWidth',2);
plot(linspace(0,51,size(x1,1)),Sim_X(:,2),'Color','blue','LineStyle',':','LineWidth',2);
plot(linspace(0,51,size(x1,1)),Sim_X(:,3),'Color','green','LineStyle',':','LineWidth',2);
plot(linspace(0,51,size(x1,1)),Sim_X(:,4),'Color','magenta','LineStyle',':','LineWidth',2);
hold off;
legend('x1(t)','x2(t)','x3(t)','x4(t)','x1(t)\_est','x2(t)\_est','x3(t)\_est','x4(t)\_est')

% The loss funcition values are:
disp("squared losses are in order(x1, x2, x3, x4): ")
disp(sq_loss)

% Squared Losses counting only 10:end elements are:
sq_loss = zeros(4,1); 
for j = 10:size(Sim_X,1)
    
    sq_loss(1) = sq_loss(1) + (Sim_X(j,1) - x1(j))^2;
    sq_loss(2) = sq_loss(2) + (Sim_X(j,2) - x1(j))^2;
    sq_loss(3) = sq_loss(3) + (Sim_X(j,3) - x1(j))^2;
    sq_loss(4) = sq_loss(4) + (Sim_X(j,4) - x1(j))^2;
    
end

disp("Squared losses are in order(x1, x2, x3, x4) when k h goes from 10: ")
disp(sq_loss)




