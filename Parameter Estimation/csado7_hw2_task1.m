%% Task 1:
function p_hat = csado7_hw2_task1(X,y)
    p_hat = pinv(X' * X) * (X' * y); 
end