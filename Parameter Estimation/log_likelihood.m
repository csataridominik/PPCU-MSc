function loglikelihood = log_likelihood(u,y,a,b)
    loglikelihood = 0;
    sigma1 = sqrt(1.2);
    sigma2 = sqrt(3);


    for k = 2:size(y,1)
        if k < 500
            loglikelihood = loglikelihood + (1/(sigma1*sqrt(2*pi))) + ...
            log(exp(((y(k)-a*y(k-1)-b*(u(k-1)))^2/(2*sigma1^2))));
        else
            loglikelihood = loglikelihood + (1/(sigma2*sqrt(2*pi))) + ...
            log(exp(((y(k)-a*y(k-1)-b*(u(k-1)))^2/(2*sigma2^2))));
        end


    end
end
