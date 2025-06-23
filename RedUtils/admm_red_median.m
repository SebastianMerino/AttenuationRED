function [r, loss , out ] = admm_red_median(A,y,mu,error,N,max_iter,inner_iters,inner_iters2,median,m,ni,rho)

x_est= ones([N 1]);
v_est=ones([N 1]);
u_est=ones([N 1]);

xPrev = x_est;
alpha = 1;
for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*beta||x - v + u||_2^2
    for j = 1:1:inner_iters
        b = (A'*y) + rho*(v_est - u_est);
        grad_e = b - (A'*A*x_est + rho*x_est);
        if norm(grad_e) == 0
            break;
        end
        res = (A'*A*grad_e) + rho*grad_e;
        beta_opt = mean(grad_e(:).*grad_e(:))/(mean(grad_e(:).*res(:)));
        x_est = x_est + beta_opt*grad_e;
    end

    % relaxation
    x_hat = alpha*x_est + (1-alpha)*v_est;

    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1:inner_iters2
        v_est = reshape(x_est,m,2*ni);
        v_est1=v_est((1:m),(1:ni));
        v_est2=v_est((1:m),(ni+1:2*ni));
        f_v_est1 = medfilt2(v_est1, [median median],'symmetric');
        f_v_est2 = medfilt2(v_est2, [median median],'symmetric');

        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);

        v_est = (rho*(x_hat + u_est) + mu*f_v_est)/(mu + rho);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_hat - v_est;

    % Stop criteria
    r(k) = norm(xPrev - x_est)/norm(xPrev);
    loss(k) = norm(A*x_hat-y)^2/2 + mu/2*x_hat'*(x_hat - f_v_est);
    if r(k) < error && k>2
        break
    end
    xPrev = x_est;
end
out = x_hat;

end
