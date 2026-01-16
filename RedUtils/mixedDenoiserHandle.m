function f = mixedDenoiserHandle(m,n,kernel,lambda)
    function fx = mixedDenoiser(x)
    x = reshape(x,[m,n,2]);
    x1 = x(:,:,1);
    x2 = x(:,:,2);
    x1 = medfilt2(x1, [kernel(1),kernel(1)],'symmetric');
    x2 = sign(x2).*max(0,abs(x2)-lambda); % L1-norm denoising
    fx= cat(3,x1,x2);
    fx = fx(:);
    end
    f = @mixedDenoiser;
end