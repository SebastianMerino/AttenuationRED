function f = medianDenoiserHandle(m,n,kernel)
    function fx = medianDenoiser(x)
    x = reshape(x,[m,n,2]);
    x1 = x(:,:,1);
    x2 = x(:,:,2);
    x1 = medfilt2(x1, [kernel(1),kernel(1)],'symmetric');
    x2 = medfilt2(x2, [kernel(2),kernel(2)],'symmetric');
    fx= cat(3,x1,x2);
    fx = fx(:);
    end
    f = @medianDenoiser;
end