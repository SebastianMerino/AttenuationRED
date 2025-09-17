function [xDiff, fpError, xk] = redPhaseGradient(A,y,lambda,L,denoiser,tol,maxIter, x0)
xDiff = nan(maxIter,1);
fpError = nan(maxIter,1);

AtA = A'*A; Aty = A'*y;
N = length(x0);
xk = x0;
fxk = denoiser(xk);
vk = fxk/L - (1-L)/L*xk;
xPrev = x0;

for k = 1:maxIter
    % First step: argmin 1/2 || Ax - y ||^2 + lambda*L/2 || x - v ||^2
    [xk,~] = cgs(AtA + lambda*L*speye(N), Aty+lambda*L*vk, [], [], [] ,[] ,xk);

    % Denoiser step
    fxk = denoiser(xk);
    vk = fxk/L - (1-L)/L*xk;

    xDiff(k) = norm(xPrev - xk)/norm(xk);
    fpError(k) = norm( AtA*xk - Aty + lambda*(xk - fxk) )/N;
    if xDiff(k) < tol && k>2
        break
    end
    xPrev = xk;
end
end