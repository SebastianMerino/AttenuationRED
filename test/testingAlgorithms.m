%% ========================= Processing loop ========================= %%
% Script to test optimzation algorithms on the layered phantom data. 

startup,

dataDir = 'Q:\NonClinicalAvendanoData\May24_Focused_Phantoms_ComposedPhantoms\SavedData2405-SMA\bf';
refDir = 'Q:\NonClinicalAvendanoData\May24_Focused_Phantoms_ComposedPhantoms\SavedData2405-SMA\ref';
sampleFiles = dir(fullfile(dataDir,'*.mat'));
refFiles = dir(fullfile(refDir,'*Ref*.mat'));

resultsDir = 'Q:\smerino\REDjournalResults\layeredPhantom\test';
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

big = true;

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 4e6; freqH = 9e6; % wide bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.2; % [cm]

% Blocksize parameters
if big
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 1.5;
    blockParams.zSup = 4;
else
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 2.2;
    blockParams.zSup = 4;
end
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Measurement ROI
c1x = 0; c1z = 3;
roiL = 2.5; roiLz = 1.2;
groundTruth = 0.53;

% Plotting constants
dynRange = [-60,0];
attRange = [0,1];
bsRange = [-10,10];
yLimits = [deadband,4];

NptodB = log10(exp(1))*20;

iAcq = 4;
sampleName = sampleFiles(iAcq).name(1:end-4);
if big 
    sampleName = sampleName + "_Big";
else
    sampleName = sampleName + "_Small";
end
disp(sampleName)
out = matfile(fullfile(dataDir,sampleFiles(iAcq).name));
xBm = out.x*1e2; % [cm]
zBm = out.z'*1e2; % [cm]
sam1 = out.rf(:,:,1);
fs = out.fs;

% Plot region of interest B-mode image
bMode = db(hilbert(sam1));
bMode = bMode - max(bMode(zBm>deadband,:),[],"all");

% get Spectra
[Sp,Sd,xAcs,zAcs,f] = getSpectrum(sam1,xBm,zBm,fs,blockParams);

% Plotting spectra
spectrumSamzf = db(squeeze(mean(Sp/2+Sd/2, 2)));
spectrumSamzf = spectrumSamzf - max(spectrumSamzf,[],2);
figure,
imagesc(f,zAcs,spectrumSamzf, [-50 0]),
ax = gca; ax.YDir = 'reverse';
hold on
xline(freqL/1e6, 'w--', 'LineWidth',2), xline(freqH/1e6, 'w--', 'LineWidth',2)
hold off
colorbar
xlim([0 12])
xlabel('f [MHz]')
ylabel('z [cm]')
title('Sample power spectrum by depth')

%% Generating Diffraction compensation
% Generating references
clear att_ref_map rfRef

for ff = 1:length(refFiles)
    out = load(fullfile(refDir,refFiles(ff).name));
    rfRef(:,:,ff) = out.rf(:,:,end);
end
att_ref_map(1,1,:) = alpha0Ref*f.^gammaRef/NptodB;

[SpRef,SdRef,~,~,~] = getSpectrum(rfRef,xBm,zBm,fs,blockParams);

% Plotting spectra
spectrumRefzf = db(squeeze(mean(SpRef/2+SdRef/2, 2)));
spectrumRefzf = spectrumRefzf - max(spectrumRefzf,[],2);
figure,
imagesc(f,zAcs,spectrumRefzf, [-50 0]),
ax = gca; ax.YDir = 'reverse';
hold on
xline(freqL/1e6, 'w--', 'LineWidth',2), xline(freqH/1e6, 'w--', 'LineWidth',2)
hold off
colorbar
xlim([0 12])
xlabel('f [MHz]')
ylabel('z [cm]')
title('Reference power spectrum by depth')
%%
save_all_figures_to_directory(resultsDir,sampleName+"_spec");
pause(0.1)
close all,

%% Setting up system
L = (zAcs(2) - zAcs(1))/(1 - blockParams.overlap)/2;   % (cm)
sld = log(Sp) - log(Sd);
sldRef = log(SpRef) - log(SdRef);
compensation = sldRef - 4*L*att_ref_map;

range = f>freqL/1e6 & f<freqH/1e6;
b = sld(:,:,range) - compensation(:,:,range);
ufr = f(range);

[m,n,p] = size(b);
A1 = kron( 4*L*ufr , speye(m*n) );
A2 = kron( ones(size(ufr)) , speye(m*n) );
A = [A1 A2];
mask = ones(m,n,p);
tol = 1e-3;

[~,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,1,roiLz);


%% Optimal mu plot
denoiser = medianDenoiserHandle(m,n,[7 7]);
mu = 10^4;
tol = 1e-5;
tic
[res,loss,u1] = admm_red_median(A'*A,A'*b(:),mu,tol,size(A'*b(:),1),1500,4,1,7,m,n,mu);
toc,
BRED = reshape(u1(1:end/2)*NptodB,m,n);
CRED = reshape(u1(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 6]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
title('B-mode')
colormap(t1,gray)
c = colorbar;
c.Label.String = 'dB';
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)

% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
    "omitnan");
disp(r)
disp(['Cost: ',num2str(norm(A*u1-b(:))^2/2 + mu/2*u1'*(u1 - denoiser(u1)))])

% Convergence
figure('Position',[100 100 600 800]),
tiledlayout(3,1)
nexttile,
semilogy(res(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('||\Delta x||/||x||')
grid on

nexttile,
semilogy(loss(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('Objective function')
grid on

nexttile,
semilogy(abs(diff(loss)), 'LineWidth',2)
xlabel('Iterations')
ylabel('\Delta Objective function')
grid on


%% ============================== MY METHOD =============================
mu = 10^4; tol = 1e-4;
tic
[res,loss,u2] = admmRedMedianv2(A,b(:),mu,tol,2*m*n,1500,7,m,n,mu);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 6]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
title('B-mode')
colormap(t1,gray)
c = colorbar;
c.Label.String = 'dB';
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)


% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
    "omitnan");
disp(r)
disp(['Cost: ',num2str(norm(A*u2-b(:))^2/2 + mu/2*u2'*(u2 - denoiser(u2)))])

% Convergence
figure('Position',[100 100 600 800]),
tiledlayout(3,1)
nexttile,
semilogy(res(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('||\Delta x||')
grid on

nexttile,
semilogy(loss(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('Objective function')
grid on

nexttile,
semilogy(abs(diff(loss)), 'LineWidth',2)
xlabel('Iterations')
ylabel('\Delta Objective function')
grid on

%% ================== NEW OPTIMIZER ==========================
denoiser = medianDenoiserHandle(m,n,[7 7]);
mu = 10^4; tol = 1e-4; L = 1.01;
tic
[xDiff, fpError, xEst] = redPhaseGradient(A,b(:),mu,L,denoiser,tol,1500);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 6]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
title('B-mode')
colormap(t1,gray)
c = colorbar;
c.Label.String = 'dB';
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)


% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
    "omitnan");
disp(r)
disp(['Cost: ',num2str(norm(A*xEst-b(:))^2/2 + mu/2*xEst'*(xEst - denoiser(xEst)))])

% Convergence
figure('Position',[100 100 600 800]),
tiledlayout(2,1)
nexttile,
semilogy(xDiff(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('||\Delta x||')
grid on

nexttile,
semilogy(fpError(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('FP error')
grid on

%% ================== NEW FASTER OPTIMIZER ==========================
denoiser = medianDenoiserHandle(m,n,[7 7]);
mu = 10^7; tol = 1e-4; L = 1.5;
tic
[xDiff, fpError, xEst] = redAcceleratedPG(A,b(:),mu,L,denoiser,tol,1500);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 24 6]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
title('B-mode')
colormap(t1,gray)
c = colorbar;
c.Label.String = 'dB';
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
ylim(yLimits)


% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
    "omitnan");
disp(r)
disp(['Cost: ',num2str(norm(A*xEst-b(:))^2/2 + mu/2*xEst'*(xEst - denoiser(xEst)))])

% Convergence
figure('Position',[100 100 600 800]),
tiledlayout(2,1)
nexttile,
semilogy(xDiff(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('||\Delta x||')
grid on

nexttile,
semilogy(fpError(2:end), 'LineWidth',2)
xlabel('Iterations')
ylabel('FP error')
grid on

%%
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



function [res_admm, loss, out] = admmRedMedianv2(A,y,lambda,tol,N,max_iter,median,m,ni,rho)

[x_est,~] = cgs(A'*A, A'*y, tol, 100);
v_est = x_est;
u_est=zeros([N 1]);

AtArho = A'*A+rho*speye(N);
Aty = A'*y;

xPrev = x_est;

for k = 1:1:max_iter

    % Part1 of the ADMM, approximates the solution of:
    % x = argmin_z 1/(2sigma^2)||Ax-y||_2^2 + 0.5*rho||x - v + u||_2^2
    [x_est,~] = cgs(AtArho, Aty + rho*(v_est-u_est), 1e-6, 5,[],[],x_est);
    
    % Part2 of the ADMM, approximates the solution of
    % v = argmin_z lambda*z'*(z-denoiser(z)) +  0.5*beta||z - x - u||_2^2
    % using gradient descent
    for j = 1:1
        v_reshap = reshape(v_est,m,2*ni);
        v_est1=v_reshap((1:m),(1:ni));
        v_est2=v_reshap((1:m),(ni+1:2*ni));
        f_v_est1 = medfilt2(v_est1, [7,7],'symmetric');
        f_v_est2 = medfilt2(v_est2, [3,3],'symmetric');

        f_v_est = [f_v_est1 f_v_est2];
        f_v_est = f_v_est(:);

        v_est = (rho*(x_est + u_est) + lambda*f_v_est)/(lambda + rho);
    end

    % Part3 of the ADMM, update the dual variable
    u_est = u_est + x_est - v_est;

    r(k) = norm(xPrev - x_est)/norm(xPrev);
    loss(k) = norm(A*x_est-y)^2/2 + lambda/2*x_est'*(x_est - f_v_est);
    if norm(xPrev - x_est)/norm(xPrev) < tol && k>2
        break
    end
    xPrev = x_est;

end
out = x_est;
res_admm = r;

end


function [xDiff, fpError, xk] = redPhaseGradient(A,y,lambda,L,denoiser,tol,maxIter)
xDiff = nan(maxIter,1);
fpError = nan(maxIter,1);

AtA = A'*A; Aty = A'*y;
[xk,~] = cgs(AtA,Aty);
N = length(xk);
fxk = denoiser(xk);
vk = fxk/L - (1-L)/L*xk;
xPrev = xk;
for k = 1:maxIter
    % First step: argmin 1/2 || Ax - y ||^2 + lambda*L/2 || x - v ||^2
    [xk,~] = cgs(AtA + lambda*L*speye(N), Aty+lambda*L*vk);

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

function [xDiff, fpError, xk] = redAcceleratedPG(A,y,lambda,L,denoiser,tol,maxIter)
xDiff = nan(maxIter,1);
fpError = nan(maxIter,1);

AtA = A'*A; Aty = A'*y;
[xk,~] = cgs(AtA,Aty);
N = length(xk);
fxk = denoiser(xk);
vk = fxk/L - (1-L)/L*xk;
xPrev = xk;
tPrev = 1;
for k = 1:maxIter
    % First step: argmin 1/2 || Ax - y ||^2 + lambda*L/2 || x - v ||^2
    [xk,~] = cgs(AtA + lambda*L*speye(N), Aty+lambda*L*vk);

    % Acceleration
    tk = (1 + sqrt(1 + 4*tPrev^2))/2;
    zk = xk + (tPrev-1)/tk*(xk - xPrev);
    
    % Denoiser step
    vk = denoiser(zk)/L - (1-L)/L*zk;

    xDiff(k) = norm(xPrev - xk)/norm(xk);
    fpError(k) = norm( AtA*xk - Aty + lambda*(xk - denoiser(xk)) )/N;
    if xDiff(k) < tol && k>2
        break
    end
    xPrev = xk;
    tPrev = tk;
end
end


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