%% ========================= Processing loop ========================= %%
startup,

baseDir = 'Q:\smerino\phantoms\myAcquisition\processed';
dataDir = fullfile(baseDir,'bf');
sampleFiles = dir(fullfile(dataDir,'T7*.mat'));

refDir = fullfile(baseDir,'ref');
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = 'Q:\smerino\REDjournalResults\newPhantom\test';
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 3.5e6; freqH = 8e6; % wide bandwidth
% freqL = 4e6; freqH = 7.5e6; % narrow bandwidth
% freqL = 3.5e6; freqH = 7.5e6; % homog bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.25; % [cm]

% Blocksize parameters
blockParams.xInf = -2;
blockParams.xSup = 2;
blockParams.zInf = 0.25;
blockParams.zSup = 5;
blockParams.blocksize = [8 12]*wl;
blockParams.overlap = 0.8;

% Plotting constants
dynRange = [-60,0];
attRange = [0.3,1.2];
bsRange = [-10,10];
yLimits = [deadband,5];

NptodB = log10(exp(1))*20;

iAcq = 1;

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
% %%
% save_all_figures_to_directory(resultsDir,sampleName+"_spec");
% pause(0.1)
% close all,

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

%% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
if iAcq == 1
    c1x = 0.1;
else
    c1x = 0;
end
c1z = 2.1;
c2z = 4;
rInc = 0.95;
roiL = 1; roiD = 0.6;
roiLz = 1;

inc = (Xq-c1x).^2 + (Zq-c1z).^2 < (rInc-0.2).^2;
back = (Xq-c1x).^2 + (Zq-c1z).^2 > (rInc+0.2).^2;
groundTruthTargets = [0.97,0.95,0.95,0.55];


%% ============================== NEW OPT =============================
denoiser = medianDenoiserHandle(m,n,[7 3]);
mu = 10^4; tol = 1e-4; L = 1.01; nIte=2000;
tic
[xDiff, fpError, xEst] = redPhaseGradient(A,b(:),mu,L,denoiser,tol,nIte);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 18 6]);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,parula)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = '\Delta BSC [dB/cm]';
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
disp(r)

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

%% ===================== Weight map from BSC changes =====================
% Initial iteration
[~, ~, xEst] = redPhaseGradient(A,b(:),mu,L,denoiser,tol,nIte);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

wBSC = 0.9*(abs(CRED) < 5)+0.1;

figure,imagesc(xAcs,zAcs, wBSC, [0 1])
axis image
colorbar
title('Weight map')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

% Weight matrices and new system
W = repmat(wBSC,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;
Aw = [A1w,A2w];

%% Weighted fidelity with BSC changes
denoiser = medianDenoiserHandle(m,n,[7 3]);
mu = 10^4; tol = 1e-4; L = 1.01; nIte=2000;
tic
[xDiff, fpError, xEst] = redPhaseGradient(Aw,bw,mu,L,denoiser,tol,nIte);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 18 6]);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,parula)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = '\Delta BSC [dB/cm]';
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
disp(r)

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


%% ===================== Weight map from residuals =====================
% Initial iteration
[~, ~, xEst] = redPhaseGradient(A,b(:),mu,L,denoiser,tol,nIte);

residuals = b(:)-A*xEst;
residuals = reshape(residuals,m,n,p);
wRes = 1/(1+(residuals/1).^2);
figure,imagesc(xAcs,zAcs, median(wRes,3))
axis image
colorbar
title('Weight map')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

% Weight matrices and new system
W = spdiags(wRes(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;
Aw = [A1w,A2w];

%% Weighted fidelity
denoiser = medianDenoiserHandle(m,n,[7 3]);
mu = 10^4; tol = 1e-4; L = 1.01; nIte=2000;
tic
[xDiff, fpError, xEst] = redPhaseGradient(Aw,bw,mu,L,denoiser,tol,nIte);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 18 6]);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,parula)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = '\Delta BSC [dB/cm]';
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
disp(r)

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

%% ===================== Weight map from SNR =====================
SNR = getSnr(sam1,xBm,zBm,blockParams);

% Calculating weights
desvMin = 30;
SNRopt = sqrt(1/(4/pi - 1));
desvSNR = abs(SNR-SNRopt)/SNRopt*100;
% wSNR = 1./(1 + (desvSNR/desvMin).^2);
wSNR = 0.9*(desvSNR < desvMin)+0.1;

figure,imagesc(xAcs,zAcs, wSNR, [0 1])
axis image
colorbar
title('Weight map')
ylabel('Axial [cm]')
xlabel('Lateral [cm]')

% Weight matrices and new system
W = repmat(wSNR,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;
Aw = [A1w,A2w];

%% Weighted fidelity
denoiser = medianDenoiserHandle(m,n,[7 3]);
mu = 10^4; tol = 1e-4; L = 1.01; nIte=2000;
tic
[xDiff, fpError, xEst] = redPhaseGradient(Aw,bw,mu,L,denoiser,tol,nIte);
toc,
BRED = reshape(xEst(1:end/2)*NptodB,m,n);
CRED = reshape(xEst(end/2+1:end)*NptodB,m,n);

figure('Units','centimeters', 'Position',[5 5 18 6]);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, CRED,bsRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,parula)
axis image
title("\mu=10^{"+log10(mu)+"}")
c = colorbar;
c.Label.String = '\Delta BSC [dB/cm]';
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off
ylim(yLimits)

% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);
AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
disp(r)

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

%% Function
function [xDiff, fpError, xk] = redPhaseGradient(A,y,lambda,L,denoiser,tol,maxIter)
xDiff = nan(maxIter,1);
fpError = nan(maxIter,1);

AtA = A'*A; Aty = A'*y;
[xk,~] = cgs(AtA,Aty);
N = length(xk);
fxk = denoiser(xk);
vk = fxk/L - (1-L)/L*xk;
% vk = zeros(N,1);
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