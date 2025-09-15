%% ========================= Processing loop ========================= %%
% Script to test some denoisers on the phantom data. Optimization algorithm
% was chosen from script testingAlgorithms.m
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

%%
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

save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_spec");
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

%% ================== median denoiser ==========================
% Initialization
[x0,~] = cgs(A'*A,A'*b(:));

denoiser = medianDenoiserHandle(m,n,[7 3]);
optimMuRed = 10^4; tol = 1e-4; L = 1;
tic
[xDiff, fpError, xEst] = redPhaseGradient(A,b(:),optimMuRed,L,denoiser,tol,1500,x0);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(optimMuRed)+"}")
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
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(optimMuRed)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
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


%% ================== Mixed denoiser ==========================
% Initialization
[x0,~] = cgs(A'*A,A'*b(:));
optimMuRed = 10^3.5; tol = 1e-4; L = 1.01;
denoiser = mixedDenoiserHandle(m,n,7,1/optimMuRed);
tic
[xDiff, fpError, xEst] = redPhaseGradient(A,b(:),optimMuRed,L,denoiser,tol,1500,x0);
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
% hold on
% contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
% hold off

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
colormap(t3,turbo)
axis image
title("\mu=10^{"+log10(optimMuRed)+"}")
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
colormap(t3,"parula")
axis image
title("\mu=10^{"+log10(optimMuRed)+"}")
c = colorbar;
c.Label.String = 'AZC [dB/cm]';
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

%% ================== NEW denoiser ==========================
[x0,~] = cgs(A'*A,A'*b(:));
x0 = reshape(x0,m,n,2);
x0(:,:,1) = mean(x0(:,:,1),'all');
x0(:,:,2) = mean(x0(:,:,2),'all');
x0 = x0(:);

lambda = 1e5;
% denoiser = bm3dDenoiserHandle(m,n,cat(3,1/lambda,1/lambda));
denoiser = dncnnDenoiserHandle(m,n);
tol = 1e-4; L = 1;
tic
[xDiff, fpError, xEst] = redPhaseGradient(A,b(:),lambda,L,denoiser,tol,1500,x0);
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
title("\mu=10^{"+log10(optimMuRed)+"}")
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
title("\mu=10^{"+log10(optimMuRed)+"}")
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

%%
function f = bm3dDenoiserHandle(m,n,noiseStd)
    function y = applyDenoiser(x)
        x = reshape(x,[m,n,2]);
        y = BM3D(x, noiseStd);
        y = double(y(:));
    end
    f = @applyDenoiser;
end


function f = dncnnDenoiserHandle(m,n)
NptodB = db(exp(1));
net = denoisingNetwork("DnCNN");
    function y = applyDenoiser(x)
        x = reshape(x,[m,n,2])*NptodB;
        acs = x(:,:,1); 
        min1 = min(acs(:)); range = max(acs(:))-min1;
        % min1 = 0; range = 2;
        acs = (acs-min1)/range;
        acs = denoiseImage(acs,net)*range+min1;
        bsc = x(:,:,2);
        min1 = min(bsc(:)); range = max(bsc(:))-min1;
        % min1 = -15; range = 30;
        bsc = (bsc-min1)/range;
        bsc = denoiseImage(bsc,net)*range+min1;
        y = cat(3,acs,bsc)/NptodB;
        y = y(:);
    end
    f = @applyDenoiser;
end


function [xDiff, fpError, xk] = redPhaseGradient(A,y,lambda,L,denoiser,tol,maxIter, x0)
xDiff = nan(maxIter,1);
fpError = nan(maxIter,1);

AtA = A'*A; Aty = A'*y;
N = length(x0);
xk = x0;
fxk = denoiser(xk);
vk = fxk/L - (1-L)/L*xk;
xPrev = x0;

% figure, tiledlayout(1,2),
% t1 = nexttile;
% t2 = nexttile;
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

    % if true% mod(k,round(maxIter/20))==1
    %     img = reshape(xk,69,121,2)*db(exp(1));
    %     imagesc(t1,img(:,:,1), [0 2])
    %     colormap(t1,"turbo")
    % 
    %     imagesc(t2,img(:,:,2), [-20 20])
    %     colormap(t2,"parula")
    %     pause(1)
    % end
end
end



function f = mixedDenoiserHandle(m,n,kernel,lambda)
    function fx = mixedDenoiser(x)
    x = reshape(x,[m,n,2]);
    x1 = x(:,:,1);
    x2 = x(:,:,2);
    x1 = medfilt2(x1, [kernel(1),kernel(1)],'symmetric');
    x2 = sign(x2).*max(0,abs(x2)-lambda); % L1-norm denoising
    % x2 = x2/(1+lambda); % L2-norm denoising
    fx= cat(3,x1,x2);
    fx = fx(:);
    end
    f = @mixedDenoiser;
end