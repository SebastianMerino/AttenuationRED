%% ========================= Processing loop ========================= %%
startup,

dataDir = "Q:\smerino\REDjournalResults\rf";

sampleName = "simuThyroid";
resultsDir = "Q:\smerino\REDjournalResults\rf\"+sampleName;
if ~exist("resultsDir","dir"); mkdir(resultsDir); end


load(fullfile(dataDir,sampleName+".mat"))
zRf = zRf';
xBm = xBm*100; zBm = zBm'*100;

big = false;
if big 
    sampleName = sampleName + "Big";
else
    sampleName = sampleName + "Small";
end
%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 3e6; freqH = 8.5e6; % wide bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.4; gammaRef = 1;  % 0.4 for simulations, 0.53 for in vivo
iAcq = 1;
groundTruthTargets = 1.19;

% Blocksize parameters
if big
    blockParams.xInf = xRf(1); % 0.8
    blockParams.xSup = xRf(end);
    blockParams.zInf = zRf(1);
    blockParams.zSup = zRf(end);
else
    blockParams.xInf = -1.1; % 0.8
    blockParams.xSup = 1.7;
    blockParams.zInf = 0.9;
    blockParams.zSup = 2.2;
end
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Measurement ROI
c1x = 0.3; c1z = 1.5;
roiL = 2; roiLz = 0.6;

% Plotting constants
dynRange = [-60,0];
attRange = [0.2,1.6];
bsRange = [-10,10];
yLimits = [zBm(1),zBm(end)];

NptodB = log10(exp(1))*20;
%%
% get Spectra
[Sp,Sd,xAcs,zAcs,f] = getSpectrum(rf,xRf,zRf,fs,blockParams);

% Plotting spectraxRf = xRf*100; zRf = zRf*100;

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
clear att_ref_map 
att_ref_map(1,1,:) = alpha0Ref*f/NptodB;

[SpRef,SdRef,~,~,~] = getSpectrum(ref,xRf,zRf,fs,blockParams);

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

%% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);

[~,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,1,roiLz);

%% For looping
muVec = 10.^(0:0.5:10);
iMu = 8;
%%
for iMu = 1:length(muVec)
%% RSLD-TV
muRsld = muVec(iMu);
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

AttInterp = interp2(X,Z,BR,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.method = 'RSLD';
r.mu = muRsld;
Metrics(iMu) = r;

%% RED no weigths
muRed = muVec(iMu);
tic
% [~ ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,7,m,n,muRed);
[~ ,~,u2] = admm_red_median(A'*A,A'*b(:),muRed,0.001,size(A'*b(:),1),1500,4,1,7,m,n,muRed);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.method = 'RED-MED';
r.mu = muRed;
Metrics(length(muVec)+iMu) = r;

%%
figure('Units','centimeters', 'Position',[5 5 24 8]);
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

t2 = nexttile;
myOverlay(t2, bMode,dynRange,xBm,zBm, BR,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
colormap(t2,turbo)
axis image
title("RSLD, \mu=10^{"+log10(muRsld)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlay(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
axis image
title("RED-MED, \mu=10^{"+log10(muRed)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)

pause(0.1)
saveas(gcf,fullfile(resultsDir,sampleName+"_mu"+iMu+".png"))
close,
end

%% Metrics plots
lw = 1.5;

T = struct2table(Metrics);
T.method = categorical(T.method);
muVec = T(T.method=='RSLD',:).mu;
tabRsld = T(T.method=='RSLD',:);
tabRed = T(T.method=='RED-MED',:);

figure,
hold on
plot(log10(muVec),tabRsld.maeInc, 'LineWidth',lw)
plot(log10(muVec),tabRed.maeInc, 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('MAE [dB/cm/MHz]')
grid on
legend('RSLD','RED')
title('MAE')
ylim([0 0.9])

[~,iMu] = min(tabRsld.maeInc);
optimMuRsld = muVec(iMu);
[~,iMu] = min(tabRed.maeInc);
optimMuRed = muVec(iMu);

colors = lines(8);
figure,
hold on
errorbar(log10(muVec),tabRsld.meanInc,tabRsld.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(1,:))
errorbar(log10(muVec),tabRed.meanInc,tabRed.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(5,:))
yline(groundTruthTargets(iAcq), '--', 'Color',colors(1,:))
yline(groundTruthTargets(end), '--', 'Color',colors(5,:))
% xline(log10(optimMuRsld), 'Color','b', 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('RSLD','RED')
ylim([0 2])



save_all_figures_to_directory(resultsDir,sampleName+"_metrics",'svg');
pause(0.1)
close all,


%% Optimal mu plot
muRsld = muVec(iMu);
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),optimMuRsld,optimMuRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));


tic
% [err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),optimMuRed,tol,2*m*n,200,7,m,n,optimMuRed);
[~ ,~,u2] = admm_red_median(A'*A,A'*b(:),optimMuRed,0.001,size(A'*b(:),1),1500,4,1,7,m,n,optimMuRed);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);

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
ylim(yLimits)

t2 = nexttile;
myOverlayInterp(t2, bMode,dynRange,xBm,zBm, BR,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
colormap(t2,turbo)
axis image
title("RSLD, \mu=10^{"+log10(optimMuRsld)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
colormap(t3,turbo)
axis image
title("RED, \mu=10^{"+log10(optimMuRed)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)


save_all_figures_to_directory(resultsDir,sampleName+"_final",'svg');
pause(0.1)
close all,

writetable(T,fullfile(resultsDir,sampleName+".xlsx"))