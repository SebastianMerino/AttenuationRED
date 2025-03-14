% ======================================================================
% ======================================================================
%% Simulated data
% Test with inclusion, only half region of interest, 
% top and bottom metric regions
startup,

dataDir = 'P:\smerino\simulation_acs\rf_data\24_04_04_inc';
sampleFiles = dir(fullfile(dataDir,'*.mat'));
sampleFiles = sampleFiles(end-2:end);

refDir = 'P:\smerino\simulation_acs\rf_data\24_04_25_ref';
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = fullfile(dataDir,'results','REDtest');
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 3.5e6; freqH = 9e6;
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.6; gammaRef = 1;
deadband = 0.1; % [cm]

% Blocksize parameters
blockParams.xInf = -2; % 0.8
blockParams.xSup = 2;
blockParams.zInf = 0.1;
blockParams.zSup = 2.35;
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Plotting constants
dynRange = [-60,0];
attRange = [0.4,1.1];
bsRange = [-10,10];
yLimits = [0.1,10];

NptodB = log10(exp(1))*20;
iAcq = 2;

out = matfile(fullfile(dataDir,sampleFiles(iAcq).name));
xBm = out.x*1e2; % [cm]
zBm = out.z'*1e2; % [cm]
sam1 = out.rf(:,:,1);
fs = out.fs;

% Plot region of interest B-mode image
bMode = db(hilbert(sam1));
bMode = bMode - max(bMode(:));

% get Spectra
[Sp,Sd,x_ACS,z_ACS,f] = getSpectrum(sam1,xBm,zBm,fs,blockParams);

% Plotting spectra
spectrumSamzf = db(squeeze(mean(Sp/2+Sd/2, 2)));
spectrumSamzf = spectrumSamzf - max(spectrumSamzf,[],2);
figure,
imagesc(f,z_ACS,spectrumSamzf, [-50 0]),
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
    out = matfile(fullfile(refDir,refFiles(ff).name));
    rfRef(:,:,ff) = out.rf(:,:,end);
end
att_ref_map(1,1,:) = alpha0Ref*f.^gammaRef/NptodB;

[SpRef,SdRef,~,~,~] = getSpectrum(rfRef,xBm,zBm,fs,blockParams);

% Plotting spectra
spectrumRefzf = db(squeeze(mean(SpRef/2+SdRef/2, 2)));
spectrumRefzf = spectrumRefzf - max(spectrumRefzf,[],2);
figure,
imagesc(f,z_ACS,spectrumRefzf, [-50 0]),
ax = gca; ax.YDir = 'reverse';
hold on
xline(freqL/1e6, 'w--', 'LineWidth',2), xline(freqH/1e6, 'w--', 'LineWidth',2)
hold off
colorbar
xlim([0 12])
xlabel('f [MHz]')
ylabel('z [cm]')
title('Reference power spectrum by depth')

save_all_figures_to_directory(resultsDir,"target"+iAcq+"_spec");
pause(0.1)
close all,

%% Setting up system
L = (z_ACS(2) - z_ACS(1))/(1 - blockParams.overlap)/2;   % (cm)
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
groundTruthBack = [0.5,0.5,0.5];
groundTruthInc = [1,1,1];

[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(xBm,zBm);
c1x = 0; c1z = 2.05;
rInc = 0.7;
roiL = 0.9; roiD = 0.35;
roiLz = 0.9;

% inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
% back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;
[~,inc] = getRegionMasks(xBm,zBm,c1x,c1z-0.2,roiL,roiD,roiLz/2);
[~,back] = getRegionMasks(xBm,zBm,c1x,c1z-1.2,roiL,roiD,roiLz/2);
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;

%% For looping
muVec = 10.^(0.5:0.5:7);
%Metrics = cell(2*length(muVec),1);
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
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'RSLD';
r.mu = muRsld;
Metrics(iMu) = r;
%% RED no weigths
muRed = muVec(iMu);
tic
[err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,5,m,n,muRed/10);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

AttInterp = interp2(X,Z,BRED,Xq,Zq);
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.meanInc = mean(AttInterp(inc),"omitnan");
r.stdInc = std(AttInterp(inc),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthBack(iAcq),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthInc(iAcq),"omitnan");
r.rmseBack = sqrt(mean( (AttInterp(back) - groundTruthBack(iAcq)).^2,"omitnan"));
r.rmseInc = sqrt(mean( (AttInterp(inc) - groundTruthInc(iAcq)).^2,"omitnan"));
r.cnr = abs(r.meanInc - r.meanBack)/sqrt(r.stdBack^2 + r.stdInc^2);
r.method = 'RED-MED';
r.mu = muRed;
Metrics(length(muVec)+iMu) = r;

%%
figure('Units','centimeters', 'Position',[5 5 24 8]);
tl = tiledlayout(1,3, "Padding","tight");

t1 = nexttile;
imagesc(xBm,zBm,bMode,dynRange)
xlim([x_ACS(1) x_ACS(end)]),
ylim([z_ACS(1) z_ACS(end)]),
xlabel('Lateral [cm]'),
ylabel('Axial [cm]')
axis image
title('B-mode')
colormap(t1,gray)
c = colorbar;
c.Label.String = 'dB';


t2 = nexttile;
myOverlay(t2, bMode,dynRange,xBm,zBm, BR,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
colormap(t2,turbo)
axis image
title("RSLD, \mu=10^{"+log10(muRsld)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1)
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off

t3 = nexttile;
myOverlay(t3, bMode,dynRange,xBm,zBm, BRED,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
axis image
title("RED-MED, \mu=10^{"+log10(muRed)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1)
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off
%%
pause(0.1)
saveas(gcf,fullfile(resultsDir,"target"+iAcq+"_mu"+iMu+".png"))
close,
end

%%
colors = lines(2);
lw = 1.5;

T = struct2table(Metrics);
T.method = categorical(T.method);
muVec = T(T.method=='RSLD',:).mu;
tabRsld = T(T.method=='RSLD',:);
tabRed = T(T.method=='RED-MED',:);

figure,
hold on
plot(log10(muVec),tabRsld.rmseInc, 'LineWidth',lw)
plot(log10(muVec),tabRed.rmseInc, 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('RMSE [dB/cm/MHz]')
grid on
legend('RSLD','RED')
title('Inclusion')
ylim([0 0.9])

figure,
hold on
plot(log10(muVec),tabRsld.rmseBack, 'LineWidth',lw)
plot(log10(muVec),tabRed.rmseBack, 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('RMSE [dB/cm/MHz]')
grid on
legend('RSLD','RED')
title('Background')
ylim([0 0.9])


figure,
hold on
errorbar(log10(muVec),tabRsld.meanInc,tabRsld.stdInc, 'LineWidth',lw)
errorbar(log10(muVec),tabRsld.meanBack,tabRsld.stdBack, 'LineWidth',lw)
yline(groundTruthInc(iAcq), '--', 'Color',colors(1,:))
yline(groundTruthBack(end), '--', 'Color',colors(2,:))
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('Inc','Back')
title('RSLD')
ylim([0 1.5])

figure,
hold on
errorbar(log10(muVec),tabRed.meanInc,tabRed.stdInc, 'LineWidth',lw)
errorbar(log10(muVec),tabRed.meanBack,tabRed.stdBack, 'LineWidth',lw)
yline(groundTruthInc(iAcq), '--', 'Color',colors(1,:))
yline(groundTruthBack(end), '--', 'Color',colors(2,:))
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('Inc','Back')
title('RED')
ylim([0 1.5])

save_all_figures_to_directory(resultsDir,"target"+iAcq+"_metrics");
pause(0.1)
% close all,

