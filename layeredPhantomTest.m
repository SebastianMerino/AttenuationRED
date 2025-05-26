%% ========================= Processing loop ========================= %%
startup,

dataDir = 'Q:\dataAvendano_May24\SavedData2405-SMA\bf';
refDir = 'Q:\dataAvendano_May24\SavedData2405-SMA\ref';
sampleFiles = dir(fullfile(dataDir,'*.mat'));
refFiles = dir(fullfile(refDir,'*Ref*.mat'));

resultsDir = 'Q:\smerino\REDjournalResults\layeredPhantom\CUEC';
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

for iRoi = 1:3

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 4.5e6; freqH = 8.5e6; % wide bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.2; % [cm]

% Blocksize parameters
if iRoi == 1
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 0.4;
    blockParams.zSup = 1.6;
elseif iRoi == 2
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 1.6;
    blockParams.zSup = 2.8;
elseif iRoi == 3
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 2.8;
    blockParams.zSup = 4;
else
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 0.4;
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
%%
iAcq = 4;
sampleName = sampleFiles(iAcq).name(1:end-4);
sampleName = string(sampleName) +"_roi"+iRoi;
out = load(fullfile(dataDir,sampleFiles(iAcq).name));
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


[u,~] = cgs(A'*A,A'*b(:));
B = (reshape(u(1:m*n)*NptodB,m,n));
fprintf("Mean ACS: %f\n",mean(B(:)))



%% Metrics
[X,Z] = meshgrid(xAcs,zAcs);
[Xq,Zq] = meshgrid(xBm,zBm);

% [~,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,1,roiLz);
inc = true(size(Xq));

%% For looping
muVec = 10.^[3,6];
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
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
    "omitnan");
r.method = 'RSLD';
r.mu = muRsld;
Metrics(iMu) = r;

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
myOverlay(t3, bMode,dynRange,xBm,zBm, B,attRange,xAcs,zAcs, 1);
xlabel('Lateral [cm]'),
axis image
title("SLD")
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


end