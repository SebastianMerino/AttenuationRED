%% ========================= Processing loop ========================= %%
startup,

dataDir = 'Q:\dataAvendano_May24\SavedData2405-SMA\bf';
refDir = 'Q:\dataAvendano_May24\SavedData2405-SMA\ref';
sampleFiles = dir(fullfile(dataDir,'*.mat'));
refFiles = dir(fullfile(refDir,'*Ref*.mat'));

resultsDir = 'Q:\smerino\REDjournalResults\layeredPhantom\final';
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
    blockParams.zInf = 0.25;
    blockParams.zSup = 4;
else
    blockParams.xInf = -2;
    blockParams.xSup = 2;
    blockParams.zInf = 3;
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
%%
for iAcq = 4 %1:length(sampleFiles)
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
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
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
r.biasInc = mean( AttInterp(inc) - groundTruth,"omitnan");
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruth).^2,...
    "omitnan") );
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruth) ),...
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
[~,iMu] = min(tabRed.maeInc,[],'omitmissing');
optimMuRed = muVec(iMu);

colors = lines(8);
figure,
hold on
errorbar(log10(muVec),tabRsld.meanInc,tabRsld.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(1,:))
errorbar(log10(muVec),tabRed.meanInc,tabRed.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(5,:))
yline(groundTruth, 'k--')
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('RSLD','RED')
ylim([0 1])



save_all_figures_to_directory(resultsDir,sampleName+"_metrics");
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


end
%% ======================== BEAMFORMING SCRIPT ======================== %%
% startup,
% baseDir = 'Q:\dataAvendano_May24\SavedData2405-SMA';
% resultsDir = fullfile(baseDir, 'bf');
% mkdir(resultsDir)
% deadBand = 0.1e-2;
% 
% folders = dir(baseDir);
% folders = folders(3:end);
% 
% for iFolder = 1:length(folders)
%     folderStr = folders(iFolder).name;
%     subFolderStr = folderStr + "_F";
% 
%     files = dir(fullfile(baseDir, folderStr, subFolderStr,"*preSet.mat"));
%     for iFile = 1:length(files)
%         samName = files(iFile).name(1:end-11);
%         fileName = fullfile(files(iFile).folder,samName);
%         presetName = fileName + "_preSet.mat";
% 
%         disp(samName)
%         if ~exist(presetName,'file'), continue; end
%         [rf,x,z,fs] = loadAndBfLinear(fileName, presetName);
% 
%         bMode = db(hilbert(rf));
%         bMode = bMode - max(bMode(z>deadBand,:),[],'all');
%         figure,
%         imagesc(x*1e2, z*1e2, bMode);
%         axis image
%         colormap gray;
%         clim([-60 0]);
%         colorbar
%         ylabel('[mm]', 'FontSize', 10);
%         xlabel('[mm]', 'FontSize', 10);
%         ylim([deadBand*100 5])
% 
% 
%         saveas(gcf, fullfile(resultsDir,samName+".png"))
%         pause(1), close,
%         save(fullfile(resultsDir,samName),'rf','x','z','fs')
%     end
% end
