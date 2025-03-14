%% ========================= Processing loop ========================= %%
startup,

baseDir = 'C:\Users\smerino.C084288\Documents\Datasets\phantomVerasonicsSebas\processed';
dataDir = fullfile(baseDir,'bf');
sampleFiles = dir(fullfile(dataDir,'T7*.mat'));

refDir = fullfile(baseDir,'ref');
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = fullfile(baseDir,'results','T7');
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 4e6; freqH = 7e6;
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.25; % [cm]

% Blocksize parameters
blockParams.xInf = -2; % 0.8
blockParams.xSup = 2;
blockParams.zInf = 0.25;
blockParams.zSup = 5;
blockParams.blocksize = [20 20]*wl;
blockParams.overlap = 0.8;

% Plotting constants
dynRange = [-60,0];
attRange = [0.4,1.1];
bsRange = [-10,10];
yLimits = [deadband,5.5];

NptodB = log10(exp(1))*20;
iAcq = 1;
%%
out = matfile(fullfile(dataDir,sampleFiles(iAcq).name));
xBm = out.x*1e2; % [cm]
zBm = out.z'*1e2; % [cm]
sam1 = out.rf(:,:,1);
fs = out.fs;

% Plot region of interest B-mode image
bMode = db(hilbert(sam1));
bMode = bMode - max(bMode(zBm>deadband,:),[],"all");

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
    out = load(fullfile(refDir,refFiles(ff).name));
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

save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_spec");
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
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(xBm,zBm);
c1x = 0; 
c1z = 2.0;
rInc = 0.95;
roiL = 1; roiD = 0.6;
roiLz = 1;

% inc = ((Xq-c1x).^2 + (Zq-c1z).^2)<= (rInc-0.1)^2;
% back = ((Xq-c1x).^2 + (Zq-c1z).^2) >= (rInc+0.1)^2;
[back,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,roiD,roiLz);
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
groundTruthTargets = [0.97,0.95,0.95,0.55];

%% For looping
muVec = 10.^(0.5:0.5:7);
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
r.meanBack = mean(AttInterp(back),"omitnan");
r.stdBack = std(AttInterp(back),"omitnan");
r.biasBack = mean( AttInterp(back) - groundTruthTargets(end),"omitnan");
r.biasInc = mean( AttInterp(inc) - groundTruthTargets(iAcq),"omitnan");
r.rmseBack = sqrt( mean( (AttInterp(back) - groundTruthTargets(end)).^2,...
    "omitnan") );
r.rmseInc = sqrt( mean( (AttInterp(inc) - groundTruthTargets(iAcq)).^2,...
    "omitnan") );
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
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
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
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
ylim(yLimits)

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
ylim(yLimits)

pause(0.1)
saveas(gcf,fullfile(resultsDir,"sample"+iAcq+"_mu"+iMu+".png"))
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
yline(groundTruthTargets(iAcq), '--', 'Color',colors(1,:))
yline(groundTruthTargets(end), '--', 'Color',colors(2,:))
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
yline(groundTruthTargets(iAcq), '--', 'Color',colors(1,:))
yline(groundTruthTargets(end), '--', 'Color',colors(2,:))
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('Inc','Back')
title('RED')
ylim([0 1.5])

save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_metrics");
pause(0.1)
close all,

%% ========================= Beamforming loop ========================= %%
% startup,
% baseDir = 'C:\Users\smerino.C084288\Documents\Datasets\SavedDataQUSPhantom';
% acqNames = ["","_2","_3"];
% deadBand = 0.1e-2;
% 
% folders = dir(baseDir);
% folders = folders(3:end);
% 
% for iFolder = 1:length(folders)
%     folderStr = folders(iFolder).name;
%     subFolderStr = folderStr + "_F";
% 
%     for iAcq = 1:length(acqNames)
%         samName = subFolderStr + acqNames(iAcq);
%         fileName = fullfile(baseDir, folderStr, subFolderStr, samName);
%         presetName = fileName + "_preSet.mat";
% 
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
%         resultsDir = fullfile(baseDir, 'bf');
%         mkdir(resultsDir)
%         saveas(gcf, fullfile(resultsDir,samName+'.png'))
%         pause(1), close,
%         save(fullfile(resultsDir,samName),'rf','x','z','fs')
%     end
% end