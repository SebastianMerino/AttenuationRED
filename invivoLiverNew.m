%% ========================= Processing loop ========================= %%
startup,

dataDir = "Q:\smerino\REDjournalResults\newLiver";

sampleName = "invivoLiver";
resultsDir = "Q:\smerino\REDjournalResults\newLiver\"+sampleName+"_med7";
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
freqL = 2e6; freqH = 6e6; % wide bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.54; gammaRef = 1;  % 0.4 for simulations, 0.53 for in vivo
iAcq = 1;
groundTruthTargets = 0.5;

% Blocksize parameters
if big
    blockParams.xInf = xRf(1); % 0.8
    blockParams.xSup = xRf(end);
    blockParams.zInf = 1.5;
    blockParams.zSup = 9;
else
    blockParams.xInf = xRf(1); % 0.8
    blockParams.xSup = xRf(end);
    blockParams.zInf = 5.5;
    blockParams.zSup = 8.5;
end
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Measurement ROI
c1x = 0; c1z = 7;
roiL = 3; roiLz = 1.5;

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.2];
bsRange = [-10,10];
yLimits = [zBm(1),9];

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
muVec = 10.^(1:0.5:9);
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
[~ ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,7,m,n,muRed);
% [~,~,u2]  =  admm_red_median(A'*A,A'*b(:),muRed,tol,2*m*n,1500,4,1,7,m,n,muRed);
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

optimMuRsld = min(muVec(tabRsld.stdInc./tabRsld.meanInc<0.1));
optimMuRed = min(muVec(tabRed.stdInc./tabRed.meanInc<0.12));

colors = lines(8);
figure,
hold on
errorbar(log10(muVec),tabRsld.meanInc,tabRsld.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(1,:))
errorbar(log10(muVec),tabRed.meanInc,tabRed.stdInc, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(5,:))
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('RSLD','RED')
ylim([-0.5 1.5])



save_all_figures_to_directory(resultsDir,sampleName+"_metrics",'svg');
pause(0.1)
close all,


%% Optimal mu plot
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),optimMuRsld,optimMuRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));


tic
[err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),optimMuRed,tol,2*m*n,200,7,m,n,optimMuRed);
% [~,~,u2]  =  admm_red_median(A'*A,A'*b(:),optimMuRed,tol,2*m*n,1500,4,1,7,m,n,optimMuRed);
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
%%

% clear, clc
% 
% baseDir = "Q:\smerino\REDjournalResults\newLiver";
% refFiles = dir(fullfile(baseDir,"RF_544_F*"));
% for iRef = 1:length(refFiles)
%     out = load(fullfile(baseDir,refFiles(iRef).name));
%     ref(:,:,iRef) = out.rf_data;
% end
% 
% out = load(fullfile(baseDir,"RF_LiverEAIntercostal_F.mat")); 
% rf = out.rf_data;
% 
% load(fullfile(baseDir,"preSet_544_F.mat")); 
% fs = preSet.Receive(1).decimSampleRate*1e6; % According to "NS200BW" Acquisition Mode
% elem_pitch = preSet.Trans.spacingMm*1e-3;
% 
% bMode = db(hilbert(rf));
% bMode = bMode - max(bMode(:));
% xBm = (1:size(rf,2))*elem_pitch; 
% xBm = xBm - mean(xBm);
% zBm = (1:size(rf,1))/fs*1540/2;
% xRf = xBm*100;
% zRf = zBm*100;
% 
% 
% figure,
% imagesc(xRf,zRf,bMode,[-70 0])
% xlabel('Lateral [cm]'),
% ylabel('Axial [cm]')
% axis image
% title('B-mode')
% colormap(gray)
% c = colorbar;
% c.Label.String = 'dB';
% ylim([0.1 8])
% 
% save("invivoLiver.mat","ref","rf","bMode","fs","xBm","zBm","xRf","zRf")
% 
% %%
% clear, clc
% 
% baseDir = "Q:\smerino\REDjournalResults\newLiver";
% refFiles = dir(fullfile(baseDir,"*ref*"));
% for iRef = 1:length(refFiles)
%     out = load(fullfile(baseDir,refFiles(iRef).name));
%     ref(:,:,iRef) = out.rf;
% end
% 
% out = load(fullfile(baseDir,"rf-liver2_cf0p8_acs0p5.mat")); 
% rf = out.beamformed_channel_data;
% bMode = db(hilbert(rf));
% bMode = bMode - max(bMode(:));
% fs = out.channel_data2.fs;
% xBm = out.channel_data2.x;
% zBm = out.channel_data2.z;
% xRf = xBm*100;
% zRf = zBm*100;
% 
% figure,
% imagesc(xRf,zRf,bMode,[-70 0])
% xlabel('Lateral [cm]'),
% ylabel('Axial [cm]')
% axis image
% title('B-mode')
% colormap(gray)
% c = colorbar;
% c.Label.String = 'dB';
% ylim([0.1 8])
% 
% save("simuLiver.mat","ref","rf","bMode","fs","xBm","zBm","xRf","zRf")