% ======================================================================
% ======================================================================
%% PHANTOMSSS
startup,

dataDir = 'Q:\smerino\phantoms\ID316V2\06-08-2023-Generic';
sampleFiles = dir(fullfile(dataDir,'*.mat'));
sampleFiles = sampleFiles(end-2:end);

refDir = 'Q:\smerino\phantoms\ID544V2\06-08-2023-Generic';
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = 'Q:\smerino\REDjournalResults\phantoms\anthony_med5';
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 2.5e6; freqH = 7.5e6;
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.1; % [cm]

% Blocksize parameters
blockParams.xInf = 0; % 0.8
blockParams.xSup = 4;
blockParams.zInf = 0.2;
blockParams.zSup = 4;
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Plotting constants
dynRange = [-60,0];
attRange = [0.4,1.1];
bsRange = [-10,10];
yLimits = [0.1 3.5];

NptodB = log10(exp(1))*20;
iAcq = 2;
%% For looping each phantom
for iAcq = 2:2
out = matfile(fullfile(dataDir,sampleFiles(iAcq).name));
xBm = out.x*1e2; % [cm]
zBm = out.z'*1e2; % [cm]
sam1 = out.RF(:,:,1);
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
    out = matfile(fullfile(refDir,refFiles(ff).name));
    rfRef(:,:,ff) = out.RF(:,:,end);
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
c1x = 1.9; 
c1z = 1.9;
% c2z = 4;
rInc = 0.95;
roiL = 0.9; roiD = 0.6;
roiLz = 1.1;

% [back,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,roiD,roiLz);
% [~,back] = getRegionMasks(xBm,zBm,c1x,c2z,roiL,roiD,roiLz);

inc = (Xq-c1x).^2 + (Zq-c1z).^2 < (rInc-0.2).^2;
back = (Xq-c1x).^2 + (Zq-c1z).^2 > (rInc+0.2).^2;

groundTruthTargets = [0.97,0.95,0.95,0.55];

% [~,backAcs] = getRegionMasks(x_ACS,z_ACS,c1x,c2z,roiL,roiD,roiLz);
% 
% sldComp = sld(:,:,:) - compensation(:,:,:);
% sldComp = sldComp.*backAcs;
% sldLine = squeeze(sum(sldComp,[1,2]))/sum(backAcs(:));
% figure,plot(f,sldLine)
% grid on
% hold on
% xline(freqL/1e6, 'k--')
% xline(freqH/1e6, 'k--')
% xlim([0 freqH*1.5/1e6])

%% For looping
muVec = 10.^(0.5:0.5:10);
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
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
r.cnr = abs(r.meanBack - r.meanInc)/sqrt(r.stdInc^2 + r.stdBack^2);
r.method = 'RSLD';
r.mu = muRsld;
Metrics(iMu) = r;

%% RED no weigths
muRed = muVec(iMu);
tic
% [err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,5,m,n,muRed);
[~ ,~,u2] = admm_red_median(A'*A,A'*b(:),muRed,0.001,size(A'*b(:),1),1500,4,1,5,m,n,muRed);
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
r.maeBack = mean(  abs( (AttInterp(back) - groundTruthTargets(end)) ),...
    "omitnan");
r.maeInc = mean(  abs( (AttInterp(inc) - groundTruthTargets(iAcq)) ),...
    "omitnan");
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
hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off

t2 = nexttile;
myOverlay(t2, bMode,dynRange,xBm,zBm, BR,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
colormap(t2,turbo)
axis image
title("RSLD, \mu=10^{"+log10(muRsld)+"}")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
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
% rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
%     'LineWidth',1, 'Curvature',1)
contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)
%%
pause(0.1)
saveas(gcf,fullfile(resultsDir,"sample"+iAcq+"_mu"+iMu+".png"))
close,
end

%% Metrics plots
lw = 2;

T = struct2table(Metrics);
T.method = categorical(T.method);
muVec = T(T.method=='RSLD',:).mu;
tabRsld = T(T.method=='RSLD',:);
tabRed = T(T.method=='RED-MED',:);

figure('Units','centimeters', 'Position',[5 5 9 6]);
hold on
plot(log10(muVec),tabRsld.maeInc/2 + tabRsld.maeBack/2, 'LineWidth',lw)
plot(log10(muVec),tabRed.maeInc/2 + tabRed.maeBack/2, 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('MAE [dB/cm/MHz]')
grid on
legend('RSLD','RED')
title('MAE')
ylim([0 0.9])

[~,iMu] = min(tabRsld.maeInc/2 + tabRsld.maeBack/2);
optimMuRsld = muVec(iMu);
[~,iMu] = min(tabRed.maeInc/2 + tabRed.maeBack/2);
optimMuRed = muVec(iMu);

colors = [0    0.4470    0.7410; 0.3010    0.7450    0.9330];

figure('Units','centimeters', 'Position',[5 5 9 6]);
hold on
errorbar(log10(muVec),tabRsld.meanInc,tabRsld.stdInc/2, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(1,:))
errorbar(log10(muVec),tabRsld.meanBack,tabRsld.stdBack/2, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(2,:))
yline(groundTruthTargets(iAcq), '--', 'Color',colors(1,:), 'LineWidth',lw*0.8)
yline(groundTruthTargets(end), '--', 'Color',colors(2,:), 'LineWidth',lw*0.8)
xline(log10(optimMuRsld), 'Color','b', 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('Inc','Bgnd')
title('RSLD')
ylim([0 1.5])

colors = [0   0.7410  0.4470    ; 0.3010  0.9330  0.7450]/1.2;
figure('Units','centimeters', 'Position',[5 5 9 6]);
hold on
errorbar(log10(muVec),tabRed.meanInc,tabRed.stdInc/2, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(1,:))
errorbar(log10(muVec),tabRed.meanBack,tabRed.stdBack/2, ...
    'd-.', 'LineWidth',lw, 'MarkerFaceColor','auto', 'Color',colors(2,:))
yline(groundTruthTargets(iAcq), '--', 'Color',colors(1,:), 'LineWidth',lw*0.8)
yline(groundTruthTargets(end), '--', 'Color',colors(2,:), 'LineWidth',lw*0.8)
xline(log10(optimMuRed), 'Color','g', 'LineWidth',lw)
hold off
xlabel('log_{10}\mu')
ylabel('ACS [dB/cm/MHz]')
grid on
legend('Inc','Bgnd')
title('RED')
ylim([0 1.5])

save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_metrics",'svg');
pause(0.1)
close all,


%% Optimal mu plot
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),optimMuRsld,optimMuRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));


tic
[~ ,~,u2] = admm_red_median(A'*A,A'*b(:),muRed,0.001,size(A'*b(:),1),1500,4,1,5,m,n,muRed);
% [err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,5,m,n,muRed);
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
myOverlayInterp(t2, bMode,dynRange,xBm,zBm, BR,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
colormap(t2,turbo)
axis image
title("RSLD")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
% contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
% contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)

t3 = nexttile;
myOverlayInterp(t3, bMode,dynRange,xBm,zBm, BRED,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
colormap(t3,turbo)
axis image
title("RED")
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
hold on
rectangle('Position',[c1x-rInc c1z-rInc 2*rInc 2*rInc], 'LineStyle','--', ...
    'LineWidth',1, 'Curvature',1, 'EdgeColor','w')
% contour(xBm,zBm,inc, [0 1], 'w--', 'LineWidth',1.5)
% contour(xBm,zBm,back, [0 1], 'w--', 'LineWidth',1.5)
hold off
ylim(yLimits)


save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_final",'svg');
pause(0.1)
close all,
T = [tabRsld;tabRed];
writetable(T,fullfile(resultsDir,'phantom.xlsx'))

end