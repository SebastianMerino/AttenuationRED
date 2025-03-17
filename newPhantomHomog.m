%% ========================= Processing loop ========================= %%
startup,

baseDir = 'C:\Users\smerino.C084288\Documents\Datasets\phantomVerasonicsSebas\processed';
dataDir = fullfile(baseDir,'bf');
sampleFiles = dir(fullfile(dataDir,'T7*.mat'));

refDir = fullfile(baseDir,'ref');
refFiles = dir(fullfile(refDir,'*.mat'));

resultsDir = fullfile(baseDir,'results','T7bottom3');
if ~exist("resultsDir","dir"); mkdir(resultsDir); end

%% Hyperparameters
% General parameters
c0 = 1540;
%freqL = 3.5e6; freqH = 7e6; % Depth :5.5cm, Arom acquisition
%freqL = 3.5e6; freqH = 8e6; % Inclusion, Arom acquisition, 
 freqL = 3.5e6; freqH = 7.5e6; % Depth : 5 cm, my acquisition
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.53; gammaRef = 1;
deadband = 0.25; % [cm]

% Blocksize parameters
blockParams.xInf = -1; % 0.8
blockParams.xSup = 1;
blockParams.zInf = 3.1;
blockParams.zSup = 5;
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Plotting constants
dynRange = [-60,0];
attRange = [0.3,1.1];
bsRange = [-10,10];
yLimits = [deadband,5.5];

NptodB = log10(exp(1))*20;
iAcq = 2;
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

figure('Position',[200 200 600 400]),
hold on
plot(f,mean(spectrumSamzf), 'LineWidth',2)
plot(f,mean(spectrumRefzf), 'LineWidth',2)
hold off
xlabel('f [MHz]')
title('Power spectrum')
grid on
hold on
xline(freqL/1e6, 'k--', 'LineWidth',2), xline(freqH/1e6, 'k--', 'LineWidth',2)
hold off
xlim([0 fs/2e6])
legend('Sample','Reference')


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


line = squeeze(mean(b,[1 2]))/4/L*NptodB;
fit = [ufr ones(length(ufr),1)]\line;

line2 = squeeze(mean(sld-compensation,[1 2]))/4/L*NptodB;

figure('Position',[200 200 600 400]),
plot(f,line2, 'LineWidth',2)
grid on
hold on
plot(f,fit(1)*f + fit(2), 'k--')
xline(freqL/1e6,'k--')
xline(freqH/1e6,'k--')
hold off
title("SLD, "+sprintf('ACS = %.2f f + %.2f',fit(1),fit(2)))
xlabel('f [MHz]')
ylabel('Attenuation [dB/cm]')
xlim([0 ufr(end)*1.5])
ylim([0 9])

save_all_figures_to_directory(resultsDir,"sample"+iAcq+"_spec");
pause(0.1)
close all,


%% Metrics
[X,Z] = meshgrid(x_ACS,z_ACS);
[Xq,Zq] = meshgrid(xBm,zBm);
c1x = 0.1; 
c1z = 2.1;
c2z = 4.0;
rInc = 0.95;
roiL = 1; roiD = 0.6;
roiLz = 1;

[~,inc] = getRegionMasks(xBm,zBm,c1x,c1z,roiL,roiD,roiLz);
[~,back] = getRegionMasks(xBm,zBm,c1x,c2z,roiL,roiD,roiLz);
x0mask = c1x - roiL/2; 
z0mask = c1z - roiLz/2;
groundTruthTargets = [0.97,0.95,0.95,0.55];

%% RSLD-TV
muRsld = 1e4;
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),muRsld,muRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));
CR = (reshape(Cn,m,n));

%% RED no weigths
muRed = 1e4;
tic
[err_fp2 ,u2]  =  admmRedMedianv2(A,b(:),muRed,tol,2*m*n,200,5,m,n,muRed/10);
toc,
BRED = reshape(u2(1:end/2)*NptodB,m,n);
CRED = reshape(u2(end/2+1:end)*NptodB,m,n);

%% SWIFT
muBswift = 1e3; muCswift = 10^0.5;
ratioCutOff = 10;
reject = 0.1;
extension = 3;

% First iteration
[~,Cn] = optimAdmmTvTikhonov(A1,A2,b(:),muBswift,muCswift,m,n,tol,mask(:));
bscMap = reshape(Cn*NptodB,m,n);

% Weight map
w = (1-reject)*(abs(bscMap)<ratioCutOff)+reject;
wExt = movmin(w,extension);

% Weight matrices and new system
W = repmat(wExt,[1 1 p]);
W = spdiags(W(:),0,m*n*p,m*n*p);
bw = W*b(:);        
A1w = W*A1;
A2w = W*A2;

% Second iteration
[Bn,cN] = optimAdmmWeightedTvTikhonov(A1w,A2w,bw,muBswift,muCswift,m,n,tol,mask(:),w);
BSWIFT = reshape(Bn*NptodB,m,n);
CSWIFT = reshape(Cn*NptodB,m,n);


%%
figure('Units','centimeters', 'Position',[5 5 32 8]);
tl = tiledlayout(1,4, "Padding","tight");

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
title(sprintf("RSLD, ACS=%.2f",mean(BR(:))))
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
ylim(yLimits)

t3 = nexttile;
myOverlay(t3, bMode,dynRange,xBm,zBm, BRED,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
axis image
title(sprintf("RED-MED, ACS=%.2f",mean(BRED(:))))
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
ylim(yLimits)

t3 = nexttile;
myOverlay(t3, bMode,dynRange,xBm,zBm, BSWIFT,attRange,x_ACS,z_ACS, 1);
xlabel('Lateral [cm]'),
axis image
title(sprintf("SWIFT, ACS=%.2f",mean(BSWIFT(:))))
c = colorbar;
c.Label.String = 'ACS [dB/cm/MHz]';
ylim(yLimits)

pause(0.1)
saveas(gcf,fullfile(resultsDir,"sample"+iAcq+".png"))
close,