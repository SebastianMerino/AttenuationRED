startup,

resultsDir = 'Q:\smerino\REDjournalResults\finalImages';
samplesDir = 'Q:\smerino\REDjournalResults\rf';

mkdir(resultsDir)

%%
sample = "simuLiver";
big = true;
if big 
    roi = "Big";
else
    roi = "Small";
end
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:);

optimMuRsld = min(Trsld1.mu(Trsld1.stdInc./Trsld1.meanInc<0.1));
optimMuRed = min(Tred1.mu(Tred1.stdInc./Tred1.meanInc<0.1));

%%
sampleName = sample;
load(fullfile(samplesDir,sampleName+".mat"))
sampleName = sampleName + roi;

zRf = zRf';
zRef = zRef';
xBm = xBm*100; zBm = zBm'*100;


%% Hyperparameters
% General parameters
c0 = 1540;
freqL = 2e6; freqH = 9e6; % wide bandwidth
wl = 2*c0/(freqL + freqH);
alpha0Ref = 0.4; gammaRef = 1;  % 0.4 for simulations, 0.53 for in vivo
iAcq = 1;
groundTruthTargets = 0.49;

% Blocksize parameters
if big
    blockParams.xInf = xRf(1); % 0.8
    blockParams.xSup = xRf(end);
    blockParams.zInf = zRf(1);
    blockParams.zSup = 5.5;
else
    blockParams.xInf = xRf(1); % 0.8
    blockParams.xSup = xRf(end);
    blockParams.zInf = 3;
    blockParams.zSup = 5.5;
end
blockParams.blocksize = [15 15]*wl;
blockParams.overlap = 0.8;

% Measurement ROI
c1x = 0; c1z = 4.2;
roiL = 2.5; roiLz = 1.2;

% Plotting constants
dynRange = [-60,0];
attRange = [0,1.2];
bsRange = [-10,10];
yLimits = [zBm(1),5.5];

NptodB = log10(exp(1))*20;


%% Sample power spectrum
% get Spectra
[Sp,Sd,xAcs,zAcs,f] = getSpectrum(rf,xRf,zRf,fs,blockParams);


%% Generating Diffraction compensation
% Generating references
clear att_ref_map 
att_ref_map(1,1,:) = alpha0Ref*f/NptodB;

[SpRef,SdRef,~,~,~] = getSpectrum(ref,xRef,zRef,fs,blockParams);

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

%% Optimal mu plot
tic
[Bn,Cn] = AlterOpti_ADMM(A1,A2,b(:),optimMuRsld,optimMuRsld,m,n,tol,mask(:));
toc
BR = (reshape(Bn*NptodB,m,n));


tic
[~,~,u2]  =  admm_red_median(A'*A,A'*b(:),optimMuRed,tol,2*m*n,1500,4,1,7,m,n,optimMuRed);
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
ylabel('Axial [cm]')
colormap(t2,turbo)
axis image
title("\mu=10^{"+log10(optimMuRsld)+"}")
hold on
contour(xBm,zBm,inc, [0 1], 'w', 'LineWidth',2)
hold off
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

pause(0.1)
saveas(gcf,fullfile(resultsDir,sampleName+".png"));
saveas(gcf,fullfile(resultsDir,sampleName+".svg"));
close all,
