startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';
samplesDir = 'Q:\smerino\REDjournalResults\rf';

colors = lines(8);
lineWidth = 1.5;
gt = 0.49;
yLimits = [-0.2 1.2];
xLimits = [0 11];
muL = 10^0.9; muH = 10^10.1;

%% Simulated liver
sample = "simuLiver"; roi = "Small";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:);

sample = "simuLiverHomo"; roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:);

sample = "simuLiver"; roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred3 = T(T.method=='RED-MED',:);
Trsld3 = T(T.method=='RSLD',:);
rangeRed = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld = Trsld1.mu> muL & Trsld1.mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld)),Trsld1.meanInc(rangeRsld), ...
    Trsld1.stdInc(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld)),Trsld2.meanInc(rangeRsld), ...
    Trsld2.stdInc(rangeRsld)/2,'vertical','d-.', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
errorbar(log10(Trsld3.mu(rangeRsld)),Trsld3.meanInc(rangeRsld), ...
    Trsld3.stdInc(rangeRsld)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('A', 'B', 'C', 'Location','southeast')
title('RSLD')

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Tred1.mu(rangeRed)),Tred1.meanInc(rangeRed), ...
    Tred1.stdInc(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed)),Tred2.meanInc(rangeRed), ...
    Tred2.stdInc(rangeRed)/2,'vertical','d-.',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
errorbar(log10(Tred3.mu(rangeRsld)),Tred3.meanInc(rangeRsld), ...
    Tred3.stdInc(rangeRsld)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('A', 'B', 'C')
title('RED')


%% In vivo liver
sample = "invivoLiver";
roi = "Small";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:); 
rangeRed = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld = Trsld1.mu> muL & Trsld1.mu<=muH;

roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:); 

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld)),Trsld1.meanInc(rangeRsld), ...
    Trsld1.stdInc(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld)),Trsld2.meanInc(rangeRsld), ...
    Trsld2.stdInc(rangeRsld)/2,'vertical','d-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('iROI 1', 'iROI 2')
title('RSLD')

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Tred1.mu(rangeRed)),Tred1.meanInc(rangeRed), ...
    Tred1.stdInc(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed)),Tred2.meanInc(rangeRed), ...
    Tred2.stdInc(rangeRed)/2,'vertical','s-',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('iROI 1', 'iROI 2')
title('RED')

%% Simulated thyroid
yLimits = [0,2.5];
gt = 1.21;

sample = "simuThyroid"; roi = "Small";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:);

sample = "simuThyroidHomo"; roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:);

sample = "simuThyroid"; roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred3 = T(T.method=='RED-MED',:);
Trsld3 = T(T.method=='RSLD',:);
rangeRed = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld = Trsld1.mu> muL & Trsld1.mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld)),Trsld1.meanInc(rangeRsld), ...
    Trsld1.stdInc(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld)),Trsld2.meanInc(rangeRsld), ...
    Trsld2.stdInc(rangeRsld)/2,'vertical','d-.', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
errorbar(log10(Trsld3.mu(rangeRsld)),Trsld3.meanInc(rangeRsld), ...
    Trsld3.stdInc(rangeRsld)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('A', 'B', 'C')
title('RSLD')

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Tred1.mu(rangeRed)),Tred1.meanInc(rangeRed), ...
    Tred1.stdInc(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed)),Tred2.meanInc(rangeRed), ...
    Tred2.stdInc(rangeRed)/2,'vertical','d-.',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
errorbar(log10(Tred3.mu(rangeRsld)),Tred3.meanInc(rangeRsld), ...
    Tred3.stdInc(rangeRsld)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('A', 'B', 'C')
title('RED')

%% In vivo thyroid
sample = "invivoThyroid";
roi = "Small";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:); 
rangeRed = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld = Trsld1.mu> muL & Trsld1.mu<=muH;

roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:); 
rangeRed = Tred2.mu> muL & Tred2.mu<=muH;
rangeRsld = Trsld2.mu> muL & Trsld2.mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld)),Trsld1.meanInc(rangeRsld), ...
    Trsld1.stdInc(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld)),Trsld2.meanInc(rangeRsld), ...
    Trsld2.stdInc(rangeRsld)/2,'vertical','d-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('iROI 1', 'iROI 2')
title('RSLD')

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Tred1.mu(rangeRed)),Tred1.meanInc(rangeRed), ...
    Tred1.stdInc(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed)),Tred2.meanInc(rangeRed), ...
    Tred2.stdInc(rangeRed)/2,'vertical','s-',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('iROI 1', 'iROI 2')
title('RED')

save_all_figures_to_directory(resultsDir,'plot','svg')