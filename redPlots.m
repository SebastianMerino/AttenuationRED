startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';
samplesDir = 'Q:\smerino\REDjournalResults\rf';

colors = lines(8);
lineWidth = 1.5;
gt = 0.49;
yLimits = [-0.2 1.2];
xLimitsRsld = [-1 8];
muLrsld = 10^0; muHrsld = 10^7;
xLimitsRed = [2 11];
muLred = 10^3; muHred = 10^10;


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
rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

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
xlim(xLimitsRsld)
ylim(yLimits)
legend('A1', 'A2', 'A3', 'Location','southeast')
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
errorbar(log10(Tred3.mu(rangeRed)),Tred3.meanInc(rangeRed), ...
    Tred3.stdInc(rangeRed)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRed)
ylim(yLimits)
legend('A1', 'A2', 'A3')
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
rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

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
xlim(xLimitsRsld)
ylim(yLimits)
legend('A1', 'A2')
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
xlim(xLimitsRed)
ylim(yLimits)
legend('A1', 'A2')
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
rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

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
xlim(xLimitsRsld)
ylim(yLimits)
legend('B1', 'B2', 'B3')
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
errorbar(log10(Tred3.mu(rangeRed)),Tred3.meanInc(rangeRed), ...
    Tred3.stdInc(rangeRed)/2,'vertical','s-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRed)
ylim(yLimits)
legend('B1', 'B2', 'B3')
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
rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

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
xlim(xLimitsRsld)
ylim(yLimits)
legend('B1', 'B2')
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
xlim(xLimitsRed)
ylim(yLimits)
legend('B1', 'B2')
title('RED')

save_all_figures_to_directory(resultsDir,'plot','svg')

%% ======================= PHANTOM DATA ======================= %%
startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';
samplesDir = 'Q:\smerino\REDjournalResults\layeredPhantom\final2';

colors = lines(8);
lineWidth = 1.5;
gt = 0.53;
yLimits = [-0.2 1.2];
xLimitsRsld = [0 11];
muLrsld = 10^0.9; muHrsld = 10^10.1;
xLimitsRed = [2 11];
muLred = 10^3; muHred = 10^10;

%% In vivo liver
sample = "8544Comp_F_3_Small";
excelFile = fullfile(samplesDir,sample+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:); 
rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

sample = "8544Comp_F_3_Big";
excelFile = fullfile(samplesDir,sample+".xlsx");
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
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRsld)
ylim(yLimits)
legend('A1', 'A2')
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
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRed)
ylim(yLimits)
legend('A1', 'A2')
title('RED')
%%
save_all_figures_to_directory(resultsDir,'phantoms','svg')

%% ================ NEW SIMULATED LIVER =====================%%
startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';
samplesDir = 'Q:\smerino\REDjournalResults\simuLiverNew';

colors = lines(8);
lineWidth = 1.5;
gt = 0.55;
yLimits = [-0.2 1.2];
xLimitsRsld = [-1 8];
muLrsld = 10^0; muHrsld = 10^7;
xLimitsRed = [2 11];
muLred = 10^3; muHred = 10^10;

excelFile = fullfile(samplesDir,"rf_prebf_liver2_cf0p4_acs0p6small.xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred1 = T(T.method=='RED-MED',:);
Trsld1 = T(T.method=='RSLD',:);

excelFile = fullfile(samplesDir,"rf_prebf_liver2_cf0p4_acs0p6big.xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:);

% sample = "simuLiver"; roi = "Big";
% excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
% opts = detectImportOptions(excelFile);
% opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
%     'double','categorical','double'});
% T = readtable(excelFile, opts);
% Tred3 = T(T.method=='RED-MED',:);
% Trsld3 = T(T.method=='RSLD',:);

rangeRed = Tred1.mu>=muLred & Tred1.mu<=muHred;
rangeRsld = Trsld1.mu>=muLrsld & Trsld1.mu<=muHrsld;

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
% errorbar(log10(Trsld3.mu(rangeRsld)),Trsld3.meanInc(rangeRsld), ...
%     Trsld3.stdInc(rangeRsld)/2,'vertical','s-', ...
%     'LineWidth',lineWidth, 'CapSize',3, ...
%     'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRsld)
ylim(yLimits)
legend('A1', 'A2', 'GT')
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
% errorbar(log10(Tred3.mu(rangeRed)),Tred3.meanInc(rangeRed), ...
%     Tred3.stdInc(rangeRed)/2,'vertical','s-', ...
%     'LineWidth',lineWidth, 'CapSize',3, ...
%     'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:)*0.75)
yline(gt, 'k--', 'LineWidth',lineWidth*0.5)
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimitsRed)
ylim(yLimits)
legend('A1', 'A2', 'GT')
title('RED')


save_all_figures_to_directory(resultsDir,'new','svg')
