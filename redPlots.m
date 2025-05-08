startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';

colors = lines(8);
lineWidth = 1.5;
gt = 0.49;
yLimits = [-0.2 1.2];
xLimits = [0 11];
muL = 10^1; muH = 10^10;

%% Liver cases
samples = ["simliverbigroi","simliversmallroi",...
    "simliverbigroihomo","simliversmallroihomo"];
for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred = T(T.Method=='RED',:);
Trsld = T(T.Method=='RSLD',:); 
rangeRed = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
% title(sample)
hold on
yline(gt, 'k--', 'LineWidth',lineWidth)
errorbar(log10(Trsld.Mu(rangeRsld)),Trsld.R1Mean(rangeRsld), ...
    Trsld.R1Std(rangeRsld)/2,'vertical','d-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
errorbar(log10(Tred.Mu(rangeRed)),Tred.R1Mean(rangeRed), ...
    Tred.R1Std(rangeRed)/2,'vertical','s-',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('', 'RSLD', 'RED')

end

%% In vivo liver
samples = ["invivoliversmallroi","invivoliverbigroi"];

sample = samples(1);
opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred1 = T(T.Method=='RED',:);
Trsld1 = T(T.Method=='RSLD',:); 
rangeRed1 = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld1 = Trsld.Mu> muL & Trsld.Mu<=muH;

sample = samples(2);
opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred2 = T(T.Method=='RED',:);
Trsld2 = T(T.Method=='RSLD',:); 
rangeRed2 = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld2 = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.Mu(rangeRsld)),Trsld1.R1Mean(rangeRsld), ...
    Trsld1.R1Std(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.Mu(rangeRsld)),Trsld2.R1Mean(rangeRsld), ...
    Trsld2.R1Std(rangeRsld)/2,'vertical','d-', ...
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
errorbar(log10(Tred1.Mu(rangeRed)),Tred1.R1Mean(rangeRed), ...
    Tred1.R1Std(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.Mu(rangeRed)),Tred2.R1Mean(rangeRed), ...
    Tred2.R1Std(rangeRed)/2,'vertical','s-',...
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

fprintf("RSLD\n")
[muMin] = min(Trsld1.Mu(Trsld1.R1Std./Trsld1.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(1),log10(muMin));
[muMin] = min(Trsld2.Mu(Trsld2.R1Std./Trsld2.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(2),log10(muMin));
fprintf("RED\n")
[muMin] = min(Tred1.Mu(Tred1.R1Std./Tred1.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(1),log10(muMin));
[muMin] = min(Tred2.Mu(Tred2.R1Std./Tred2.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(2),log10(muMin));


%% Thyroid cases
samples = ["simthyrbigroi","simthyrsmallroi",...
    "simthyrbigroihomo","simthyrsmallroihomo"];
yLimits = [0,2.5];
gt = 1.21;

for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred = T(T.Method=='RED',:);
Trsld = T(T.Method=='RSLD',:); 
rangeRed = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
title(sample)
hold on
yline(gt, 'k--', 'LineWidth',lineWidth)
errorbar(log10(Trsld.Mu(rangeRsld)),Trsld.R1Mean(rangeRsld), ...
    Trsld.R1Std(rangeRsld)/2,'vertical','d-', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
errorbar(log10(Tred.Mu(rangeRed)),Tred.R1Mean(rangeRed), ...
    Tred.R1Std(rangeRed)/2,'vertical','s-',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:))
hold off
xlabel('log_{10}(\mu)')
ylabel('ACS [dB/cm/MHz]')
grid on
xlim(xLimits)
ylim(yLimits)
legend('', 'RSLD', 'RED')


end

%% In vivo thyroid
samples = ["invivothyrsmallroi","invivothyrbigroi"];

sample = samples(1);
opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred1 = T(T.Method=='RED',:);
Trsld1 = T(T.Method=='RSLD',:); 
rangeRed1 = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld1 = Trsld.Mu> muL & Trsld.Mu<=muH;

sample = samples(2);
opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical','double'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred2 = T(T.Method=='RED',:);
Trsld2 = T(T.Method=='RSLD',:); 
rangeRed2 = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld2 = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.Mu(rangeRsld)),Trsld1.R1Mean(rangeRsld), ...
    Trsld1.R1Std(rangeRsld)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.Mu(rangeRsld)),Trsld2.R1Mean(rangeRsld), ...
    Trsld2.R1Std(rangeRsld)/2,'vertical','d-', ...
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
errorbar(log10(Tred1.Mu(rangeRed)),Tred1.R1Mean(rangeRed), ...
    Tred1.R1Std(rangeRed)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.Mu(rangeRed)),Tred2.R1Mean(rangeRed), ...
    Tred2.R1Std(rangeRed)/2,'vertical','s-',...
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

fprintf("RSLD\n")
[muMin] = min(Trsld1.Mu(Trsld1.R1Std./Trsld1.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(1),log10(muMin));
[muMin] = min(Trsld2.Mu(Trsld2.R1Std./Trsld2.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(2),log10(muMin));
fprintf("RED\n")
[muMin] = min(Tred1.Mu(Tred1.R1Std./Tred1.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(1),log10(muMin));
[muMin] = min(Tred2.Mu(Tred2.R1Std./Tred2.R1Mean < 0.1));
fprintf("Sample %s, Mu: %.1f\n",samples(2),log10(muMin));

save_all_figures_to_directory(resultsDir,'plot','svg')