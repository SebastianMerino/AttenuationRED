startup,

resultsDir = 'Q:\smerino\REDjournalResults';

colors = lines(8);
lineWidth = 1.5;
gt = 0.49;
yLimits = [-0.2 1.2];
xLimits = [0 11];
muL = 10^0; muH = 10^10;

%% Liver cases
samples = ["simliverbigroi","simliversmallroi",...
    "simliverbigroihomo","simliversmallroihomo"];
for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical'});
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
samples = ["invivoliverbigroi","invivoliversmallroi"];

for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred = T(T.Method=='RED',:);
Trsld = T(T.Method=='RSLD',:); 
rangeRed = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
% title(sample)
hold on
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
legend('RSLD', 'RED')

end
%% Thyroid cases
samples = ["simthyrbigroi","simthyrsmallroi",...
    "simthyrbigroihomo","simthyrsmallroihomo"];
yLimits = [0,2];
gt = 1.21;

for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical'});
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
samples = ["invivothyrbigroi","invivothyrsmallroi"];

for sample = samples

opts = detectImportOptions(fullfile(resultsDir,"resultstesis.xlsx"));
opts.SelectedVariableNames = opts.VariableNames(3:end);  % Keeps columns 3 to end
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','double','double','double','double','double','categorical'});
opts.Sheet = sample;
T = readtable(fullfile(resultsDir,"resultstesis.xlsx"), opts);
Tred = T(T.Method=='RED',:);
Trsld = T(T.Method=='RSLD',:); 
rangeRed = Tred.Mu> muL & Tred.Mu<=muH;
rangeRsld = Trsld.Mu> muL & Trsld.Mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
% title(sample)
hold on
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
legend('RSLD', 'RED')

end


save_all_figures_to_directory(resultsDir,'plot','svg')