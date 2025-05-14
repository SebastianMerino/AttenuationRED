startup,

resultsDir = 'Q:\smerino\REDjournalResults\plots';
samplesDir = 'Q:\smerino\REDjournalResults\rf';

colors = lines(8);
lineWidth = 1.5;
gt = 0.49;
yLimits = [-0.2 1.2];
xLimits = [0 11];
muL = 10^0.9; muH = 10^10.1;

%% Liver cases
samples = ["simuLiver","simuLiverHomo"];
rois = ["Small","Big"];
for sample = samples
    for roi = rois
        excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
        opts = detectImportOptions(excelFile);
        opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
            'double','categorical','double'});
        T = readtable(excelFile, opts);
        Tred = T(T.method=='RED-MED',:);
        Trsld = T(T.method=='RSLD',:);
        rangeRed = Tred.mu> muL & Tred.mu<=muH;
        rangeRsld = Trsld.mu> muL & Trsld.mu<=muH;

        figure('Units','centimeters', 'Position',[5 5 12 6]),
        % title(sample)
        hold on
        yline(gt, 'k--', 'LineWidth',lineWidth)
        errorbar(log10(Trsld.mu(rangeRsld)),Trsld.meanInc(rangeRsld), ...
            Trsld.stdInc(rangeRsld)/2,'vertical','d-', ...
            'LineWidth',lineWidth, 'CapSize',3, ...
            'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
        errorbar(log10(Tred.mu(rangeRed)),Tred.meanInc(rangeRed), ...
            Tred.stdInc(rangeRed)/2,'vertical','s-',...
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
end
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
rangeRed1 = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld1 = Trsld1.mu> muL & Trsld1.mu<=muH;

roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:); 
rangeRed2 = Tred2.mu> muL & Tred2.mu<=muH;
rangeRsld2 = Trsld2.mu> muL & Trsld2.mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld1)),Trsld1.meanInc(rangeRsld1), ...
    Trsld1.stdInc(rangeRsld1)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld2)),Trsld2.meanInc(rangeRsld2), ...
    Trsld2.stdInc(rangeRsld2)/2,'vertical','d-', ...
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
errorbar(log10(Tred1.mu(rangeRed1)),Tred1.meanInc(rangeRed1), ...
    Tred1.stdInc(rangeRed1)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed2)),Tred2.meanInc(rangeRed2), ...
    Tred2.stdInc(rangeRed2)/2,'vertical','s-',...
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

%% Thyroid cases
yLimits = [0,2.5];
gt = 1.21;

samples = ["simuThyroid","simuThyroidHomo"];
rois = ["Small","Big"];
for sample = samples
    for roi = rois
        excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
        opts = detectImportOptions(excelFile);
        opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
            'double','categorical','double'});
        T = readtable(excelFile, opts);
        Tred = T(T.method=='RED-MED',:);
        Trsld = T(T.method=='RSLD',:);
        rangeRed = Tred.mu> muL & Tred.mu<=muH;
        rangeRsld = Trsld.mu> muL & Trsld.mu<=muH;

        figure('Units','centimeters', 'Position',[5 5 12 6]),
        % title(sample)
        hold on
        yline(gt, 'k--', 'LineWidth',lineWidth)
        errorbar(log10(Trsld.mu(rangeRsld)),Trsld.meanInc(rangeRsld), ...
            Trsld.stdInc(rangeRsld)/2,'vertical','d-', ...
            'LineWidth',lineWidth, 'CapSize',3, ...
            'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:))
        errorbar(log10(Tred.mu(rangeRed)),Tred.meanInc(rangeRed), ...
            Tred.stdInc(rangeRed)/2,'vertical','s-',...
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
end

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
rangeRed1 = Tred1.mu> muL & Tred1.mu<=muH;
rangeRsld1 = Trsld1.mu> muL & Trsld1.mu<=muH;

roi = "Big";
excelFile = fullfile(samplesDir,sample,sample+roi+".xlsx");
opts = detectImportOptions(excelFile);
opts = setvartype(opts, {'double', 'double', 'double', 'double', ...
    'double','categorical','double'});
T = readtable(excelFile, opts);
Tred2 = T(T.method=='RED-MED',:);
Trsld2 = T(T.method=='RSLD',:); 
rangeRed2 = Tred2.mu> muL & Tred2.mu<=muH;
rangeRsld2 = Trsld2.mu> muL & Trsld2.mu<=muH;

figure('Units','centimeters', 'Position',[5 5 12 6]),
hold on
errorbar(log10(Trsld1.mu(rangeRsld1)),Trsld1.meanInc(rangeRsld1), ...
    Trsld1.stdInc(rangeRsld1)/2,'vertical','o:', ...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(1,:) + 0.2)
errorbar(log10(Trsld2.mu(rangeRsld2)),Trsld2.meanInc(rangeRsld2), ...
    Trsld2.stdInc(rangeRsld2)/2,'vertical','d-', ...
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
errorbar(log10(Tred1.mu(rangeRed1)),Tred1.meanInc(rangeRed1), ...
    Tred1.stdInc(rangeRed1)/2,'vertical','o:',...
    'LineWidth',lineWidth, 'CapSize',3, ...
    'MarkerFaceColor','auto', 'MarkerSize',4, 'Color',colors(5,:) *1.2)
errorbar(log10(Tred2.mu(rangeRed2)),Tred2.meanInc(rangeRed2), ...
    Tred2.stdInc(rangeRed2)/2,'vertical','s-',...
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