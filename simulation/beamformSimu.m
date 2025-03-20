startup,
baseDir = 'P:\smerino\simulation_acs\rf_data\25_03_14_inc';
targetFiles = dir(fullfile(baseDir,'rf_prebf_ref*.mat'));
resultsDir = fullfile(baseDir,'ref');
mkdir(resultsDir);
% targetFiles(5)

%%
% Setting
fNumber = 2;
deadband = 0.1e-2; % [cm]
source_f0 = 6.66e6;
c0 = 1540;
ii = 1;

%%
for ii = 6:length(targetFiles)
    file = targetFiles(ii);
    disp(file.name)
    tic
    load(fullfile(baseDir,file.name))
    toc

    tic
    t0 = max(time_delays) + 3.5/source_f0/2;
    rf = bfLineAcq(rf_prebf,x,z,t0,c0,fs,fNumber);
    izValid = z>deadband;
    z = z(izValid);
    rf = rf(izValid,:);
    toc
    %%
    bm = db(hilbert(rf));
    bm = bm - max(bm(:));
    figure('Position',[300 100 400 600])
    imagesc(x*100,z*100,bm, [-60 0])
    colormap gray
    axis image
    ylim([0.1 4.9]);

    pause(0.1)
    %%
    save(fullfile(resultsDir,targetFiles(ii).name([1:3,10:end])),...
        'rf','x','z','fs')
    % save(fullfile(resultsDir,"rf_ref"+(ii)+".mat"),...
    %     'rf','x','z','fs')
end

