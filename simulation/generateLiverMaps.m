%% ===================== NEW SIMULATION BATCH ===================== %%
startup,
baseDir = 'Q:\smerino\simulation_acs\rf_data\25_05_28_liver';
im = imread(fullfile(baseDir, 'liver2.png'));
im = im(:,:,1);
dx = 8.25e-5;
[m,n] = size(im);
xMask = (0:n-1)*dx; zMask = (0:m-1)*dx;
xMask = xMask- mean(xMask);
[X,Z] = meshgrid(xMask,zMask);

load(fullfile(baseDir,"refCords.mat"))

% [homog, fat, muscle, connective, liver]
imValues = [0,32,64,96,128];

acsLiver = 0.55;
for compressionFactor = 0.3:0.1:0.5

maskMap = interp2(X,Z*compressionFactor,im,rx,rz, 'nearest');
maskMap(isnan(maskMap)) = 0;
maskMap(rz>5e-2*compressionFactor) = imValues(end);

% Iterating for the rest of maps
rho0 = 1000; c0 = 1540; acs0 = 0.5; ba0 = 6;
% Sources: Douglas  
% [homog, fat, muscle, connective, liver]
rhoLayers = [rho0,950,1050,1120,rho0];
cLayers = [c0,c0,c0,c0,c0];
stdLayers = [4,2,4,8,4]/1000;
acsLayers = [acs0,0.48,1.09,1.57,acsLiver];
baLayers = [ba0,9.6,8,8,7.6];

soundSpeedMap = c0*ones(size(maskMap));
densityMap = rho0*ones(size(maskMap));
acsMap = acs0*ones(size(maskMap));
baMap = ba0*ones(size(maskMap));
stdMap = ones(size(maskMap));
for ii = 1:length(imValues)
    px = imValues(ii);
    soundSpeedMap(maskMap == px) = cLayers(ii);
    acsMap(maskMap == px) = acsLayers(ii);
    baMap(maskMap == px) = baLayers(ii);
    stdMap(maskMap == px) = stdLayers(ii);
    densityMap(maskMap == px) = rhoLayers(ii);
end
% Echogenicity map

stdMap(rz<=0) = 0;
% Plot
x = rx(1,:)*100;
z = rz(:,1)*100;
figure('Units','centimeters', 'Position',[5 5 30 8]),
tiledlayout(1,5, 'TileSpacing','compact', 'Padding','compact')
nexttile, imagesc(x,z,densityMap, [rho0*0.9,rho0*1.1])
axis image
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
title('Density')
colorbar

nexttile, imagesc(x,z,stdMap)
axis image
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
title('STD of SOS')
colorbar

nexttile, imagesc(x,z,soundSpeedMap, [c0*0.9,c0*1.1])
axis image
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
title('Sound speed')
colorbar

tacs = nexttile; imagesc(x,z,acsMap, [0 1.5])
axis image
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
title('ACS')
colorbar
colormap(tacs,"turbo")

tba = nexttile; imagesc(x,z,baMap, [5 12])
axis image
xlabel('Lateral [cm]')
ylabel('Axial [cm]')
title('B/A')
colorbar
colormap(tba,"turbo")
%%
name = sprintf("liver2_cf%.1f_acs%.1f",compressionFactor,acsLiver);
name = strrep(name,'.','p');
disp(name)
save_all_figures_to_directory(baseDir,name+"_fig")
pause(0.1)
close,
save(fullfile(baseDir,name),'densityMap','stdMap',...
    'soundSpeedMap','acsMap','baMap');

end