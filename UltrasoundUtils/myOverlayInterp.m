function [hF,hB,hColor] = myOverlayInterp(ax, B,climB,xBm,zBm, overlay,clim,x,z, alpha)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm

[X,Z] = meshgrid(x,z);
[Xq,Zq] = meshgrid(xBm,zBm);
imgInterp = interp2(X,Z,overlay,Xq,Zq);
emptyRegion = isnan(imgInterp);
newRoi = ~emptyRegion;

B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(ax, xBm,zBm,B);%
axis image on;
colormap(ax, gray)

hColor = colorbar; colormap(ax, turbo);
hold on;
hF = imagesc(ax, xBm,zBm, imgInterp,clim);

% Make the foreground image transparent
alphadata = alpha.*(newRoi);
set(hF,'AlphaData',alphadata);
hold off
axis image