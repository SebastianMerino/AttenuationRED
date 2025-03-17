function [hF,hB,hColor] = myOverlayInterp(ax, B,climB,xBm,zBm, overlay,clim,x,z, alpha)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm

if isempty(ax)
    figure, ax = gca;
end

[X,Z] = meshgrid(x,z);

dx = 0.02;
xNew = xBm(1):dx:xBm(end);
zNew = zBm(1):dx:zBm(end);
[Xq,Zq] = meshgrid(xNew,zNew);
overlayInterp = interp2(X,Z,overlay,Xq,Zq, 'cubic');

B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(ax, xBm,zBm,B);%
axis image on;
colormap(ax, gray)

hColor = colorbar; colormap(ax, turbo);
hold on;
hF = imagesc(ax, xNew,zNew,overlayInterp,clim);

% Make the foreground image transparent
alphadata = alpha.*(~isnan(overlayInterp));
set(hF,'AlphaData',alphadata);
hold off
axis image
