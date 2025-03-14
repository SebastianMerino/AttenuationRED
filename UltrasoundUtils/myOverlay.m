function [hF,hB,hColor] = myOverlay(ax, B,climB,xBm,zBm, overlay,clim,x,z, alpha)
% IMOVERLAY(B,F) displays the image SWS transparently over the image B.
%   alpha:  transparency
%   x:      lateral coordinate in mm
%   z:      depth in mm

if isempty(ax)
    figure, ax = gca;
end

B = repmat(mat2gray(double(B),double(climB)),[1,1,3]);

hB = imagesc(ax, xBm,zBm,B);%
axis image on;
colormap(ax, gray)

hColor = colorbar; colormap(ax, turbo);
hold on;
hF = imagesc(ax, x,z,overlay,clim);

% Make the foreground image transparent
alphadata = alpha.*(ones(size(overlay)));
set(hF,'AlphaData',alphadata);
hold off
axis image
