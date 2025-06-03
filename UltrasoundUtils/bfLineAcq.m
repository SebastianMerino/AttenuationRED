function rf = bfLineAcq(rf_prebf,x,z,t0,c0,fs,fNumber)
% Beamforms line by line acquisition

pitch = x(2) - x(1);
xElem = (1:size(rf_prebf,2))*pitch;
xElem = xElem - mean(xElem);
rf = zeros(size(rf_prebf,[1 3]));
delays = calcDelaysLine(xElem,z,c0);

for iLine = 1:size(rf_prebf,3)
    rfLineArr = rf_prebf(:,:,iLine);
    rf(:,iLine) = getLine(xElem,z,delays + t0, fs, fNumber,rfLineArr);
end

end

%% Functions
function rfLine = getLine(x,z,delays,fs,fn,rfLineArr)
% x:    lateral coordinates of pitch elements
% z:    axial coordinates of points

pitch = x(2) - x(1);
rfLine = zeros(length(z),1);

for iz = 1:length(z)
    aperture = z(iz)/fn;
    aperture = max(aperture,2*pitch);
    rxElementMask = abs(x)<=aperture/2;
    rxElements = 1:length(x);
    rxElements = rxElements(rxElementMask);

    delaySam = round((delays(iz,rxElements))*fs + 1);
    for ix = 1:length(rxElements)
        if delaySam(ix)>size(rfLineArr,1), continue; end
        rfLine(iz) = rfLine(iz) + rfLineArr(delaySam(ix),rxElements(ix));
    end
end
end

function delays = calcDelaysLine(x,z,c0)
% x:    lateral coordinates of pitch elements
% z:    axial coordinates of points
% c0:   sound speed
  
if iscolumn(x); x = x'; end
if ~iscolumn(z); z = z'; end

delays = z/c0 + sqrt(x.^2+z.^2)/c0;
end