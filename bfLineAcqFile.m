function [rf,x,z,fs] = bfLineAcqFile(fileName,c0,fNumber)
tic
% Beamforms line by line acquisition
source_f0 = 6.66e6;
source_cycles = 3.5;

out = load(fileName);
t0 = max(out.time_delays) + source_cycles/source_f0/2;
x = out.x; z = out.z; fs = out.fs;

[nSamples,nElems,nLines] = size(out.rf_prebf);

pitch = x(2) - x(1);
xElem = (1:nElems)*pitch;
xElem = xElem - mean(xElem);
rf = zeros(nSamples,nLines);
delays = calcDelaysLine(xElem,z,c0);
toc
tic
for iLine = 1:nLines
    rfLineArr = out.rf_prebf(1:nSamples,1:nElems,iLine);
    rf(:,iLine) = getLine(xElem,z,delays + t0, fs, fNumber,rfLineArr);
end
toc
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