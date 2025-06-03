function [rf,x,z,sampleFreq] = loadAndBfLinear(fileName, presetName)
    [~, ~, ext] = fileparts(fileName);

    % Whether is a '.vrs' or '.mat' recover the RF Channel Data
    if (ext == ".vrs") % Data stream format
        channelData = vsv.file.readVRSFile(fileName);
    else %(ext == ".mat") % Classical Matlab format
        channelData = load(fileName);
        channelData = cell2mat(channelData.RcvData);
    end

    % Frame Selector
    channelData = channelData(:,:,end); % Last frame from Buffer Saved
    
    presetData = load(presetName);
    presetData = presetData.preSet; 

    Trans = presetData.Trans; 
    Receive = presetData.Receive;
    P = presetData.P;

    centralFreq = Receive(1).demodFrequency*1e6; % Central freq. [Hz]
    sampleFreq = Receive(1).decimSampleRate*1e6; % Sample freq. [Hz]

    %nPulses = P.numRays;
    nPulses = 128;
    endDepth = P.endDepth;
    nElements = Trans.numelements;
    nSamples = Receive(1).endSample - Receive(1).startSample + 1;
    %soundSpeed = presetData.speedOfSound; % [m/s]
    soundSpeed = 1540; % [m/s]
    wvl = soundSpeed/centralFreq;

    % Aux variables
    t = (0:(nSamples-1))/sampleFreq; % [sec.] time domain 0:T:(N_sample-1)*T
    z = soundSpeed*t/2;
    elementPosX = Trans.ElementPos(:, 1)*wvl; % x [m] X center of Trans.element
    elementPosZ = Trans.ElementPos(:, 3)*wvl; % z [m] Z center of Trans.element
    focus = soundSpeed*t(:)/2;
    elementPitch = Trans.spacingMm*1e-3;

    % Get reception delays
    rxDelays = getRXDelaysLinear(elementPosX, elementPosZ, focus, nElements, nPulses, soundSpeed);
    rxSamples = round(rxDelays*sampleFreq);
    rxSamples(rxSamples<1) = 1;
    rxSamples(rxSamples>nSamples) = nSamples;

    % Get apertures
    fNum = 2;
    maxAperSize = 64;
    dynAperture = getDynAperMat(fNum, z, elementPitch, maxAperSize, nElements, nPulses);

    % Organize and BMode
    rfChannel = getOrderedData(channelData, Receive, nSamples, nElements, nPulses);
    rf = getBModeLinear(rfChannel, rxSamples, dynAperture, nSamples, nPulses, nElements, endDepth);

    x = elementPosX'; % Each scanline centered at each element of the array [m]
end

%% Aux Functions
% Dynamic Focusing Delays
function rx_delays = getRXDelaysLinear(element_pos_x, element_pos_z, focus, n_elements, n_pulses, sound_speed)
    rx_delays = zeros(length(focus), n_elements, n_pulses); % length(focus) is same as length(t)
    for n = 1:n_pulses
        xfocus = element_pos_x(n);
        zfocus = element_pos_z(n) + focus;
        for e = 1:n_elements
            rx_delays(:,e,n) = (focus + ...
                sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
        end
    end
end

% Dynamic aperture matrix
function dyn_aperture = getDynAperMat(f_num, z, element_pitch, maxAperSize, n_elements, n_pulses)
    dyn_aperture = zeros(length(z), n_elements, n_pulses);
    for n = 1:n_pulses
        for z_i = 1:length(z)
            a = z(z_i)/(2*f_num);
            halfAperSize = floor(a / element_pitch);
            if (halfAperSize > maxAperSize/2)
                halfAperSize = floor(maxAperSize / 2);
            end
            a_i = -halfAperSize:halfAperSize;
            
            aper_center = n;
            aper = aper_center + a_i;
            
            aper = aper(aper>=1);
            aper = aper(aper<=128);

            dyn_aperture(z_i, aper, n) = 1;
        end
    end
end

% Organize data
function rf_channel = getOrderedData(BuffData, Receive, n_samples, n_elements, n_pulses)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(n).startSample:Receive(n).endSample, :); % RF Data from Buffer
    end
end

% Get BMode
function rf_data = getBModeLinear(rf_channel, rx_samples, dyn_aperture, n_samples, n_pulses, n_elements, endDepth)
    rf_data = zeros(n_samples, n_pulses);
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            rf_channel(:, e, n) = rf_channel(rx_samples(:, e, n), e, n);
        end
        % Dynamic Aperture
        rf_channel(:, :, n) = rf_channel(:, :, n) .* dyn_aperture(:, :, n);
        % Summing
        rf_data(:, n) = sum(rf_channel(:, :, n), 2);
    end
end