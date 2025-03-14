function [rf,xPolar,zPolar,z0Polar,th,r,fs] = loadAndBf(fileName, presetName, attCoef, soundSpeed)

% genBModeMonofocal.m
% Función para generar imágenes en modo-B a partir de datos de RF de ultrasonido
% utilizando un enfoque monofocal con apodización dinámica.
%
% Esta función procesa datos crudos de RF para generar imágenes de ultrasonido
% en modo-B, aplicando técnicas de enfoque dinámico y apodización para mejorar
% la calidad de imagen.
%
% Sintaxis:
%   genBModeMonofocal(fileName, presetName, ax, attCoef, soundSpeed)
%
% Parámetros de Entrada:
%   fileName    - String. Ruta del archivo que contiene los datos de RF.
%                 Soporta archivos '.vrs' o '.mat'
%   presetName  - String. Ruta del archivo de preset que contiene la 
%                 configuración del transductor y parámetros de adquisición
%   attCoef     - Double. Coeficiente de atenuación para compensación TGC
%   soundSpeed  - Double. Velocidad del sonido en el medio [m/s]

    [~, ~, ext] = fileparts(fileName);
    % Whether is a '.vrs' or '.mat' recover the RF Channel Data
    if (ext == ".vrs")
        channelData = vsv.file.readVRSFile(fileName);
        
    elseif (ext == ".mat")
        channelData = load(fileName);
        channelData = cell2mat(channelData.RcvData);
    end

    presetData = load(presetName);
    presetData = presetData.preSet;
    
    Trans = presetData.Trans;
    Receive = presetData.Receive;
    P = presetData.P;
    txNum = 1;

    centralFreq = Receive(1).demodFrequency*1e6; % Central freq. [Hz]
    sampleFreq = Receive(1).decimSampleRate*1e6; % Sample freq. [Hz]

    nPulses = P.numRays;
    endDepth = P.endDepth;
    nElements = Trans.numelements;
    nSamples = Receive(1).endSample - Receive(1).startSample + 1;
    wvl = soundSpeed/centralFreq;


    % Aux variables
    t = (0:(nSamples-1))/sampleFreq; % [sec.] time domain 0:T:(N_sample-1)*T
    z = soundSpeed*t/2*1e3;
    elementPosX = Trans.ElementPos(:, 1)*wvl; % x [m] X center of Trans.element
    elementPosZ = Trans.ElementPos(:, 3)*wvl; % z [m] Z center of Trans.element
    phi = Trans.ElementPos(:, 4); % phi [rad] angle for every Trans.element
    focus = soundSpeed*t(:)/2;
    tgcDbCm = attCoef * (centralFreq * 1e-6) ^ 1;
    tgcNpM = tgcDbCm / 8.686 * 100;
    r = exp(tgcNpM * z*1e-3); % [m]

    % Get reception delays
    rxDelays = getRXDelays(elementPosX, elementPosZ, phi, focus, nElements, nPulses, soundSpeed);
    rxSamples = round(rxDelays*sampleFreq);
    rxSamples(rxSamples<1) = 1;
    rxSamples(rxSamples>nSamples) = nSamples;

    % Apodization
    maxAperSize = 65;
    apodAperture = getApodCoeff(maxAperSize, nElements, nPulses);

    % Organize and BMode (Frame 1)
    rfChannel = getOrderedData(channelData(:, :, end), Receive, nSamples, nElements, nPulses, txNum);
    rf = getRfC52v(rfChannel, rxSamples, apodAperture, nSamples, nPulses, nElements, endDepth, r);

    % For polar grid
    param = getparam('C5-2v');
    siz = size(rf);
    z = soundSpeed*t/2;
    zmax = z(end);
    R = param.radius;
    p = param.pitch;
    N = param.Nelements;
    L = 2*R*sin(asin(p/2/R)*(N-1)); % chord length
    d = sqrt(R^2-L^2/4); % apothem
    z0 = -d;
    
    th = -(linspace(atan2(L/2,d),atan2(-L/2,d),siz(2)))*180/pi;
    r = linspace(R+p,-z0+zmax,siz(1));
    
    % To Polar Coordinates
    [xPolar,zPolar, z0Polar] = impolgrid(size(rf), z(end),param);
    fs = sampleFreq;
end

%% Aux Functions

% Dynamic Focusing Delays
function rx_delays = getRXDelays(element_pos_x, element_pos_z, phi, focus, n_elements, n_pulses, sound_speed)
    rx_delays = zeros(length(focus), n_elements, n_pulses); % length(focus) is same as length(t)
    for n = 1:n_pulses
        xfocus = element_pos_x(n) + sin(phi(n)) * focus;
        zfocus = element_pos_z(n) + cos(phi(n)) * focus;
        for e = 1:n_elements
            rx_delays(:,e,n) = (focus + ...
                sqrt((zfocus- element_pos_z(e)).^2 + (xfocus - element_pos_x(e)).^2))/sound_speed;
        end
    end
end

% Dynamic aperture matrix
function  apodAperture = getApodCoeff(maxAperSize, nElements, nPulses)
    apodAperture = zeros(nElements, nPulses);
    for i = 1:nPulses
        aperCenter = i;
        halfAperSize = floor(maxAperSize/2);
        aIndex = -halfAperSize:halfAperSize;
        aperture = aperCenter + aIndex;
        aperture = aperture(aperture>=1);
        aperture = aperture(aperture<=nElements);
        apodAperture(aperture, i) = 1;
    end
end

% Organize data
function rf_channel = getOrderedData(BuffData, Receive, n_samples, n_elements, n_pulses, k)
    rf_channel = zeros(n_samples, n_elements, n_pulses);
    for n = 1:n_pulses  
        % Select RF Channel Data from every Buffer
        rf_channel(:, :, n) = BuffData(Receive(k*(n-1)+1).startSample:Receive(k*(n-1)+1).endSample, :); % RF Data from Buffer
    end
end

% Get beamformed rf
function rf = getRfC52v(rf_channel, rx_samples, apodAperture, n_samples, n_pulses, n_elements, endDepth, tgcr)
    rf_data = zeros(n_samples, n_pulses);
    for n = 1:n_pulses
        % Delaying
        for e = 1:n_elements
            rf_channel(:, e, n) = rf_channel(rx_samples(:, e, n), e, n);
        end
        % Apodization
        rf_channel(:, :, n) = rf_channel(:, :, n) .* apodAperture(n, :);
        % Summing
        rf_data(:, n) = sum(rf_channel(:, :, n), 2);
    end

    % Depth segmentation
    plotDepth = floor(endDepth*2*4);
    rf_data = rf_data(1:plotDepth, :);

    % TGC
    tgcr = tgcr(1:plotDepth);
    rf_data = bsxfun(@times, tgcr', rf_data);

    rf = rf_data;
end

