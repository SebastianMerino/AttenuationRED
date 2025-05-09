%
% File name: SetUpC5_2v_128RyLns.m - Example of curved array imaging with
%                                    focus transmits
% Description:
%   Sequence programming file for C5-2v curved array using 128 ray lines
%   (focus transmits) and receive acquisitions. Of the 128 transmit channels
%   only 65 are used, with 32 transmitters on each side of the center
%   element (where possible). All 128 receive channels are active for each
%   acquisition. The receive acquisitions use 100% bandwidth to improve DMA
%   transfers. Processing is asynchronous with respect to acquisition.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
% 
% Copyright © 2013-2023 Verasonics, Inc.

clear all;

P.numRays = 128;  % no. of raylines to program
P.numTx = 14;
P.txFocus = 100;  % focal point in wavelength
P.startDepth = 0; % startDepth and endDepth in wavelength % def 5
P.endDepth = 400; % 192 - default
%->% WO modifying buffer size (RowsperFrame) tops at 418 wl ~~ 18 cm at 3.70 MHz (418)
%->% At 3.70MHz 602 ~~ to 25cm, no simulation. Default: 192

% Specify system parameters.
Resource.Parameters.numTransmit = 128;  % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 1;

% Specify Trans structure array.
Trans.name = 'C5-2v';
Trans.units = 'wavelengths';  
Trans = computeTrans(Trans);  % C5-2v transducer is a 'known' transducer so we can use computeTrans.
%Trans.frequency = 5.10; % Uncomment to modify central frequency (default: 3.57 MHz)
Trans.maxHighVoltage = 80;  % set maximum high voltage limit for pulser supply.
% Limits for voltage [1.6 - 96] V (Not recommended to use max. hardware voltage [DAMAGE RISK])
radius = Trans.radius;
scanangle = Trans.numelements*Trans.spacing/radius;
dtheta = scanangle/P.numRays;
theta = -(scanangle/2) + 0.5*dtheta; % angle to left edge from centerline
Angle = theta:dtheta:(-theta);

% Specify PData structure array.
PData(1).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
sizeRows = 10 + ceil((P.endDepth + radius - (radius * cos(scanangle/2)))/PData(1).PDelta(3));
sizeCols = 10 + ceil(2*(P.endDepth + radius)*sin(scanangle/2)/PData(1).PDelta(1));
PData(1).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
PData(1).Origin(1,1) = (P.endDepth+radius)*sin(-scanangle/2) - 5;
PData(1).Origin(1,2) = 0;
PData(1).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
% Define PData Regions for numRays scanlines
for j = 1:P.numRays
    PData(1).Region(j) = struct(...
        'Shape',struct('Name','Sector',...
                       'Position',[0,0,-radius],...
                       'r1',radius+P.startDepth,...
                       'r2',radius+P.endDepth,...
                       'angle',dtheta,...
                       'steer',Angle(j)));
end
PData(1).Region = computeRegions(PData(1));

%  Media points for curved array.
% - Uncomment for speckle
% Media.numPoints = (20000);
% Media.MP = rand(Media.numPoints,4);
% Media.MP(:,2) = 0;
% Media.MP(:,4) = 0.01 + 0.04*Media.MP(:,4);  % Random low amplitude
% RandR = endDepthMm *Media.MP(:,1)+radius;
% RandTheta = scanangle*(Media.MP(:,3)-0.5);
% Media.MP(:,1) = RandR.*sin(RandTheta);
% Media.MP(:,3) = RandR.*cos(RandTheta)-radius;
wvl = Resource.Parameters.speedOfSound/(Trans.frequency*1e6); % [m]
scalemm2wvl = 1/wvl;

% - Define points
%Media.MP(1,:) = [0,0,70,1.0];

% Media.MP(1,:) = [0,0,30,1.0]; % 30  
% Media.MP(2,:) = [0,0,60,1.0]; % 60 
% Media.MP(3,:) = [0,0,90,1.0]; % 90 

Media.MP(1,:) = [0,0,10*scalemm2wvl*1e-3,1.0]; % [mm] 
Media.MP(2,:) = [0,0,20*scalemm2wvl*1e-3,1.0]; % [mm] 
Media.MP(3,:) = [0,0,30*scalemm2wvl*1e-3,1.0]; % [mm] 

%Media.MP(1,:) = [0,0,50,1.0];
% Media.MP(2,:) = [(radius+10)*sin(-0.2608),0,(radius+10)*cos(-0.2608)-radius,1.0];
% Media.MP(3,:) = [(radius+10)*sin(0.2608),0,(radius+10)*cos(0.2608)-radius,1.0];
% Media.MP(4,:) = [(radius+10)*sin(-0.5267),0,(radius+10)*cos(-0.5267)-radius,1.0];
% Media.MP(5,:) = [(radius+10)*sin(0.5267),0,(radius+10)*cos(0.5267)-radius,1.0];
% Media.MP(6,:) = [0,0,40,1.0];
% Media.MP(7,:) = [0,0,70,1.0];
% Media.MP(8,:) = [(radius+70)*sin(-0.2608),0,(radius+70)*cos(-0.2608)-radius,1.0];
% Media.MP(9,:) = [(radius+70)*sin(0.2608),0,(radius+70)*cos(0.2608)-radius,1.0];
% Media.MP(10,:) = [(radius+70)*sin(-0.5267),0,(radius+70)*cos(-0.5267)-radius,1.0];
% Media.MP(11,:) = [(radius+70)*sin(0.5267),0,(radius+70)*cos(0.5267)-radius,1.0];
% Media.MP(12,:) = [0,0,100,1.0];
% Media.MP(13,:) = [0,0,130,1.0];
% Media.MP(14,:) = [(radius+130)*sin(-0.2608),0,(radius+130)*cos(-0.2608)-radius,1.0];
% Media.MP(15,:) = [(radius+130)*sin(0.2608),0,(radius+130)*cos(0.2608)-radius,1.0];
% Media.MP(16,:) = [(radius+130)*sin(-0.5267),0,(radius+130)*cos(-0.5267)-radius,1.0];
% Media.MP(17,:) = [(radius+130)*sin(0.5267),0,(radius+130)*cos(0.5267)-radius,1.0];
% Media.MP(18,:) = [0,0,160,1.0];
% Media.MP(19,:) = [0,0,190,1.0];
% Media.MP(20,:) = [(radius+190)*sin(-0.2608),0,(radius+190)*cos(-0.2608)-radius,1.0];
% Media.MP(21,:) = [(radius+190)*sin(0.2608),0,(radius+190)*cos(0.2608)-radius,1.0];
% Media.MP(22,:) = [(radius+190)*sin(-0.5267),0,(radius+190)*cos(-0.5267)-radius,1.0];
% Media.MP(23,:) = [(radius+190)*sin(0.5267),0,(radius+190)*cos(0.5267)-radius,1.0];
% Media.MP(24,:) = [0, 0, 380, 1.0];
% Media.MP(25,:) = [(radius+380)*sin(-0.2608),0,(radius+380)*cos(-0.2608)-radius,1.0];
% Media.MP(26,:) = [(radius+380)*sin(0.2608),0,(radius+380)*cos(0.2608)-radius,1.0];
% Media.MP(27,:) = [(radius+380)*sin(-0.5267),0,(radius+380)*cos(-0.5267)-radius,1.0];
% Media.MP(28,:) = [(radius+380)*sin(0.5267),0,(radius+380)*cos(0.5267)-radius,1.0];
% Media.MP(29,:) = [0, 0, 580, 1.0];
% Media.MP(30,:) = [(radius+580)*sin(-0.2608),0,(radius+580)*cos(-0.2608)-radius,1.0];
% Media.MP(31,:) = [(radius+580)*sin(0.2608),0,(radius+580)*cos(0.2608)-radius,1.0];
% Media.MP(32,:) = [(radius+580)*sin(-0.5267),0,(radius+580)*cos(-0.5267)-radius,1.0];
% Media.MP(33,:) = [(radius+580)*sin(0.5267),0,(radius+580)*cos(0.5267)-radius,1.0];
Media.attenuation = -0.5;
Media.function = 'movePoints';

% Specify Resources.
% Continuous acquisition
num_samples = 4096;
%->% Modified to allow more than 512 wavelenghts -> 4096/2/4
%->for 25 cm (4096+2048);
numFrames = 10;
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*num_samples;  % 128 segments of 4096 samples
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = numFrames;
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = numFrames;
Resource.DisplayWindow(1).Title = 'C5-2v_128RyLns';
Resource.DisplayWindow(1).pdelta = 0.45;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = numFrames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow.Colormap = gray(256);


% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

% Specify TX structure array.
% - We need 128 transmit specifications.
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements)), 1, P.numRays);
% - Set event specific TX attributes.
%    numTx/2 is the number of elements to include on each side of the
%    center element, for the specified focus and sensitivity cutoff.
%    Thus the full transmit aperture will be numTx + 1 elements.
for n = 1:P.numRays   % numRays transmit events
    TX(n).waveform = 1;  % Set transmit waveform
    % Set transmit Origins.
    TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
    ce = round(1+127*(Angle(n) - theta)/(-2*theta));
    % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > Trans.numelements, rt = Trans.numelements; end
    TX(n).Apod(lft:rt) = 1.0;
    TX(n).Delay = computeTXDelays(TX(n));
end

% Specify Receive structure arrays.
% - We need numRays Receives for each frame.
% -- Compute the maximum receive path length, using the law of cosines.
maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, P.numRays*Resource.RcvBuffer(1).numFrames);
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    Receive(P.numRays*(i-1)+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(P.numRays*(i-1)+j).framenum = i;
        Receive(P.numRays*(i-1)+j).acqNum = j;
    end
end

% Specify TGC Waveform structure.
TGC.CntrlPts = [153,308,410,520,605,665,705,760];
TGC.rangeMax = P.endDepth;
TGC.Waveform = computeTGCWaveform(TGC);

% Specify Recon structure arrays.
Recon = struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame',-1, ...
               'ImgBufDest', [1,-1], ...
               'RINums', 1:P.numRays);

% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'replaceIntensity', ...  % replace data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, P.numRays);
% - Set specific ReconInfo attributes.
for i = 1:P.numRays
    ReconInfo(i).txnum = i;
    ReconInfo(i).rcvnum = i;
    ReconInfo(i).regionnum = i;
end

% Specify Process structure array.
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',1.0,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};

% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; %  jump back to start.
SeqControl(1).argument = 1;
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 250; % 250 usec/ray * 128 = 32msec
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 1000; % 1 msec
SeqControl(4).command = 'returnToMatlab';
nsc = 5;

record_check = 0;
% Specify Event structure arrays.
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numRays  % Acquire rays
        Event(n).info = 'Acquire ray line';
        Event(n).tx = j;
        Event(n).rcv = P.numRays*(i-1)+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = 3; % time between frames
    Event(n).info = 'Transfer frame to host.';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 0;
    Event(n).seqControl = nsc;
       SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
       nsc = nsc+1;
    n = n+1;

    Event(n).info = 'recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if (floor(i/5) == i/5)&&(i ~= Resource.RcvBuffer(1).numFrames)  % Exit to Matlab every 5th frame
        Event(n).seqControl = 4;
    end
    n = n+1;
end

Event(n).info = 'Jump back';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;

% User specified UI Control Elements

import vsv.seq.uicontrol.VsSliderControl
import vsv.seq.uicontrol.VsButtonControl
import vsv.seq.uicontrol.VsButtonGroupControl

% - Sensitivity Cutoff
UI(1).Control = VsSliderControl('LocationCode','UserB7',...
    'Label','Sens. Cutoff',...
    'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
    'SliderStep',[0.025,0.1],'ValueFormat','%1.3f',...
    'Callback',@SensCutoffCallback);

% - Range Change
wls2mm = 1;
AxesUnit = 'wls';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'mm';
        wls2mm = Resource.Parameters.speedOfSound/1000/Trans.frequency;
    end
end
UI(2).Control = VsSliderControl('LocationCode','UserA1',...
    'Label',['Range (',AxesUnit,')'],...
    'SliderMinMaxVal',[64,P.endDepth,P.endDepth]*wls2mm,...
    'SliderStep',[0.1,0.2],'ValueFormat','%3.0f',...
    'Callback',@RangeChangeCallback);

% - Transmit focus change
UI(3).Control = VsSliderControl('LocationCode','UserB4',...
    'Label',['TX Focus (',AxesUnit,')'],...
    'SliderMinMaxVal',[50,P.endDepth,P.txFocus]*wls2mm,...
    'SliderStep',[10/250,20/250],'ValueFormat','%3.0f',...
    'Callback',@TxFocusCallback);

% - F number change
UI(4).Control = VsSliderControl('LocationCode','UserB3',...
    'Label','F Number',...
    'SliderMinMaxVal',[0.8,20,round(P.txFocus/(P.numTx*Trans.spacing))],...
    'SliderStep',[0.05,0.1],'ValueFormat','%2.1f',...
    'Callback',@FNumCallback);

% - Recorded frame change
n_select = Resource.ImageBuffer(1).numFrames;
UI(5).Control = VsSliderControl('LocationCode', 'UserB2', ...
    'Label', 'Recorded Frame', ...
    'SliderMinMaxVal',[1, n_select, n_select],...
    'SliderStep',[1/(n_select-1), 1/(n_select-1)],'ValueFormat','%2.0f',...
    'Callback', @RecordedFrameCallback);

% - Save Frame RF
UI(6).Control = VsButtonControl('LocationCode', 'UserC2', ...
    'Label', 'Save Frame', ...
    'Callback', @SaveFrameRFCallback);

% - Save Buffer RF
UI(7).Control = VsButtonControl('LocationCode', 'UserC1', ...
    'Label', 'Save Buffer', ...
    'Callback', @SaveBufferRFCallback);

% - Plane Acquisition Selector
acq_plane = 'LP';
UI(8).Control = VsButtonGroupControl('LocationCode', 'UserC4', ...
    'Title', 'Acquisition Plane', ...
    'PossibleCases', {'Longitudinal plane', 'Transverse plane', 'Subcostal oblique plane', ...
                    'Free Plane 1', 'Free Plane 2'},...
    'Callback', @PlaneAcquisitionCallback);

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 5;

% Save all the structures to a .mat file.
save('MatFiles/C5-2v_128RyLns_mod_csg_patient_data');

%% **** Callback routines used by UIControls (UI) ****

function SensCutoffCallback(~, ~, UIValue)
    ReconL = evalin('base', 'Recon');
    for i = 1:size(ReconL,2)
        ReconL(i).senscutoff = UIValue;
    end
    assignin('base','Recon',ReconL);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'Recon'};
    assignin('base','Control', Control);
end

function RangeChangeCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No range change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.endDepth'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

    P = evalin('base','P');
    P.endDepth = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.endDepth = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    scanangle = evalin('base','scanangle');
    radius = evalin('base','radius');
    theta = evalin('base','theta');
    dtheta = -2*theta/P.numRays;
    Angle = theta:dtheta:(-theta);
    height = P.endDepth + radius - (radius*cos(scanangle/2));
    % Modify PData for new range
    PData = evalin('base','PData');
    PData(1).Size(1) = 10+ceil(height/PData(1).PDelta(3));
    for j = 1:P.numRays
        PData(1).Region(j).Shape = struct( ...
                           'Name','Sector',...
                           'Position',[0,0,-radius],...
                           'r1',radius+P.startDepth,...
                           'r2',radius+P.endDepth,...
                           'angle',dtheta,...
                           'steer',Angle(j));
        PData(1).Region(j).PixelsLA = [];
        PData(1).Region(j).numPixels = 0;
    end
    assignin('base','PData',PData);
    evalin('base','PData(1).Region = computeRegions(PData(1));');
    evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
    Receive = evalin('base', 'Receive');
    maxAcqLength = ceil(sqrt((P.endDepth+radius)^2 + radius^2 - 2*(P.endDepth+radius)*radius*cos(scanangle)));
    for i = 1:size(Receive,2)
        Receive(i).endDepth = maxAcqLength;
    end
    assignin('base','Receive',Receive);
    evalin('base','TGC.rangeMax = P.endDepth;');
    evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'PData','InterBuffer','ImageBuffer','Receive','TGC','Recon','DisplayWindow'};
    assignin('base','Control', Control);
    assignin('base', 'action', 'displayChange');
end

function TxFocusCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    % No focus change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',evalin('base','P.txFocus'));
        return
    end
    Trans = evalin('base','Trans');
    Resource = evalin('base','Resource');
    scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

    P = evalin('base','P');
    P.txFocus = UIValue;
    if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
        if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
            P.txFocus = UIValue*scaleToWvl;
        end
    end
    assignin('base','P',P);

    TX = evalin('base', 'TX');
    for n = 1:128   % 128 transmit events
        TX(n).focus = P.txFocus;
        TX(n).Delay = computeTXDelays(TX(n));
    end
    assignin('base','TX', TX);
    % Update Fnumber based on new txFocus
    evalin('base','set(UI(4).handle(2),''Value'',round(P.txFocus/(P.numTx*Trans.spacing)));');
    evalin('base','set(UI(4).handle(3),''String'',num2str(round(P.txFocus/(P.numTx*Trans.spacing)),''%2.1f''));');
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function FNumCallback(hObject, ~, UIValue)
    simMode = evalin('base','Resource.Parameters.simulateMode');
    P = evalin('base','P');
    Trans = evalin('base','Trans');
    % No F number change if in simulate mode 2.
    if simMode == 2
        set(hObject,'Value',round(P.txFocus/(P.numTx*Trans.spacing)));
        return
    end
    txFNum = UIValue;
    P.numTx = round(P.txFocus/(txFNum*Trans.spacing));
    radius = Trans.radius;
    scanangle = Trans.numelements*Trans.spacing/radius;
    dtheta = scanangle/P.numRays;
    theta = -(scanangle/2) + 0.5*dtheta; % angle to left edge from centerline
    Angle = theta:dtheta:(-theta);
    % - Redefine event specific TX attributes for the new numTx.
    TX = evalin('base', 'TX');
    for n = 1:P.numRays  % numRays transmit events
        TX(n).waveform = 1;
        % Set transmit Origins.
        TX(n).Origin = [radius*sin(Angle(n)), 0.0, radius*cos(Angle(n))-radius];
        ce = round(1+127*(Angle(n) - theta)/(-2*theta));
        % Set transmit Apodization so that a maximum of numTx + 1 transmitters are active.
        lft = round(ce - P.numTx/2);
        if lft < 1, lft = 1; end
        rt = round(ce + P.numTx/2);
        if rt > Trans.numelements, rt = Trans.numelements; end
        TX(n).Apod = zeros(1,128);
        TX(n).Apod(lft:rt) = 1.0;
        TX(n).Delay = computeTXDelays(TX(n));
    end
    assignin('base','TX', TX);
    Control = evalin('base','Control');
    Control.Command = 'update&Run';
    Control.Parameters = {'TX'};
    assignin('base','Control', Control);
end

function RecordedFrameCallback(hObject, ~, UIValue)

    if evalin('base','freeze')==0   % no action if not in freeze
        n_default = evalin('base', 'n_select');
        set(hObject, 'Value', n_default);  % reset slider to default
        msgbox('Please freeze VSX');
        return
    end

    frame_n = round(UIValue);
    assignin('base', 'frame_n', frame_n);

    return
end

function SaveFrameRFCallback(hObject, ~)
    if ~isempty(findobj('tag','UI')) % running VSX
        if evalin('base','freeze')==0   % no action if not in freeze
            msgbox('Please freeze VSX');
            return
        else
            Control.Command = 'copyBuffers';
            runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
        end
    else % not running VSX
        if evalin('base','exist(''RcvData'',''var'');')
            RcvData = evalin('base','RcvData');
        else
            disp('RcvData does not exist!');
            return
        end
    end

    if evalin('base', 'exist(''frame_n'', ''var'')')
        frame_n = evalin('base', 'frame_n');
    else
        frame_n = evalin('base', 'n_select');
    end
    RcvData = cell2mat(RcvData);
    RcvData = RcvData(:,:,frame_n);
    RcvData = {RcvData};

    acq_plane = evalin('base', 'acq_plane');
    patientID = evalin('base', 'patientID');
    
    RFfilename = [patientID, '_', acq_plane, '_F', ...
        num2str(frame_n)];
        
    [fn,pn,~] = uiputfile('*.mat','Save RF data as',RFfilename);
    if ~isequal(fn,0) % fn will be zero if user hits cancel
        fn = strrep(fullfile(pn,fn), '''', '''''');
        save(fn,'RcvData','-v7.3');
        fprintf('The RF data has been saved at %s \n',fn);
        cmd = ['save ' fn ' Trans ' 'TX ' 'P ' 'Resource ' 'Receive ' 'TGC ' ...
            'PData ' 'TPC ' '-append'];
        evalin('base',cmd)
        disp(['Additional data saved to ' fn]);
    else
        disp('The RF data is not saved.');
    end

    return

end

function SaveBufferRFCallback(hObject, ~)
    if ~isempty(findobj('tag','UI')) % running VSX
        if evalin('base','freeze')==0   % no action if not in freeze
            msgbox('Please freeze VSX');
            return
        else
            Control.Command = 'copyBuffers';
            runAcq(Control); % NOTE:  If runAcq() has an error, it reports it then exits MATLAB.
        end
    else % not running VSX
        if evalin('base','exist(''RcvData'',''var'');')
            RcvData = evalin('base','RcvData');
        else
            disp('RcvData does not exist!');
            return
        end
    end

    acq_plane = evalin('base', 'acq_plane');
    patientID = evalin('base', 'patientID');
    
    RFfilename = [patientID, '_', acq_plane, '_B'];
    
    RcvLastFrame = size(RcvData,3);
    if (~evalin('base','simButton'))
        RcvLastFrame = Resource.RcvBuffer(1).lastFrame;
    end
    
    [fn,pn,~] = uiputfile('*.mat','Save RF data as',RFfilename);
    if ~isequal(fn,0) % fn will be zero if user hits cancel
        fn = strrep(fullfile(pn,fn), '''', '''''');
        save(fn,'RcvData','RcvLastFrame','-v7.3');
        fprintf('The RF data has been saved at %s \n',fn);
        cmd = ['save ' fn ' Trans ' 'TX ' 'P ' 'Resource ' 'Receive ' 'TGC ' ...
            'PData ' 'TPC ' '-append'];
        evalin('base',cmd)
        disp(['Additional data saved to ' fn]);
    else
        disp('The RF data is not saved.');
    end

    return

end

function PlaneAcquisitionCallback(hObject, ~, UIValue)
    if evalin('base', 'exist(''acq_plane'', ''var'')')
        acq_plane = evalin('base', 'acq_plane');
    else
        acq_plane = 'LP';
    end

    plane_n = UIValue;
    switch plane_n
        case 1
            acq_plane = 'LP';
        case 2
            acq_plane = 'TP';
        case 3
            acq_plane = 'SOP';
        case 4
            acq_plane = 'FP1';
        case 5
            acq_plane = 'FP2';
        otherwise
            warning('Unexpected plane selected. No plane.');
    end
    assignin('base', 'acq_plane', acq_plane);
end