
% This example simulates liver and abdominall wall
% Created on April 24th, 2024
% Author: Sebastian Merino

simuNames = {'inc4'}; 
rng(23)
addpath(genpath(pwd))

% medium parameters
c0              = 1540;     % sound speed [m/s]
rho0            = 1000;     % density [kg/m^3]

% source parameters
source_f0       = 6.66e6;   % source frequency [Hz]
source_amp      = 1e6;      % source pressure [Pa]
source_cycles   = 3.5;      % number of toneburst cycles
source_focus    = 4e-2;   % focal length [m]
element_pitch   = 0.3e-3;   % pitch [m]
element_width   = 0.25e-3;  % width [m]
focal_number_Tx = 4;
focal_number_Rx = 2;
nLines          = 128;       % Number of beams;

% grid parameters
grid_size_x     = 5.6e-2;    % [m]
grid_size_y     = 4e-2;     % [m]

% transducer position
translation     = [-2.7e-2, 0];
rotation        = 0;

% computational parameters
DATA_CAST       = 'gpuArray-single'; % set to 'single' or 'gpuArray-single'
ppw             = 6;        % number of points per wavelength, 4 to 8
depth           = 5.5e-2;   % imaging depth [m]
cfl             = 0.3;      % CFL number, could be 0.3 or 0.5
PMLSize         = [41,41];   % 

%% For looping simulations

for iSim = 1:length(simuNames)
    %% GRID

    % calculate the grid spacing based on the PPW and F0
    dx = c0 / (ppw * source_f0);   % [m]

    % compute the size of the grid
    Nx = roundEven(grid_size_x / dx);
    Ny = roundEven(grid_size_y / dx);

    % create the computational grid
    kgrid = kWaveGrid(Nx, dx, Ny, dx);

    % create the time array
    t_end           = depth*2/c0;     % [s];    % total compute time [s]
    kgrid.makeTime(c0, cfl, t_end);

    %% MEDIUM
    rz = kgrid.x - translation(1);
    rx = kgrid.y;

    background_std = 0.002;
    inc_std = 0.008;
    background_alpha = 0.5;       % [dB/(MHz^y cm)]
    inc_alpha = 1;

    % Background
    medium = addRegionSimu([],c0,rho0,background_std,...
        background_alpha,ones(Nx,Ny));

    % Inclusion
    cz = 20e-3; cx = 0;
    if iSim < 3
        r = 0.9e-2;
    else
        r = 1e-2;
    end
    maskNodule = (rz-cz).^2 + (rx-cx).^2 < r^2;
    medium = addRegionSimu(medium,c0,rho0,inc_std,...
        inc_alpha,maskNodule);

    medium.alpha_power = 1;
    medium.alpha_mode = 'no_dispersion';
    medium.sound_speed_ref = 1540;

    figure('Units','centimeters', 'Position',[5 5 25 10]),
    tiledlayout(1,3),
    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.sound_speed)
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Sound speed')
    c = colorbar; ylabel(c,'m/s')
    axis image

    nexttile,
    imagesc(100*rx(1,:),100*rz(:,1),medium.density)
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Density')
    c = colorbar; ylabel(c,'kg/m^3')
    colorbar,
    axis image

    t3 = nexttile;
    imagesc(100*rx(1,:),100*rz(:,1),medium.alpha_coeff, [0.4 1.7])
    xlabel('x [cm]'), ylabel('z [cm]')
    title('Absorption')
    c = colorbar; ylabel(c,'dB/cm/MHz')
    axis image
    colormap(t3,"turbo")

    %% SOURCE
    aperture_Rx = source_focus/focal_number_Rx;
    aperture_Tx = source_focus/focal_number_Tx;
    element_num_Tx = floor(aperture_Tx/element_pitch);
    element_num = floor(aperture_Rx/element_pitch);
    
    amp_vector = source_amp*ones(element_num,1);
    no_Tx_elements = floor((element_num - element_num_Tx)/2);
    amp_vector(1:no_Tx_elements) = 0;
    amp_vector(end-no_Tx_elements+1:end) = 0;

    % set indices for each element
    ids = (0:element_num-1) - (element_num-1)/2;

    % set time delays for each element to focus at source_focus
    time_delays = -(sqrt((ids .* element_pitch).^2 + source_focus.^2) - source_focus) ./ c0;
    time_delays = time_delays - min(time_delays);

    % create time varying source signals (one for each physical element)
    source_sig = amp_vector .* toneBurst(1/kgrid.dt, source_f0, ...
        source_cycles, 'SignalOffset', round(time_delays / kgrid.dt));

    % create empty kWaveArray
    karray = kWaveArray('BLITolerance', 0.05, 'UpsamplingRate', 10);

    % add rectangular elements
    for ind = 1:element_num

        % set element y position
        y_pos = 0 - (element_num * element_pitch/2 - element_pitch/2) + (ind-1) * element_pitch;

        % add element (set rotation angle to match the global rotation angle)
        karray.addRectElement([0, y_pos], element_width/4, element_width, rotation);
    end

    %% For Looping
    yCords = ( (0:nLines-1) - (nLines-1)/2 )* element_pitch; % Lateral cord of each element
    rf_prebf = zeros(kgrid.Nt ,element_num, nLines);

    for iLine = 1:nLines
        %% SOURCE
        % move the array
        translation(2) = yCords(iLine);
        karray.setArrayPosition(translation, rotation)

        % assign binary mask
        source.p_mask = karray.getArrayBinaryMask(kgrid);

        % assign source signals
        source.p = karray.getDistributedSourceSignal(kgrid, source_sig);
        clc
        disp(['Line: ',num2str(iLine),' of ',num2str(nLines)]);

        %% SENSOR

        % set sensor mask to record central plane
        sensor.mask = karray.getArrayBinaryMask(kgrid);
        sensor.directivity_size = 10*kgrid.dx;
        sensor.directivity_angle = zeros(size(sensor.mask));

        %% SIMULATION

        % set input options
        input_args = {...
            'PMLInside', false, ...
            'PMLSize', PMLSize, ...
            'DataCast', DATA_CAST, ...
            'DataRecast', true, ...
            'PlotSim', true};

        % MATLAB CPU
        sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor, ...
            input_args{:});

        % combine sensor data
        combined_sensor_data = karray.combineSensorData(kgrid, sensor_data);

        % beamforming
        fs = 1/kgrid.dt;

        rf_prebf(:,:,iLine) = combined_sensor_data';
    end

    %% VISUALISATION
    offset = 1;

    axAxis = 0:kgrid.Nt-1; axAxis = axAxis*kgrid.dt*c0/2;
    latAxis = yCords;

    z = axAxis(offset:end);
    x = latAxis;
    
    simuName = simuNames{iSim};
    save("rf_prebf_"+simuName+".mat",...
        'rf_prebf','x','z','fs','time_delays');

end
