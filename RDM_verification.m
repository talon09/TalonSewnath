%% Read in the data and store 
file_path = 'F:\Radar\Experimental Data\Experiment_stationary_target_1m_1_data.hdf5';
%file_path = 'E:\Radar\Experimental Data\Experiment_stationary_target_2m_1_data.hdf5';
%file_path = 'E:\Radar\Experimental Data\Experiment_stationary_target_1m_pose_data.hdf5';
data = h5read(file_path,"/Sensors/TI_Radar/Data/Frame_90/frame_data");
dims = size(data.r);

%Extract first channel chirp and sample data
chirps_samples_2D_r = data.r;
chirps_samples_2D_i = data.i;
chans = dims(1);
rows = dims(2);
cols = dims(3);

chirps_samples_2D = zeros(chans,rows,cols);
for i = 1:chans
    for j = 1:rows
        for k = 1:cols
            chirps_samples_2D(i,j,k) = chirps_samples_2D_r(i,j,k) + chirps_samples_2D_i(i,j,k)*1j;
        end
    end
end


%% Read in parameters

nSamples = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/numAdcSamples'));
nChirps = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/frameCfg/numChirps'));
nChannels = 12;
c = 3e8;
f_start = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/startFreq')) * 1e9;
frequency_slope = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/freqSlopeConst')) * 1e12;
sampling_rate = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/digOutSampleRate')) * 1e3;
time_sweep = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/rampEndTime')) * 1e-6 ;  
time_idle = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/idleTime')) *1e-6 ;  

T_chirp = time_idle + time_sweep;
bandwidth = (nSamples/sampling_rate)*frequency_slope;

wave_length = c/(f_start + 0.5*bandwidth);

T_data = (1/sampling_rate) * nSamples;

max_range = (c*sampling_rate)/(2*frequency_slope);
range_resolution = c/(2*bandwidth);
max_vel = wave_length / (4 * T_chirp);
vel_res = wave_length / (2 * T_chirp * nChirps);

%% Range Doppler Map

windowed_data = zeros(nChannels, nChirps, nSamples);
range_fft = zeros(nChannels,nChirps,nSamples);
doppler_fft = zeros(nChannels,nChirps,nSamples);
window = hanning(nChirps); 

for i = 1:chans
    windowed_data(i, :, :) = squeeze(chirps_samples_2D(i, :, :)) .* window; 
    range_fft(i, :, :) = fft(squeeze(windowed_data(i, :, :)), nSamples, 2); % FFT along range dimension
    doppler_fft(i, :, :) = (fft(squeeze(range_fft(i, :, :)), nChirps, 1)); % FFT along Doppler dimension

    doppler_fft(i, :, :) = fftshift(squeeze(doppler_fft(i, :, :)), 1);
    doppler_fft(i, :, :) = fliplr(squeeze(doppler_fft(i, :, :)));
end

range_bins = (0:nSamples) * range_resolution;
doppler_bins = (floor(-nChirps/2):floor(nChirps/2)) * vel_res;

% Plot the Range-Doppler Map
figure;
imagesc(range_bins, doppler_bins,20*log10(abs(squeeze(doppler_fft(1, :, :)))));
xlabel('Range (m)');
ylabel('Velocity (m/s)');
title('Range-Doppler Map');
colorbar;
axis xy;  % Flip the y-axis to have positive velocities upwards
colormap(jet);  % Apply the jet colormap

%% CFAR Implementation
cfar_range_guard = 3;
cfar_range_training = 8;
cfar_doppler_guard = 4;
cfar_doppler_training = 10;
PFA = 10e-4;
cfar_threshold = 0.5;
detections = CA_CFAR(squeeze(doppler_fft(1, :, :)), cfar_range_guard, cfar_range_training, cfar_doppler_guard, cfar_doppler_training, PFA, cfar_threshold);
n = sum(sum(detections));

%% Angle of Arrival

azimuth_array = zeros(1,n);
elevation_array = zeros(1,n);

duplicates = []; % Create an array that the azimuth angles related to the same range as each other

count = 1;
for i = 1:nChirps
    for j = 1:nSamples
        if detections(i,j) == 1
            channel = get_channel(i,j,chirps_samples_2D);
            elevation_array(count) = elevation_estimation(channel);
            
            azimuth_vector = azimuth_estimation(channel); % Assign all the azimuth angles in this temporary vector
            p = length(azimuth_vector);
            if p > 1
                for k = 2:p % This loop stores every angle after the first into the matrix of dupplicate angles 
                    new_row = [range_bins(j) , azimuth_vector(p), elevation_array(count)];
                    duplicates = [duplicates; new_row]; 
                end
            end

            azimuth_array(count) = azimuth_vector(1);
            count = count + 1;
        end
    end
end


%% Spherical to Cartesian 
range_array = zeros(1,n);
count_2 = 1;

for i = 1:nChirps
    for j = 1:nSamples
        if detections(i,j) == 1
            range_array(count_2) = range_bins(j);
            count_2 = count_2 + 1;
        end
    end
end
xyz = zeros(n,3);

for j = 1:n
    xyz(j,1) = range_array(j)*sind(elevation_array(j))*cosd(azimuth_array(j));
    xyz(j,2) = range_array(j)*sind(elevation_array(j))*sind(azimuth_array(j));
    xyz(j,3) = range_array(j)*cosd(elevation_array(j));
end

%% Adding the duplicate ranged angles to the point cloud
dims2 = size(duplicates);

for i = 1:dims2(1)
    x = duplicates(i,1)*sind(duplicates(i,3))*cosd(duplicates(i,2));
    y = duplicates(i,1)*sind(duplicates(i,3))*sind(duplicates(i,2));
    z = duplicates(i,1)*cosd(duplicates(i,3));
    new = [x, y, z];
    xyz = [xyz; new];
end

%% Plot Point Cloud
figure; % Create a new figure window
plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'o', 'LineWidth', 2);
grid on; % Add a grid to the plot
xlabel('X'); % Label the X axis
ylabel('Y'); % Label the Y axis
zlabel('Z'); % Label the Z axis
title('3D Coordinates Plot');
axis equal; % Set axis scaling to be equal for better visual representation
view(3); % Set default 3D view


%% Helper Functions

function [channel] = get_channel(row,col,rdm)
    channel = zeros(1,12);
    channel = transpose(squeeze(rdm(:,row,col)));
end

function [azimuth] = azimuth_estimation(channel_array)
    
    % extraction of azimuth channels from virtual antennas
    azimuth_array = zeros(1,8); 
    for i = 1:4
        azimuth_array(i) = channel_array(i);
    end
    for i = 9:12
        azimuth_array(i-4) = channel_array(i);
    end
    
    num_fft_points = 4096;  % Increase FFT points to improve resolution (Zero-padding)
    fft_result = fftshift(fft(azimuth_array, num_fft_points)); 

    angle_axis = linspace(-90, 90, length(fft_result));  % AoA axis from -90 to 90 degrees

    [~, peak_indices] = findpeaks(abs(fft_result), 'MinPeakHeight', max(abs(fft_result)) * 0.7);

    % Get the corresponding angles for the peaks
    azimuth = angle_axis(peak_indices);

end

function [elevation] = elevation_estimation(channel_array)

    elevation_array = zeros(1,2);
    
    elevation_array(1) = channel_array(8);
    elevation_array(2) = channel_array(12);

    num_fft_points = 4096;  % Increase FFT points to improve resolution (Zero-padding)
    fft_result = fftshift(fft(elevation_array, num_fft_points)); 

    angle_axis = linspace(-90, 90, length(fft_result));  % AoA axis from -90 to 90 degrees

    % Find the peak in the FFT result to estimate the AoA
    [~, max_idx] = max(abs(fft_result));  % Index of the maximum value
    elevation = angle_axis(max_idx);  % Corresponding AoA value

end