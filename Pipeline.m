filepath = string(pwd);
files = dir;



for k = 4:99
    disp(files(k).name);
    target = files(k).name;
    %[~,ind]=sort({files.name});
    file_path_2 = strcat(filepath, '\', target); % Accesses the exact file
    [~,name, ~] = fileparts(target);

    %% Read In Parameters
    nChannels = 12;
    c = 3e8;
    
    [nSamples, nChirps, f_start, frequency_slope, sampling_rate, time_sweep, time_idle] = read_params(file_path_2);
    
    T_chirp = time_idle + time_sweep; % Chirp period
    bandwidth = (nSamples/sampling_rate)*frequency_slope; % Bandwidth of the system
    wave_length = c/(f_start + 0.5*bandwidth);  % Wavelength of the signal
    T_data = (1/sampling_rate) * nSamples;  % Total sampling time
    
    max_range = (c*sampling_rate)/(2*frequency_slope);
    range_resolution = c/(2*bandwidth);
    max_vel = wave_length / (4 * T_chirp);
    vel_res = wave_length / (2 * T_chirp * nChirps);
    
    rdm = zeros(9, nChannels, nChirps, nSamples);
    RDM_sum = zeros(nChirps, nSamples);
    
    range_bins = (0:nSamples-1) * range_resolution;
    doppler_bins = (floor(-nChirps/2):floor(nChirps/2)) * (1/T_chirp);
    
    
    for i = 11:19
        %% Read In Data One Frame Per Loop
        location = strcat("/Sensors/TI_Radar/Data/Frame_", string(i),"/frame_data");
        data = h5read(file_path_2, location);
        dims = size(data.r);
        
        chirps_samples_2D_r = data.r;
        chirps_samples_2D_i = data.i;
        chans = dims(1);
        rows = dims(2);
        cols = dims(3);
        chirps_samples_2D = chirps_samples_2D_r + 1j * chirps_samples_2D_i;
        
        %% Range Doppler
        windowed_data = complex(zeros(nChannels, nChirps, nSamples));
        range_fft = complex(zeros(nChannels, nChirps, nSamples));
        doppler_fft = complex(zeros(nChannels, nChirps, nSamples));
        window = hanning(nChirps); 
        
        for j = 1:chans
            windowed_data(j, :, :) = squeeze(chirps_samples_2D(j, :, :)) .* window; 
            range_fft(j, :, :) = fft(squeeze(windowed_data(j, :, :)), nSamples, 2); % FFT along range dimension
            doppler_fft(j, :, :) = (fft(squeeze(range_fft(j, :, :)), nChirps, 1)); % FFT along Doppler dimension
        
            doppler_fft(j, :, :) = fftshift(squeeze(doppler_fft(j, :, :)), 1);
            doppler_fft(j, :, :) = fliplr(squeeze(doppler_fft(j, :, :)));
        end
    
        rdm(i-10, : , : , : ) = doppler_fft; 
    
    
        for a = 1:nChannels
            RDM_sum = abs(squeeze(squeeze(doppler_fft(a, :, :)))) + RDM_sum;
        end
    
    end
    
    
    %{
    % Plot the Range-Doppler Map
    figure;
    imagesc(range_bins, doppler_bins, 20*log10(abs(squeeze(squeeze(rdm(2, 1, :, :))))));
    xlabel('Range (m)');
    ylabel('Doppler Frequency (Hz)');
    title('Range-Doppler Map');
    colorbar;
    axis xy;  % Flip the y-axis to have positive Doppler upwards
    colormap(jet);  % Apply the jet colormap
    
    figure;
    imagesc(range_bins, doppler_bins, 20*log10(RDM_sum));
    xlabel('Range (m)');
    ylabel('Doppler Frequency (Hz)');
    title('Range-Doppler Map');
    colorbar;
    axis xy;  % Flip the y-axis to have positive Doppler upwards
    colormap(jet);  % Apply the jet colormap
    %}
    
    
    cfar_range_guard = 1;
    cfar_range_training = 1;
    cfar_doppler_guard = 2;
    cfar_doppler_training = 3;
    PFA = 10e-4;
    cfar_threshold = 1;
    detections = TR_CA_CFAR(RDM_sum, cfar_range_guard, cfar_range_training, cfar_doppler_guard, cfar_doppler_training, PFA, cfar_threshold);
    n = sum(sum(detections));
    
    %% Angle of Arrival
    
    azimuth_array = zeros(1,n);
    elevation_array = zeros(1,n);
    
    duplicates = []; % Create an array that the azimuth angles related to the same range as each other
    
    rho_theta_phi_I = [];
    
    for i = 1:nChirps
        for j = 1:nSamples
            if detections(i,j) == 1
                el_samples = get_channel_el(i,j,rdm); % 2 channels and 9 samples are retrieved
                el_angles = ph_beam_former_el(el_samples);
                az_samples = get_channel_az(i,j,rdm); % 8 channels and 9 samples are retrieved
                az_angles = ph_beam_former_az(az_samples);
                
                
                for z = 1:10
                    rho = range_bins(j);
                    theta = az_angles(z);
                    phi = el_angles(z);
                    rtp = [rho, theta, phi, RDM_sum(i,j)];
                    rho_theta_phi_I = [rho_theta_phi_I; rtp];
                end
            end
        end
    end
    
    stp = n*10;
    xyz = zeros(stp,3);
    count = 0;
    xyz2 = [];
    
    for t = 1:stp
        xyz(t,1) = rho_theta_phi_I(t,1) * sind(rho_theta_phi_I(t,3)) * cosd(rho_theta_phi_I(t,2));
        xyz(t,2) = rho_theta_phi_I(t,1) * sind(rho_theta_phi_I(t,3)) * sind(rho_theta_phi_I(t,2));
        xyz(t,3) = rho_theta_phi_I(t,1) * cosd(rho_theta_phi_I(t,3));
        
        if xyz(t,3) > 1 && xyz(t,3) < 3
            count = count + 1;
            temp = [xyz(t,2), xyz(t,1), xyz(t,3), rho_theta_phi_I(t,4)]; % width by height by depth
            xyz2 = [xyz2;temp];
        end
    end
    
    %% Plot Point Cloud
    figure; % Create a new figure window
    plot3(xyz2(:,3), xyz2(:,1), xyz2(:,2), 'o', 'LineWidth', 2); % depth by height by width 
    grid on; % Add a grid to the plot
    xlabel('Z'); 
    ylabel('X'); 
    zlabel('Y'); 
    title('3D Coordinates Plot');
    axis equal; % Set axis scaling to be equal for better visual representation
    view(3); % Set default 3D view
    
    reshaped_map = map_reshape(xyz2);
    save_path = strcat(filepath,'/Data_Files/',name,'.mat');
    save(save_path,'reshaped_map');
end

%% Heper Functions
function [el_samples] = get_channel_el(row,col,rdm)
    el_samples = zeros(2,9);
    el_samples(1,:) = transpose(squeeze(rdm(:,8,row,col)));
    el_samples(2,:) = transpose(squeeze(rdm(:,12,row,col)));
end

function [az_samples] = get_channel_az(row,col,rdm)
    az_samples = zeros(8,9);
    az_samples(1,:) = transpose(squeeze(rdm(:,1,row,col)));
    az_samples(2,:) = transpose(squeeze(rdm(:,2,row,col)));
    az_samples(3,:) = transpose(squeeze(rdm(:,3,row,col)));
    az_samples(4,:) = transpose(squeeze(rdm(:,4,row,col)));
    az_samples(5,:) = transpose(squeeze(rdm(:,9,row,col)));
    az_samples(6,:) = transpose(squeeze(rdm(:,10,row,col)));
    az_samples(7,:) = transpose(squeeze(rdm(:,11,row,col)));
    az_samples(8,:) = transpose(squeeze(rdm(:,12,row,col)));
end



