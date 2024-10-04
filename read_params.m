function [nSamples, nChirps, f_start, frequency_slope, sampling_rate, time_sweep, time_idle] = read_params(file_path)
%   [nSamples, nChirps, f_start, frequency_slope, sampling_rate, time_sweep, time_idle] = READ_PARAMS(file_path)
%
%   Inputs:
%       file_path - String specifying the path to the HDF5 file.
%
%   Outputs:
%       nSamples        - Number of ADC samples (double).
%       nChirps         - Number of chirps (double).
%       f_start         - Start frequency in Hz (double).
%       frequency_slope - Frequency slope constant in Hz/s (double).
%       sampling_rate   - Sampling rate in Hz (double).
%       time_sweep      - Sweep time in seconds (double).
%       time_idle       - Idle time in seconds (double).
%
    % Read parameters directly without loops
    nSamples = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/numAdcSamples'));

    nChirps = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/frameCfg/numChirps'));

    f_start = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/startFreq')) * 1e9; % Convert GHz to Hz

    frequency_slope = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/freqSlopeConst')) * 1e12; % Convert units to Hz/s

    sampling_rate = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/digOutSampleRate')) * 1e3; % Convert kHz to Hz

    time_sweep = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/rampEndTime')) * 1e-6; % Convert microseconds to seconds

    time_idle = double(h5read(file_path, '/Sensors/TI_Radar/Parameters/profileCfg/idleTime')) * 1e-6; % Convert microseconds to seconds

end
