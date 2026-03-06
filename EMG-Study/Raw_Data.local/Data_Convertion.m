
%% --- Initialization ---
clear all          
clc                

% Files : table with file names/paths
% Data  : cell array, each cell holds a struct loaded from one .mat file
Raw_Data = fullfile('Raw_Data');
addpath(fullfile(pwd, 'Functions'))

[Files, Data] = load_all_mat_structs(Raw_Data);  

%% --- Build EMG_recordings from Data ---
j = 0;  % Subject/record counter (for shared ID across SP/NT pairs)

% Loop over odd indices from 1 to 63 (1:2:64), pairing i and i+1
for i = 1:2:64

    j = j + 1;    

    % -------- SP entry (index i) --------
    EMG_recordings(i).ID                  = j;                              % Shared ID for the pair
    EMG_recordings(i).sex                 = Data{i}.info.sex;              % Sex from source
    EMG_recordings(i).age                 = '--';                          % Placeholder (fill later)
    EMG_recordings(i).posture             = 'SP';                          % Posture/SP condition
    EMG_recordings(i).amplitude_units     = 'uV';                          % EMG amplitude units
    EMG_recordings(i).frequency_units     = 'Hz';                          % Frequency units
    EMG_recordings(i).time_units          = 's';                           % Time units
    EMG_recordings(i).sampling_frequency  = 2000;                         % Sampling rate (Hz)
    EMG_recordings(i).data.Brachioradialis = Data{i}.movements.sources.signals.signal_1.data;
    EMG_recordings(i).data.BicepsBrachi    = Data{i}.movements.sources.signals.signal_2.data;

    % -------- NT entry (index i+1) --------
    EMG_recordings(i+1).ID                  = j;                           % Same ID as SP
    EMG_recordings(i+1).sex                 = Data{i+1}.info.sex;          % Sex from source
    EMG_recordings(i+1).age                 = '--';                        % Placeholder
    EMG_recordings(i+1).posture             = 'NT';                        % Posture/NT condition
    EMG_recordings(i+1).amplitude_units     = 'uV';                        % EMG amplitude units
    EMG_recordings(i+1).frequency_units     = 'Hz';                        % Frequency units
    EMG_recordings(i+1).time_units          = 's';                         % Time units
    EMG_recordings(i+1).sampling_frequency  = 2000;                       % Sampling rate (Hz)
    EMG_recordings(i+1).data.Brachioradialis = Data{i+1}.movements.sources.signals.signal_1.data;
    EMG_recordings(i+1).data.BicepsBrachi    = Data{i+1}.movements.sources.signals.signal_2.data;
end

%% --- Auto-save up to NT-32 with distinct variable names per file ---
% Each file will contain exactly one variable, named EMG_recording_SPxx / EMG_recording_NTxx.
% Requires EMG_recordings to have at least 64 entries (SP-01/NT-01 ... SP-32/NT-32).

assert(numel(EMG_recordings) >= 64, ...
    'EMG_recordings must contain at least 64 entries.');

idx = 1;   % Points to the current EMG_recordings element to save

for k = 1:32
    % -------- Save SP-k --------
    spFile = sprintf('SP-%02d.mat', k);           % Output file name for SP
    spVar  = sprintf('EMG_recording_SP%02d', k);  % Variable name inside the .mat
    payload = EMG_recordings(idx);                % If cell array, use EMG_recordings{idx}
    S = struct();                                  % Wrapper
    S.(spVar) = payload;                           % Dynamic field -> becomes variable name in .mat
    save(spFile, '-struct', 'S', '-v7.3');         % Save as v7.3 (HDF5); remove if not needed
    idx = idx + 1;                                 % Advance to next element

    % -------- Save NT-k --------
    ntFile = sprintf('NT-%02d.mat', k);
    ntVar  = sprintf('EMG_recording_NT%02d', k);
    payload = EMG_recordings(idx);
    S = struct();
    S.(ntVar) = payload;
    save(ntFile, '-struct', 'S', '-v7.3');
    idx = idx + 1;
end

% Also keep a single file with the whole collection for convenience
save("EMG_recordings.mat","EMG_recordings")

