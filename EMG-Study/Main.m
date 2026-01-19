%% ========================================================================
%  EMG Spectral Analysis and Preprocessing Pipeline
%
%  Description:
%  ------------
%  This script performs a complete preprocessing and spectral analysis
%  pipeline for surface EMG (sEMG) signals recorded from two forearm muscles:
%    - Biceps brachii
%    - Brachioradialis
%
%  Signals are analyzed under two experimental conditions:
%    - Neutral forearm position (NT)
%    - Supinated forearm position (SP)
%
%  The pipeline includes:
%    1. Loading EMG datasets for multiple subjects
%    2. Visual inspection of raw EMG signals
%    3. Automatic segmentation of muscle activation events
%    4. Removal of invalid subjects
%    5. Spectral analysis of each EMG segment (FFT-based)
%    6. Spectrum smoothing using MTF filtering
%    7. Spectrum normalization (max and energy normalization)
%    8. Separation of data by muscle
%    9. Saving preprocessed datasets for further analysis
%
%  Output:
%  -------
%  The script generates normalized and smoothed EMG spectra ready for
%  statistical analysis and comparison across subjects, muscles, and
%  experimental conditions.
%
%  Authors: Piotr Wawryka and Ludwin Molina Arias
%  Afiliation: AGH University of Krakow
%  ========================================================================

clear; close all; clc;

%% ========================================================================
%  1. Load EMG data
%  ========================================================================

NT_path = fullfile('Data', 'Neutral Forearm');
SP_path = fullfile('Data', 'Supinated Forearm');
addpath(fullfile(pwd, 'Functions'))

% Load data only if not already present in workspace
if ~exist('NT_Data','var') || isempty(NT_Data)
    [~, NT_Data] = load_all_mat_structs(NT_path);
end

if ~exist('SP_Data','var') || isempty(SP_Data)
    [~, SP_Data] = load_all_mat_structs(SP_path);
end

%% ========================================================================
%  2. General information
%  ========================================================================

nSubjects = 32;   % Initial number of subjects

%% ========================================================================
%  3. Raw EMG signal visualization (example subject)
%  ========================================================================

% Example subject
iSubj = 14;

% EMG signals
x1 = SP_Data{iSubj}.data.BicepsBrachi;
x2 = SP_Data{iSubj}.data.Brachioradialis;
x3 = NT_Data{iSubj}.data.BicepsBrachi;
x4 = NT_Data{iSubj}.data.Brachioradialis;

% Time vectors
N1  = length(x1);
dt = 1 / SP_Data{iSubj}.sampling_frequency;
t1  = linspace(0, dt*N1, N1);

N2  = length(x3);
t2  = linspace(0, dt*N2, N2);

% Subplot letters
subplotLetters = {'a','b'};

% Create figure with white background
figure('Name','Raw sEMG data','NumberTitle','off','Color','w');

% First subplot: SP condition
ax1 = subplot(2,1,1);
plot(t1, x1, 'LineWidth',1.5); hold on
plot(t1, x2, 'LineWidth',1.5); hold off
ylim([-8500 8500])
xlim([0 60])
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
ylabel('EMG amplitude (µV)','FontSize',12,'FontWeight','bold')
legend({'Biceps brachii','Brachioradialis'},'Location','northeast','FontSize',12)
box off
set(ax1,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{1},')'],'FontSize',12,'FontWeight','bold','Interpreter','none')


% Second subplot: NT condition
ax2 = subplot(2,1,2);
plot(t2, x3, 'LineWidth',1.5); hold on
plot(t2, x4, 'LineWidth',1.5); hold off
ylim([-8000 8000])
xlabel('Time (s)','FontSize',12,'FontWeight','bold')
xlim([0 60])
ylabel('EMG amplitude (µV)','FontSize',12,'FontWeight','bold')
legend({'Biceps brachii','Brachioradialis'},'Location','northeast','FontSize',12)
box off
set(ax2,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{2},')'],'FontSize',14,'FontWeight','bold','Interpreter','none')


% Save figure in Nature style
FigArticle('Raw_data','vector',14);  % vector PDF in ./Results folder

%% ========================================================================
%  4. EMG segmentation parameters
%  ========================================================================

fs = SP_Data{1}.sampling_frequency;

params.minAmplitude    = 500;   % Minimum peak amplitude
params.minPeakDistance = 3.2;   % Minimum distance between peaks [s]
params.minActions      = 8;     % Minimum number of detected actions
params.treshold        = 20;    % Detection threshold

%% ========================================================================
%  5. Automatic segmentation of muscle activations
%  ========================================================================

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    % Concatenate muscles (columns)
    x = [NT_Data{i}.data.BicepsBrachi, NT_Data{i}.data.Brachioradialis];
    y = [SP_Data{i}.data.BicepsBrachi, SP_Data{i}.data.Brachioradialis];

    % Segment EMG activity
    [y_segments, ~] = segmentMuscleActions(x, fs, params);
    [x_segments, ~] = segmentMuscleActions(y, fs, params);

    EMG_activity_segments_SP.(fieldName) = y_segments;
    EMG_activity_segments_NT.(fieldName) = x_segments;
end

nSeg = 10;   % Number of repetitions per muscle

%---------------Plot---------------

px = [SP_Data{iSubj}.data.BicepsBrachi, SP_Data{iSubj}.data.Brachioradialis];
plotSegmentation(px,fs,params)
FigArticle('Segmentation','vector',14)

%% ========================================================================
% 5b. Apply notch filter to remove power-line interference (50 Hz)
% ========================================================================

f0 = 50;      % Notch frequency (Hz)
Q  = 35;      % Quality factor (narrow notch)

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    % SP condition
    nSeg_SP = numel(EMG_activity_segments_SP.(fieldName));  % actual number of segments
    for j = 1:nSeg_SP
        raw_signal = EMG_activity_segments_SP.(fieldName)(j).data;
        EMG_activity_segments_SP.(fieldName)(j).data = notchFilterEMG(raw_signal, fs, f0, Q);
    end

    % NT condition
    nSeg_NT = numel(EMG_activity_segments_NT.(fieldName));  % actual number of segments
    for j = 1:nSeg_NT
        raw_signal = EMG_activity_segments_NT.(fieldName)(j).data;
        EMG_activity_segments_NT.(fieldName)(j).data = notchFilterEMG(raw_signal, fs, f0, Q);
    end
end

%% ========================================================================
%  6. Remove invalid subject (subject 06)
%  ========================================================================

disp('➡ Removing subject 06 from segment structures...');

segmentStructs = {'EMG_activity_segments_SP','EMG_activity_segments_NT'};

for s = 1:length(segmentStructs)
    S = eval(segmentStructs{s});

    if isfield(S, 'subject06')
        S = rmfield(S, 'subject06');
    end

    % Renumber remaining subjects
    for i = 7:32
        oldField = sprintf('subject%02d', i);
        newField = sprintf('subject%02d', i-1);
        if isfield(S, oldField)
            S.(newField) = S.(oldField);
            S = rmfield(S, oldField);
        end
    end

    assignin('base', segmentStructs{s}, S);
end

nSubjects = nSubjects - 1;
disp(['Subject removed. Updated nSubjects = ', num2str(nSubjects)]);

%% ========================================================================
%  7. Determine global maximum segment length
%  ========================================================================

maxLength = 0;

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    maxLength = max([maxLength, ...
        max(cellfun(@(s) size(s.data,1), num2cell(EMG_activity_segments_SP.(fieldName)))), ...
        max(cellfun(@(s) size(s.data,1), num2cell(EMG_activity_segments_NT.(fieldName))))]);
end

fprintf('Global maximum segment length: %d samples\n', maxLength);

%% ========================================================================
%  8. Spectral analysis of EMG segments
%  ========================================================================

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    for j = 1:nSeg*2
        EMG_spectrum_SP.(fieldName)(j) = ...
            computeEMGSpectrum(EMG_activity_segments_SP.(fieldName)(j).data, fs, maxLength);

        EMG_spectrum_NT.(fieldName)(j) = ...
            computeEMGSpectrum(EMG_activity_segments_NT.(fieldName)(j).data, fs, maxLength);
    end
end

%% ========================================================================
%  9. Spectrum smoothing (MTF)
%  ========================================================================

Tu = 40;
multiplier = 4.5;

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    for j = 1:nSeg*2
        EMG_SmoothSpectrum_SP.(fieldName)(j) = ...
            mtfSmoothSpectrum(EMG_spectrum_SP.(fieldName)(j).power, Tu, 5, multiplier, ...
                              EMG_spectrum_SP.(fieldName)(j).freq);

        EMG_SmoothSpectrum_NT.(fieldName)(j) = ...
            mtfSmoothSpectrum(EMG_spectrum_NT.(fieldName)(j).power, Tu, 5, multiplier, ...
                              EMG_spectrum_NT.(fieldName)(j).freq);
    end
end

%% ========================================================================
% 10. Spectrum normalization (max & energy)
% ========================================================================
% This section normalizes the smoothed EMG power spectra in two different ways:
%
% 1) Max normalization:
%    - Each spectrum is divided by its maximum value
%    - Preserves the spectral shape
%    - Useful for comparing frequency distribution independently of amplitude
%
% 2) Energy normalization:
%    - Each spectrum is normalized by its total spectral energy
%      (computed using numerical integration)
%    - Preserves relative power contribution across frequencies
%    - Useful for comparing power distribution between conditions/muscles
%
% Both normalizations are stored separately for further analysis.
% ========================================================================

EMG_SmoothSpectrum_NT_maxNorm    = struct();
EMG_SmoothSpectrum_SP_maxNorm    = struct();
EMG_SmoothSpectrum_NT_energyNorm = struct();
EMG_SmoothSpectrum_SP_energyNorm = struct();

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);
    nSegments = numel(EMG_SmoothSpectrum_NT.(fieldName));

    for j = 1:nSegments

        % ---------------- MAX NORMALIZATION ----------------
        smoothNT = EMG_SmoothSpectrum_NT.(fieldName)(j).smoothSpectrum;
        smoothSP = EMG_SmoothSpectrum_SP.(fieldName)(j).smoothSpectrum;

        EMG_SmoothSpectrum_NT_maxNorm.(fieldName)(j).frequency = ...
            EMG_SmoothSpectrum_NT.(fieldName)(j).frequency;
        EMG_SmoothSpectrum_SP_maxNorm.(fieldName)(j).frequency = ...
            EMG_SmoothSpectrum_SP.(fieldName)(j).frequency;

        EMG_SmoothSpectrum_NT_maxNorm.(fieldName)(j).smoothSpectrum = ...
            smoothNT ./ max(max(smoothNT), eps);

        EMG_SmoothSpectrum_SP_maxNorm.(fieldName)(j).smoothSpectrum = ...
            smoothSP ./ max(max(smoothSP), eps);

        % ---------------- ENERGY NORMALIZATION ----------------
        freqNT = EMG_SmoothSpectrum_NT.(fieldName)(j).frequency;
        freqSP = EMG_SmoothSpectrum_SP.(fieldName)(j).frequency;

        energyNT = trapz(freqNT, smoothNT);
        energySP = trapz(freqSP, smoothSP);

        EMG_SmoothSpectrum_NT_energyNorm.(fieldName)(j).frequency = freqNT;
        EMG_SmoothSpectrum_SP_energyNorm.(fieldName)(j).frequency = freqSP;

        EMG_SmoothSpectrum_NT_energyNorm.(fieldName)(j).smoothSpectrum = ...
            smoothNT ./ max(energyNT, eps);

        EMG_SmoothSpectrum_SP_energyNorm.(fieldName)(j).smoothSpectrum = ...
            smoothSP ./ max(energySP, eps);
    end
end

%% ========================================================================
% 11. Split spectra by muscle
% ========================================================================
% Each segmented EMG spectrum contains interleaved muscle data:
%
%   Odd indices  (1,3,5,...) → Biceps brachii
%   Even indices (2,4,6,...) → Brachioradialis
%
% This section separates spectra by muscle while preserving:
%   - Subject
%   - Forearm condition (NT / SP)
%   - Normalization method (max / energy)
%
% The resulting structures are organized and ready for statistical analysis.
% ========================================================================

for i = 1:nSubjects
    fieldName = sprintf('subject%02d', i);

    % ---------------- Brachioradialis ----------------
    idx = 1;
    for j = 2:2:2*nSeg
        Brachioradialis_SP_energy.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_SP_energyNorm.(fieldName)(j);
        Brachioradialis_NT_energy.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_NT_energyNorm.(fieldName)(j);

        Brachioradialis_SP_max.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_SP_maxNorm.(fieldName)(j);
        Brachioradialis_NT_max.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_NT_maxNorm.(fieldName)(j);

        idx = idx + 1;
    end

    % ---------------- Biceps brachii ----------------
    idx = 1;
    for j = 1:2:2*nSeg
        BicepsBrachi_SP_energy.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_SP_energyNorm.(fieldName)(j);
        BicepsBrachi_NT_energy.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_NT_energyNorm.(fieldName)(j);

        BicepsBrachi_SP_max.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_SP_maxNorm.(fieldName)(j);
        BicepsBrachi_NT_max.(fieldName)(idx) = ...
            EMG_SmoothSpectrum_NT_maxNorm.(fieldName)(j);

        idx = idx + 1;
    end
end

%% ========================================================================
% 12. Save final preprocessed dataset
%  ========================================================================

save("Preprocessed_EMG.mat", ...
     "BicepsBrachi_SP_energy","BicepsBrachi_NT_energy", ...
     "Brachioradialis_SP_energy","Brachioradialis_NT_energy", ...
     "BicepsBrachi_SP_max","BicepsBrachi_NT_max", ...
     "Brachioradialis_SP_max","Brachioradialis_NT_max")

disp('EMG preprocessing completed successfully.');

%% Preprocessing (Figures)

f=EMG_SmoothSpectrum_SP.subject01(1).frequency;
x1=EMG_SmoothSpectrum_SP.subject01(1).smoothSpectrum;
x2=EMG_SmoothSpectrum_SP.subject01(2).smoothSpectrum;

x01=EMG_spectrum_SP.subject01(1).power;
x02=EMG_spectrum_SP.subject01(2).power;

x3=EMG_SmoothSpectrum_NT.subject01(1).smoothSpectrum;
x4=EMG_SmoothSpectrum_NT.subject01(2).smoothSpectrum;

x03=EMG_spectrum_NT.subject01(1).power;
x04=EMG_spectrum_NT.subject01(2).power;

% Subplot letters
subplotLetters = {'a','b','c','d'};

% Nature-style colors
colors     = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.4940 0.1840 0.5560]};

% Line widths
lw_primary = 2;  % thicker line for first signal
lw_secondary = 1.2;  % thinner line for second signal

% Transparency for second line
alpha_secondary = 0.3;

% Create figure with white background
fig = figure('Name','Processed data','NumberTitle','off','Color','w');

% -------------------------
% Subplot 1: SP condition, Biceps
ax1 = subplot(2,2,1);
plot(f, x1, 'Color', colors{2}, 'LineWidth', lw_primary); hold on
plot(f, x01, 'Color', [colors{1} alpha_secondary], 'LineWidth', lw_secondary); hold off
xlim([0 200])
xlabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
ylabel('PSD (a.u.)','FontSize',12,'FontWeight','bold')
box off
set(ax1,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{1},')'],'FontSize',12,'FontWeight','bold','Interpreter','none')

% Subplot 2: SP condition, Brachioradialis
ax2 = subplot(2,2,2);
plot(f, x2, 'Color', colors{2}, 'LineWidth', lw_primary); hold on
plot(f, x02, 'Color', [colors{1} alpha_secondary], 'LineWidth', lw_secondary); hold off
xlim([0 200])
xlabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
ylabel('PSD (a.u.)','FontSize',12,'FontWeight','bold')
box off
set(ax2,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{2},')'],'FontSize',12,'FontWeight','bold','Interpreter','none')

% Subplot 3: NT condition, Biceps
ax3 = subplot(2,2,3);
plot(f, x3, 'Color', colors{2}, 'LineWidth', lw_primary); hold on
plot(f, x03, 'Color', [colors{1} alpha_secondary], 'LineWidth', lw_secondary); hold off
xlim([0 200])
xlabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
ylabel('PSD (a.u.)','FontSize',12,'FontWeight','bold')
box off
set(ax3,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{3},')'],'FontSize',12,'FontWeight','bold','Interpreter','none')

% Subplot 4: NT condition, Brachioradialis
ax4 = subplot(2,2,4);
plot(f, x4, 'Color', colors{2}, 'LineWidth', lw_primary); hold on
plot(f, x04, 'Color', [colors{1} alpha_secondary], 'LineWidth', lw_secondary); hold off
xlim([0 200])
xlabel('Frequency (Hz)','FontSize',12,'FontWeight','bold')
ylabel('PSD (a.u.)','FontSize',12,'FontWeight','bold')
box off
set(ax4,'FontSize',12,'LineWidth',1,'TickDir','out')
title(['(', subplotLetters{4},')'],'FontSize',12,'FontWeight','bold','Interpreter','none')

% Save figure in Nature style
FigArticle('Processed_data','vector',14);  
