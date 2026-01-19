function emg_filtered = notchFilterEMG(emg_signal, fs, f0, Q)
%NOTCHFILTEREMG Applies a notch filter to remove power-line interference from sEMG
%
%   emg_filtered = notchFilterEMG(emg_signal, fs, f0, Q)
%
%   Inputs:
%       emg_signal : vector of sEMG data
%       fs         : sampling frequency (Hz)
%       f0         : notch frequency (Hz), e.g., 50 or 60 Hz
%       Q          : quality factor (higher Q = narrower notch)
%
%   Output:
%       emg_filtered : filtered sEMG signal
%
%   This version uses MATLAB's designNotchPeakIIR function for robust
%   coefficient calculation and filtfilt for zero-phase filtering.

    % --- Normalized center frequency (0-1, relative to Nyquist) ---
    w0 = f0 / (fs/2);

    % --- Design notch filter using designNotchPeakIIR ---
    [B, A] = designNotchPeakIIR( ...
        Response="notch", ...
        CenterFrequency=w0, ...
        QualityFactor=Q, ...
        FilterOrder=2 ...
    );

    % --- Apply zero-phase filtering ---
    emg_filtered = filtfilt(B, A, emg_signal);

end
