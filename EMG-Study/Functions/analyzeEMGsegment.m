function out = FFT_Transformation(segment, fs, Nfft)
% Computes single-sided amplitude and power spectrum of EMG
% segment [µV] with optional zero-padding
%
% Inputs:
%   segment : column vector [µV]
%   fs      : sampling frequency [Hz]
%   zpFactor: zero-padding factor (default 1)
%
% Outputs:
%   out.freq      : frequency vector [Hz]
%   out.amplitude : single-sided amplitude spectrum [µV]
%   out.power     : single-sided power spectrum [µV^2]

segment = double(segment(:));
N = length(segment);           % original signal length

% Window
win = hann(N);
xWin = segment .* win;

% FFT with zero-padding
Y = fft(xWin, Nfft);

% Single-sided amplitude spectrum
P2 = abs(Y / N);               % divide by N (original length) to keep units
nyqIdx = floor(Nfft/2)+1;
P1 = P2(1:nyqIdx);
P1(2:end-1) = 2*P1(2:end-1);

% Power spectrum
Power = P1.^2;

% Frequency vector
freq = (0:nyqIdx-1)' * fs / Nfft;

% Output
out = struct();
out.freq = freq;
out.amplitude = P1;
out.power = Power;
out.Nfft = Nfft;
out.N = N;
end
