function out = mtfSmoothSpectrum(rawSpectrum, smoothingWindow, filterType, multiplier, freqVec)
% mtfSmoothSpectrum
% -----------------
% Designs an MTF smoothing filter and applies it to a spectrum.
% Ensures smoothed spectrum is never negative.
%
% INPUTS
%   rawSpectrum     : column vector (magnitude or power spectrum)
%   smoothingWindow : smoothing window (Tu equivalent, in samples)
%   filterType      : MTF type (integer, e.g., 5 = z3 trend filter)
%   multiplier      : multiplier to adjust smoothing window (0 < multiplier <= 1)
%   freqVec         : frequency vector corresponding to rawSpectrum [Hz]
%
% OUTPUTS
%   out : struct with fields
%       .smoothSpectrum  : smoothed spectrum (same length as input, >=0)
%       .filterKernel    : MTF kernel used for smoothing
%       .frequency       : frequency vector [Hz]

%% Ensure column vector
rawSpectrum = rawSpectrum(:);
N = length(rawSpectrum);

if N < 10
    error('Spectrum too short for MTF smoothing.');
end

%% Check frequency vector
if nargin < 5 || isempty(freqVec)
    freqVec = (0:N-1)'; % fallback: normalized 0..1
elseif length(freqVec) ~= N
    error('Frequency vector length must match length of rawSpectrum.');
end

%% Adjust smoothing window
adjustedWindow = smoothingWindow * multiplier;

%% MTF filter parameters
Nh1 = [1 1 1.38 2 2*1.112 1.38 2 2 1.38 1.78 1.38];
nh1 = Nh1(filterType);

filterSize = max(round(adjustedWindow * nh1), 3);
halfFilter = filterSize - 1;
filterLength = 2*filterSize - 1;

%% MTF filter design (z-trend)
t = (-filterSize+1:0)' / filterSize;
designMatrix = ones(filterSize,1);
if filterType >= 2, designMatrix = [designMatrix t]; end
if filterType >= 3, designMatrix = [designMatrix t.^2]; end
if filterType >= 4, designMatrix = [designMatrix t.^3]; end

projection = designMatrix * ((designMatrix' * designMatrix) \ designMatrix');

filterKernel = zeros(filterLength,1);
for k = 1:filterLength
    idx = max(1,k-filterSize+1):min(filterSize,k);
    for i = idx
        filterKernel(k) = filterKernel(k) + projection(i, filterSize-k+i);
    end
end
filterKernel = filterKernel / filterSize;

%% Extend spectrum safely (left reflection)
extendedSpectrum = [rawSpectrum(max(halfFilter:-1:1,1)); rawSpectrum];

%% Apply filter
smoothSpectrum = zeros(N,1);
for n = 1:N
    startIndex = n;
    indexRange = startIndex:(startIndex+filterLength-1);
    
    % avoid exceeding extendedSpectrum
    indexRange = indexRange(indexRange <= length(extendedSpectrum));
    smoothSpectrum(n) = sum(extendedSpectrum(indexRange) .* filterKernel(1:length(indexRange)));
end

%% Force non-negative values
smoothSpectrum(smoothSpectrum < 0) = 0;

%% Return struct
out.smoothSpectrum = smoothSpectrum;
out.filterKernel   = filterKernel;
out.frequency      = freqVec;
end

