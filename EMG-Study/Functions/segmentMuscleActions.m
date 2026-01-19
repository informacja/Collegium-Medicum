function [segments, info] = segmentMuscleActions(signal, fs, params)
% SEGMENTMUSCLEACTIONS
% Segmentation of muscle actions from 2-channel EMG signal
%
% signal : [Nx2] EMG data (channel 1, channel 2)
% fs     : sampling frequency [Hz]
% params : struct with fields:
%          .minAmplitude
%          .minPeakDistance (seconds)
%          .minActions
%          .rmsWindow

% -------------------- defaults --------------------
if ~isfield(params,'minAmplitude'),     params.minAmplitude = 500; end
if ~isfield(params,'minPeakDistance'),  params.minPeakDistance = 0.35; end
if ~isfield(params,'minActions'),       params.minActions = 10; end
if ~isfield(params,'rmsWindow'),        params.rmsWindow = 1111; end
if ~isfield(params,'params.treshold'),  params.rmsWindow = 15; end

% -------------------- preparation --------------------
y = signal(:,1) + signal(:,2);   % combined activity
segments = struct('data',{},'channel',{},'index',{});
nrs = 1;

% -------------------- peak detection --------------------
[pks, locs] = findpeaks( ...
    y, ...
    'MinPeakProminence', params.minAmplitude, ...
    'MinPeakDistance',  params.minPeakDistance * fs, ...
    'Threshold', params.treshold);

numActions = numel(locs);

if numActions < params.minActions
    warning('Not enough muscle actions detected.');
    info.valid = false;
    return
end

% add last boundary
locs(end+1) = length(y);

% -------------------- segmentation --------------------
n1 = 1;
for k = 1:numActions

    segRange = locs(k):locs(k+1);

    % RMS envelope to find split point
    e = envelope(y(segRange), params.rmsWindow, 'rms');
    [~, idxMin] = min(e);
    Nbf = segRange(1) + idxMin - 1;

    for ch = 1:2
        s = signal(n1:Nbf, ch);

        segments(nrs).data    = s;
        segments(nrs).channel = ch;
        segments(nrs).index   = k;

        nrs = nrs + 1;
    end

    n1 = Nbf + 1;
end

% -------------------- output info --------------------
info.valid          = true;
info.numActions     = numActions;
info.numSegments    = numel(segments);
info.segmentLengths = arrayfun(@(x) length(x.data), segments);

end
