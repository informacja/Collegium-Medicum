function plotSegmentation(signal, fs, params)

% -------------------- defaults --------------------
if ~isfield(params,'minAmplitude'),     params.minAmplitude = 500; end
if ~isfield(params,'minPeakDistance'),  params.minPeakDistance = 0.35; end
if ~isfield(params,'rmsWindow'),        params.rmsWindow = 1200; end
if ~isfield(params,'treshold'),         params.treshold = 15; end

% -------------------- Nature colors --------------------
blueLightNature = [0.3 0.55 0.85];  % light blue signal
orangeNature    = [0.85 0.33 0.1];  % RMS envelope

% -------------------- preparation --------------------
y = signal(:,1) + signal(:,2);
t = (0:length(y)-1)/fs;

% -------------------- global RMS envelope --------------------
envRMS = envelope(y, params.rmsWindow, 'rms');

% -------------------- peak detection --------------------
[pks, locs] = findpeaks( ...
    y, ...
    'MinPeakProminence', params.minAmplitude, ...
    'MinPeakDistance',  params.minPeakDistance * fs, ...
    'Threshold', params.treshold);

locs(end+1) = length(y);

% -------------------- figure --------------------
figure('Name','sEMG Segmentation','NumberTitle','off','Color','w'); hold on; box off;

% Signal (light blue, thin)
hSignal = plot(t, y, ...
    'Color', blueLightNature, ...
    'LineWidth', 0.8);

% Continuous RMS envelope (Nature orange)
hEnv = plot(t, envRMS, ...
    'Color', orangeNature, ...
    'LineWidth', 2.0);

% -------------------- segmentation lines --------------------
for k = 1:numel(locs)-1

    segRange = locs(k):locs(k+1);

    % find split point using global envelope
    [~, idxMin] = min(envRMS(segRange));
    Nbf = segRange(1) + idxMin - 1;

    xline(t(Nbf), '--', ...
        'Color', [0.25 0.25 0.25], ...
        'LineWidth', 1.5, ...
        'HandleVisibility','off');
end

% Peaks
hPeaks = plot(t(locs(1:end-1)), pks, 'kv', ...
    'MarkerFaceColor','k', ...
    'MarkerSize', 8);

% -------------------- Nature formatting --------------------
set(gca, ...
    'FontName','Helvetica', ...
    'FontSize', 14, ...
    'LineWidth', 0.75, ...
    'TickDir','out');

xlabel('Time (s)','FontWeight','bold')
ylabel('Combined amplitude [uV]','FontWeight','bold')

legend([hSignal hEnv hPeaks], ...
       {'Signal','RMS envelope','Peaks'}, ...
       'Location','northoutside', ...
       'Orientation','horizontal', ...
       'Box','off',FontSize=14);

xlim([t(1) 60])
grid off

end
