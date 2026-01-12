%  Purpose:
%    - Load EMG recordings from folder structure (MAT files)
%    - Segment repetitions for selected exercises
%    - Apply Hann windowing and zero-padding
%    - Compute spectra and spectral centroids
%    - Compare conditions using Minkowski distances (several variants)
%    - Compose & export publication-ready figures (IEEE inches)

%% --------------------------- Configuration ------------------------------
DEBUG         = 1;   % 1 = verbose, figures left open; 0 = quiet (e.g., compiled)
forArticle    = 1;   % 1 = clean caches + export figures for publication
experimental  = 0;   % 1 = enable extra two-band plotting (experimental)
% Parseval flag is set later if not defined; see DEBUG section

% Language codes (used for labels on plots/strings)
PL = 1; EN = 2;
lang = EN;

% In deployed mode, reduce verbosity and set Parseval behavior
if (isdeployed)
    DEBUG    = 0;
    % Parseval semantics:
    %  -1 = negative (relative Twygl vs. segment length)
    %   0 = do not recompute Twygl, segments without zeros
    %   1 = add zeros
    Parseval = -1;
end

% If preparing figures for an article, remove caches and previous figure exports
if (forArticle)
    % WARNING: destructive actions
    % delete segments.mat    % Uncomment only if you want to force re-segmentation
    % delete signals.mat     % Uncomment only if you want to recompute signals
    %delete spectrums.mat     % force recomputation of spectra
    %delete centroids.mat     % force recomputation of centroids
    %delete figBase/*         % remove previously saved figures
    % TODO: add a safety check or backup before deleting artifacts
end

% Debug-mode housekeeping
if (DEBUG)
    close all;               % keep variables, close figures
    % delete segments.mat    % optional: will force segmentation (~10s)
    % delete signals.mat     % optional: will force signal prep (Hann + zeros)
    % delete spectrums.mat
    % delete centroids.mat

    % If Parseval was not defined earlier, default to adding zeros
    if (~exist("Parseval","var"))
        Parseval = 1;        % 1 = add zeros
    end
    % Parseval = 0;          % example: 0 = do not add zeros
end

% Display cutoff for centroid plots (Hz)
fWyswieltCentroidow = 555;

% Globals used elsewhere in helper functions
global cntCol;
global txBR txBB;

% Stopwatch for entire pipeline
allElapsedTime = tic;

% Processing flags
windowing           = 1;   % 1 = apply Hann windowing
printCentroids      = 0;   % 1 = verbose centroid tables in console
plotAllFigures      = 0;   % 1 = show all figures (can be heavy)
compareExampleData  = 0;   % 1 = load demo Neurosoft data (outside Noraxon/Qualisys)
deprecated          = 0;   % code that is not used

% Root data dir:
dirname = 'data'; % NOTE: author recorded twice (sum patients ≈ 34)
% dirname = '../../Ortheo3D/data'; % example path

% Pre-allocate recording container
v = [];                   % array of structs; each entry stores a single muscle stream
files     = dir(fullfile(dirname,'**','*.mat'));  % recursive search for MAT files
datafiles = fullfile({files.folder},{files.name});

% Exercise text labels (Polish/English)
txPr    = "Pośredni";         % neutral grip (forearm)
txPc    = "Podchwyt";         % supinated position
txInLab = "wysiłek w pracowni EMG";
if (lang == EN)
    txPr    = "Neutral grip";
    txPc    = "Supinated position";
    txInLab = "in EMG Lab";
end
txBR = "Brachioradialis";
txBB = "Biceps brachii";

% Segmentation thresholds
minActions                 = 10;   % minimal number of repetitions for selection
skipedTrainingReppetinons  = 0;    % tracking skipped repetitions

Qualisys = 0;   % 1 = Qualisys/Delsys input; 0 = Noraxon Ultium input

% Index bookkeeping
k = 0; 
global fromSzU;   % used when loading Neurosoft demo data

% Fail early if no data
if (isempty(datafiles))
    error("No datafiles found\n");
end

addpath(fullfile('additionalFiles'));
addpath(fullfile('processData'));

%% -------------------- Optional: Load demo Neurosoft data ----------------
if (compareExampleData)
    % NOTE: This branch loads Neurosoft data (not Noraxon/Qualisys)
    % dirnameN = '../../PhD/data/SzU/...';  % set a valid path on your machine
    dirnameN = '../../PhD/data/SzU/Neurosoft/EDIpublication/neurogenic (tibialis anterior)/T.K.27.02.2025K72';
    [v, k, fromSzU] = LoadNeurosoftData(dirnameN, k);
end
skipSegmentationIndex = k;

%% -------------------------- Load all MAT files --------------------------
% Each MAT file may contain either Qualisys/Delsys or Noraxon Ultium structure
for f = datafiles
    k  = k + 1;
    df = load(string(f));            % load MAT structure
    vars = fieldnames(df);           % top-level variables in the MAT

    for i = 1:length(vars)
        assignin('base', 'tmp', df.(vars{i}));

        if (isfield(tmp,"Analog"))   % ---------- Qualisys / Delsys Trigno ----------
            Qualisys = 1;

            % Remove DC offset per channel
            tmpData = tmp.Analog.Data;
            tmpData = tmpData - mean(tmpData, 2);

            % Board name check (only implemented for Delsys Trigno API)
            tmpName = tmp.Analog.BoardName;
            if (tmpName == 'Delsys Trigno API')
                for (c = tmp.Analog.ChannelNumbers)
                    switch (tmp.Analog.Labels{c})
                        case 'L_Rectus Femoris'
                            % TODO: Verify mapping. "Rectus Femoris" is a quadriceps muscle,
                            % not a forearm muscle (BR). If labels are from leg recordings,
                            % this mapping is incorrect. Ensure channel labels match BR/BB.
                            v(k).infoRDisp = txBR;
                            v(k).dataR     = tmpData(c,:);         % assign to BR
                            v(k).infoRecord= tmp.File(end-4);      % NOTE: fragile index
                            v(k).infoRName = tmp.Analog.Labels{c}; % raw label
                        case 'L_Vastus Lateralis'
                            % Same caution as above: validate mapping.
                            v(k).dataB     = tmpData(c,:);         % assign to BB
                            v(k).infoRecord= tmp.File(end-4);
                            v(k).infoBName = tmp.Analog.Labels{c};
                            v(k).infoBDisp = txBB;
                    end
                end
            else
                error("Unimplemented BoardName for Qualisys");
            end

        else                          % ---------- Noraxon Ultium ----------
            % Noraxon stores channels as signal_1 / signal_2 with names
            tmpData = tmp.movements.sources.signals.signal_1.data;
            tmpName = tmp.movements.sources.signals.signal_1.name;

            if (tmpName == "Ultium EMG-BRACHIORAD. RT")
                v(k).dataR      = tmpData;
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoRName  = tmpName;
                v(k).infoRDisp  = txBR;
            end
            if (tmpName == "Ultium EMG-BICEPS BR. RT")
                v(k).dataB      = tmpData;
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoBName  = tmpName;
                v(k).infoBDisp  = txBB;
            end

            % Some sessions store the second muscle as signal_2
            tmpData = tmp.movements.sources.signals.signal_2.data;
            tmpName = tmp.movements.sources.signals.signal_2.name;

            if (tmpName == "Ultium EMG-BRACHIORAD. RT")
                v(k).dataR      = tmpData;
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoRName  = tmpName;
                v(k).infoRDisp  = txBR;
            end
            if (tmpName == "Ultium EMG-BICEPS BR. RT")
                v(k).dataB      = tmpData;
                v(k).infoRecord = tmp.info.record_name;
                v(k).infoBName  = tmpName;
                v(k).infoBDisp  = txBB;
            end
            if(tmpName == "SN204 16 kanałowe EMG-BRACHIORAD. RT") k=k-1; end;
        end % source type
    end % vars loop
end % files loop

% Summary: count of subjects (heuristic: each file contains two muscles)
fprintf(1,"Number of subjects (approx.): %d, lang=%d\n", int32(length(v)/2), lang); 
toc(allElapsedTime);

fprintf(1,"Matrix v size: %d (measurements * number of exercises)\n", length(v));

%% ------------------------ Sampling & Units Setup ------------------------
% Infer Y units and sampling frequency from the last loaded sample
if (Qualisys)
    Yunits   = 'uV';
    fpom     = tmp.Analog.Frequency;
    minActions = 2;     % often fewer repetitions recorded via Qualisys
else
    Yunits   = tmp.movements.sources.signals.signal_1.units;
    fpom     = tmp.movements.sources.signals.signal_1.frequency;
    % Else: Neurosoft branch (if enabled) would set via LoadNeurosoftData
end

dtpom = 1 / fpom;   % sampling period [s]

%% ----------------------------- Segmentation -----------------------------
% If cached segmentation exists, load it; otherwise compute
if (exist("segments.mat", "file"))
    load segments.mat
else
    nrF = 20;             % figure number for segmentation visualizations
    segmentActions;       % external function: populates 'segment', 'segMio', 'segTraining'
end

% Optional reporting of distribution of repetitions per grip
if (skipedTrainingReppetinons)
    fprintf(1, "Distribution: %d (neutral) %d (supinated). Skipped reps total: %d\n", ...
            countPosredni, countPodchwyt, skipedTrainingReppetinons);
end

fprintf(1,"Number of segments: %d (files * muscles * repetitions)\n", length(segment));
fprintf(1,"Estimated segments: (%d)\n", length(v) * 2 * 10);

% Time quantization for segments (ensure even length)
nrFw  = 201;         % figure number for spectra
lfrow = 4; lc = 2;   % subplot grid defaults
Tsyg  = ceil(ceil(m*dtpom)/2) * 2;   % max action duration [s], rounded to even
lSyg  = round(Tsyg/dtpom);           % samples per action (target length)

clear Syg;
ksyg = 0;     
kol  = 'kbrm';        % colors
kf   = 0;

%% --------------------- Windowing & Zero-padding -------------------------
% If cached, load; otherwise compute windowed/zero-padded signals
if (exist("signals.mat", "file"))
    load signals.mat
else
    tic

    % Determine min/max segment lengths (for diagnostics)
    sLmin = 1e10; sLmax = 0;
    for i = 1:length(segment)
        sLmin = min(sLmin, length(segment(i).data));
        sLmax = max(sLmax, length(segment(i).data));
    end

    % Enforce consistent target length for zero-padding
    sLmax = lSyg;

    if (windowing), fprintf(1, "Applying Hann window...\n"); end

    for i = 1:length(segment)
        SygRawLen(i) = length(segment(i).data');

        % Apply Hann window per segment (column/row safe)
        if (windowing)
            win = hann(SygRawLen(i));
            if (iscolumn(segment(i).data))
                segment(i).data = segment(i).data .* win;
            else
                segment(i).data = segment(i).data' .* win;
            end
        end

        % Zero-pad to target length lSyg
        Syg(i,1:lSyg) = [segment(i).data' zeros(1, lSyg - length(segment(i).data))];

        % Category coding:
        % 'SygKat' indexes by patient/training; 'sygKat' indexes by muscle/grip
        SygKat(i)       = segMio(i) + (segTraining(i)-1) * 2;
        sygKat(i)       = segment(i).miesien + (segment(i).gest-1) * 2;
        segment(i).kat  = sygKat(i);
        segment(i).len  = SygRawLen(i);
        % Mapping:
        % 1 - Neutral BR, 2 - Neutral BB, 3 - Supinated BR, 4 - Supinated BB
    end

    toc;
    save signals.mat Syg SygKat SygRawLen segment sLmax sLmin
end

%% ---------------------------- Spectral Analysis -------------------------
nrF = 500; 
selectTraining;                 % external function: select subset of segments ('j')

% Trend containers (used by helper functions)
MTF(1).Tu = []; MTF(2).Tu = []; MTF(3).Tu = [];

% If cached spectra exist, load; otherwise compute them
if (exist("spectrums.mat", "file"))
    load spectrums.mat
else
    nrF = 1000; fprintf(1, "Computing spectra...\n");
    spectrumTrend;               % external function: fills 'Widma' etc.
end

fprintf(1, "Spectra size: %dx%d\n", size(Widma));

%% --------------------------- Spectral Centroids -------------------------
if (exist("centroids.mat", "file"))
    load centroids.mat
else
    nrF = 2000; fprintf(1, "Computing centroids...\n");
    centroid;                    % external function: fills 'CentrWidm', 'lpacj'
end
fprintf(1, "Centroids size: %dx%d\n", size(CentrWidm));
fprintf(1, "Averaged centroids size: %dx%d\n", size(lpacj));

% Frequency limit for centroid-based comparisons (75 Hz)
flim12 = 0.075;                            % kHz
nlim   = find(xf * dtpom <= flim12);       % harmonic index threshold
nflim  = nlim(end);

% Compute centroid deltas/differences
dCentr;

%% --------- Minkowski Distances — spectra shape vs. centroid features ----
nrF = 4000;
for (jakieDist = 1:4)
    if (mod(jakieDist, 2) == 1)
        nrF = nrF + 1; bf = 0;   % base figure increment
    else
        bf = 4;                  % band offset (experimental)
    end

    % Flag to use maxima or full spectra in metric
    if (jakieDist > 2), flagaMaxima = 0; else flagaMaxima = 1; end

    nrFig = jakieDist * 2 + i + 500;
    wybrJakieDist = jakieDist;

    minkowskiDist;        % external function: populates Sb(1:2).dCM etc.
    figure(nrF); 
    disppolt;             % external plotting of distance results

    if (experimental)
        % Show multi-band experimental plots
        for (nband = bf+1)
            figure(nrF + nband);
            disppolt;
        end

        nrF = nrF + 3;
        figure(nrF);

        kol = 'rgbk';
        for (i = 1:size(Sb(2).dCM,1))
            plot(Sb(1).dCM(i,:), Sb(2).dCM(1,:), [kol(v(i).infoTraining) 'o']); hold on;
        end
        hold off; nrF = nrF + 3; figure(nrF);
        for (i = 1:size(Sb(2).dCM,2))
            plot(Sb(1).dCM(:,i), Sb(2).dCM(:,1), [kol(v(i).infoTraining) 'o']); hold on;
        end
        hold off; nrF = nrF + 3; figure(nrF);

        for (i = 1:28)
            plot(Sb(1).dCM(:,i), Sb(2).dCM(:,i), [kol(v(i).infoTraining),'o']); hold on;
        end
        hold off
    end
end

% Additional Minkowski variants on centroids only
for (jakieDist = 5:6)
    if (mod(jakieDist,2) == 1), nrF = nrF + 1; bf = 0; figure(nrF); else bf = 4; end
    if (jakieDist > 5), flagaMaxima = 0; else flagaMaxima = 1; end
    nrFig = jakieDist * 2 + 500;

    minkowskiDist;
    figure(nrF); 
    ddplot;          % external plotting function for centroid-only distances
end

fprintf(1, "main = "); toc(allElapsedTime);

%% ---------------------- Figure Composition & Export ---------------------
nrF = 5e3; 
clear s;

% Language constants (if used downstream)
PL = 1; EN = 2;

% Figure selection for article (IDs must exist before copy)
s.figNrList = [  
        25
        41
        42
        503
        509
        512
        523
        527
        1001
        1013
        1023
        2005
        2008
        2069
        2078
        4003 ];

% Grid specification per figure {figure_number, rows, cols}
s.figNrGridSize = [
        {25,   2,2}
        {41,   3,1}
        {42,   3,1}
        {503,  2,2}
        {509,  2,2}
        {512,  2,2}
        {523,  2,2}
        {527,  2,2}
        {1001, 2,2}
        {1013, 1,1}
        {1023, 1,1}
        {2005, 2,4}
        {2008, 1,1}
        {2069, 1,1}
        {2078, 2,2}
        {4003, 3,2} ];

% Compose figures: copy axes from base figures to target layouts
% NOTE: This depends on Children indices — fragile; prefer tagging axes.
ax1Chil = figure(4001).Children;
f2 = figure(41);
copyobj(ax1Chil(1), f2)
copyobj(ax1Chil(2), f2)
copyobj(ax1Chil(3), f2)
f2.Children(1).Subtitle.String = "a)";
f2.Children(2).Subtitle.String = "b)";
f2.Children(3).Subtitle.String = "c)";
copyobj([ax1Chil(8) ax1Chil(9)], f2)
f2.Children(2).Visible = "off";
for (i = 1:numel(f2.Children(2).Children)), f2.Children(2).Children(i).Visible = "off"; end
set(f2.Children(1),'Position',[0.3093 0.3579 0.0769 0.1038]);

ax1Chil = figure(4002).Children;
f3 = figure(42);
copyobj(ax1Chil(1), f3)
copyobj(ax1Chil(2), f3)
copyobj(ax1Chil(3), f3)
f3.Children(1).Subtitle.String = "a)";
f3.Children(2).Subtitle.String = "b)";
f3.Children(3).Subtitle.String = "c)";
copyobj([ax1Chil(8) ax1Chil(9)], f3);
f3.Children(2).Visible = "off";
for (i = 1:numel(f3.Children(2).Children)), f3.Children(2).Children(i).Visible = "off"; end
set(f3.Children(1),'Position',[0.3046 0.3429 0.0864 0.1338]);

ax1Chil = figure(2005).Children;
f4 = figure(25);
copyobj(ax1Chil(5), f4)
copyobj(ax1Chil(3), f4)
copyobj(ax1Chil(2), f4)
copyobj(ax1Chil(4), f4)
% Optional: set subtitles "a) ... d)" here

% Export block: adjust sizes and save figures if forArticle is ON
if (forArticle)
    % s.exportPath defines destination folder
    s.exportPath = "article/fi/"; 
    s.inch3dot25 = -1.618;      % NOTE: custom scaling factor

    % Adjust figure layout to article specs and save
    adjust4article(s);
    save4article(s);
    fprintf(1, "main + forArticle(1) = "); toc(allElapsedTime);
end
return