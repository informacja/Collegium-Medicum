%% ========================================================================
%  Inter-subject EMG Spectral Variability Analysis
%
%  Description:
%  ------------
%  This script quantifies inter-subject variability of EMG spectral
%  representations by computing pairwise distances between subjects.
%
%  The analysis is performed on subject-specific mean EMG spectra,
%  previously obtained from:
%    - Biceps brachii
%    - Brachioradialis
%
%  and under two forearm postures:
%    - Supinated (SP)
%    - Neutral (NT)
%
%  Two normalization strategies are considered:
%    - Energy normalization
%    - Max normalization
%
%  For each muscle–posture–normalization combination, the script:
%    1. Computes pairwise inter-subject spectral distances
%    2. Extracts descriptive statistics (mean, STD, min, max)
%    3. Performs paired statistical testing (SP vs NT)
%    4. Visualizes inter-subject variability using boxplots, scatter plots,
%       and distribution fitting
%
%  Distances are computed in the spectral domain using the Minkowski metrics,
%  providing a frequency-based characterization of inter-subject variability.
%
%  Author: Piotr Wawryka and Ludwin Molina Arias
%  Afiliation: AGH University of Krakow
%  ========================================================================
clear; close; clc; 

distType = 1;

%% ========================================================================
% 1. Configuration
% ========================================================================

nSubjects = 31;                         % Number of subjects
normTypes = {'energy','max'};           % Normalization methods
if(distType == 1) distanceMetric = 'Euclidean'; end           % Inter-subject distance metric
if(distType == 2) distanceMetric = 'Manhattan'; end
if(distType == 3) distanceMetric = 'Chebyshev'; end 
muscleLabelsI = {'BB/SP','BB/NT','BR/SP','BR/NT'};

muscleFields = { ...
    'BicepsBrachi_SP', ...
    'BicepsBrachi_NT', ...
    'Brachioradialis_SP', ...
    'Brachioradialis_NT'};

muscleLabels = { ...
    'Biceps brachii / SP', ...
    'Biceps brachii / NT', ...
    'Brachioradialis / SP', ...
    'Brachioradialis / NT'};

%% ========================================================================
% 2. Load subject-level mean spectra
% ========================================================================
addpath("Results/")
% Mean spectra obtained from the intra-subject analysis
load("Distances_Energy.mat","MeanSpectraEnergy");
load("Distances_Max.mat","MeanSpectraMax");

%% ========================================================================
% 3. Pairwise inter-subject distance computation
% ========================================================================
%
% For each normalization and muscle–posture condition, a symmetric
% distance matrix (subjects × subjects) is computed using the
% Euclidean distance between subject-specific mean spectra.
%


DistancesInter = struct();

for n = 1:numel(normTypes)

    normType = normTypes{n};

    if strcmp(normType,'energy')
        MeanSpectra = MeanSpectraEnergy;
    else
        MeanSpectra = MeanSpectraMax;
    end

    for m = 1:numel(muscleFields)

        distMat = zeros(nSubjects);

        for i = 1:nSubjects
            spec_i = MeanSpectra(i).(muscleFields{m});
            for j = i+1:nSubjects
                spec_j = MeanSpectra(j).(muscleFields{m});

                % --- ÚNICO CAMBIO: mapeo explícito del distType ---
                diff_ij = spec_i - spec_j;
                switch distType
                    case 1      % Euclidean (2-norm)
                        d = norm(diff_ij, 2);
                    case 2      % Manhattan (1-norm)
                        d = norm(diff_ij, 1);
                    case 3      % Chebyshev (Inf-norm)
                        d = norm(diff_ij, Inf);
                    otherwise
                        error('Unsupported distType: %d (use 1=Euclidean, 2=Manhattan, 3=Chebyshev).', distType);
                end
                % ---------------------------------------------------

                distMat(i,j) = d;
                distMat(j,i) = d;
            end
        end

        DistancesInter.(normType).(muscleFields{m}) = distMat;
    end
end

%% ========================================================================
% 4. Inter-subject descriptive statistics
% ========================================================================
%
% Pairwise distances are vectorized using the upper triangular part
% of each distance matrix (excluding the diagonal).
%

StatsInter = table( ...
    'Size',[numel(normTypes)*numel(muscleFields) 7], ...
    'VariableTypes',{'string','string','double','double','double','double','double'}, ...
    'VariableNames',{'Muscle','Normalization','Mean','STD','CV','Min','Max'});

rowIdx = 1;

for n = 1:numel(normTypes)
    normType = normTypes{n};

    for m = 1:numel(muscleFields)

        distMat = DistancesInter.(normType).(muscleFields{m});
        d = distMat(triu(true(nSubjects),1));  

        StatsInter.Muscle(rowIdx)        = muscleLabels{m};
        StatsInter.Normalization(rowIdx)= normType;
        StatsInter.Mean(rowIdx)         = mean(d);
        StatsInter.STD(rowIdx)          = std(d);
        StatsInter.Min(rowIdx)          = min(d);
        StatsInter.CV(rowIdx)           = std(d)/mean(d);
        StatsInter.Max(rowIdx)          = max(d);

        rowIdx = rowIdx + 1;
    end
end

disp('Inter-subject descriptive statistics:');
disp(StatsInter);

%% ========================================================================
% 5. Statistical comparison between postures (SP vs NT)
% ========================================================================
%
% A Wilcoxon signed-rank test is applied to paired inter-subject
% distance distributions (SP vs NT) for each muscle and normalization.
%

musclesToCompare = {'BicepsBrachi','Brachioradialis'};

for n = 1:numel(normTypes)

    normType = normTypes{n};
    fprintf('\nNormalization: %s\n', normType);

    for i = 1:numel(musclesToCompare)

        muscle = musclesToCompare{i};

        dSP = DistancesInter.(normType).([muscle '_SP']);
        dNT = DistancesInter.(normType).([muscle '_NT']);

        vecSP = dSP(triu(true(nSubjects),1));
        vecNT = dNT(triu(true(nSubjects),1));

        [p,~,stats] = signrank(vecSP, vecNT);

        fprintf('%s | p = %.4f (W = %.2f)\n', muscle, p, stats.signedrank);

        if p < 0.05
            fprintf('  → Significant inter-subject difference between postures\n');
        else
            fprintf('  → No significant inter-subject difference\n');
        end
    end
end

%% ========================================================================
% 7. Inter-person variability (Energy & Max)
% ========================================================================
%
% This section visualizes inter-subject variability for each muscle and
% normalization type using boxplots. Each box represents the distribution
% of pairwise Euclidean distances between subjects for a given muscle-condition.
% Two subplots are created in one figure: Left = Max, Right = Energy.
%

normTypes = {'max','energy'};  % Order for subplot arrangement
boxColorsNorm = [0 0 0.7;       % Dark blue for Max
                 0.85 0.33 0.1]; % Orange for Energy
muscleColors  = lines(4);       % Distinct colors for 4 muscles

figure('Name','Inter-person variability (Max & Energy)', 'NumberTitle','off', 'Color','w');

subplotLetters = {'(a)','(b)'};  % Subplot labels

for n = 1:numel(normTypes)
    
    normType = normTypes{n};
    
    % Initialize data and labels
    allData = [];
    groupLabels = {};
    
    for m = 1:4
        distMat = DistancesInter.(normType).(muscleFields{m});
        d = distMat(triu(true(nSubjects),1));  % vectorize upper triangle
        allData = [allData d];
        groupLabels = [groupLabels repmat({muscleLabelsI{m}},1,numel(d))];
    end
    
    % Create subplot
    subplot(1,2,n)
    
    % === Boxplot with notch ===
    h = boxplot(allData, groupLabels, 'Notch','on', 'Widths',0.5, ...
                'MedianStyle','line', 'Symbol','o');
    
    % Customize boxes: color for normalization (Max=blue, Energy=orange)
    boxes = findobj(gca,'Tag','Box'); boxes = flipud(boxes); % correct order
    for k = 1:length(boxes)
        patch(get(boxes(k),'XData'), get(boxes(k),'YData'), ...
              boxColorsNorm(n,:), 'FaceAlpha',0.3, 'EdgeColor',boxColorsNorm(n,:), 'LineWidth',1.2);
    end
    
    % Whiskers
    whiskers = findobj(gca,'Tag','Whisker');
    set(whiskers,'Color','k','LineWidth',1.5);
    
    % Caps
    caps = findobj(gca,'Tag','Upper Whisker'); set(caps,'Color','k','LineWidth',1.5);
    caps = findobj(gca,'Tag','Lower Whisker'); set(caps,'Color','k','LineWidth',1.5);
    
    % Medians
    medians = findobj(gca,'Tag','Median'); set(medians,'Color','k','LineWidth',1.8);
    
    % Outliers
    outliers = findobj(gca,'Tag','Outliers'); set(outliers,'Color','k','MarkerSize',6);
    
    % Labels and formatting
    ylabel('Mean spectral distance between subjects','FontSize',12,'FontWeight','bold');
    xlabel('Muscle/Posture','FontWeight','bold')
    title(subplotLetters{n},'FontWeight','bold');  
    box off;
    set(gca,'FontSize',12,'LineWidth',1);
end

% Save the figure
FigArticle('InterSubj_Boxplot','vector',14);

%% ========================================================================
% 8. Inter-person distances – Scatter plots with vertical offset
% ========================================================================
%
% Visualize pairwise inter-subject distances for each muscle and normalization.
% Horizontal jitter and vertical offset help separate points by category.
%

subnames ={'(a)','(b)'};
figure('Name','Inter-person distances by category','NumberTitle','off','Color','w');

colors = lines(4);       % colors for the 4 muscles
nNorms = length(normTypes);
gap = 0.85;              % separation between categories on X-axis

for n = 1:nNorms
    normType = normTypes{n};
    subplot(1,nNorms,n)
    hold on
    xlabel('Muscle / Posture','FontWeight','bold','FontSize',12)
    ylabel([ distanceMetric ' spectral distances'],'FontWeight','bold','FontSize',12)
    title(subnames(n),'FontWeight','bold','FontSize',12)
    
    xOffset = 0;  % initial x-offset for first block
    
    for m = 1:4
        distMat = DistancesInter.(normType).(muscleFields{m});
        dVec = distMat(triu(true(nSubjects),1));  % upper-triangle distances
        nPairs = length(dVec);
        
        % Horizontal jitter
        jitter = (rand(1,nPairs)-0.5)*0.6;
        xVals = xOffset + jitter;
        
        % Vertical offset for visibility
        yOffset = strcmp(normType,'energy')*0.1*(m-1) + strcmp(normType,'max')*2.0*(m-1);
        yVals = dVec + yOffset;
        
        scatter(xVals, yVals, 20, 'filled', 'MarkerFaceColor', colors(m,:))
        
        % Advance offset for next category
        xOffset = xOffset + gap;
    end
    
    % X-axis ticks
    xticks(gap*(0:3))
    xticklabels(muscleLabelsI)
    xlim([-gap/2 xOffset-0.5])
    box off
end

FigArticle('ScatterPlot','vector',14)
%% ================= Histograms with model fits (separate figures for Max & Energy) =================
nBins = 30;  % number of bins
boxColors = [0 0 0.7;       % dark blue for Max
             0.85 0.33 0.1]; % orange for Energy

subplotLetters = {'(a)','(b)','(c)','(d)'};
for n = 1:length(normTypes)
    normType = normTypes{n};
    
    % Create figure for each normalization
    fig = figure('Name',sprintf('Histograms - %s',normType), 'NumberTitle','off', 'Color','w');
    
    for m = 1:4
        distMat = DistancesInter.(normType).(muscleFields{m});
        distVec = distMat(triu(true(nSubjects),1));
        distVec = distVec(~isnan(distVec));  % remove NaNs
        
        subplot(2,2,m)
        hold on
        
        % Fixed color per normalization
        if strcmp(normType,'max')
            col = boxColors(1,:);
        else
            col = boxColors(2,:);
        end
        
        % ===== Histogram =====
        [counts, edges] = histcounts(distVec, nBins, 'Normalization','pdf');
        centers = edges(1:end-1) + diff(edges)/2;
        hBar = bar(centers, counts, 'FaceAlpha',0.4, 'FaceColor', col, 'EdgeColor','none');
        set(get(get(hBar,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % bar not in legend
        
        % ===== Normal fit =====
        pd_norm = fitdist(distVec,'Normal');
        y_norm = pdf(pd_norm, centers);
        plot(centers, y_norm, 'k-', 'LineWidth',1.5, 'DisplayName','Normal')
        
        % ===== Double Exponential (Laplace) fit =====
        mu_hat = median(distVec);
        b_hat  = mean(abs(distVec - mu_hat));
        y_lap = (1/(2*b_hat)) * exp(-abs(centers - mu_hat)/b_hat);
        plot(centers, y_lap, 'r-', 'LineWidth',1.5, 'DisplayName','Double Exp')
        
        % ===== Maxwell-Boltzmann (GEV) fit =====
        try
            pd_gev = fitdist(distVec,'GeneralizedExtremeValue');
            y_gev = pdf(pd_gev, centers);
            plot(centers, y_gev, 'b-', 'LineWidth',1.5, 'DisplayName','Maxwell-Boltzmann')
        catch
            warning('GEV/Maxwell fit failed for %s (%s)', normType, muscleFields{m});
        end
        
        xlabel( [distanceMetric ' Spectral Distance'],'FontWeight','bold')
        ylabel('PDF','FontWeight','bold')
        title(subplotLetters(m))
        legend('Location','best')
        grid on
        
        % ===== Chi-squared tests =====
        % Normal
        try
            [~,p_norm] = chi2gof(distVec, 'CDF', pd_norm, 'NParams',2, 'Alpha',0.05);
        catch
            p_norm = NaN;
        end
        
        % Double Exponential
        cdf_lap = @(x) 0.5 + 0.5*sign(x-mu_hat).*(1 - exp(-abs(x-mu_hat)/b_hat));
        try
            [~,p_lap] = chi2gof(distVec, 'CDF', cdf_lap, 'Alpha',0.05);
        catch
            p_lap = NaN;
        end
        
        % Maxwell-Boltzmann / GEV
        try
            cdf_gev = @(x) cdf(pd_gev, x);
            [~,p_gev] = chi2gof(distVec, 'CDF', cdf_gev, 'Alpha',0.05);
        catch
            p_gev = NaN;
        end
        
        fprintf('\n%s (%s) - Chi-squared test results:\n', muscleLabels{m}, normType);
        fprintf('Normal: p = %.4f -> %s\n', p_norm, ternary(p_norm<0.05,'Reject H0','Fail to reject H0'));
        fprintf('Double Exp: p = %.4f -> %s\n', p_lap, ternary(p_lap<0.05,'Reject H0','Fail to reject H0'));
        fprintf('Maxwell-Boltzmann: p = %.4f -> %s\n', p_gev, ternary(p_gev<0.05,'Reject H0','Fail to reject H0'));
    end
    
    % Save figure
    FigArticle(sprintf('Histograms_%s',normType),'vector',14);
end

%% ===== Auxiliary ternary function =====
function out = ternary(cond, valTrue, valFalse)
    if cond
        out = valTrue;
    else
        out = valFalse;
    end
end
