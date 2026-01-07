function [pvg,pve,pvM,chi2g,chi2e,chi2M,bins,Nemp] = mhistMGEfromChat(y, lbins)
% mhistMGE_merge - histogram + chi-square GOF for Normal, Laplace, Maxwell-Boltzmann
% with automatic merging of adjacent bins so that expected counts >= minExp.
%
% USAGE:
%  [pvg,pve,pvM,chi2g,chi2e,chi2M,bins,Nemp] = mhistMGE_merge(y, lbins)
%
% Default lbins = 30 if not provided.
% Returns p-values and chi2 statistics for each model as well as the
% original bin centers and empirical counts.
%
% Notes:
%  - Merging is done separately for each model when computing chi2 and p.
%  - Degrees of freedom corrected: df = nGroups - 1 - nParamsEstimated.
%  - Minimum expected count for each merged group is minExp (5 by default).

    if nargin < 2, lbins = 30; end
    minExp = 5;                 % minimum expected count per merged group
    N = length(y);

    % Basic stats
    ymin = min(y); ymax = max(y);
    mu = mean(y); sigma = std(y);

    % Histogram (counts and bin centers)
    [Nemp, edges] = histcounts(y, lbins);
    bins = edges(1:end-1) + diff(edges)/2;
    Dy = edges(2) - edges(1);

    % Plot histogram and hold
    figure; bar(bins, Nemp, 'FaceAlpha', 0.3); hold on;

    % Normalization constant for Gaussian
    CN = 1/(sqrt(2*pi)*sigma);

    % Fit Maxwell-Boltzmann scale 'a' by simple grid search (like before)
    % initial guess from modal bin
    [~, idxMax] = max(Nemp);
    modeGuess = abs(bins(idxMax))/sqrt(2);
    amin = 0.5 * modeGuess; amax = 1.5 * modeGuess;
    La = 1000;
    aGrid = linspace(amin, amax, La);
    bestJ = inf;
    for ii = 1:La
        aTry = aGrid(ii);
        % compute expected per bin under MB using current aTry
        fMvals = arrayfun(@(z) pdfMBn(z, aTry), bins);
        EM = fMvals * N * Dy;
        % safe chi2 (avoid division by zero)
        J = sum( (Nemp - EM).^2 ./ max(EM, eps) );
        if J < bestJ
            bestJ = J; aOpt = aTry;
        end
    end
    a = aOpt;

    % Create smooth theoretical curves for plotting
    ldx = 500;
    xx = linspace(ymin, ymax, ldx);
    fg_x = CN .* exp(-((xx - mu)./(sqrt(2)*sigma)).^2);      % normal pdf
    fe_x = (0.5/sigma) .* exp(-abs(xx - mu)/sigma);         % Laplace pdf (location=mu, scale=sigma)
    fM_x = arrayfun(@(z) pdfMBn(z, a), xx);                 % MB pdf
    plot(xx, N .* fg_x * Dy, 'r', xx, N .* fe_x * Dy, 'g', xx, N .* fM_x * Dy, 'b');

    % For each model: compute expected per original bin, merge adjacent bins until each E >= minExp,
    % then compute chi2 and p-value with corrected df.

    % Helper: merge adjacent bins (observed and expected) until all expected >= minExp
    function [ObsGroup, ExpGroup] = merge_bins_for_expected(Obs, Exp)
        ObsGroup = Obs(:);
        ExpGroup = Exp(:);
        % If any expected < minExp, merge with the neighbor that gives a larger expected after merge.
        while any(ExpGroup < minExp) && length(ExpGroup) > 1
            idx = find(ExpGroup < minExp, 1, 'first'); % take first small expected group
            % consider merging with left or right neighbor
            if idx == 1
                % must merge with right
                ExpGroup(2) = ExpGroup(2) + ExpGroup(1);
                ObsGroup(2) = ObsGroup(2) + ObsGroup(1);
                ExpGroup(1) = []; ObsGroup(1) = [];
            elseif idx == length(ExpGroup)
                % merge with left
                ExpGroup(end-1) = ExpGroup(end-1) + ExpGroup(end);
                ObsGroup(end-1) = ObsGroup(end-1) + ObsGroup(end);
                ExpGroup(end) = []; ObsGroup(end) = [];
            else
                % merge with neighbor that yields larger resulting expected
                if ExpGroup(idx-1) >= ExpGroup(idx+1)
                    ExpGroup(idx-1) = ExpGroup(idx-1) + ExpGroup(idx);
                    ObsGroup(idx-1) = ObsGroup(idx-1) + ObsGroup(idx);
                    ExpGroup(idx) = []; ObsGroup(idx) = [];
                else
                    ExpGroup(idx+1) = ExpGroup(idx+1) + ExpGroup(idx);
                    ObsGroup(idx+1) = ObsGroup(idx+1) + ObsGroup(idx);
                    ExpGroup(idx) = []; ObsGroup(idx) = [];
                end
            end
        end
        % If after repeated merging we still have small expected (rare), merge from ends until okay
        while any(ExpGroup < minExp) && length(ExpGroup) > 1
            if ExpGroup(1) < minExp
                ExpGroup(2) = ExpGroup(2) + ExpGroup(1);
                ObsGroup(2) = ObsGroup(2) + ObsGroup(1);
                ExpGroup(1) = []; ObsGroup(1) = [];
            elseif ExpGroup(end) < minExp
                ExpGroup(end-1) = ExpGroup(end-1) + ExpGroup(end);
                ObsGroup(end-1) = ObsGroup(end-1) + ObsGroup(end);
                ExpGroup(end) = []; ObsGroup(end) = [];
            else
                break;
            end
        end
    end

    % Normal model
    % Compute expected counts under Normal distribution using CDF differences
% edges = histogram bin edges (length = lbins+1)
edges = [bins - Dy/2, bins(end) + Dy/2];  % reconstruct edges from bin centers
Eg = N .* ( normcdf(edges(2:end), mu, sigma) - normcdf(edges(1:end-1), mu, sigma) );

% Merge bins so all expected >= 5
[ObsG, ExpG] = merge_bins_for_expected(Nemp, Eg);

% Chi-square statistic
chi2g = sum( (ObsG - ExpG).^2 ./ max(ExpG, eps) );

% Degrees of freedom: groups - 1 - estimated parameters (mu, sigma)
df_g = max(length(ObsG) - 1 - 2, 1);

% p-value
pvg = 1 - chi2cdf(chi2g, df_g);


    % Laplace model (location mu, scale sigma as proxy => 2 params)
    Ee = N .* ((0.5/sigma) .* exp(-abs(bins - mu)/sigma)) * Dy;
    [ObsE, ExpE] = merge_bins_for_expected(Nemp, Ee);
    chi2e = sum( (ObsE - ExpE).^2 ./ max(ExpE, eps) );
    df_e = max(length(ObsE) - 1 - 2, 1);
    pve = 1 - chi2cdf(chi2e, df_e);

    % Maxwell-Boltzmann model (1 parameter 'a')
    EM = N .* arrayfun(@(z) pdfMBn(z, a), bins) * Dy;
    [ObsM, ExpM] = merge_bins_for_expected(Nemp, EM);
    chi2M = sum( (ObsM - ExpM).^2 ./ max(ExpM, eps) );
    df_M = max(length(ObsM) - 1 - 1, 1);
    pvM = 1 - chi2cdf(chi2M, df_M);

    hold off;
end
