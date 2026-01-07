function [pvg,pve,pvM,chi2g,chi2e,chi2M,bins,Nemp] = mhistMGE_CHAT(y,lbins)
% mhistMGE - Fits Normal, Exponential, and Maxwell-Boltzmann distributions
%            to the histogram of data and performs chi-square tests.
%
% OUTPUTS:
%   pvg   - p-value for Normal fit
%   pve   - p-value for Exponential fit
%   pvM   - p-value for Maxwell-Boltzmann fit
%   chi2g - chi-square statistic for Normal
%   chi2e - chi-square statistic for Exponential
%   chi2M - chi-square statistic for Maxwell-Boltzmann
%   bins  - histogram bin centers
%   Nemp  - observed bin counts

if nargin < 2, lbins = 30; end
N = length(y);

% Basic stats
mu    = mean(y);
sigma = std(y);

% Histogram
[Nemp, bins] = hist(y, lbins);
Dy = (max(y) - min(y)) / lbins;
bar(bins, Nemp); hold on;

% --- Compute expected counts using CDF differences (more accurate than midpoints) ---
edges = [bins - Dy/2, bins(end) + Dy/2];  % reconstruct edges from centers

% 1) Normal distribution
Eg = N .* ( normcdf(edges(2:end), mu, sigma) - normcdf(edges(1:end-1), mu, sigma) );

% 2) Laplace(ysr, sigma) (two-sided exponential, as in your code)
Fe = @(x) 0.5*(1 + sign(x-mu).*(1 - exp(-abs(x-mu)/sigma)));
Ee = N .* ( Fe(edges(2:end)) - Fe(edges(1:end-1)) );

% 3) Maxwell-Boltzmann
% First optimize scale parameter 'a'
nm = find(Nemp == max(Nemp), 1);
a0 = bins(nm)/sqrt(2);
amin = a0/2; amax = 1.5*a0;
a_grid = linspace(amin, amax, 200);
J = inf(size(a_grid));
for k = 1:numel(a_grid)
    fM = @(x) pdfMBn(x, a_grid(k));
    EM = N .* ( integral(fM, edges(1:end-1), edges(2:end)) ); % expected by bin
    J(k) = sum((Nemp - EM).^2 ./ max(EM, eps));
end
[~, idx] = min(J);
aOpt = a_grid(idx);
fM = @(x) pdfMBn(x, aOpt);
EM = N .* ( arrayfun(@(i) integral(fM, edges(i), edges(i+1)), 1:lbins) );

% --- Chi-square tests ---
% Merge bins if needed so all expected >= 5
[ObsG, ExpG] = merge_bins_for_expected(Nemp, Eg);
[ObsE, ExpE] = merge_bins_for_expected(Nemp, Ee);
[ObsM, ExpM] = merge_bins_for_expected(Nemp, EM);

chi2g = sum((ObsG - ExpG).^2 ./ max(ExpG, eps));
chi2e = sum((ObsE - ExpE).^2 ./ max(ExpE, eps));
chi2M = sum((ObsM - ExpM).^2 ./ max(ExpM, eps));

% Degrees of freedom
df_g = max(length(ObsG) - 1 - 2, 1);  % Normal: 2 params (mu,sigma)
df_e = max(length(ObsE) - 1 - 2, 1);  % Laplace: 2 params (mu,sigma)
df_M = max(length(ObsM) - 1 - 1, 1);  % MB: 1 param (a)

% p-values
pvg = 1 - chi2cdf(chi2g, df_g);
pve = 1 - chi2cdf(chi2e, df_e);
pvM = 1 - chi2cdf(chi2M, df_M);

% --- Plot theoretical fits ---
xgrid = linspace(min(y), max(y), 500);
fg = normpdf(xgrid, mu, sigma);
fe = (1/(2*sigma)) * exp(-abs(xgrid-mu)/sigma);
fMvals = arrayfun(fM, xgrid);

plot(xgrid, N*fg*Dy, 'r', 'LineWidth', 1.5);
plot(xgrid, N*fe*Dy, 'g', 'LineWidth', 1.5);
plot(xgrid, N*fMvals*Dy, 'b', 'LineWidth', 1.5);

hold off;
end
