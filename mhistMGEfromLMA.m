function results = mhistMGEfromLMA(x, nbins, doplot)
% GOODNESS_OF_FIT - Ocena dopasowania danych do rozkładów:
% Normalnego, Double Exponential (Laplace), wykładniczego oraz Maxwella-Boltzmanna.
%
% Wejście:
%   x      : dane eksperymentalne (wektor)
%   nbins  : liczba przedziałów histogramu (np. 30)
%   doplot : true/false – czy rysować wykres porównawczy
%
% Wyjście:
%   results : struktura zawierająca oszacowane parametry oraz
%             wartości p (p-values) testów zgodności:
%             - test chi-kwadrat (p_chi2)
%
% Autor: [Twoje imię/nazwisko]

    if nargin < 2
        nbins = 30;
    end
    if nargin < 3
        doplot = true;
    end

    % --- Histogram do testu chi2 ---
    [counts, edges] = histcounts(x, nbins);
    mids = (edges(1:end-1) + edges(2:end))/2;
    n = numel(x);
    bin_width = edges(2) - edges(1);

    % --- Rozkład normalny ---
    mu = mean(x);
    sigma = std(x);
    pd_norm = makedist('Normal', 'mu', mu, 'sigma', sigma);

    expected = n * pdf(pd_norm, mids) * bin_width;
    [~, p_chi2_norm] = chi2gof(x, 'Ctrs', mids, 'Expected', expected, 'NParams', 2);
    results.Normal.mu = mu;
    results.Normal.sigma = sigma;
    results.Normal.p_chi2 = p_chi2_norm;

    % --- Rozkład Double Exponential (Laplace) ---
    mu_lap = median(x);
    b_lap = mean(abs(x - mu_lap));
    laplace_pdf = @(z) (1/(2*b_lap)) * exp(-abs(z - mu_lap)/b_lap);

    expected = n * laplace_pdf(mids) * bin_width;
    % Zapobiega zerowym wartościom expected
    expected(expected==0) = eps;
    [~, p_chi2_lap] = chi2gof(x, 'Ctrs', mids, 'Expected', expected, 'NParams', 2);

    results.DoubleExponential.mu = mu_lap;
    results.DoubleExponential.b = b_lap;
    results.DoubleExponential.p_chi2 = p_chi2_lap;

    % --- Rozkład Maxwella-Boltzmanna ---
    maxwell_pdf = @(a, z) sqrt(2/pi) * (z.^2 ./ a.^3) .* exp(-z.^2 ./ (2*a.^2));
    mle_fun = @(a) -sum(log(maxwell_pdf(a, x(x>0))));
    scale = fminsearch(mle_fun, std(x));

    expected = n * maxwell_pdf(scale, mids) * bin_width;
    [~, p_chi2_max] = chi2gof(x, 'Ctrs', mids, 'Expected', expected, 'NParams', 1);

    results.Maxwell.scale = scale;
    results.Maxwell.p_chi2 = p_chi2_max;

    % --- Wykres porównawczy ---
    if doplot
        figure; hold on;
        histogram(x, nbins, 'Normalization', 'pdf', 'FaceAlpha', 0.5);

        xx = linspace(min(x), max(x), 500);
        plot(xx, pdf(pd_norm, xx), 'r--', 'LineWidth', 2);
        plot(xx, maxwell_pdf(scale, xx), 'b--', 'LineWidth', 2);
        plot(xx, laplace_pdf(xx), 'm--', 'LineWidth', 2);

        legend_entries = {'Dane', 'Normalny'};
        legend_entries = [legend_entries, {'Maxwell-Boltzmann','Double Exponential'}];
        legend(legend_entries{:});

        title('Goodness of Fit - porównanie rozkładów');
        hold off;
    end
end