function results = mhistMGEfromLMA(x, nbins, doplot)
    % GOODNESS_OF_FIT - Ocena dopasowania danych do rozkładów:
    % Normalnego, wykładniczego oraz Maxwella-Boltzmanna.
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
        nbins = 30; % domyślna liczba przedziałów histogramu
    end
    if nargin < 3
        doplot = true; % domyślnie rysujemy wykres
    end

    % --- Histogram do testu chi2 ---
    [counts, edges] = histcounts(x, nbins); % zliczenia w koszykach
    mids = (edges(1:end-1) + edges(2:end)) / 2; % środki przedziałów
    n = numel(x); % liczba obserwacji
    bin_width = edges(2) - edges(1); % szerokość przedziału

    % --- Rozkład normalny ---
    mu = mean(x); % średnia
    sigma = std(x); % odchylenie standardowe
    pd_norm = makedist('Normal', 'mu', mu, 'sigma', sigma);

    expected = n * pdf(pd_norm, mids) * bin_width;
    [~, p_chi2_norm] = chi2gof(x, 'Ctrs', mids, 'Expected', expected, 'NParams', 2);

    results.Normal.mu = mu;
    results.Normal.sigma = sigma;
    results.Normal.p_chi2 = p_chi2_norm;  

    % --- Rozkład wykładniczy ---
    pd_exp = fitdist(x(:), 'Exponential'); % dopasowanie MLE
    expected = n * pdf(pd_exp, mids) * bin_width;
    [~, p_chi2_exp] = chi2gof(x, 'Ctrs', mids, 'Expected', expected, 'NParams', 1);

    results.Exponential.lambda = 1 / pd_exp.mu;
    results.Exponential.p_chi2 = p_chi2_exp;

    % --- Rozkład Maxwella-Boltzmanna ---
    % PDF: f(x;a) = sqrt(2/pi) * (x^2 / a^3) * exp(-x^2 / (2a^2)), x>0
    maxwell_pdf = @(a, z) sqrt(2/pi) * (z.^2 ./ a.^3) .* exp(-z.^2 ./ (2*a.^2));
    mle_fun = @(a) -sum(log(maxwell_pdf(a, x(x>0)))); % log-likelihood
    scale = fminsearch(mle_fun, std(x)); % estymacja parametru a

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
        plot(xx, pdf(pd_exp, xx), 'g--', 'LineWidth', 2);
        plot(xx, maxwell_pdf(scale, xx), 'b--', 'LineWidth', 2);

        legend('Dane', 'Normalny', 'Wykładniczy', 'Maxwell-Boltzmann');
        title('Goodness of Fit - porównanie rozkładów');
        hold off;
    end
end