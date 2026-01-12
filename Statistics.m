
clear all
clc

%% Example: Statistical comparison of two groups with normality check and visualization
rng(1); % For reproducibility
groupA = randn(50,1)*0.5 + 2;  % Simulated distances for Group A
groupB = randn(50,1)*0.5 + 2.8;  % Simulated distances for Group B

%% 1. Normality check using Kolmogorov-Smirnov
% Report the normality test used (e.g., KS test), p-values for each group,
% and whether data met normality assumptions.
[hA, pA] = kstest((groupA - mean(groupA))/std(groupA));
[hB, pB] = kstest((groupB - mean(groupB))/std(groupB));

fprintf('Group A Normality: %s (p=%.4f)\n', ternary(hA==0,'Normal','Not normal'), pA);
fprintf('Group B Normality: %s (p=%.4f)\n', ternary(hB==0,'Normal','Not normal'), pB);

%% 2. Normality plots
% Include histograms and QQ-plots as figures or supplementary material
% to visually support the normality assessment.
figure('Name','Normality Check','NumberTitle','off');
subplot(2,2,1);
histogram(groupA);
title('Histogram Group A');

subplot(2,2,2);
qqplot(groupA);
title('QQ-plot Group A');

subplot(2,2,3);
histogram(groupB);
title('Histogram Group B');

subplot(2,2,4);
qqplot(groupB);
title('QQ-plot Group B');


if hA==0 && hB==0
    % Normal data: mean ± SD
    fprintf('Group A: %.2f ± %.2f\n', mean(groupA), std(groupA));
    fprintf('Group B: %.2f ± %.2f\n', mean(groupB), std(groupB));
else
    % Non-normal data: median and IQR
    fprintf('Group A: median = %.2f [IQR: %.2f–%.2f]\n', median(groupA), prctile(groupA,25), prctile(groupA,75));
    fprintf('Group B: median = %.2f [IQR: %.2f–%.2f]\n', median(groupB), prctile(groupB,25), prctile(groupB,75));
end

%% 3. Hypothesis
% H0: The two groups have equal distributions (or means).
% H1: The two groups differ.
% Clearly state these hypotheses in the Methods section.

%% 4. Select appropriate test based on normality
% Report which test was chosen (t-test or Mann-Whitney),
% the test statistic, p-value, and confidence intervals.
if hA==0 && hB==0
    fprintf('\nBoth groups are normal -> Using independent t-test.\n');
    [h,p,ci,stats] = ttest2(groupA, groupB);
    fprintf('t-test result: t=%.4f, df=%d, p=%.4f\n', stats.tstat, stats.df, p);
    % In article: Include t-statistic, degrees of freedom, p-value, and CI.
else
    fprintf('\nAt least one group is not normal -> Using Mann-Whitney (ranksum).\n');
    [p,h,stats] = ranksum(groupA, groupB);
    fprintf('Mann-Whitney result: U=%.4f, p=%.4f\n', stats.ranksum, p);
    % In article: Include U statistic and p-value.
end

%% 5. Interpretation
% State whether H0 was rejected and interpret the result.
if p < 0.05
    fprintf('Conclusion: Reject H0. Groups are significantly different.\n');
else
    fprintf('Conclusion: Fail to reject H0. No significant difference.\n');
end

%%
if hA==0 && hB==0
    % Cohen's d for t-test
    meanA = mean(groupA);
    meanB = mean(groupB);
    sdA = std(groupA);
    sdB = std(groupB);
    nA = length(groupA);
    nB = length(groupB);
    pooledSD = sqrt(((nA-1)*sdA^2 + (nB-1)*sdB^2)/(nA+nB-2));
    cohens_d = (meanA - meanB)/pooledSD;
    fprintf('Effect size (Cohen''s d): %.4f\n', cohens_d);
else
    % Approximate r for Mann-Whitney
    N = length(groupA) + length(groupB);
    % MATLAB does not give Z directly, but you can compute it if needed
    % For simplicity, report ranksum statistic and N
    fprintf('Effect size r can be computed from Z-score manually.\n');
end

%% What to include:
% - Normality test results (p-values for each group).
% - Chosen statistical test and justification.
% - Test statistic, p-value, and confidence intervals.
% - Effect size (Cohen's d for t-test or r for Mann-Whitney).
% - Descriptive statistics: mean ± SD (if normal) or median [IQR] (if not normal).
% - Figures: histograms, QQ-plots, and boxplots for group comparison.

%% Helper function
function out = ternary(cond,trueVal,falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end
