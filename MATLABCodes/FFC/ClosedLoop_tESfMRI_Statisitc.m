%% ===================== Connectivity: Statistical Analysis Pipeline =====================
% This script assumes your CSV has columns: [ScanID, TR1, TR2, Test] with 15 rows per ScanID.
% It reproduces a behavioral-style pipeline for connectivity:
% - Subject-level mean across 15 connections at TR1/TR2/Test
% - Training = mean(TR1,TR2); DeltaTest = Test - Training
% - Normality, within-group tests, between-group tests, ANCOVA
% - Effect sizes (Cohen's d, Hedges' g) + 95% CI
% - Optional permutation test for between-group DeltaTest difference
% - Save summary table to CSV

clearvars -except tbl TR1_exp TR2_exp Test_exp TR1_ctrl TR2_ctrl Test_ctrl expScans ctrlScans
close all; clc

%% -------------------- User Inputs --------------------
filePath = '/Users/ghazaleh/Downloads/Connectivity.csv';

% Group assignments (as you used)
expScans = {'E16426','E16502','E16520','E16537','E16561','E16591','E16616','E16844','E16979','E17009'};
ctrlScans = {'E16449','E16488','E16518','E16531','E16599','E16628','E16682','E16909','E16942'};

% Analysis options
doPermutation   = true;     % run permutation test for ΔTest between-groups
Nperm           = 20000;    % permutation iterations
alpha           = 0.05;     % significance threshold
outCSV          = '/Users/ghazaleh/Downloads/Connectivity_Stats_Summary.csv';

%% -------------------- Load / Parse (robust to re-runs) --------------------
if ~exist('tbl','var') || isempty(tbl)
    tbl = readtable(filePath);
end

allScans = unique(string(tbl{:,1}));
nConn = 15;

TR1_exp = []; TR2_exp = []; Test_exp = [];
TR1_ctrl = []; TR2_ctrl = []; Test_ctrl = [];
okExp = {}; okCtrl = {};

for i = 1:numel(allScans)
    scanID = string(allScans(i));
    idx = strcmp(tbl{:,1}, scanID);  % rows for this scan

    this_TR1  = tbl{idx, 2};
    this_TR2  = tbl{idx, 3};
    this_Test = tbl{idx, 4};

    if numel(this_TR1) ~= nConn
        warning('Skipping %s: expected %d rows, got %d', scanID, nConn, numel(this_TR1));
        continue;
    end

    if ismember(scanID, expScans)
        TR1_exp(end+1,:)  = this_TR1';
        TR2_exp(end+1,:)  = this_TR2';
        Test_exp(end+1,:) = this_Test';
        okExp{end+1,1} = char(scanID);
    elseif ismember(scanID, ctrlScans)
        TR1_ctrl(end+1,:)  = this_TR1';
        TR2_ctrl(end+1,:)  = this_TR2';
        Test_ctrl(end+1,:) = this_Test';
        okCtrl{end+1,1} = char(scanID);
    else
        warning('Scan ID %s not in group lists', scanID);
    end
end

okExp  = string(okExp);
okCtrl = string(okCtrl);

%% -------------------- Subject-level means (across 15 connections) --------------------
mTR1_exp  = mean(TR1_exp,  2, 'omitnan');
mTR2_exp  = mean(TR2_exp,  2, 'omitnan');
mTest_exp = mean(Test_exp, 2, 'omitnan');

mTR1_ctrl  = mean(TR1_ctrl,  2, 'omitnan');
mTR2_ctrl  = mean(TR2_ctrl,  2, 'omitnan');
mTest_ctrl = mean(Test_ctrl, 2, 'omitnan');

% Pack as matrices [nSubj x 3] columns = [TR1 TR2 Test]
exp_raw  = [mTR1_exp,  mTR2_exp,  mTest_exp ];
ctrl_raw = [mTR1_ctrl, mTR2_ctrl, mTest_ctrl];

% Drop subjects with NaN across any timepoint
valid_exp  = all(~isnan(exp_raw),  2);
valid_ctrl = all(~isnan(ctrl_raw), 2);
exp_raw  = exp_raw(valid_exp,  :);   okExp  = okExp(valid_exp);
ctrl_raw = ctrl_raw(valid_ctrl, :);   okCtrl = okCtrl(valid_ctrl);

if isempty(exp_raw) || isempty(ctrl_raw)
    error('No valid subjects after NaN filtering.');
end

% Training averages & ΔTest
train_exp  = mean(exp_raw(:,1:2), 2);
train_ctrl = mean(ctrl_raw(:,1:2), 2);
delta_exp  = exp_raw(:,3)  - train_exp;    % Test - Training
delta_ctrl = ctrl_raw(:,3) - train_ctrl;

%% -------------------- Normality checks (Lilliefors) --------------------
fprintf('\n=== NORMALITY (Lilliefors) ===\n');
phaseNames = {'TR1','TR2','Test','Training','DeltaTest'};
toCheck_exp  = {exp_raw(:,1),exp_raw(:,2),exp_raw(:,3),train_exp,delta_exp};
toCheck_ctrl = {ctrl_raw(:,1),ctrl_raw(:,2),ctrl_raw(:,3),train_ctrl,delta_ctrl};

for k = 1:numel(phaseNames)
    pE = lillietest(toCheck_exp{k});
    pC = lillietest(toCheck_ctrl{k});
    fprintf('%-9s | Exp p=%.4f (%s)  | Ctrl p=%.4f (%s)\n', phaseNames{k}, ...
        pE, ternary(pE>alpha,'normal','non-normal'), ...
        pC, ternary(pC>alpha,'normal','non-normal'));
end

%% -------------------- Within-group tests --------------------
% 1) Paired TR1 vs Test, TR2 vs Test (paired t + nonparametric)
fprintf('\n=== WITHIN-GROUP TESTS ===\n');

[~,p_t_exp_TR1Test,~,stats_t_exp_TR1Test] = ttest(exp_raw(:,1), exp_raw(:,3));
[~,p_t_exp_TR2Test,~,stats_t_exp_TR2Test] = ttest(exp_raw(:,2), exp_raw(:,3));
[~,p_t_ctrl_TR1Test,~,stats_t_ctrl_TR1Test] = ttest(ctrl_raw(:,1), ctrl_raw(:,3));
[~,p_t_ctrl_TR2Test,~,stats_t_ctrl_TR2Test] = ttest(ctrl_raw(:,2), ctrl_raw(:,3));

p_w_exp_TR1Test   = signrank(exp_raw(:,1), exp_raw(:,3));
p_w_exp_TR2Test   = signrank(exp_raw(:,2), exp_raw(:,3));
p_w_ctrl_TR1Test  = signrank(ctrl_raw(:,1), ctrl_raw(:,3));
p_w_ctrl_TR2Test  = signrank(ctrl_raw(:,2), ctrl_raw(:,3));

% 2) One-sample test on ΔTest vs 0
[~,p_t_exp_delta,~,stats_t_exp_delta] = ttest(delta_exp, 0);
[~,p_t_ctrl_delta,~,stats_t_ctrl_delta] = ttest(delta_ctrl, 0);
p_w_exp_delta  = signrank(delta_exp, 0);
p_w_ctrl_delta = signrank(delta_ctrl, 0);

% Effect sizes for paired differences (Cohen's dz; Hedges' g for paired)
dz_exp_TR1Test  = cohen_d_paired(exp_raw(:,3), exp_raw(:,1));  % Test - TR1
dz_exp_TR2Test  = cohen_d_paired(exp_raw(:,3), exp_raw(:,2));  % Test - TR2
dz_ctrl_TR1Test = cohen_d_paired(ctrl_raw(:,3), ctrl_raw(:,1));
dz_ctrl_TR2Test = cohen_d_paired(ctrl_raw(:,3), ctrl_raw(:,2));
g_exp_delta     = hedges_g_onesample(delta_exp, 0);
g_ctrl_delta    = hedges_g_onesample(delta_ctrl, 0);

fprintf('Exp: TR1 vs Test  t(%d)=%.3f, p=%.4f | Wilcoxon p=%.4f | dz=%.3f\n', ...
    stats_t_exp_TR1Test.df, stats_t_exp_TR1Test.tstat, p_t_exp_TR1Test, p_w_exp_TR1Test, dz_exp_TR1Test);
fprintf('Exp: TR2 vs Test  t(%d)=%.3f, p=%.4f | Wilcoxon p=%.4f | dz=%.3f\n', ...
    stats_t_exp_TR2Test.df, stats_t_exp_TR2Test.tstat, p_t_exp_TR2Test, p_w_exp_TR2Test, dz_exp_TR2Test);
fprintf('Ctrl: TR1 vs Test t(%d)=%.3f, p=%.4f | Wilcoxon p=%.4f | dz=%.3f\n', ...
    stats_t_ctrl_TR1Test.df, stats_t_ctrl_TR1Test.tstat, p_t_ctrl_TR1Test, p_w_ctrl_TR1Test, dz_ctrl_TR1Test);
fprintf('Ctrl: TR2 vs Test t(%d)=%.3f, p=%.4f | Wilcoxon p=%.4f | dz=%.3f\n', ...
    stats_t_ctrl_TR2Test.df, stats_t_ctrl_TR2Test.tstat, p_t_ctrl_TR2Test, p_w_ctrl_TR2Test, dz_ctrl_TR2Test);

fprintf('Exp: ΔTest one-sample t(%d)=%.3f, p=%.4f | signrank p=%.4f | Hedges g=%.3f\n', ...
    stats_t_exp_delta.df, stats_t_exp_delta.tstat, p_t_exp_delta, p_w_exp_delta, g_exp_delta.g);
fprintf('Ctrl: ΔTest one-sample t(%d)=%.3f, p=%.4f | signrank p=%.4f | Hedges g=%.3f\n', ...
    stats_t_ctrl_delta.df, stats_t_ctrl_delta.tstat, p_t_ctrl_delta, p_w_ctrl_delta, g_ctrl_delta.g);

%% -------------------- Between-group tests on ΔTest --------------------
fprintf('\n=== BETWEEN-GROUP: ΔTest (Test - Training) ===\n');

% t-test
[~,p_t_2,~,stats_t_2] = ttest2(delta_exp, delta_ctrl, 'Vartype','unequal');
d_ind = cohen_d_independent(delta_exp, delta_ctrl);
g_ind = hedges_g_independent(delta_exp, delta_ctrl);

% rank-sum
p_rs = ranksum(delta_exp, delta_ctrl);

% 95% CI for mean difference
diffMean = mean(delta_exp) - mean(delta_ctrl);
se_diff  = sqrt(var(delta_exp)/numel(delta_exp) + var(delta_ctrl)/numel(delta_ctrl));
df_satter = satterthwaite_df(var(delta_exp), numel(delta_exp), var(delta_ctrl), numel(delta_ctrl));
tcrit = tinv(1 - alpha/2, df_satter);
ci_diff = diffMean + tcrit*[-1 1]*se_diff;

fprintf('ΔTest: t(%0.1f)=%.3f, p=%.4f  | rank-sum p=%.4f | mean diff=%.4f, 95%% CI [%.4f, %.4f]\n', ...
    df_satter, stats_t_2.tstat, p_t_2, p_rs, diffMean, ci_diff(1), ci_diff(2));
fprintf('Effect sizes: Cohen d=%.3f, Hedges g=%.3f\n', d_ind, g_ind.g);

% Optional permutation
if doPermutation
    p_perm = permtest_meandiff(delta_exp, delta_ctrl, Nperm);
    fprintf('Permutation (N=%d) p=%.4f for |mean(Exp) - mean(Ctrl)|\n', Nperm, p_perm);
else
    p_perm = NaN;
end

%% -------------------- ANCOVA: Test ~ Group + Training --------------------
% Stack rows: [exp; ctrl]
Test_all   = [exp_raw(:,3); ctrl_raw(:,3)];
Train_all  = [train_exp;     train_ctrl];
Group_all  = [ones(size(exp_raw,1),1); zeros(size(ctrl_raw,1),1)];  % 1=Exp, 0=Ctrl

T = table(Test_all, Train_all, Group_all, 'VariableNames', {'Test','Training','Group'});
lm = fitlm(T, 'Test ~ Group + Training'); % OLS ANCOVA

coefs = lm.Coefficients;
p_group = coefs.pValue(strcmp(coefs.Properties.RowNames,'Group'));
beta_group = coefs.Estimate(strcmp(coefs.Properties.RowNames,'Group'));
ci_group   = coefCI(lm, alpha);
ci_idx     = find(strcmp(coefs.Properties.RowNames,'Group'));
ci_grp     = ci_group(ci_idx,:);

fprintf('\n=== ANCOVA (Test ~ Group + Training) ===\n');
disp(lm)
fprintf('Group coefficient: %.4f, 95%% CI [%.4f, %.4f], p=%.4f\n', beta_group, ci_grp(1), ci_grp(2), p_group);

%% -------------------- Summarize & Save --------------------
semfun = @(x) std(x)/sqrt(numel(x));
summary_rows = {
    'Training mean (Exp)', mean(train_exp), semfun(train_exp), numel(train_exp);
    'Training mean (Ctrl)', mean(train_ctrl), semfun(train_ctrl), numel(train_ctrl);
    'Test mean (Exp)', mean(exp_raw(:,3)), semfun(exp_raw(:,3)), numel(exp_raw(:,3));
    'Test mean (Ctrl)', mean(ctrl_raw(:,3)), semfun(ctrl_raw(:,3)), numel(ctrl_raw(:,3));
    'ΔTest mean (Exp)', mean(delta_exp), semfun(delta_exp), numel(delta_exp);
    'ΔTest mean (Ctrl)', mean(delta_ctrl), semfun(delta_ctrl), numel(delta_ctrl);
    'ΔTest mean diff (Exp-Ctrl)', diffMean, se_diff, NaN;
    };

stats_tbl = cell2table(summary_rows, ...
    'VariableNames', {'Metric','Mean','SEM','N'});

% Add key p-values and effects as a separate table to save alongside
keyvals = table( ...
    p_t_2, p_rs, p_perm, d_ind, g_ind.g, ...
    'VariableNames', {'p_t_betweenDelta','p_ranksum_betweenDelta','p_perm_betweenDelta','Cohen_d_between','Hedges_g_between'});

% Write both (two sheets if Excel, or two CSVs if preferred)
try
    writetable(stats_tbl, outCSV);
    [pth,base,~] = fileparts(outCSV);
    outCSV2 = fullfile(pth, [base '_KeyTests.csv']);
    writetable(keyvals, outCSV2);
    fprintf('\nSaved: %s and %s\n', outCSV, outCSV2);
catch ME
    warning('Could not write CSV: %s', ME.message);
end

%% ===================== Helper Functions =====================
function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end

function d = cohen_d_paired(x2, x1)
% Cohen's dz for paired samples: mean(diff)/std(diff)
d = (mean(x2 - x1, 'omitnan')) / std(x2 - x1, 0, 'omitnan');
end

function d = cohen_d_independent(x, y)
% Pooled SD Cohen's d (independent)
nx = numel(x); ny = numel(y);
sx = var(x, 'omitnan'); sy = var(y, 'omitnan');
sp2 = ((nx-1)*sx + (ny-1)*sy) / (nx+ny-2);
d = (mean(x,'omitnan') - mean(y,'omitnan')) / sqrt(sp2);
end

function df = satterthwaite_df(vx, nx, vy, ny)
% Satterthwaite approximate df for unequal variances
num = (vx/nx + vy/ny)^2;
den = (vx^2/((nx^2)*(nx-1))) + (vy^2/((ny^2)*(ny-1)));
df = num / den;
end

function g = hedges_g_onesample(x, mu0)
% Hedges g for one-sample difference (x - mu0)
n  = numel(x);
dx = mean(x) - mu0;
sx = std(x, 0);
d  = dx / sx;
J  = 1 - (3/(4*n - 9)); % small sample correction
g.g = J * d;
g.d = d;
end

function g = hedges_g_independent(x, y)
% Hedges g for independent groups
nx = numel(x); ny = numel(y);
sx2 = var(x, 'omitnan'); sy2 = var(y, 'omitnan');
sp2 = ((nx-1)*sx2 + (ny-1)*sy2) / (nx+ny-2);
d = (mean(x,'omitnan') - mean(y,'omitnan')) / sqrt(sp2);
J = 1 - (3/(4*(nx+ny) - 9));
g.g = J * d;
g.d = d;
end

function p = permtest_meandiff(x, y, N)
% Two-sided permutation test for difference in means |mean(x)-mean(y)|
x = x(:); y = y(:);
obs = abs(mean(x) - mean(y));
XY = [x; y];
nx = numel(x);
count = 0;
for i = 1:N
    idx = randperm(numel(XY));
    x_perm = XY(idx(1:nx));
    y_perm = XY(idx(nx+1:end));
    if abs(mean(x_perm) - mean(y_perm)) >= obs
        count = count + 1;
    end
end
p = (count + 1) / (N + 1); % add-one correction
end
