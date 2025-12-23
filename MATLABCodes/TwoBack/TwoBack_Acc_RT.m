% Load data and treat empty cells as NaN
clear all
close all
clc 
filePath = '/Users/ghazaleh/Downloads/2-BackResults.xlsx';
data = readmatrix(filePath, 'TreatAsMissing', '');

% Validate size
[nRows, nCols] = size(data);
if nRows ~= 120 || mod(nCols, 6) ~= 0
    error('Expected data size of 120 rows and a multiple of 6 columns. Got %d x %d.', nRows, nCols);
end

nSubjects = nCols / 6;

% Preallocate
reshaped_all = [];

for i = 1:nSubjects
    cols = (i-1)*6 + (1:6);  % Select 6 columns for this subject
    subject_block = data(:, cols);  % 120x6
    reshaped_all = [reshaped_all; subject_block];  % Stack downward
end

% Final matrix: 120 rows × 6 columns × 19 subjects → 2280 × 6
disp(size(reshaped_all));  % Should show [2280, 6]


% List of Scan IDs (order should match the column order in your Excel file)
scanIDs = {'E16426','E16449','E16488','E16502','E16518','E16520','E16531','E16537','E16561','E16591',...
           'E16599','E16616','E16628','E16682','E16844','E16909','E16942','E16979','E17009'};

nSubjects = numel(scanIDs);
rowsPerSubject = 120;

% Create repeated ScanID labels
ScanID_col = repelem(scanIDs, rowsPerSubject)';  % 2280 × 1 cell array

% Combine with your reshaped matrix
finalTable = table(ScanID_col, ...
                   reshaped_all(:,1), reshaped_all(:,2), reshaped_all(:,3), ...
                   reshaped_all(:,4), reshaped_all(:,5), reshaped_all(:,6), ...
                   'VariableNames', {'ScanID','Col1','Col2','Col3','Col4','Col5','Col6'});

finalTable.Properties.VariableNames = { ...
    'ScanID', ...
    'TR1', ...
    'TR1_timing', ...
    'TR2', ...
    'TR2_timing', ...
    'Test', ...
    'Test_timing' ...
};


%% Plot Accuracy and Reaction Time 


% Assumes finalTable has: ScanID, TR1, TR1_timing, TR2, TR2_timing, Test, Test_timing

% Define subject groups
expScans = {'E16426','E16502','E16520','E16537','E16561','E16591','E16616','E16844','E16979','E17009'};
ctrlScans = {'E16449','E16488','E16518','E16531','E16599','E16628','E16682','E16909','E16942'};

phases = {'TR1', 'TR2', 'Test'};
timing_phases = {'TR1_timing', 'TR2_timing', 'Test_timing'};
nPhases = numel(phases);

% Preallocate
acc_exp = nan(nPhases, 1); rt_exp = nan(nPhases, 1);
acc_ctrl = nan(nPhases, 1); rt_ctrl = nan(nPhases, 1);
sem_acc_exp = nan(nPhases, 1); sem_rt_exp = nan(nPhases, 1);
sem_acc_ctrl = nan(nPhases, 1); sem_rt_ctrl = nan(nPhases, 1);

for i = 1:nPhases
    % Experimental
    exp_idx = ismember(finalTable.ScanID, expScans);
    acc_vals_exp = finalTable{exp_idx, phases{i}};
    rt_vals_exp = finalTable{exp_idx, timing_phases{i}};
    acc_exp(i) = mean(acc_vals_exp, 'omitnan');
    rt_exp(i) = mean(rt_vals_exp, 'omitnan');
    sem_acc_exp(i) = std(acc_vals_exp, 'omitnan') / sqrt(sum(~isnan(acc_vals_exp)));
    sem_rt_exp(i) = std(rt_vals_exp, 'omitnan') / sqrt(sum(~isnan(rt_vals_exp)));

    % Control
    ctrl_idx = ismember(finalTable.ScanID, ctrlScans);
    acc_vals_ctrl = finalTable{ctrl_idx, phases{i}};
    rt_vals_ctrl = finalTable{ctrl_idx, timing_phases{i}};
    acc_ctrl(i) = mean(acc_vals_ctrl, 'omitnan');
    rt_ctrl(i) = mean(rt_vals_ctrl, 'omitnan');
    sem_acc_ctrl(i) = std(acc_vals_ctrl, 'omitnan') / sqrt(sum(~isnan(acc_vals_ctrl)));
    sem_rt_ctrl(i) = std(rt_vals_ctrl, 'omitnan') / sqrt(sum(~isnan(rt_vals_ctrl)));
end



% Define group colors
expColor = [0.1 0.7 0.7];
ctrlColor = [0.8 0.2 0.2];

% Plot
x = 1:nPhases;
labels = {'TR1', 'TR2', 'Test'};

figure;

% Accuracy Plot
subplot(1,2,1); hold on;
errorbar(x, acc_ctrl, sem_acc_ctrl, '-o', ...
    'Color', ctrlColor, 'LineWidth', 2, 'MarkerFaceColor', ctrlColor, 'MarkerSize', 8);
errorbar(x, acc_exp, sem_acc_exp, '-^', ...
    'Color', expColor, 'LineWidth', 2, 'MarkerFaceColor', expColor, 'MarkerSize', 8);
xlim([0.8 3.2]); ylim([0 1]);
xticks(x); xticklabels(labels);
ylabel('Accuracy');
title('2-Back Accuracy');
legend('Down-Regulation','Up-Regulation','Location','southwest');

% Reaction Time Plot
subplot(1,2,2); hold on;
errorbar(x, rt_ctrl, sem_rt_ctrl, '-o', ...
    'Color', ctrlColor, 'LineWidth', 2, 'MarkerFaceColor', ctrlColor, 'MarkerSize', 8);
errorbar(x, rt_exp, sem_rt_exp, '-^', ...
    'Color', expColor, 'LineWidth', 2, 'MarkerFaceColor', expColor, 'MarkerSize', 8);
xlim([0.8 3.2]);
xticks(x); xticklabels(labels);
ylabel('Reaction Time (s)');
title('2-Back Reaction Time');
legend('Down-Regulation','Up-Regulation','Location','northwest');

sgtitle('2-Back Task Performance by Group');


%% ===================== Statistical analyses (subject-level) =====================
% Goal:
%   1) Phase (TR1, TR2, Test) × Group (Down vs Up) for Accuracy and RT
%   2) Planned contrasts:
%       a) within each group: Test vs Baseline(mean(TR1,TR2))
%       b) between groups: delta = Test - Baseline

% ---- 1) Make subject-level means (one row per ScanID) ----
scanList = unique(finalTable.ScanID, 'stable');
nS = numel(scanList);

Subj = table('Size',[nS 0],'VariableTypes',{},'VariableNames',{});
Subj.ScanID = scanList;

% define group label by ScanID membership
Subj.Group = repmat("NA", nS, 1);
Subj.Group(ismember(Subj.ScanID, ctrlScans)) = "Down-Regulation";
Subj.Group(ismember(Subj.ScanID, expScans))  = "Up-Regulation";

% subject-level means per phase
Subj.TR1_acc  = nan(nS,1); Subj.TR2_acc  = nan(nS,1); Subj.Test_acc  = nan(nS,1);
Subj.TR1_rt   = nan(nS,1); Subj.TR2_rt   = nan(nS,1); Subj.Test_rt   = nan(nS,1);

for s = 1:nS
    idx = strcmp(finalTable.ScanID, Subj.ScanID{s});
    Subj.TR1_acc(s)  = mean(finalTable.TR1(idx),         'omitnan');
    Subj.TR2_acc(s)  = mean(finalTable.TR2(idx),         'omitnan');
    Subj.Test_acc(s) = mean(finalTable.Test(idx),        'omitnan');

    Subj.TR1_rt(s)   = mean(finalTable.TR1_timing(idx),  'omitnan');
    Subj.TR2_rt(s)   = mean(finalTable.TR2_timing(idx),  'omitnan');
    Subj.Test_rt(s)  = mean(finalTable.Test_timing(idx), 'omitnan');
end

% Baseline = mean(TR1,TR2)
Subj.Base_acc = mean([Subj.TR1_acc Subj.TR2_acc], 2, 'omitnan');
Subj.Base_rt  = mean([Subj.TR1_rt  Subj.TR2_rt ], 2, 'omitnan');

% Deltas
Subj.Delta_acc = Subj.Test_acc - Subj.Base_acc;
Subj.Delta_rt  = Subj.Test_rt  - Subj.Base_rt;

% Keep only participants with known group
keep = Subj.Group ~= "NA";
Subj = Subj(keep,:);

fprintf('\n===== SUBJECT-LEVEL SUMMARY (used for stats) =====\n');
disp(groupsummary(Subj, "Group", "mean", ["TR1_acc","TR2_acc","Test_acc","Base_acc","Delta_acc","TR1_rt","TR2_rt","Test_rt","Base_rt","Delta_rt"]));
vars = ["TR1_acc","TR2_acc","Test_acc","Base_acc","Delta_acc", ...
        "TR1_rt","TR2_rt","Test_rt","Base_rt","Delta_rt"];

% Mean and SD (safe)
G_mean = groupsummary(Subj, "Group", "mean", vars);
G_sd   = groupsummary(Subj, "Group", "std",  vars);

% Rename columns
G_mean.Properties.VariableNames = strrep(G_mean.Properties.VariableNames,"mean_","Mean_");
G_sd.Properties.VariableNames   = strrep(G_sd.Properties.VariableNames,  "std_", "SD_");

% ---- Compute N manually (robust) ----
Groups = categories(categorical(Subj.Group));
Ntbl = table(Groups, 'VariableNames', {'Group'});

for v = vars
    Nv = zeros(numel(Groups),1);
    for g = 1:numel(Groups)
        idx = Subj.Group == Groups{g};
        Nv(g) = sum(~isnan(Subj{idx, v}));
    end
    Ntbl.("N_"+v) = Nv;
end

% ---- Merge everything ----
SummaryStats = G_mean(:,["Group","Mean_"+vars]);
SummaryStats = [SummaryStats, G_sd(:, "SD_"+vars)];
SummaryStats = [SummaryStats, Ntbl(:, "N_"+vars)];

fprintf('\n===== SUBJECT-LEVEL SUMMARY (used for stats) =====\n');
disp(SummaryStats);

%% ---- 2) Phase × Group model (Accuracy) ----
% Use repeated-measures ANOVA (fitrm/ranova) if available; otherwise fall back to mixed model.
hasFitrm = exist('fitrm','file')==2 && exist('ranova','file')==2;

if hasFitrm
    % Accuracy RM-ANOVA
    Tacc = Subj(:, {'Group','TR1_acc','TR2_acc','Test_acc'});
    Tacc.Group = categorical(Tacc.Group);

    within = table(categorical(["TR1";"TR2";"Test"]), 'VariableNames', {'Phase'});
    rm_acc = fitrm(Tacc, 'TR1_acc- Test_acc ~ Group', 'WithinDesign', within);
    fprintf('\n===== Accuracy: Phase × Group (RM-ANOVA) =====\n');
    disp(ranova(rm_acc, 'WithinModel','Phase'));

    % RT RM-ANOVA
    Trt = Subj(:, {'Group','TR1_rt','TR2_rt','Test_rt'});
    Trt.Group = categorical(Trt.Group);

    rm_rt = fitrm(Trt, 'TR1_rt- Test_rt ~ Group', 'WithinDesign', within);
    fprintf('\n===== RT: Phase × Group (RM-ANOVA) =====\n');
    disp(ranova(rm_rt, 'WithinModel','Phase'));

else
    % Mixed-effects fallback (requires Statistics Toolbox; more robust than nothing)
    % Long format
    ScanID = string(Subj.ScanID);
    Group  = categorical(Subj.Group);

    longAcc = table;
    longAcc.ScanID = [ScanID; ScanID; ScanID];
    longAcc.Group  = [Group;  Group;  Group];
    longAcc.Phase  = categorical([repmat("TR1",nS,1); repmat("TR2",nS,1); repmat("Test",nS,1)]);
    longAcc.Acc    = [Subj.TR1_acc; Subj.TR2_acc; Subj.Test_acc];

    longRT = table;
    longRT.ScanID = [ScanID; ScanID; ScanID];
    longRT.Group  = [Group;  Group;  Group];
    longRT.Phase  = categorical([repmat("TR1",nS,1); repmat("TR2",nS,1); repmat("Test",nS,1)]);
    longRT.RT     = [Subj.TR1_rt; Subj.TR2_rt; Subj.Test_rt];

    fprintf('\n===== Accuracy: Phase × Group (LME fallback) =====\n');
    mAcc = fitlme(longAcc, 'Acc ~ Group*Phase + (1|ScanID)');
    disp(anova(mAcc));

    fprintf('\n===== RT: Phase × Group (LME fallback) =====\n');
    mRT = fitlme(longRT, 'RT ~ Group*Phase + (1|ScanID)');
    disp(anova(mRT));
end

%% ---- 3) Planned contrasts ----
% Within-group: Test vs Baseline (paired t-test; if non-normal, switch to signrank)
fprintf('\n===== Planned contrasts: within-group Test vs Baseline =====\n');

grpNames = categories(categorical(Subj.Group));
for g = 1:numel(grpNames)
    gg = grpNames{g};
    idx = (Subj.Group == string(gg));

    % Accuracy
    x = Subj.Test_acc(idx);
    y = Subj.Base_acc(idx);
    ok = ~isnan(x) & ~isnan(y);

    if sum(ok) >= 3
        [~,p,~,st] = ttest(x(ok), y(ok));  % paired
        fprintf('Accuracy (%s): Test vs Base: n=%d, t=%.3f, df=%d, p=%.4g\n', gg, sum(ok), st.tstat, st.df, p);
    else
        fprintf('Accuracy (%s): not enough data for paired test\n', gg);
    end

    % RT
    x = Subj.Test_rt(idx);
    y = Subj.Base_rt(idx);
    ok = ~isnan(x) & ~isnan(y);

    if sum(ok) >= 3
        [~,p,~,st] = ttest(x(ok), y(ok));  % paired
        fprintf('RT (%s): Test vs Base: n=%d, t=%.3f, df=%d, p=%.4g\n', gg, sum(ok), st.tstat, st.df, p);
    else
        fprintf('RT (%s): not enough data for paired test\n', gg);
    end
end

% Between-group: compare Delta(Test-Base)
fprintf('\n===== Planned contrasts: between-group Δ(Test-Base) =====\n');

isDown = Subj.Group == "Down-Regulation";
isUp   = Subj.Group == "Up-Regulation";

% Accuracy delta (Welch)
x = Subj.Delta_acc(isDown); y = Subj.Delta_acc(isUp);
x = x(~isnan(x)); y = y(~isnan(y));
if numel(x) >= 2 && numel(y) >= 2
    [~,p,~,st] = ttest2(x,y,'Vartype','unequal');
    fprintf('ΔAccuracy: Down n=%d vs Up n=%d: t=%.3f, df=%.1f, p=%.4g\n', numel(x), numel(y), st.tstat, st.df, p);
else
    fprintf('ΔAccuracy: not enough data\n');
end

% RT delta (Welch)
x = Subj.Delta_rt(isDown); y = Subj.Delta_rt(isUp);
x = x(~isnan(x)); y = y(~isnan(y));
if numel(x) >= 2 && numel(y) >= 2
    [~,p,~,st] = ttest2(x,y,'Vartype','unequal');
    fprintf('ΔRT: Down n=%d vs Up n=%d: t=%.3f, df=%.1f, p=%.4g\n', numel(x), numel(y), st.tstat, st.df, p);
else
    fprintf('ΔRT: not enough data\n');
end
