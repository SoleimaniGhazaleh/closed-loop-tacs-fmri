%% ExploreModeration_TwoBack_Generalizability.m
% ------------------------------------------------------------
% What this script does (end-to-end):
% 1) Reads your 2-BackResults.xlsx (trial-level) and builds subject-level outcomes:
%       - PreAcc (mean of TR1/TR2), TestAcc, ΔAcc = Test - Pre
%       - PreRT  (mean of TR1/TR2 timing), TestRT, ΔRT = Test - Pre
% 2) Merges OPTIONAL covariates (age, sex, motion, impedance, E-field) from files you point to
% 3) Runs exploratory moderation models:
%       ΔAcc ~ Group + Moderator + Group×Moderator
% 4) Builds subgroup (stratified) effect sizes and forest plots for generalizability:
%       - Continuous moderators split at median (Low vs High)
%       - Sex shown as Female vs Male (if available)
% 5) Saves:
%       - Moderation table (interaction betas + p-values)
%       - Subgroup effect-size table (Hedges g + 95% CI)
%       - Forest plot PNG(s)
%
% Notes:
% - This is EXPLORATORY (report effect sizes + CIs; don’t over-interpret p-values).
% - Script is robust: if a covariate file is missing, it skips that moderator gracefully.
% ------------------------------------------------------------

clear; clc;

%% ===================== USER SETTINGS =====================
% ---- Behavioral Excel (2-back) ----
filePath = '/Users/ghazaleh/Downloads/2-BackResults.xlsx';

% ---- Output folder ----
outDir = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/TwoBack_ModerationReports';
if ~exist(outDir,'dir'); mkdir(outDir); end

% ---- Group colors (KEEP EXACTLY AS YOU REQUESTED) ----
COLORS.Up   = [0.1 0.7 0.7];  % Up-regulation
COLORS.Down = [0.8 0.2 0.2];  % Down-regulation

% ---- OPTIONAL: Covariate files (leave "" if you don't have them yet) ----
% These files should contain a participant identifier column that matches your SubjectID
% or matches "sub-<SubjectID>" (script handles both).
demogFile   = "";  % e.g., '/.../phenotype/demographics.tsv'   columns: participant_id, age, sex
motionFile  = "";  % e.g., '/.../qc/motion_summary.tsv'        columns: participant_id, mean_fd
impFile     = "";  % e.g., '/.../stimulation/impedance.tsv'    columns: participant_id, impedance_mean
efieldFile  = "";  % e.g., '/.../efields/efield_roi.tsv'       columns: participant_id, mean_efield

% Column names in those files (edit to match your actual headers if needed)
DEMOG_ID_COL = "participant_id";  DEMOG_AGE_COL = "age";   DEMOG_SEX_COL = "sex";
MOTION_ID_COL= "participant_id";  MOTION_FD_COL = "mean_fd";
IMP_ID_COL   = "participant_id";  IMP_COL       = "impedance_mean";
EF_ID_COL    = "participant_id";  EF_COL        = "mean_efield";

% Which outcome(s) to analyze
DO_ACC = true;   % ΔAccuracy
DO_RT  = true;   % ΔRT

%% ===================== SUBJECT MAPPING (ScanID -> SubjectID -> Group) =====================
% IMPORTANT: Update SubjectID mapping if needed; Group must be Control/Experimental initially.
map = table( ...
    ["E16426";"E16449";"E16488";"E16502";"E16518";"E16520";"E16531";"E16537";"E16561";"E16591"; ...
     "E16599";"E16616";"E16628";"E16682";"E16844";"E16909";"E16942";"E16979";"E17009"], ...
    ["BK144";"AZ290";"AA115";"AV531";"BK798";"BK832";"AR165";"BL003";"BA977";"AZ275"; ...
     "AV221";"AP913";"AQ975";"AU393";"AB018";"AV438";"AV503";"BK060";"BH126"], ...
    ["Experimental";"Control";"Control";"Experimental";"Control";"Experimental";"Control";"Experimental";"Experimental";"Experimental"; ...
     "Control";"Experimental";"Control";"Control";"Experimental";"Control";"Control";"Experimental";"Experimental"], ...
    'VariableNames', {'ScanID','SubjectID','Group'} ...
);

map.Group = string(map.Group);
map.Group(map.Group=="Control")      = "Down-Regulation";
map.Group(map.Group=="Experimental") = "Up-Regulation";

%% ===================== READ 2-BACK EXCEL =====================
data = readmatrix(filePath, 'TreatAsMissing','');

[nRows, nCols] = size(data);
if nRows ~= 120 || mod(nCols, 6) ~= 0
    error('Expected 120 rows and columns multiple of 6. Got %d x %d.', nRows, nCols);
end

nSub = nCols/6;
if height(map) ~= nSub
    error('Mismatch: Excel has %d subjects (nCols/6), but mapping table has %d rows.', nSub, height(map));
end

%% ===================== RESHAPE TO (120*nSub) x 6 =====================
reshaped_all = nan(nRows*nSub, 6);
for i = 1:nSub
    cols = (i-1)*6 + (1:6);
    block = data(:, cols);
    r0 = (i-1)*nRows + 1;
    r1 = i*nRows;
    reshaped_all(r0:r1, :) = block;
end

%% ===================== BUILD TRIAL-LEVEL TABLE =====================
rowsPerSubject = nRows;
scanIDs = string(map.ScanID);
ScanID_col = repelem(scanIDs, rowsPerSubject, 1);

if size(reshaped_all,1) ~= numel(ScanID_col)
    error('Row mismatch: reshaped_all=%d rows, ScanID_col=%d rows.', size(reshaped_all,1), numel(ScanID_col));
end

finalTable = table(ScanID_col, ...
    reshaped_all(:,1), reshaped_all(:,2), reshaped_all(:,3), ...
    reshaped_all(:,4), reshaped_all(:,5), reshaped_all(:,6), ...
    'VariableNames', {'ScanID','TR1','TR1_timing','TR2','TR2_timing','Test','Test_timing'});

finalTable = outerjoin(finalTable, map, 'Keys','ScanID', 'MergeKeys',true);

%% ===================== SUBJECT-LEVEL OUTCOMES =====================
phasesAcc = ["TR1","TR2","Test"];
phasesRT  = ["TR1_timing","TR2_timing","Test_timing"];

Subj = unique(finalTable(:, {'ScanID','SubjectID','Group'}), 'rows');
Subj = sortrows(Subj, 'ScanID');

for p = 1:numel(phasesAcc)
    accVar = phasesAcc(p);
    rtVar  = phasesRT(p);

    mAcc = nan(height(Subj),1);  sdAcc = nan(height(Subj),1);
    mRT  = nan(height(Subj),1);  sdRT  = nan(height(Subj),1);

    for s = 1:height(Subj)
        idx = finalTable.ScanID == Subj.ScanID(s);

        a = finalTable{idx, accVar};
        r = finalTable{idx, rtVar};

        mAcc(s) = mean(a, 'omitnan');
        sdAcc(s)= std(a,  'omitnan');

        mRT(s)  = mean(r, 'omitnan');
        sdRT(s) = std(r,  'omitnan');
    end

    Subj.(accVar + "_Mean") = mAcc;
    Subj.(accVar + "_SD")   = sdAcc;
    Subj.(rtVar  + "_Mean") = mRT;
    Subj.(rtVar  + "_SD")   = sdRT;
end

Subj.PreAcc_Mean = mean([Subj.TR1_Mean, Subj.TR2_Mean], 2, 'omitnan');
Subj.PreRT_Mean  = mean([Subj.TR1_timing_Mean, Subj.TR2_timing_Mean], 2, 'omitnan');

Subj.DeltaAcc_TestMinusPre = Subj.Test_Mean - Subj.PreAcc_Mean;
Subj.DeltaRT_TestMinusPre  = Subj.Test_timing_Mean - Subj.PreRT_Mean;

% Group coding (0=Down, 1=Up) for models
Subj.GroupBin = double(Subj.Group=="Up-Regulation");

%% ===================== MERGE OPTIONAL COVARIATES =====================
Subj = addCovariate(Subj, demogFile,  DEMOG_ID_COL,  {DEMOG_AGE_COL, DEMOG_SEX_COL});
Subj = addCovariate(Subj, motionFile, MOTION_ID_COL, {MOTION_FD_COL});
Subj = addCovariate(Subj, impFile,    IMP_ID_COL,    {IMP_COL});
Subj = addCovariate(Subj, efieldFile, EF_ID_COL,     {EF_COL});

% Standardize expected names (so the rest of the script is stable)
Subj = normalizeCovariateNames(Subj, DEMOG_AGE_COL, DEMOG_SEX_COL, MOTION_FD_COL, IMP_COL, EF_COL);

% Z-score continuous moderators if present
Subj = addZ(Subj, "Age",        "zAge");
Subj = addZ(Subj, "PreAcc_Mean","zBaseline");
Subj = addZ(Subj, "MeanFD",     "zMotion");
Subj = addZ(Subj, "Impedance",  "zImp");
Subj = addZ(Subj, "Efield",     "zEfield");

%% ===================== DEFINE MODERATORS TO TEST =====================
% Each row: {variableNameInSubj, displayLabel, type}
% type: "cont" uses median split for subgroup forest + interaction term in model
%       "bin"  uses levels directly (e.g., Sex)
mods = {};
mods = addMod(mods, "zAge",      "Age (z)",                   "cont");
mods = addMod(mods, "Sex",       "Sex",                       "bin");
mods = addMod(mods, "zBaseline", "Baseline WM (PreAcc, z)",   "cont");
mods = addMod(mods, "zEfield",   "E-field (z)",               "cont");
mods = addMod(mods, "zImp",      "Impedance (z)",             "cont");
mods = addMod(mods, "zMotion",   "Motion (Mean FD, z)",       "cont");

% Keep only those actually present in Subj
mods = mods(cellfun(@(v) ismember(v, Subj.Properties.VariableNames), mods(:,1)), :);

if isempty(mods)
    warning('No moderators found in Subj (covariate files empty or columns not present). The script will stop after saving behavioral summary.');
    writetable(Subj, fullfile(outDir,'Subj_TwoBack_BehavioralOnly.csv'));
    return;
end

%% ===================== MODERATION MODELS =====================
% Models:
%   Outcome ~ GroupBin + M + GroupBin*M
% Saves interaction beta + p-value, plus N used.
T_mod_ACC = table();
T_mod_RT  = table();

if DO_ACC
    T_mod_ACC = runModerationLoop(Subj, "DeltaAcc_TestMinusPre", mods);
    writetable(T_mod_ACC, fullfile(outDir,'Table1_Moderation_DeltaAcc.csv'));
end

if DO_RT
    T_mod_RT  = runModerationLoop(Subj, "DeltaRT_TestMinusPre",  mods);
    writetable(T_mod_RT,  fullfile(outDir,'Table2_Moderation_DeltaRT.csv'));
end

%% ===================== SUBGROUP EFFECT SIZES + FOREST PLOTS =====================
% Forest plot outcome = ΔAcc or ΔRT (effect size = Hedges g: Up vs Down)
if DO_ACC
    [T_sub_ACC, figFile] = subgroupForest(Subj, "DeltaAcc_TestMinusPre", mods, ...
        'Δ Accuracy (Test - Pre)', outDir);
    writetable(T_sub_ACC, fullfile(outDir,'Table3_Subgroups_Forest_DeltaAcc.csv'));
    fprintf('Saved forest plot: %s\n', figFile);
end

if DO_RT
    [T_sub_RT, figFile] = subgroupForest(Subj, "DeltaRT_TestMinusPre", mods, ...
        'Δ RT (Test - Pre)', outDir);
    writetable(T_sub_RT, fullfile(outDir,'Table4_Subgroups_Forest_DeltaRT.csv'));
    fprintf('Saved forest plot: %s\n', figFile);
end

fprintf('\nDONE. All outputs saved to:\n%s\n', outDir);

%% ===================== LOCAL FUNCTIONS =====================

function T = addMod(T, varName, label, type)
T = [T; {varName, label, type}]; %#ok<AGROW>
end

function Subj = addCovariate(Subj, filePath, idCol, covCols)
% Robust join by SubjectID; handles participant_id like "sub-BK144" or "BK144"
if strlength(string(filePath))==0
    return;
end
if ~exist(filePath,'file')
    warning('Covariate file not found: %s (skipping)', filePath);
    return;
end

[~,~,ext] = fileparts(filePath);
if strcmpi(ext,'.tsv')
    C = readtable(filePath, 'FileType','text', 'Delimiter','\t');
else
    C = readtable(filePath);
end

% Ensure ID column exists
if ~ismember(idCol, string(C.Properties.VariableNames))
    warning('ID column "%s" not found in %s (skipping)', idCol, filePath);
    return;
end

% Clean IDs in covariate file
C.ParticipantID_clean = cleanID(C.(idCol));

% Clean SubjectID in Subj
Subj.SubjectID_clean = cleanID(Subj.SubjectID);

% Keep only requested covariate cols that exist
keepCols = ["ParticipantID_clean"];
for k = 1:numel(covCols)
    col = string(covCols{k});
    if ismember(col, string(C.Properties.VariableNames))
        keepCols(end+1) = col; %#ok<AGROW>
    else
        warning('Covariate column "%s" not found in %s (skipping that column)', col, filePath);
    end
end

C2 = C(:, cellstr(keepCols));

% Join
Subj = outerjoin(Subj, C2, ...
    'LeftKeys', 'SubjectID_clean', ...
    'RightKeys','ParticipantID_clean', ...
    'MergeKeys', true);

% Cleanup
if ismember('ParticipantID_clean', Subj.Properties.VariableNames)
    Subj.ParticipantID_clean = []; % redundant after merge
end
end

function ids = cleanID(x)
% Convert anything to string, remove "sub-" prefix if present
ids = string(x);
ids = strip(ids);
ids = replace(ids, "sub-", "");
ids = replace(ids, "SUB-", "");
end

function Subj = normalizeCovariateNames(Subj, ageCol, sexCol, fdCol, impCol, efCol)
% Map whatever the joined columns are named into stable names:
% Age, Sex, MeanFD, Impedance, Efield
vn = string(Subj.Properties.VariableNames);

if ismember(ageCol, vn) && ~ismember("Age", vn)
    Subj.Age = Subj.(ageCol);
end
if ismember(sexCol, vn) && ~ismember("Sex", vn)
    % Convert to categorical where possible
    s = Subj.(sexCol);
    if iscellstr(s) || isstring(s)
        Subj.Sex = categorical(string(s));
    else
        % If numeric: assume 0/1 or 1/2, keep as categorical but label unknown
        Subj.Sex = categorical(s);
    end
end
if ismember(fdCol, vn) && ~ismember("MeanFD", vn)
    Subj.MeanFD = Subj.(fdCol);
end
if ismember(impCol, vn) && ~ismember("Impedance", vn)
    Subj.Impedance = Subj.(impCol);
end
if ismember(efCol, vn) && ~ismember("Efield", vn)
    Subj.Efield = Subj.(efCol);
end
end

function Subj = addZ(Subj, src, dst)
if ~ismember(src, Subj.Properties.VariableNames); return; end
x = Subj.(src);
if ~isnumeric(x); return; end
Subj.(dst) = nan(size(x));
m = ~isnan(x);
if sum(m)>=2
    Subj.(dst)(m) = (x(m) - mean(x(m))) / std(x(m));
end
end

function T = runModerationLoop(Subj, outcomeVar, mods)
% Outcome ~ GroupBin + M + GroupBin*M
rows = {};
for i = 1:size(mods,1)
    M = string(mods{i,1});
    label = string(mods{i,2});
    type = string(mods{i,3});

    % Build model table
    if strcmpi(type,"bin")
        % Ensure categorical for bin moderators
        mcol = Subj.(M);
        if ~iscategorical(mcol)
            mcol = categorical(mcol);
        end
        tmp = table(Subj.(outcomeVar), Subj.GroupBin, mcol, ...
            'VariableNames', {'Y','GroupBin','M'});
    else
        tmp = table(Subj.(outcomeVar), Subj.GroupBin, Subj.(M), ...
            'VariableNames', {'Y','GroupBin','M'});
    end

    tmp = rmmissing(tmp);
    N = height(tmp);

    if N < 8
        rows(end+1,:) = {label, N, NaN, NaN, NaN, NaN, NaN, NaN}; %#ok<AGROW>
        continue;
    end

    % Fit
    % Use formula with interaction
    lm = fitlm(tmp, 'Y ~ GroupBin*M');

    % Extract interaction term
    coefNames = string(lm.CoefficientNames);
    % Interaction name in MATLAB will look like: 'GroupBin:M' (numeric) or 'GroupBin:M_<level>' for categorical
    ix = find(contains(coefNames, "GroupBin:M"), 1, 'first');

    if isempty(ix)
        beta = NaN; se = NaN; t = NaN; p = NaN;
    else
        beta = lm.Coefficients.Estimate(ix);
        se   = lm.Coefficients.SE(ix);
        t    = lm.Coefficients.tStat(ix);
        p    = lm.Coefficients.pValue(ix);
    end

    Rsq = lm.Rsquared.Ordinary;
    rows(end+1,:) = {label, N, beta, se, t, p, Rsq, lm.ModelCriterion.AIC}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'Moderator','N','InteractionBeta','SE','tStat','pValue','R2','AIC'});
end

function [T_sub, figFile] = subgroupForest(Subj, outcomeVar, mods, xLabel, outDir)
% Builds subgroup table (Hedges g Up vs Down) and makes a forest plot.

rows = {};
for i = 1:size(mods,1)
    M = string(mods{i,1});
    label = string(mods{i,2});
    type = string(mods{i,3});

    if strcmpi(type,"bin")
        % Use each level as subgroup
        g = Subj.(M);
        if ~iscategorical(g); g = categorical(g); end
        lvls = categories(g);
        for k = 1:numel(lvls)
            lvl = string(lvls{k});
            idx = (g == lvls{k});
            [gH, lo, hi, nD, nU] = hedgesG(Subj.(outcomeVar)(idx), Subj.Group(idx));
            rows(end+1,:) = {label, lvl, nD, nU, gH, lo, hi}; %#ok<AGROW>
        end
    else
        % Median split for continuous moderators
        x = Subj.(M);
        if ~isnumeric(x)
            continue;
        end
        med = median(x, 'omitnan');
        idxLow  = x <= med;
        idxHigh = x >  med;

        [gL, loL, hiL, nD1, nU1] = hedgesG(Subj.(outcomeVar)(idxLow),  Subj.Group(idxLow));
        [gH, loH, hiH, nD2, nU2] = hedgesG(Subj.(outcomeVar)(idxHigh), Subj.Group(idxHigh));

        rows(end+1,:) = {label, "Low (≤ median)",  nD1, nU1, gL, loL, hiL}; %#ok<AGROW>
        rows(end+1,:) = {label, "High (> median)", nD2, nU2, gH, loH, hiH}; %#ok<AGROW>
    end
end

T_sub = cell2table(rows, 'VariableNames', ...
    {'Moderator','Level','N_down','N_up','Hedges_g','CI_low','CI_high'});

% Order rows nicely
T_sub.Moderator = string(T_sub.Moderator);
T_sub.Level     = string(T_sub.Level);

% ---------- Forest plot ----------
fig = figure('Color','w'); 
ax = axes(fig); hold(ax,'on');

n = height(T_sub);
y = 1:n;                 % MUST be increasing for yticks

% Plot CI lines + points
for i = 1:n
    if isnan(T_sub.Hedges_g(i))
        continue;
    end
    plot(ax, [T_sub.CI_low(i), T_sub.CI_high(i)], [y(i) y(i)], 'k-', 'LineWidth', 1.5);
    plot(ax, T_sub.Hedges_g(i), y(i), 'ko', 'MarkerFaceColor','k', 'MarkerSize', 6);
end

% Reference line at 0
plot(ax, [0 0], [0 n+1], 'k--', 'LineWidth', 1);

% Y labels
ylabs = strings(n,1);
for i = 1:n
    ylabs(i) = sprintf('%s: %s  (nD=%d, nU=%d)', ...
        T_sub.Moderator(i), T_sub.Level(i), T_sub.N_down(i), T_sub.N_up(i));
end

yticks(ax, y);
yticklabels(ax, ylabs);

% Put first row at top (visual “forest plot” convention)
set(ax, 'YDir','reverse');

xlabel(ax, 'Effect size (Hedges g): Up vs Down');
title(ax, sprintf('Subgroup generalizability: %s', xLabel));
grid(ax, 'on'); box(ax, 'on');

% x-limits
lo = min(T_sub.CI_low, [], 'omitnan');
hi = max(T_sub.CI_high, [], 'omitnan');
if ~isnan(lo) && ~isnan(hi)
    xlim(ax, [lo-0.1, hi+0.1]);
end

figFile = fullfile(outDir, sprintf('Forest_%s.png', char(strrep(outcomeVar,":","_"))));
saveas(fig, figFile);
end

function [g, ciLo, ciHi, nDown, nUp] = hedgesG(values, groupStr)
% Hedges g for difference between Up and Down on "values"
% groupStr must be string array with "Down-Regulation"/"Up-Regulation"
g = NaN; ciLo = NaN; ciHi = NaN;

if isempty(values) || isempty(groupStr); nDown=0; nUp=0; return; end
values = values(:);
groupStr = string(groupStr(:));

xD = values(groupStr=="Down-Regulation"); xD = xD(~isnan(xD));
xU = values(groupStr=="Up-Regulation");   xU = xU(~isnan(xU));

nDown = numel(xD); nUp = numel(xU);
if nDown < 2 || nUp < 2
    return;
end

mD = mean(xD); sD = std(xD);
mU = mean(xU); sU = std(xU);

% Pooled SD
sp = sqrt(((nDown-1)*sD^2 + (nUp-1)*sU^2) / (nDown+nUp-2));
if sp==0
    return;
end

d  = (mU - mD) / sp;

% Hedges correction
J = 1 - (3 / (4*(nDown+nUp) - 9));
g  = J * d;

% Approx SE and 95% CI
% Var(g) approx:
var_g = (nDown+nUp)/(nDown*nUp) + (g^2)/(2*(nDown+nUp-2));
se_g  = sqrt(var_g);

ciLo = g - 1.96*se_g;
ciHi = g + 1.96*se_g;
end
