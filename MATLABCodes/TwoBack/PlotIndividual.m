%% TwoBack_Analysis_WithStatsAndPlots.m
% Fixes:
% - makeGroupSummary row-size bug (1x6 vs 1x7)
% - keeps everything else the same

clear; clc; close all;
COLORS.Up   = [0.1 0.7 0.7];   % Up-regulation (teal)
COLORS.Down = [0.8 0.2 0.2];   % Down-regulation (red)


%% ========== INPUTS ==========
filePath = '/Users/ghazaleh/Downloads/2-BackResults.xlsx';
CLOSE_FIGS = false;   % true = batch mode, false = interactive


outDir = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/TwoBack_Reports';
if ~exist(outDir, 'dir'); mkdir(outDir); end

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

%% ========== READ EXCEL ==========
data = readmatrix(filePath, 'TreatAsMissing','');

[nRows, nCols] = size(data);
if nRows ~= 120 || mod(nCols, 6) ~= 0
    error('Expected 120 rows and columns multiple of 6. Got %d x %d.', nRows, nCols);
end

nSub = nCols/6;

if height(map) ~= nSub
    error('Mismatch: Excel has %d subjects (nCols/6), but mapping table has %d rows.', nSub, height(map));
end

%% ========== RESHAPE TO 2280×6 ==========
reshaped_all = nan(nRows*nSub, 6);

for i = 1:nSub
    cols = (i-1)*6 + (1:6);
    block = data(:, cols);
    r0 = (i-1)*nRows + 1;
    r1 = i*nRows;
    reshaped_all(r0:r1, :) = block;
end

%% ========== BUILD TRIAL-LEVEL TABLE ==========
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

%% ========== SUBJECT-LEVEL SUMMARY ==========
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

%% ========== GROUP-LEVEL SUMMARY TABLES (WITH SD + SEM) ==========
phaseRows = table();
for p = 1:numel(phasesAcc)
    accM = Subj.(phasesAcc(p)+"_Mean");
    rtM  = Subj.(phasesRT(p)+"_Mean");

    phaseRows = [phaseRows;
        makeGroupSummary(Subj, accM, "Accuracy",      phasesAcc(p));
        makeGroupSummary(Subj, rtM,  "ReactionTime",  phasesAcc(p))];
end
T_phase = phaseRows;

T_delta = [ ...
    makeGroupSummary(Subj, Subj.DeltaAcc_TestMinusPre, "DeltaAccuracy", "Test-Pre"); ...
    makeGroupSummary(Subj, Subj.DeltaRT_TestMinusPre,  "DeltaRT",       "Test-Pre") ...
];

writetable(T_phase, fullfile(outDir, 'Table1_TwoBack_GroupSummary_ByPhase.csv'));
writetable(T_delta, fullfile(outDir, 'Table2_TwoBack_GroupSummary_Deltas.csv'));

%% ========== STATS ==========
T_stats_between = betweenGroupStats(Subj);
T_stats_within  = withinGroupPrePostStats(Subj);

writetable(T_stats_between, fullfile(outDir, 'Table3_TwoBack_BetweenGroup_Stats.csv'));
writetable(T_stats_within,  fullfile(outDir, 'Table4_TwoBack_WithinGroup_PreVsPost.csv'));

%% ========== PLOTS ==========
groups = ["Down-Regulation","Up-Regulation"];
phaseLabels = ["TR1","TR2","Test"];

plotGroupMeans(Subj, ["TR1_Mean","TR2_Mean","Test_Mean"], phaseLabels, groups, ...
    'Accuracy', [0 1], fullfile(outDir, 'Fig1_Accuracy_GroupMeans_SD.png'), COLORS);

plotGroupMeans(Subj, ["TR1_timing_Mean","TR2_timing_Mean","Test_timing_Mean"], phaseLabels, groups, ...
    'Reaction Time (s)', [], fullfile(outDir, 'Fig2_RT_GroupMeans_SD.png'), COLORS);

plotSpaghetti(Subj, ["TR1_Mean","TR2_Mean","Test_Mean"], phaseLabels, groups, ...
    'Accuracy', [0 1], fullfile(outDir, 'Fig3_Accuracy_Spaghetti.png'), COLORS);

plotSpaghetti(Subj, ["TR1_timing_Mean","TR2_timing_Mean","Test_timing_Mean"], phaseLabels, groups, ...
    'Reaction Time (s)', [], fullfile(outDir, 'Fig4_RT_Spaghetti.png'), COLORS);

plotDelta(Subj, "DeltaAcc_TestMinusPre", groups, 'Δ Accuracy (Test - Pre)', [], ...
    fullfile(outDir, 'Fig5_DeltaAcc_TestMinusPre.png'), COLORS);

plotDelta(Subj, "DeltaRT_TestMinusPre", groups, 'Δ RT (Test - Pre)', [], ...
    fullfile(outDir, 'Fig6_DeltaRT_TestMinusPre.png'), COLORS);

fprintf('\nSaved outputs to:\n%s\n', outDir);

%% ===================== LOCAL FUNCTIONS =====================

function T = makeGroupSummary(Subj, values, measureName, phaseName)
% Returns 2 rows (one per group) with N, Mean, SD, SEM
groups = ["Down-Regulation","Up-Regulation"];

Phase   = strings(0,1);
Group   = strings(0,1);
Measure = strings(0,1);
Ncol    = zeros(0,1);
Meancol = nan(0,1);
SDcol   = nan(0,1);
SEMcol  = nan(0,1);

for g = 1:numel(groups)
    grp = groups(g);
    idx = Subj.Group == grp;

    x = values(idx);
    x = x(~isnan(x));

    N = numel(x);
    if N==0
        mu = NaN; sd = NaN; sem = NaN;
    else
        mu = mean(x);
        sd = std(x);
        sem = sd/sqrt(N);
    end

    Phase(end+1,1)   = string(phaseName);
    Group(end+1,1)   = string(grp);
    Measure(end+1,1) = string(measureName);
    Ncol(end+1,1)    = N;
    Meancol(end+1,1) = mu;
    SDcol(end+1,1)   = sd;
    SEMcol(end+1,1)  = sem;
end

T = table(Phase, Group, Measure, Ncol, Meancol, SDcol, SEMcol, ...
    'VariableNames', {'Phase','Group','Measure','N','Mean','SD','SEM'});
end

function T = betweenGroupStats(Subj)
rows = {};

pairs = { ...
    'Test_Mean',            'Accuracy (Test)'; ...
    'Test_timing_Mean',     'RT (Test)'; ...
    'DeltaAcc_TestMinusPre','Δ Accuracy (Test-Pre)'; ...
    'DeltaRT_TestMinusPre', 'Δ RT (Test-Pre)' ...
};

for i = 1:size(pairs,1)
    var = pairs{i,1};
    label = pairs{i,2};

    xD = Subj.(var)(Subj.Group=="Down-Regulation"); xD = xD(~isnan(xD));
    xU = Subj.(var)(Subj.Group=="Up-Regulation");   xU = xU(~isnan(xU));

    ND = numel(xD); NU = numel(xU);
    mD = mean(xD, 'omitnan'); sD = std(xD, 'omitnan');
    mU = mean(xU, 'omitnan'); sU = std(xU, 'omitnan');

    if ND>=2 && NU>=2 && ~(sD==0 && sU==0)
        [~,p,~,st] = ttest2(xD, xU, 'Vartype','unequal');
        t = st.tstat; df = st.df;
    else
        p = NaN; t = NaN; df = NaN;
    end

    rows(end+1,:) = {label, ND, mD, sD, NU, mU, sU, t, df, p}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'Measure','N_down','Mean_down','SD_down','N_up','Mean_up','SD_up','t','df','p_value'});
end

function T = withinGroupPrePostStats(Subj)
rows = {};
groups = ["Down-Regulation","Up-Regulation"];

for g = 1:numel(groups)
    grp = groups(g);
    idx = Subj.Group==grp;

    preA  = Subj.PreAcc_Mean(idx);
    postA = Subj.Test_Mean(idx);

    preR  = Subj.PreRT_Mean(idx);
    postR = Subj.Test_timing_Mean(idx);

    [nA, mPreA, mPostA, tA, dfA, pA] = pairedStats(preA, postA);
    rows(end+1,:) = {char(grp), 'Accuracy', nA, mPreA, mPostA, tA, dfA, pA}; %#ok<AGROW>

    [nR, mPreR, mPostR, tR, dfR, pR] = pairedStats(preR, postR);
    rows(end+1,:) = {char(grp), 'ReactionTime', nR, mPreR, mPostR, tR, dfR, pR}; %#ok<AGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'Group','Measure','N_pairs','Mean_Pre','Mean_Post','t','df','p_value'});
end

function [nPairs, mPre, mPost, t, df, p] = pairedStats(pre, post)
mask = ~isnan(pre) & ~isnan(post);
pre  = pre(mask);
post = post(mask);

nPairs = numel(pre);
mPre  = mean(pre, 'omitnan');
mPost = mean(post,'omitnan');

if nPairs>=2 && ~(std(pre)==0 && std(post)==0)
    [~,p,~,st] = ttest(post, pre);
    t  = st.tstat;
    df = st.df;
else
    p = NaN; t = NaN; df = NaN;
end
end


function plotSpaghetti(Subj, varNames, xLabels, groups, ylab, ylims, outPng, COLORS)

figure('Color','w'); hold on;
x = 1:numel(varNames);

for g = 1:numel(groups)
    grp = groups(g);
    idx = find(Subj.Group==grp);

    thisColor = getGroupColor(grp, COLORS);

    for k = 1:numel(idx)
        s = idx(k);
        y = nan(1,numel(varNames));
        for i = 1:numel(varNames)
            y(i) = Subj.(varNames(i))(s);
        end
        plot(x, y, '-o', ...
            'Color', thisColor, ...
            'MarkerFaceColor', thisColor, ...
            'LineWidth', 1);
    end
end

xticks(x); xticklabels(xLabels);
xlabel('Phase'); ylabel(ylab);
title([ylab ' (Subject trajectories)']);
if ~isempty(ylims); ylim(ylims); end
grid on;

% saveas(gcf, outPng); close(gcf);
end


function plotGroupMeans(Subj, varNames, xLabels, groups, ylab, ylims, outPng, COLORS)
% Bar (mean) + errorbar (SEM) + subject scatter overlay with group-specific markers
% Down-Regulation: red circles
% Up-Regulation: teal triangles

figure('Color','w'); hold on;
x = 1:numel(varNames);
barW = 0.36;           % bar width
jitW = 0.10;           % jitter width for scatter
alphaBar = 0.9;

for i = 1:numel(varNames)
    % Collect group values
    for g = 1:numel(groups)
        grp = groups(g);
        idx = Subj.Group == grp;

        v = Subj.(varNames(i))(idx);
        v = v(~isnan(v));

        mu(g)  = mean(v,'omitnan'); %#ok<AGROW>
        sd(g)  = std(v,'omitnan');  %#ok<AGROW>
        sem(g) = sd(g) / sqrt(max(1,numel(v))); %#ok<AGROW>

        % x position for this group at this phase
        if contains(lower(grp),"down")
            xg = x(i) - barW/2;
        else
            xg = x(i) + barW/2;
        end

        % ---- bar ----
        c = getGroupColor(grp, COLORS);
        bh = bar(xg, mu(g), barW, 'FaceColor', c, 'EdgeColor', 'none');
        bh.FaceAlpha = alphaBar;

        % ---- error bar (SEM) ----
        errorbar(xg, mu(g), sem(g), 'k', 'LineWidth', 1.2, 'CapSize', 8);

        % ---- scatter (subjects) ----
        [mk, mfc, mec] = getGroupMarker(grp, COLORS);
        if ~isempty(v)
            jitter = (rand(size(v)) - 0.5) * jitW;
            scatter(xg + jitter, v, 42, ...
                'Marker', mk, ...
                'MarkerFaceColor', mfc, ...
                'MarkerEdgeColor', mec, ...
                'LineWidth', 0.8);
        end
    end
end

xticks(x); xticklabels(xLabels);
xlabel('Phase'); ylabel(ylab);
title([ylab ' (Mean ± SEM + subjects)']);
grid on;
if ~isempty(ylims); ylim(ylims); end

% Legend (matching example)
h1 = plot(nan,nan,'o','MarkerFaceColor',COLORS.Down,'MarkerEdgeColor','k','LineWidth',0.8);
h2 = plot(nan,nan,'^','MarkerFaceColor',COLORS.Up,'MarkerEdgeColor','k','LineWidth',0.8);
legend([h1 h2], {'Down-Regulation','Up-Regulation'}, 'Location','best');

% saveas(gcf, outPng); close(gcf);
end


function plotDelta(Subj, varName, groups, ylab, ylims, outPng, COLORS)
% Bar (mean) + errorbar (SEM) + subject scatter overlay with group-specific markers

figure('Color','w'); hold on;

xpos = 1:numel(groups);
barW = 0.55;
jitW = 0.15;

mu  = nan(1,numel(groups));
sd  = nan(1,numel(groups));
sem = nan(1,numel(groups));

for g = 1:numel(groups)
    grp = groups(g);
    idx = Subj.Group == grp;

    v = Subj.(varName)(idx);
    v = v(~isnan(v));

    mu(g)  = mean(v,'omitnan');
    sd(g)  = std(v,'omitnan');
    sem(g) = sd(g) / sqrt(max(1,numel(v)));

    c = getGroupColor(grp, COLORS);

    % bar
    bh = bar(xpos(g), mu(g), barW, 'FaceColor', c, 'EdgeColor','none');
    bh.FaceAlpha = 0.9;

    % errorbar (SEM)
    errorbar(xpos(g), mu(g), sem(g), 'k', 'LineWidth', 1.2, 'CapSize', 10);

    % scatter
    [mk, mfc, mec] = getGroupMarker(grp, COLORS);
    if ~isempty(v)
        jitter = (rand(size(v)) - 0.5) * jitW;
        scatter(xpos(g) + jitter, v, 55, ...
            'Marker', mk, ...
            'MarkerFaceColor', mfc, ...
            'MarkerEdgeColor', mec, ...
            'LineWidth', 0.8);
    end
end

xticks(xpos); xticklabels(cellstr(groups));
ylabel(ylab);
title([ylab ' (Mean ± SEM + subjects)']);
grid on;
if ~isempty(ylims); ylim(ylims); end

% saveas(gcf, outPng); close(gcf);
end


function c = getGroupColor(grp, COLORS)
grp = string(grp);
if contains(lower(grp), "up")
    c = COLORS.Up;
elseif contains(lower(grp), "down")
    c = COLORS.Down;
else
    c = [0 0 0];
end
end

function [mk, mfc, mec] = getGroupMarker(grp, COLORS)
% Down: red circles, Up: teal triangles
grp = string(grp);
if contains(lower(grp), "down")
    mk  = 'o';
    mfc = COLORS.Down;
    mec = 'k';
elseif contains(lower(grp), "up")
    mk  = '^';
    mfc = COLORS.Up;
    mec = 'k';
else
    mk  = 'o';
    mfc = [0.5 0.5 0.5];
    mec = 'k';
end
end
