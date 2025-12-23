%% ============================================================
%  BAECL - POST ONLY (after_fps_arm_1)
%  Proper Likert reporting: n (%) per response level (by group)
%  + Between-group tests appropriate for categorical/ordinal data
%
%  INCLUDED measures: ONLY core items baecl_<number> (NO a/b, NO *_other)
%  Groups:
%    Control       -> down-regulation
%    Experimental  -> up-regulation
%
%  OUTPUT CSVs (saved in same folder as the TSV file):
%   1) Table1_BAECL_POSTONLY_CountsPerc_ByGroup.csv
%   2) Table2_BAECL_POSTONLY_BetweenGroup_Tests.csv
%   3) Table3_BAECL_POSTONLY_Missingness.csv
% ============================================================

clear; clc;

%% -------------------- PATH --------------------
baeclFile = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/tacs_blindness_and_aversive_effects_check_list.tsv';
outDir    = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/MATLABCode/Blindness';

%% -------------------- SUBJECT LIST + GROUP --------------------
SubjectID = { ...
    'BK144'; 'AZ290'; 'AA115'; 'AV531'; 'BK798'; 'BK832'; ...
    'AR165'; 'BL003'; 'BA977'; 'AZ275'; 'AV221'; 'AP913'; ...
    'AQ975'; 'AU393'; 'AB018'; 'AV438'; 'AV503'; 'BK060'; 'BH126'};

Group = { ...
    'Experimental'; 'Control'; 'Control'; 'Experimental'; 'Control'; ...
    'Experimental'; 'Control'; 'Experimental'; 'Experimental'; 'Experimental'; ...
    'Control'; 'Experimental'; 'Control'; 'Control'; 'Experimental'; ...
    'Control'; 'Control'; 'Experimental'; 'Experimental'};

participants = table(string(SubjectID(:)), string(Group(:)), ...
    'VariableNames', {'SubjectID','Group'});

%% -------------------- READ TSV --------------------
if ~isfile(baeclFile)
    error('File not found:\n%s', baeclFile);
end

opts = detectImportOptions(baeclFile, 'FileType','text', 'Delimiter','\t');
opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','NaN','nan','N/A','n/a','',' '});
T = readtable(baeclFile, opts);

mustHave = ["participant_id","session"];
for k = 1:numel(mustHave)
    if ~ismember(mustHave(k), string(T.Properties.VariableNames))
        error('Missing required column: %s', mustHave(k));
    end
end

%% -------------------- NORMALIZE IDs + KEEP POST --------------------
T.participant_id = string(T.participant_id);
T.SubjectID      = strtrim(erase(T.participant_id, "sub-"));
T.session        = strtrim(string(T.session));

postLabel = "after_fps_arm_1";
T = T(T.session == postLabel, :);
if height(T)==0
    error('No rows found for session == %s', postLabel);
end

T = T(ismember(T.SubjectID, participants.SubjectID), :);
if height(T)==0
    error('No matching participant_id values found between TSV and SubjectID list (after removing sub-).');
end

%% -------------------- MERGE GROUP + RENAME --------------------
X = innerjoin(T, participants, 'Keys','SubjectID');

X.Group = string(X.Group);
X.Group(X.Group=="Control")      = "down-regulation";
X.Group(X.Group=="Experimental") = "up-regulation";
X.Group = categorical(X.Group, ["down-regulation","up-regulation"]);

%% -------------------- MEASURES LIST (ONLY MAIN BAECL ITEMS) --------------------
vars = string(X.Properties.VariableNames);

% keep only baecl_<number> (no a, no b, no _other)
isMainBAECL = ~cellfun(@isempty, regexp(vars, '^baecl_\d+$', 'once'));
measures = vars(isMainBAECL);

% sort numerically
nums = cellfun(@(x) str2double(regexp(x,'\d+','match','once')), cellstr(measures));
[~, ord] = sort(nums);
measures = measures(ord);

%% -------------------- CONVERT TO NUMERIC ROBUSTLY --------------------
for i = 1:numel(measures)
    vn = measures(i);
    x  = X.(vn);

    if isnumeric(x)
        X.(vn) = double(x);
    else
        xs = string(x);
        tok = regexp(xs, '(-?\d+(\.\d+)?)', 'tokens', 'once');
        num = nan(size(xs));
        hasTok = ~cellfun(@isempty, tok);
        num(hasTok) = str2double(string(cellfun(@(c)c{1}, tok(hasTok), 'UniformOutput', false)));
        X.(vn) = num;
    end
end

% drop measures entirely missing in POST
keep = false(size(measures));
for i = 1:numel(measures)
    m = char(measures(i));
    keep(i) = any(~isnan(X.(m)));
end
measures = measures(keep);

%% -------------------- QUESTION LABELS (CORE ITEMS ONLY) --------------------
L = containers.Map();
L('baecl_1')  = 'Blinding guess: Sham (1) vs Real (2)';
L('baecl_2')  = 'Confidence in guess (0–10)';
L('baecl_3')  = 'Perceived N-back improvement (0–4)';
L('baecl_4')  = 'Headache?';
L('baecl_5')  = 'Neck pain?';
L('baecl_6')  = 'Scalp pain?';
L('baecl_7')  = 'Tingling?';
L('baecl_8')  = 'Itching?';
L('baecl_9')  = 'Burning sensation?';
L('baecl_10') = 'Skin redness?';
L('baecl_11') = 'Sleepiness?';
L('baecl_12') = 'Trouble concentrating?';
L('baecl_13') = 'Acute mood change?';
L('baecl_14') = 'Other adverse effects?';

getLabel = @(m) localGetLabel(L, m);

%% -------------------- SCALE DEFINITIONS (levels to report) --------------------
% You can edit these if your coding differs.
scale = struct();

scale.baecl_1  = 1:2;    % Sham/Real
scale.baecl_2  = 0:10;   % 0-10
scale.baecl_3  = 0:4;    % 0-4

% side-effects presence: typically 0/1 in your data
for k = 4:14
    scale.(sprintf('baecl_%d',k)) = 0:1;
end

%% ============================================================
% TABLE 1: n (%) per level, by group (Likert-style reporting)
% ============================================================
countTbl = table;
groups = categories(X.Group);

for mi = 1:numel(measures)
    m = char(measures(mi));
    q = string(getLabel(m));
    q = strrep(q, "‚Äì", "–"); % encoding cleanup

    if isfield(scale, m)
        levels = scale.(m);
    else
        % fallback: observed unique values
        vv = X.(m);
        levels = unique(vv(~isnan(vv)))';
    end

    for gi = 1:numel(groups)
        g = groups{gi};
        v = double(X{X.Group==g, m});
        v = v(~isnan(v));
        N = numel(v);

        for li = 1:numel(levels)
            lv = levels(li);
            n  = sum(v == lv);
            pct = (N>0) * (100*n/N);

            countTbl = [countTbl; table( ...
                "post", string(g), string(m), q, ...
                lv, n, pct, N, ...
                'VariableNames',{'Session','Group','Measure','Question','Level','n','percent','N_nonmissing'})]; %#ok<AGROW>
        end
    end
end

%% ============================================================
% TABLE 2: Between-group tests (categorical/ordinal-appropriate)
%   - Fisher exact for:
%       * binary presence items (0/1) and baecl_1 (2 levels)
%   - Rank-sum (Mann–Whitney U) for:
%       * baecl_2 (0–10), baecl_3 (0–4)
%   - If a table is degenerate (no variation / empty), p=NaN or p=1
% ============================================================
testTbl = table;

for mi = 1:numel(measures)
    m = char(measures(mi));
    q = string(getLabel(m));
    q = strrep(q, "‚Äì", "–");

    vDown = double(X{X.Group=="down-regulation", m});
    vUp   = double(X{X.Group=="up-regulation",   m});
    vDown = vDown(~isnan(vDown));
    vUp   = vUp(~isnan(vUp));

    Ndown = numel(vDown);
    Nup   = numel(vUp);

    meanDown = mean(vDown,'omitnan'); sdDown = std(vDown,'omitnan');
    meanUp   = mean(vUp,'omitnan');   sdUp   = std(vUp,'omitnan');

    testName = "NA"; stat = NaN; df = NaN; p = NaN;

    if Ndown==0 || Nup==0
        testName = "NA";
        p = NaN;
    else
        allVals = [vDown(:); vUp(:)];
        uAll = unique(allVals(~isnan(allVals)));

        if numel(uAll)==1
            testName = "NoVariation";
            p = 1;
        else
            isBinary = isequal(sort(uAll(:))', [0 1]) || (numel(uAll)==2 && all(ismember(uAll,[0 1])));

            if strcmp(m,'baecl_1') || isBinary
                % Fisher exact on 2x2 table.
                u = unique(allVals(~isnan(allVals)));
                if numel(u)~=2
                    testName = "FisherExact";
                    p = NaN;
                else
                    a = u(1); b = u(2);
                    tbl2x2 = [sum(vDown==a) sum(vDown==b); sum(vUp==a) sum(vUp==b)];
                    try
                        [~, p] = fishertest(tbl2x2);
                        testName = "FisherExact";
                    catch
                        testName = "FisherExact";
                        p = NaN; % if stats toolbox missing
                    end
                end

            else
                % Ordinal (Likert-like): rank-sum preferred
                try
                    p = ranksum(vDown, vUp);
                    testName = "RankSum";
                catch
                    testName = "RankSum";
                    p = NaN;
                end
            end
        end
    end

    testTbl = [testTbl; table( ...
        "post", string(m), q, ...
        Ndown, meanDown, sdDown, ...
        Nup,   meanUp,   sdUp, ...
        testName, stat, df, p, ...
        'VariableNames',{'Session','Measure','Question', ...
        'N_down','Mean_down','SD_down','N_up','Mean_up','SD_up', ...
        'Test','Stat','df','p_value'})]; %#ok<AGROW>
end

%% ============================================================
% TABLE 3: Missingness by measure (post only)
% ============================================================
missTbl = table;
for mi = 1:numel(measures)
    m = char(measures(mi));
    v = X.(m);
    N_total   = numel(v);
    N_missing = sum(isnan(v));
    missTbl = [missTbl; table(string(m), N_missing, N_total, ...
        'VariableNames',{'Measure','N_missing','N_total'})]; %#ok<AGROW>
end

%% -------------------- SAVE CSVs --------------------
f1 = fullfile(outDir, 'Table1_BAECL_POSTONLY_CountsPerc_ByGroup.csv');
f2 = fullfile(outDir, 'Table2_BAECL_POSTONLY_BetweenGroup_Tests.csv');
f3 = fullfile(outDir, 'Table3_BAECL_POSTONLY_Missingness.csv');

writetable(countTbl, f1);
writetable(testTbl,  f2);
writetable(missTbl,  f3);

fprintf('\n===== Saved reported tables (in %s) =====\n', outDir);
fprintf('1) %s\n', f1);
fprintf('2) %s\n', f2);
fprintf('3) %s\n', f3);

%% -------------------- DISPLAY --------------------
fprintf('\n===== BAECL POST-ONLY COUNTS/PERCENTS (by group, per level) =====\n');
disp(countTbl);

fprintf('\n===== BAECL POST-ONLY BETWEEN-GROUP (categorical/ordinal tests) =====\n');
disp(testTbl);

fprintf('\n===== Missingness by measure (post only) =====\n');
disp(missTbl);

%% =================== Local helper ===================
function q = localGetLabel(L, m)
    if isKey(L, m)
        q = L(m);
    else
        q = m;
    end
end
