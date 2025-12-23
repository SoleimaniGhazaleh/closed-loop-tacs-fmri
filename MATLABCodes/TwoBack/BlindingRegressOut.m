%% ============================================================
%  Blinding & Side-effects analysis (BAECL) + regression on ΔAccuracy
%  Requires:
%   - finalTable from your 2-back pipeline (trial-level/row-level)
%   - BAECL TSV file (post-only has data)
% ============================================================

%% ---------- USER PATHS ----------
baeclFile = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/tacs_blindness_and_aversive_effects_check_list.tsv';

postSessionLabel = "after_fps_arm_1";  % only this has values per your note

%% ---------- 1) Subject-level 2-back summary ----------
% We want one row per ScanID with TR1/TR2/Test and delta metrics
subjTbl = groupsummary(finalTable, "ScanID", "mean", ...
    ["TR1","TR2","Test","TR1_timing","TR2_timing","Test_timing"]);

% Rename groupsummary columns to simpler names
subjTbl.Properties.VariableNames = strrep(subjTbl.Properties.VariableNames, "mean_", "");
subjTbl = renamevars(subjTbl, "GroupCount", "Ntrials");

% Delta accuracy: Test - mean(TR1,TR2)
subjTbl.Acc_TRmean = mean([subjTbl.TR1 subjTbl.TR2], 2, "omitnan");
subjTbl.dAcc_TestMinusTR = subjTbl.Test - subjTbl.Acc_TRmean;

% Delta RT: Test_timing - mean(TR1_timing,TR2_timing)
subjTbl.RT_TRmean = mean([subjTbl.TR1_timing subjTbl.TR2_timing], 2, "omitnan");
subjTbl.dRT_TestMinusTR = subjTbl.Test_timing - subjTbl.RT_TRmean;

% Add Group (down/up) using your scan lists
subjTbl.Group = repmat("NA", height(subjTbl), 1);
subjTbl.Group(ismember(subjTbl.ScanID, ctrlScans)) = "down-regulation";
subjTbl.Group(ismember(subjTbl.ScanID, expScans))  = "up-regulation";
subjTbl.Group = categorical(subjTbl.Group, ["down-regulation","up-regulation"]);

%% ---------- 2) Read BAECL TSV (post only) ----------
opts = detectImportOptions(baeclFile,'FileType','text','Delimiter','\t');
opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','NaN','nan','N/A','n/a','',' '});
B = readtable(baeclFile, opts);

% Normalize IDs: participant_id like "sub-BK144" -> "BK144"
B.participant_id = string(B.participant_id);
B.SubjectID = strtrim(erase(B.participant_id, "sub-"));

% Keep post session
B.session = strtrim(string(B.session));
B = B(B.session == postSessionLabel, :);

% Keep ONLY your study subjects (from ScanID->SubjectID mapping table you gave earlier)
% Build ScanID -> SubjectID mapping from your participant list:
scan2sub = table( ...
    string({'E16426','E16449','E16488','E16502','E16518','E16520','E16531','E16537','E16561','E16591', ...
            'E16599','E16616','E16628','E16682','E16844','E16909','E16942','E16979','E17009'})', ...
    string({'BK144','AZ290','AA115','AV531','BK798','BK832','AR165','BL003','BA977','AZ275', ...
            'AV221','AP913','AQ975','AU393','AB018','AV438','AV503','BK060','BH126'})', ...
    'VariableNames', {'ScanID','SubjectID'});

% Add SubjectID to subjTbl via ScanID
subjTbl = outerjoin(subjTbl, scan2sub, 'Keys','ScanID', 'MergeKeys', true);

% Join BAECL to subject table
M = innerjoin(subjTbl, B, 'Keys','SubjectID');

if height(M)==0
    error('No matching rows after joining 2-back summary with BAECL. Check SubjectID mapping and session label.');
end

%% ---------- 3) Extract blinding variables ----------
% baecl_1: 1=Sham, 2=Real
% baecl_2: confidence 0-10

% Convert to numeric robustly
M.baecl_1 = local_toNumeric(M.baecl_1);
M.baecl_2 = local_toNumeric(M.baecl_2);

% "GuessIsReal" binary
M.GuessIsReal = double(M.baecl_1 == 2);  % 1 if guessed Real, 0 if guessed Sham

% Confidence (0-10), keep as numeric
M.Confidence = M.baecl_2;

%% ---------- 4) Define "true" assignment + correctness ----------
% IMPORTANT: edit if your mapping differs
% Here we assume: up-regulation = Real, down-regulation = Sham
M.TrueIsReal = double(M.Group == "up-regulation");

M.CorrectGuess = double(M.GuessIsReal == M.TrueIsReal);

%% ---------- 5) Blinding index (simple + confidence-weighted) ----------
% Unweighted: +1 correct, 0 missing, -1 incorrect
M.BlindIndex_simple = nan(height(M),1);
ok = ~isnan(M.GuessIsReal) & ~isnan(M.TrueIsReal);
M.BlindIndex_simple(ok) = 2*(M.CorrectGuess(ok)) - 1;   % correct=+1, incorrect=-1

% Confidence-weighted: same sign, scaled by confidence/10 (range -1..+1)
M.BlindIndex_weighted = nan(height(M),1);
ok2 = ok & ~isnan(M.Confidence);
M.BlindIndex_weighted(ok2) = (2*(M.CorrectGuess(ok2)) - 1) .* (M.Confidence(ok2)/10);

%% ---------- 6) Side-effect burden score (POST) ----------
% Use presence items only (baecl_4..baecl_14) as you requested (no a/b)
sideVars = "baecl_" + string(4:14);

for v = sideVars
    if ismember(v, string(M.Properties.VariableNames))
        M.(v) = local_toNumeric(M.(v)); % 0/1
    else
        M.(v) = nan(height(M),1);
    end
end

% Burden: sum of yes(=1) across available items (omit missing)
M.SideEffectBurden = sum(M{:, sideVars}, 2, "omitnan");

% Optional: number of non-missing side-effect items
M.SideEffectNanswered = sum(~isnan(M{:, sideVars}), 2);

%% ---------- 7) Primary regression: ΔAccuracy ~ Guess + SideEffects (+ covariates optional) ----------
% Model outcome: Test-run Δaccuracy (Test - mean(TR1,TR2))
Y = M.dAcc_TestMinusTR;

% Predictors:
%   GuessIsReal (0/1)
%   Confidence (optional)
%   SideEffectBurden
%   Group (optional; often include to separate true arm from guess)
%
% Since your question: "regress Δaccuracy on guessed assignment and side-effect scores"
% We'll run:
%   Model A: ΔAcc ~ GuessIsReal + SideEffectBurden
%   Model B: ΔAcc ~ GuessIsReal + SideEffectBurden + Group
%   Model C: ΔAcc ~ GuessIsReal + Confidence + SideEffectBurden + Group

Treg = M(:, {'ScanID','SubjectID','Group','dAcc_TestMinusTR','GuessIsReal','Confidence','SideEffectBurden','CorrectGuess'});
Treg = rmmissing(Treg, 'DataVariables', {'dAcc_TestMinusTR','GuessIsReal','SideEffectBurden'}); % keep rows needed

mdlA = fitlm(Treg, 'dAcc_TestMinusTR ~ GuessIsReal + SideEffectBurden');
mdlB = fitlm(Treg, 'dAcc_TestMinusTR ~ GuessIsReal + SideEffectBurden + Group');
% Confidence often missing for some; only fit when present
TregC = rmmissing(Treg, 'DataVariables', {'Confidence'});
mdlC = fitlm(TregC, 'dAcc_TestMinusTR ~ GuessIsReal + Confidence + SideEffectBurden + Group');

fprintf('\n===== Regression A: ΔAcc ~ Guess + SideEffects =====\n');
disp(mdlA);
fprintf('\n===== Regression B: ΔAcc ~ Guess + SideEffects + Group =====\n');
disp(mdlB);
fprintf('\n===== Regression C: ΔAcc ~ Guess + Confidence + SideEffects + Group =====\n');
disp(mdlC);

%% ---------- 8) Repeat primary models excluding CORRECT guessers ----------
T_wrong = Treg(Treg.CorrectGuess==0, :);
if height(T_wrong) >= 5
    mdlA_wrong = fitlm(T_wrong, 'dAcc_TestMinusTR ~ GuessIsReal + SideEffectBurden');
    mdlB_wrong = fitlm(T_wrong, 'dAcc_TestMinusTR ~ GuessIsReal + SideEffectBurden + Group');

    fprintf('\n===== EXCLUDING CORRECT-GUESS participants =====\n');
    fprintf('N kept = %d (wrong guess only)\n', height(T_wrong));
    fprintf('\n--- Regression A (wrong only) ---\n'); disp(mdlA_wrong);
    fprintf('\n--- Regression B (wrong only) ---\n'); disp(mdlB_wrong);
else
    fprintf('\nNot enough participants after excluding correct guessers to refit models.\n');
end

%% ---------- 9) Report blinding summary by group (counts + %) ----------
%% ---------- 9) Report blinding summary by group (counts + %) ----------
% Use groupsummary with "sum" and grab N from GroupCount (works across versions)

tmp = groupsummary(M, "Group", "sum", "CorrectGuess");   % gives sum_CorrectGuess + GroupCount
tmp = renamevars(tmp, "sum_CorrectGuess", "N_correct");
tmp.N = tmp.GroupCount;
tmp.CorrectPct = 100 * tmp.N_correct ./ tmp.N;

blSummary = tmp(:, {'Group','N','N_correct','CorrectPct'});

fprintf('\n===== Blinding accuracy by group (post) =====\n');
disp(blSummary);


%% ---------- 10) Save key tables (optional) ----------
outDir = fileparts(baeclFile);
writetable(Treg,     fullfile(outDir, 'Table_2Back_Blinding_SideEffects_RegressionData.csv'));
writetable(blSummary,fullfile(outDir, 'Table_BlindingAccuracy_ByGroup.csv'));

fprintf('\nSaved:\n  %s\n  %s\n', ...
    fullfile(outDir, 'Table_2Back_Blinding_SideEffects_RegressionData.csv'), ...
    fullfile(outDir, 'Table_BlindingAccuracy_ByGroup.csv'));

%% ===================== LOCAL HELPER =====================
function xnum = local_toNumeric(x)
    % Robust conversion of table column to numeric double
    if isnumeric(x)
        xnum = double(x);
        return;
    end
    xs = string(x);
    xs = strtrim(xs);
    xs(xs=="") = "NaN";
    tok = regexp(xs, '(-?\d+(\.\d+)?)', 'tokens', 'once');
    xnum = nan(size(xs));
    hasTok = ~cellfun(@isempty, tok);
    xnum(hasTok) = str2double(string(cellfun(@(c)c{1}, tok(hasTok), 'UniformOutput', false)));
end
