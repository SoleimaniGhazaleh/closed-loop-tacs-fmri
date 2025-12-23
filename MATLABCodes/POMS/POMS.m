%% ============================================================
%  POMS: read TSV, strip "sub-" prefix, keep only your subjects,
%  robustly convert items to numeric, compute subscales,
%  merge groups, rename groups (down/up regulation),
%  keep ONLY before/after sessions (drop during),
%  DESCRIPTIVES by Group x Session,
%  (A) Between-group comparisons within each session (Welch t-tests)
%  (B) Pre vs Post within each group (paired t-tests where possible)
%  save outputs.
% ============================================================

clear; clc;

%% ---- 1) Paths ----
pomsFile = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/profile_of_mood_states_poms.tsv';

%% ---- 2) Define your participant list (SUBJECTS + GROUP) ----
ScanID = { ...
    'E16426'; 'E16449'; 'E16488'; 'E16502'; 'E16518'; 'E16520'; ...
    'E16531'; 'E16537'; 'E16561'; 'E16591'; 'E16599'; 'E16616'; ...
    'E16628'; 'E16682'; 'E16844'; 'E16909'; 'E16942'; 'E16979'; 'E17009'};

SubjectID = { ...
    'BK144'; 'AZ290'; 'AA115'; 'AV531'; 'BK798'; 'BK832'; ...
    'AR165'; 'BL003'; 'BA977'; 'AZ275'; 'AV221'; 'AP913'; ...
    'AQ975'; 'AU393'; 'AB018'; 'AV438'; 'AV503'; 'BK060'; 'BH126'};

Age = [ ...
    44; 46; 52; 24; 20; 35; 28; 53; 41; 45; 39; 32; ...
    37; 53; 42; 34; 33; 46; 20];

Gender = { ...
    'F'; 'F'; 'F'; 'M'; 'M'; 'F'; 'F'; 'F'; 'F'; 'F'; ...
    'M'; 'M'; 'M'; 'F'; 'F'; 'F'; 'M'; 'M'; 'M'};

Group = { ...
    'Experimental'; 'Control'; 'Control'; 'Experimental'; 'Control'; ...
    'Experimental'; 'Control'; 'Experimental'; 'Experimental'; 'Experimental'; ...
    'Control'; 'Experimental'; 'Control'; 'Control'; 'Experimental'; ...
    'Control'; 'Control'; 'Experimental'; 'Experimental'};

participants = table(ScanID, SubjectID, Age, Gender, Group);
participants.SubjectID = string(participants.SubjectID);
participants.Gender = categorical(participants.Gender);
participants.Group  = categorical(participants.Group);

%% ---- 3) Read POMS TSV ----
if ~isfile(pomsFile)
    error('POMS file not found:\n%s', pomsFile);
end

opts = detectImportOptions(pomsFile, 'FileType','text', 'Delimiter','\t');
opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','NaN','nan','N/A','n/a','',' '});
poms = readtable(pomsFile, opts);

if ~ismember('participant_id', poms.Properties.VariableNames)
    error('TSV must contain a column named "participant_id".');
end
if ~ismember('session', poms.Properties.VariableNames)
    error('TSV must contain a column named "session".');
end

%% ---- 4) Normalize IDs and keep only before/after ----
poms.participant_id = string(poms.participant_id);
poms.SubjectID = erase(poms.participant_id, "sub-");
poms.SubjectID = strtrim(poms.SubjectID);

poms.session = string(poms.session);
poms.session = strtrim(poms.session);

% Keep only before/after (drop during)
sessionKeep = ["before_fps_arm_1","after_fps_arm_1"];
poms = poms(ismember(poms.session, sessionKeep), :);

%% ---- 5) Keep only subjects in your participant list ----
inList = ismember(poms.SubjectID, participants.SubjectID);
if ~any(inList)
    error('No matching IDs found after removing "sub-". Check ID formatting.');
end
poms = poms(inList,:);

%% ---- 6) Robust conversion: poms_1..poms_65 to numeric 0–4 ----
itemVars = "poms_" + string(1:65);
itemVars = itemVars(ismember(itemVars, string(poms.Properties.VariableNames)));

for i = 1:numel(itemVars)
    vn = itemVars(i);
    x  = poms.(vn);

    if isnumeric(x)
        poms.(vn) = double(x);
        continue;
    end

    xs = string(x);
    tok = regexp(xs, '([0-4])', 'tokens', 'once'); % works for "2 - moderately"
    num = nan(size(xs));
    hasTok = ~cellfun(@isempty, tok);
    num(hasTok) = str2double(string(tok(hasTok)));
    poms.(vn) = num;
end

%% ---- 7) Compute POMS subscale scores ----
poms.poms_anger_score = ...
    poms.poms_3  + poms.poms_12 + poms.poms_17 + poms.poms_24 + ...
    poms.poms_31 + poms.poms_33 + poms.poms_39 + poms.poms_42 + ...
    poms.poms_47 + poms.poms_52 + poms.poms_53 + poms.poms_57;

poms.poms_tension_score = ...
    poms.poms_2  + poms.poms_10 + poms.poms_16 + poms.poms_20 + ...
    (4 - poms.poms_22) + poms.poms_26 + poms.poms_27 + ...
    poms.poms_34 + poms.poms_41;

poms.poms_depression_score = ...
    poms.poms_5  + poms.poms_9  + poms.poms_14 + poms.poms_18 + ...
    poms.poms_21 + poms.poms_23 + poms.poms_32 + poms.poms_35 + ...
    poms.poms_36 + poms.poms_44 + poms.poms_45 + poms.poms_48 + ...
    poms.poms_58 + poms.poms_61 + poms.poms_62;

poms.poms_fatigue_score = ...
    poms.poms_4  + poms.poms_11 + poms.poms_29 + poms.poms_40 + ...
    poms.poms_46 + poms.poms_49 + poms.poms_65;

poms.poms_confusion_score = ...
    poms.poms_8  + poms.poms_28 + poms.poms_37 + poms.poms_50 + ...
    (4 - poms.poms_54) + poms.poms_59 + poms.poms_64;

poms.poms_vigour_score = ...
    poms.poms_7  + poms.poms_15 + poms.poms_19 + poms.poms_38 + ...
    poms.poms_51 + poms.poms_56 + poms.poms_60 + poms.poms_63;

poms.poms_tmd_score = ...
    poms.poms_tension_score + poms.poms_depression_score + poms.poms_anger_score + ...
    poms.poms_fatigue_score + poms.poms_confusion_score - poms.poms_vigour_score;

%% ---- 8) Merge participant info ----
poms_all = innerjoin(poms, participants, 'Keys','SubjectID');

%% ---- 9) Rename groups ----
poms_all.Group = string(poms_all.Group);
poms_all.Group(poms_all.Group == "Control")       = "down-regulation";
poms_all.Group(poms_all.Group == "Experimental") = "up-regulation";
poms_all.Group = categorical(poms_all.Group, ["down-regulation","up-regulation"]);

%% ---- 10) Clean session labels ----
poms_all.SessionRaw = poms_all.session;

poms_all.Session = poms_all.session;
poms_all.Session(poms_all.Session == "before_fps_arm_1") = "pre";
poms_all.Session(poms_all.Session == "after_fps_arm_1")  = "post";
poms_all.Session = categorical(poms_all.Session, ["pre","post"]);

%% ---- 11) Descriptives by Group x Session ----
scoreVars = ["poms_tension_score","poms_depression_score","poms_anger_score", ...
             "poms_fatigue_score","poms_confusion_score","poms_vigour_score","poms_tmd_score"];

grpLevels = categories(poms_all.Group);
sesLevels = categories(poms_all.Session);

descGS = table;

for si = 1:numel(sesLevels)
    sname = sesLevels{si};
    for gi = 1:numel(grpLevels)
        gname = grpLevels{gi};

        idx = (poms_all.Session == sname) & (poms_all.Group == gname);

        for vi = 1:numel(scoreVars)
            v = scoreVars(vi);
            x = poms_all{idx, v};
            x = x(~isnan(x));

            N   = numel(x);
            mu  = mean(x);
            sd  = std(x);
            sem = sd / sqrt(max(N,1));

            descGS = [descGS; table( ...
                string(sname), string(gname), string(v), N, mu, sd, sem, ...
                'VariableNames',{'Session','Group','Measure','N','Mean','SD','SEM'})]; %#ok<AGROW>
        end
    end
end

fprintf('\n===== POMS DESCRIPTIVES: Group x Session =====\n');
disp(descGS);

%% ---- 12A) Between-group comparisons WITHIN each session (Welch t-tests) ----
betweenTbl = table;

for si = 1:numel(sesLevels)
    sname = sesLevels{si};

    for vi = 1:numel(scoreVars)
        v = scoreVars(vi);

        xDown = poms_all{(poms_all.Session == sname) & (poms_all.Group == "down-regulation"), v};
        xUp   = poms_all{(poms_all.Session == sname) & (poms_all.Group == "up-regulation"),   v};

        xDown = xDown(~isnan(xDown));
        xUp   = xUp(~isnan(xUp));

        if numel(xDown) >= 2 && numel(xUp) >= 2
            [~, p, ~, st] = ttest2(xDown, xUp, 'Vartype','unequal'); % Welch
            tval = st.tstat;
            df   = st.df;
        else
            p = NaN; tval = NaN; df = NaN;
        end

        betweenTbl = [betweenTbl; table( ...
            string(sname), string(v), ...
            numel(xDown), mean(xDown,'omitnan'), std(xDown,'omitnan'), ...
            numel(xUp),   mean(xUp,'omitnan'),   std(xUp,'omitnan'), ...
            tval, df, p, ...
            'VariableNames',{'Session','Measure', ...
            'N_down','Mean_down','SD_down', ...
            'N_up','Mean_up','SD_up', ...
            't','df','p_value'})]; %#ok<AGROW>
    end
end

fprintf('\n===== BETWEEN-GROUP COMPARISONS (within each session; Welch t-tests) =====\n');
disp(betweenTbl);

%% ---- 12B) Pre vs Post WITHIN each group (paired t-tests on matched subjects) ----
withinTbl = table;

for gi = 1:numel(grpLevels)
    gname = grpLevels{gi};

    grpData = poms_all(poms_all.Group == gname, :);

    for vi = 1:numel(scoreVars)
        v = scoreVars(vi);

        % Build paired table: one row per SubjectID with pre and post
        tmp = grpData(:, {'SubjectID','Session', char(v)});
        tmp.Properties.VariableNames{3} = 'Score';

        wide = unstack(tmp, 'Score', 'Session');  % columns: pre, post (if present)
        if ~ismember('pre', wide.Properties.VariableNames) || ~ismember('post', wide.Properties.VariableNames)
            % If one column is missing entirely, skip
            withinTbl = [withinTbl; table(string(gname), string(v), 0, NaN, NaN, NaN, NaN, ...
                'VariableNames',{'Group','Measure','N_pairs','Mean_pre','Mean_post','t','df','p_value'})]; %#ok<AGROW>
            continue;
        end

        pre  = wide.pre;
        post = wide.post;

        % Keep only complete pairs
        ok = ~isnan(pre) & ~isnan(post);
        pre  = pre(ok);
        post = post(ok);

        Npairs = numel(pre);

        if Npairs >= 2
            [~, p, ~, st] = ttest(pre, post); % paired t-test
            tval = st.tstat;
            df   = st.df;
        else
            p = NaN; tval = NaN; df = NaN;
        end

        withinTbl = [withinTbl; table( ...
            string(gname), string(v), Npairs, ...
            mean(pre,'omitnan'), mean(post,'omitnan'), ...
            tval, df, p, ...
            'VariableNames',{'Group','Measure','N_pairs','Mean_pre','Mean_post','t','df','p_value'})]; %#ok<AGROW>
    end
end

fprintf('\n===== WITHIN-GROUP PRE vs POST (paired t-tests on matched subjects) =====\n');
disp(withinTbl);

%% ---- 13) Save outputs ----
outDir = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/MATLABCode';
writetable(poms_all,    fullfile(outDir, 'poms_scored_merged_onlySubjectList_prepost.csv'));
writetable(descGS,      fullfile(outDir, 'poms_descriptives_GroupXSession_prepost.csv'));
writetable(betweenTbl,  fullfile(outDir, 'poms_betweenGroup_withinSession_Welch.csv'));
writetable(withinTbl,   fullfile(outDir, 'poms_withinGroup_pre_vs_post_paired.csv'));

fprintf('\nSaved:\n  %s\n  %s\n  %s\n  %s\n', ...
    fullfile(outDir, 'poms_scored_merged_onlySubjectList_prepost.csv'), ...
    fullfile(outDir, 'poms_descriptives_GroupXSession_prepost.csv'), ...
    fullfile(outDir, 'poms_betweenGroup_withinSession_Welch.csv'), ...
    fullfile(outDir, 'poms_withinGroup_pre_vs_post_paired.csv'));

%% ---- 14) Quick sanity check: item range (should be ~0-4) ----
allItems = poms_all{:, cellstr(itemVars)};
fprintf('\nSanity check (items): Min = %.2f, Max = %.2f\n', ...
    min(allItems(:),[],'omitnan'), max(allItems(:),[],'omitnan'));

%% ---- Save the three reported tables ----
outDir = fileparts(pomsFile);

descFile    = fullfile(outDir, 'Table1_POMS_Descriptives_GroupXSession.csv');
betweenFile = fullfile(outDir, 'Table2_POMS_BetweenGroup_WithinSession.csv');
withinFile  = fullfile(outDir, 'Table3_POMS_WithinGroup_PrePost.csv');

writetable(descGS,     descFile);
writetable(betweenTbl, betweenFile);
writetable(withinTbl,  withinFile);

fprintf('\nSaved reported tables:\n');
fprintf('1) %s\n', descFile);
fprintf('2) %s\n', betweenFile);
fprintf('3) %s\n', withinFile);

