%% ============================================================
%  STAI-State ITEM-LEVEL ANALYSIS (named questions) + optional total score
%
%  - Reads TSV
%  - Fixes participant_id by removing "sub-"
%  - Keeps ONLY your subject list
%  - Keeps ONLY pre/post sessions (drops during)
%  - Converts stais_1..stais_20 to numeric (1–4)
%  - Computes optional total score (stais_state_score) with reverse scoring
%  - Renames groups: Control->down-regulation, Experimental->up-regulation
%
%  OUTPUTS (and saves) 3 tables:
%   (1) Descriptives: Group x Session x Measure + Question text
%   (2) Between-group: Welch t-tests within each Session x Measure
%   (3) Within-group: Pre vs Post paired t-tests within Group x Measure
% ============================================================

clear; clc;

%% ---- 1) Paths ----
staiFile = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/phenotype/statetrait_anxiety_inventory_staistate.tsv';

%% ---- 2) Your participant list (SUBJECTS + GROUP) ----
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

%% ---- 3) Read TSV ----
if ~isfile(staiFile)
    error('STAI file not found:\n%s', staiFile);
end

opts = detectImportOptions(staiFile, 'FileType','text', 'Delimiter','\t');
opts = setvaropts(opts, opts.VariableNames, 'TreatAsMissing', {'NA','NaN','nan','N/A','n/a','',' '});
stai = readtable(staiFile, opts);

if ~ismember('participant_id', stai.Properties.VariableNames)
    error('TSV must contain a column named "participant_id".');
end
if ~ismember('session', stai.Properties.VariableNames)
    error('TSV must contain a column named "session".');
end

%% ---- 4) Normalize IDs and keep only pre/post sessions ----
stai.participant_id = string(stai.participant_id);
stai.SubjectID = strtrim(erase(stai.participant_id, "sub-"));

stai.session = strtrim(string(stai.session));

sessionKeep = ["before_fps_arm_1","after_fps_arm_1"];
stai = stai(ismember(stai.session, sessionKeep), :);

%% ---- 5) Keep only subjects in your list ----
stai = stai(ismember(stai.SubjectID, participants.SubjectID), :);
if height(stai) == 0
    error('No matching IDs found after removing "sub-".');
end

%% ---- 6) Convert items stais_1..stais_20 to numeric 1–4 ----
items = "stais_" + string(1:20);
missingItems = setdiff(items, string(stai.Properties.VariableNames));
if ~isempty(missingItems)
    error('TSV missing required items:\n%s', strjoin(missingItems, ', '));
end

for i = 1:numel(items)
    vn = items(i);
    x  = stai.(vn);

    if isnumeric(x)
        stai.(vn) = double(x);
    else
        xs  = string(x);
        tok = regexp(xs, '([1-4])', 'tokens', 'once'); % handles "3 - Moderately so"
        num = nan(size(xs));
        hasTok = ~cellfun(@isempty, tok);
        num(hasTok) = str2double(string(tok(hasTok)));
        stai.(vn) = num;
    end
end

%% ---- 7) Compute optional total score (recommended for classic STAI reporting) ----
includeTotal = true; % set false if you ONLY want item-level questions

if includeTotal
    rev = [1 2 5 8 10 11 15 16 19 20];
    dir = setdiff(1:20, rev);

    Xrev = double(stai{:, cellstr("stais_" + string(rev))});
    Xdir = double(stai{:, cellstr("stais_" + string(dir))});

    stai.stais_state_score = sum(5 - Xrev, 2, 'omitnan') + sum(Xdir, 2, 'omitnan');
    measures = [items, "stais_state_score"];
else
    measures = items;
end

%% ---- 8) Human-readable question labels (Measure -> Question) ----
staiLabels = containers.Map( ...
    { ...
    'stais_1','stais_2','stais_3','stais_4','stais_5','stais_6','stais_7','stais_8','stais_9','stais_10', ...
    'stais_11','stais_12','stais_13','stais_14','stais_15','stais_16','stais_17','stais_18','stais_19','stais_20' ...
    }, ...
    { ...
    'I feel calm', ...
    'I feel secure', ...
    'I am tense', ...
    'I feel strained', ...
    'I feel at ease', ...
    'I feel upset', ...
    'I am presently worrying over possible misfortunes', ...
    'I feel satisfied', ...
    'I feel frightened', ...
    'I feel comfortable', ...
    'I feel self-confident', ...
    'I feel nervous', ...
    'I am jittery', ...
    'I feel indecisive', ...
    'I am relaxed', ...
    'I feel content', ...
    'I am worried', ...
    'I feel confused', ...
    'I feel steady', ...
    'I feel pleasant' ...
    } ...
);

if includeTotal
    staiLabels('stais_state_score') = 'STAI-State Total Score';
end

%% ---- 9) Merge participant info and rename groups ----
stai_all = innerjoin(stai, participants, 'Keys','SubjectID');

stai_all.Group = string(stai_all.Group);
stai_all.Group(stai_all.Group == "Control")       = "down-regulation";
stai_all.Group(stai_all.Group == "Experimental") = "up-regulation";
stai_all.Group = categorical(stai_all.Group, ["down-regulation","up-regulation"]);

stai_all.SessionRaw = stai_all.session;
stai_all.Session = stai_all.session;
stai_all.Session(stai_all.Session == "before_fps_arm_1") = "pre";
stai_all.Session(stai_all.Session == "after_fps_arm_1")  = "post";
stai_all.Session = categorical(stai_all.Session, ["pre","post"]);

grpLevels = categories(stai_all.Group);
sesLevels = categories(stai_all.Session);

%% ---- 10) TABLE 1: Descriptives (Group x Session x Measure + Question) ----
descTbl = table;

for mi = 1:numel(measures)
    m = measures(mi);
    mchar = char(m);

    if isKey(staiLabels, mchar)
        qtxt = staiLabels(mchar);
    else
        qtxt = mchar;
    end

    for si = 1:numel(sesLevels)
        sname = sesLevels{si};
        for gi = 1:numel(grpLevels)
            gname = grpLevels{gi};

            idx = (stai_all.Session == sname) & (stai_all.Group == gname);
            x = double(stai_all{idx, m});
            x = x(~isnan(x));

            N   = numel(x);
            mu  = mean(x);
            sd  = std(x);
            sem = sd / sqrt(max(N,1));

            descTbl = [descTbl; table( ...
                string(sname), string(gname), string(m), string(qtxt), ...
                N, mu, sd, sem, ...
                'VariableNames',{'Session','Group','Measure','Question','N','Mean','SD','SEM'})]; %#ok<AGROW>
        end
    end
end

fprintf('\n===== STAI DESCRIPTIVES (Group x Session; named questions) =====\n');
disp(descTbl);

%% ---- 11) TABLE 2: Between-group comparisons within each session (Welch t-tests) ----
betweenTbl = table;

for mi = 1:numel(measures)
    m = measures(mi);
    mchar = char(m);

    if isKey(staiLabels, mchar)
        qtxt = staiLabels(mchar);
    else
        qtxt = mchar;
    end

    for si = 1:numel(sesLevels)
        sname = sesLevels{si};

        xDown = double(stai_all{(stai_all.Session == sname) & (stai_all.Group == "down-regulation"), m});
        xUp   = double(stai_all{(stai_all.Session == sname) & (stai_all.Group == "up-regulation"),   m});

        xDown = xDown(~isnan(xDown));
        xUp   = xUp(~isnan(xUp));

        if numel(xDown) >= 2 && numel(xUp) >= 2
            [~, p, ~, st] = ttest2(xDown, xUp, 'Vartype','unequal'); % Welch
            tval = st.tstat; df = st.df;
        else
            p = NaN; tval = NaN; df = NaN;
        end

        betweenTbl = [betweenTbl; table( ...
            string(sname), string(m), string(qtxt), ...
            numel(xDown), mean(xDown,'omitnan'), std(xDown,'omitnan'), ...
            numel(xUp),   mean(xUp,'omitnan'),   std(xUp,'omitnan'), ...
            tval, df, p, ...
            'VariableNames',{'Session','Measure','Question', ...
            'N_down','Mean_down','SD_down', ...
            'N_up','Mean_up','SD_up', ...
            't','df','p_value'})]; %#ok<AGROW>
    end
end

fprintf('\n===== STAI BETWEEN-GROUP (within session; named questions) =====\n');
disp(betweenTbl);

%% ---- 12) TABLE 3: Within-group Pre vs Post (paired t-tests on matched subjects) ----
withinTbl = table;

for mi = 1:numel(measures)
    m = measures(mi);
    mchar = char(m);

    if isKey(staiLabels, mchar)
        qtxt = staiLabels(mchar);
    else
        qtxt = mchar;
    end

    for gi = 1:numel(grpLevels)
        gname = grpLevels{gi};

        grpData = stai_all(stai_all.Group == gname, :);

        tmp = grpData(:, {'SubjectID','Session', mchar});
        tmp.Properties.VariableNames{3} = 'Score';

        wide = unstack(tmp, 'Score', 'Session'); % columns: pre, post

        if ~ismember('pre', wide.Properties.VariableNames) || ~ismember('post', wide.Properties.VariableNames)
            withinTbl = [withinTbl; table(string(gname), string(m), string(qtxt), 0, NaN, NaN, NaN, NaN, ...
                'VariableNames',{'Group','Measure','Question','N_pairs','Mean_pre','Mean_post','t','df','p_value'})]; %#ok<AGROW>
            continue;
        end

        pre  = double(wide.pre);
        post = double(wide.post);

        ok = ~isnan(pre) & ~isnan(post);
        pre  = pre(ok);
        post = post(ok);

        Npairs = numel(pre);

        if Npairs >= 2
            [~, p, ~, st] = ttest(pre, post);
            tval = st.tstat; df = st.df;
        else
            p = NaN; tval = NaN; df = NaN;
        end

        withinTbl = [withinTbl; table( ...
            string(gname), string(m), string(qtxt), ...
            Npairs, mean(pre,'omitnan'), mean(post,'omitnan'), ...
            tval, df, p, ...
            'VariableNames',{'Group','Measure','Question','N_pairs','Mean_pre','Mean_post','t','df','p_value'})]; %#ok<AGROW>
    end
end

fprintf('\n===== STAI WITHIN-GROUP PRE vs POST (paired; named questions) =====\n');
disp(withinTbl);

%% ---- 13) Save the 3 reported tables ----
outDir = '/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/Data_DropBox/MATLABCode/STAI';

t1 = fullfile(outDir, 'Table1_STAI_NamedQuestions_Descriptives_GroupXSession.csv');
t2 = fullfile(outDir, 'Table2_STAI_NamedQuestions_BetweenGroup_WithinSession.csv');
t3 = fullfile(outDir, 'Table3_STAI_NamedQuestions_WithinGroup_PrePost.csv');

writetable(descTbl,    t1);
writetable(betweenTbl, t2);
writetable(withinTbl,  t3);

fprintf('\nSaved reported tables:\n1) %s\n2) %s\n3) %s\n', t1, t2, t3);

%% ---- 14) Sanity checks ----
Xitems = double(stai_all{:, cellstr(items)});
fprintf('\nSanity check (items): Min=%.2f, Max=%.2f\n', ...
    min(Xitems(:),[],'omitnan'), max(Xitems(:),[],'omitnan'));

if includeTotal
    fprintf('Sanity check (stais_state_score): Min=%.2f, Max=%.2f\n', ...
        min(stai_all.stais_state_score,[],'omitnan'), max(stai_all.stais_state_score,[],'omitnan'));
end
