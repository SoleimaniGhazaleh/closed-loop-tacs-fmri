clear all
close all
clc

% Load data
filePath = '/Users/ghazaleh/Downloads/Connectivity.csv';
tbl = readtable(filePath);

% Define group assignments
expScans = {'E16426','E16502','E16520','E16537','E16561','E16591','E16616','E16844','E16979','E17009'};
ctrlScans = {'E16449','E16488','E16518','E16531','E16599','E16628','E16682','E16909','E16942'};

allScans = unique(tbl{:,1});  % first column is Scan ID
nConn = 15;

% Initialize
TR1_exp = []; TR2_exp = []; Test_exp = [];
TR1_ctrl = []; TR2_ctrl = []; Test_ctrl = [];

for i = 1:length(allScans)
    scanID = string(allScans(i));
    idx = strcmp(tbl{:,1}, scanID);  % rows for this subject

    % Extract 15 rows of TR1, TR2, Test
    this_TR1 = tbl{idx, 2};  % 2nd column = TR1
    this_TR2 = tbl{idx, 3};  % 3rd column = TR2
    this_Test = tbl{idx, 4}; % 4th column = Test

    if numel(this_TR1) ~= nConn
        warning('Skipping %s: not 15 rows', scanID);
        continue;
    end

    % Append to correct group
    if ismember(scanID, expScans)
        TR1_exp(end+1,:) = this_TR1';
        TR2_exp(end+1,:) = this_TR2';
        Test_exp(end+1,:) = this_Test';
    elseif ismember(scanID, ctrlScans)
        TR1_ctrl(end+1,:) = this_TR1';
        TR2_ctrl(end+1,:) = this_TR2';
        Test_ctrl(end+1,:) = this_Test';
    else
        warning('Scan ID %s not in group list', scanID);
    end
end

% Combine data across timepoints
exp_all = [TR1_exp, TR2_exp, Test_exp];    % size: nExp × 45
ctrl_all = [TR1_ctrl, TR2_ctrl, Test_ctrl]; % size: nCtrl × 45

% Compute mean and SEM
mean_exp = mean(exp_all, 1);  sem_exp = std(exp_all, 0, 1) ./ sqrt(size(exp_all,1));
mean_ctrl = mean(ctrl_all, 1); sem_ctrl = std(ctrl_all, 0, 1) ./ sqrt(size(ctrl_all,1));

% Plot
x = 1:45;
expColor = [0.1 0.7 0.7];
ctrlColor = [0.8 0.2 0.2];

figure; hold on;

% Experimental: shaded error and line
fill([x fliplr(x)], [mean_exp - sem_exp fliplr(mean_exp + sem_exp)], ...
    expColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p1 = plot(x, mean_exp, 'o-', 'Color', expColor, 'LineWidth', 2);

% Control: shaded error and line
fill([x fliplr(x)], [mean_ctrl - sem_ctrl fliplr(mean_ctrl + sem_ctrl)], ...
    ctrlColor, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
p2 = plot(x, mean_ctrl, 'o-', 'Color', ctrlColor, 'LineWidth', 2);

% Timepoint markers
xline(15.5, '--k'); xline(30.5, '--k');
xticks([7.5, 22.5, 37.5]); xticklabels({'TR1', 'TR2', 'Test'});

xlabel('Timepoint');
ylabel('Connectivity');
legend([p1 p2], {'Experimental', 'Control'}, 'Location', 'best');
title('Connectivity (mean ± SEM) across Timepoints');
grid on;


