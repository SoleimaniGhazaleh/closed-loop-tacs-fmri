% ----------------------------
% Polar Plot: Phase Difference by Group + Styled Enhancements
% ----------------------------

% Frequency data (radius)
clear all
close all
clc

load('/Users/ghazaleh/Documents/1.PROJECTS/Closed-Loop tES-fMRI/StimParam.mat');
Frequency = StimParam(:,1);
Phase_deg = StimParam(:,2);
% Group assignment: 0 = Control, 1 = Experimental
Group = [0, 1, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0];



% Convert to radians
Phase_rad = deg2rad(Phase_deg);

% Group colors
expColor = [0.1 0.7 0.7];  % Teal
ctrlColor = [0.8 0.2 0.2]; % Red

% Create figure and polar axes
figure('Color', 'w');
ax = polaraxes;
hold(ax, 'on');

% Style polar axes
ax.ThetaZeroLocation = 'right';         % 0° on the right
ax.ThetaDir = 'counterclockwise';       % 90° at top
ax.ThetaTick = 0:30:330;
ax.RLim = [0 ceil(max(Frequency))+1];
ax.FontSize = 12;
ax.Color = [0.97 0.97 1];                % Light background
ax.GridColor = [0.6 0.6 0.6];
ax.GridAlpha = 0.3;
ax.ThetaColor = [0.3 0.3 0.3];
ax.RColor = [0.3 0.3 0.3];

% Optional: color each point individually (rainbow style)
use_colormap = false;
if use_colormap
    cmap = turbo(length(Frequency));  % or jet, parula
    for i = 1:length(Frequency)
        polarplot(ax, Phase_rad(i), Frequency(i), 'o', ...
            'MarkerFaceColor', cmap(i,:), ...
            'MarkerEdgeColor', 'k', 'MarkerSize', 8);
    end
    colormap(cmap);
    colorbar('Ticks', [], 'Location', 'eastoutside', 'Box', 'off');
else
    % Standard 2-color group plot
    polarplot(ax, Phase_rad(Group==0), Frequency(Group==0), 'o', ...
        'MarkerFaceColor', ctrlColor, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, 'DisplayName', 'Control');

    polarplot(ax, Phase_rad(Group==1), Frequency(Group==1), '^', ...
        'MarkerFaceColor', expColor, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 8, 'DisplayName', 'Experimental');
end

% Optional: Add subject number labels
% for i = 1:length(Frequency)
%     text(ax, Phase_rad(i), Frequency(i), num2str(i), ...
%         'FontSize', 8, 'Color', [0.2 0.2 0.8], ...
%         'VerticalAlignment', 'bottom', ...
%         'HorizontalAlignment', 'left');
% end

% Frequency labels placed at 270° (bottom)
r_ticks = ax.RTick;
theta_label_rad = deg2rad(225);  % Bottom direction
for r = r_ticks
    text(theta_label_rad, r, sprintf('%d Hz', r), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'top', ...
        'FontSize', 10, 'FontWeight', 'bold', ...
        'Color', [0.2 0.2 0.2]);
end

% Title and Legend
title(ax, 'Polar Plot: Phase Difference vs. Frequency by Group', ...
      'FontWeight', 'bold', 'FontSize', 14);

if ~use_colormap
    legend('Location', 'southoutside', 'Orientation', 'horizontal', ...
           'FontSize', 11, 'Box', 'off');
end

% Optional: Save high-res figure
% print('polar_phase_colored', '-dpng', '-r300');
