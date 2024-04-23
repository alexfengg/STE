clear
close all

res = readmatrix('result/fund_result_rotation.csv');
MEDIAN = res(:,1:2:end);
MEAN = res(:,2:2:end);

%% rotation

% Define the dataset names
data_names = {
    'Alamo', 'Ellis Island', 'Madrid Metropolis',...
    'Montreal ND', 'NYC Library', 'Notre Dame', 'Piazza del Popolo',...
    'Piccadilly', 'Roman Forum', 'Tower of London',...
    'Union Square', 'Vienna Cathedral', 'Yorkminster','Gendarmenmarkt'};

% Create a bar graph with a more modern color scheme
figure('Position', [100, 100, 1200, 400]);
bar_handle = bar(1:14,MEDIAN, 'grouped');

% Use the lines colormap for the bars
colormap(lines);

% Enhance the grid lines
grid on;
set(gca, 'GridLineStyle', '--');

% Add x-axis labels with rotation for better readability
set(gca, 'xtick', 1:length(data_names), 'xticklabel', data_names);
xtickangle(45);

% Set y-axis label with a larger font size
ylabel('Rotation Error (Median)', 'FontSize', 14);

% Add a legend with a more readable font size
legend('RANSAC','DEGENSAC','LO-RANSAC','TME','FMS','SFMS','STE', ...
    'FontSize',12,'Location','best');

% Adjust the width of the bars
for k = 1:5
    set(bar_handle(k), 'BarWidth', 0.85);
end

% Tighten the axes for a better fit
axis tight;

% Improve the figure's appearance by modifying the following properties:
set(gca, 'FontSize', 14); % General font size
set(gca, 'Box', 'off'); % Turn off the box surrounding the whole axes
set(gcf, 'Color', 'w'); % Set the background color to white

% Save the figure to a file
saveas(gcf, 'result/fund_median_r.jpg');

%%

% Define the dataset names

% Create a bar graph with a more modern color scheme
figure('Position', [100, 100, 1200, 400]);
bar_handle = bar(1:14,MEAN, 'grouped');

% Use the lines colormap for the bars
colormap(lines);

% Enhance the grid lines
grid on;
set(gca, 'GridLineStyle', '--');

% Add x-axis labels with rotation for better readability
set(gca, 'xtick', 1:length(data_names), 'xticklabel', data_names);
xtickangle(45);

% Set y-axis label with a larger font size
ylabel('Rotation Error (Mean)', 'FontSize', 14);

% Add a legend with a more readable font size
legend('RANSAC','DEGENSAC','LO-RANSAC','TME','FMS','SFMS','STE', ...
    'FontSize',12,'Location','best');

% Adjust the width of the bars
for k = 1:5
    set(bar_handle(k), 'BarWidth', 0.85);
end

% Tighten the axes for a better fit
axis tight;

% Improve the figure's appearance by modifying the following properties:
set(gca, 'FontSize', 14); % General font size
set(gca, 'Box', 'off'); % Turn off the box surrounding the whole axes
set(gcf, 'Color', 'w'); % Set the background color to white

% Save the figure to a file
saveas(gcf, 'result/fund_mean_r.jpg');

%% translation

res = readmatrix('result/fund_result_translation.csv');
MEDIAN = res(:,1:2:end);
MEAN = res(:,2:2:end);

% Create a bar graph with a more modern color scheme
figure('Position', [100, 100, 1200, 400]);
bar_handle = bar(1:14,MEDIAN, 'grouped');

% Use the lines colormap for the bars
colormap(lines);

% Enhance the grid lines
grid on;
set(gca, 'GridLineStyle', '--');

% Add x-axis labels with rotation for better readability
set(gca, 'xtick', 1:length(data_names), 'xticklabel', data_names);
xtickangle(45);

% Set y-axis label with a larger font size
ylabel('Direction Vector Error (Median)', 'FontSize', 14);

% Add a legend with a more readable font size
legend('RANSAC','DEGENSAC','LO-RANSAC','TME','FMS','SFMS','STE', ...
    'FontSize',10,'Location','best');
% Adjust the width of the bars
for k = 1:5
    set(bar_handle(k), 'BarWidth', 0.85);
end

% Tighten the axes for a better fit
axis tight;

% Improve the figure's appearance by modifying the following properties:
set(gca, 'FontSize', 14); % General font size
set(gca, 'Box', 'off'); % Turn off the box surrounding the whole axes
set(gcf, 'Color', 'w'); % Set the background color to white

% Save the figure to a file
saveas(gcf, 'result/fund_median_t.jpg');

%%

% Create a bar graph with a more modern color scheme
figure('Position', [100, 100, 1200, 400]);
bar_handle = bar(1:14,MEAN, 'grouped');

% Use the lines colormap for the bars
colormap(lines);

% Enhance the grid lines
grid on;
set(gca, 'GridLineStyle', '--');

% Add x-axis labels with rotation for better readability
set(gca, 'xtick', 1:length(data_names), 'xticklabel', data_names);
xtickangle(45);

% Set y-axis label with a larger font size
% ylabel('Directional Error (Median)', 'FontSize', 14);
ylabel('Direction Vector Error (Mean)', 'FontSize', 14);

% Add a legend with a more readable font size
legend('RANSAC','DEGENSAC','LO-RANSAC','TME','FMS','SFMS','STE', ...
    'FontSize',10,'Location','best');

% Adjust the width of the bars
for k = 1:5
    set(bar_handle(k), 'BarWidth', 0.85);
end

% Tighten the axes for a better fit
axis tight;

% Improve the figure's appearance by modifying the following properties:
set(gca, 'FontSize', 14); % General font size
set(gca, 'Box', 'off'); % Turn off the box surrounding the whole axes
set(gcf, 'Color', 'w'); % Set the background color to white

% Save the figure to a file
saveas(gcf, 'result/fund_mean_t.jpg');
