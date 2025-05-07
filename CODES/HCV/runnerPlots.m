
close all; clear all;clc;
load('resultsep.mat');
load('resultsphi.mat');
load('resultstau.mat');
%%
close all
% Define an array of markers
% customLineStyles = {
%     [1 1]   % Custom solid line
%     [2 2]   % Custom dashed line (equal dashes and spaces)
%     [4 1]   % Custom dash-dot line (dash with smaller dot)
%     [2 1 1 1] % Custom dash-dot-dot (shorter dash, two dots)
% };
customLineStyles = {
    '-'      % Solid line
    '--'     % Dashed line
    ':'      % Dotted line
    '-.'     % Dash-dot line
};
figure;
i=3;
hold on
for j=1:4
        % Select a marker based on the loop index
    lineStyle = customLineStyles{mod(j-1, length(customLineStyles)) + 1};
        %semilogy(timephi{j},resultsphi{j}(:,i),'LineStyle',lineStyle,'LineWidth',2)
        %semilogy(timeep{j},resultsep{j}(:,i),'LineStyle',lineStyle,'LineWidth',2)
        semilogy(timetau{j},resultstau{j}(:,i),'LineStyle',lineStyle,'LineWidth',2)
        set(gca, 'YScale', 'log');
        years = 5:5:ceil(max(timephi{end})/365); % Calculate number of years
        xticks(years*365); % Set ticks at the start of each year
        xticklabels(years);
         legend('Isolated','$1\%$ mixing','$5\%$ mixing','Homogeneous mixing','location','bestoutside','Orientation','horizontal')
end
hfig = gcf;
fname = 'TauIE1';