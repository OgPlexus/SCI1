%% Initialization
close all; 
clear all; 
clc;

%% Simulation parameters
numRuns = 5;                   % Number of simulation runs
num_steps = 250;               % Number of time steps in each simulation

% Storage for results across runs
I2_all = zeros(numRuns, num_steps + 1);     % Infections in Group 2
I1_all = zeros(numRuns, num_steps + 1);     % Infections in Group 1
newI2_all = zeros(numRuns, num_steps + 1);  % New infections in Group 2
newI1_all = zeros(numRuns, num_steps + 1);  % New infections in Group 1
eft_contacts_all = zeros(numRuns, 2000);    

% Parameters
num_agents = 2000;             % Total number of agents
radius = 0.02;                 % Radius for neighborhood interactions
group2_fraction = 0.5;         % 50% belong to Group 2, rest to Group 1
move_range = 0.08;             % Maximum step size for random movement

% Different movement fractions for Group 2 agents moving to Group 1 region
movement_fractionD = [0, 0.02];  

% To store epidemic growth rate and mean infected over time
EGR = zeros(length(movement_fractionD), 1);
mean_I2m = zeros(num_steps + 1, length(movement_fractionD));
mean_I1m = zeros(num_steps + 1, length(movement_fractionD));
mean_Totm = zeros(num_steps + 1, length(movement_fractionD));
%% %%%%%%%%%%
for itr1 = 1:length(movement_fractionD)
    movement_fraction = movement_fractionD(itr1);  % Proportion of Group 2 agents temporarily moving to Group 1 region
    mean_exposure_P = [];
    mean_infection_P = [];
    
    for run = 1:numRuns
        rng(1233 + run);  % Ensure reproducibility per run

        % Initialize exposure and infection timers
        time_to_exposure = [];
        exposure_time = zeros(num_agents, 1);
        infection_time = zeros(num_agents, 1);

        % Divide agents into Group 1 and Group 2
        group2_size = round(num_agents * group2_fraction);
        group1_size = num_agents - group2_size;

        % Initial agent health state: 1 = Susceptible, 2 = Exposed, 3 = Infectious, 4 = Recovered
        S = ones(num_agents, 1);
        positions = rand(num_agents, 2);  % Initialize random positions

        % Assign x-coordinates based on group regions
        positions(1:group2_size, 1) = 0.49 * rand(group2_size, 1);               % Group 2 on left
        positions(group2_size+1:end, 1) = 0.51 + 0.49 * rand(group1_size, 1);    % Group 1 on right

        % Midline used to restrict movement between groups
        midline = max(positions(:,1)) / 2;
        original_positions = positions;

        % Assign contact rates (log-normal) for each group
        mu_C2 = 5; sigma_C2 = 0.4;   % Group 2: high contact
        mu_C1 = 0.5; sigma_C1 = 0.3; % Group 1: low contact
        C = zeros(num_agents, 1);
        C(1:group2_size) = lognrnd(mu_C2, sigma_C2, group2_size, 1);
        C(group2_size+1:end) = lognrnd(mu_C1, sigma_C1, group1_size, 1);

        % Assign infection probabilities (Beta-distributed) per group
        p = zeros(num_agents, 1);
        p(1:group2_size) = betarnd(2, 5, group2_size, 1);      % Group 2: higher
        p(group2_size+1:end) = betarnd(1, 5, group1_size, 1);  % Group 1: lower

        % Assign infection durations (exponential) and compute group-wise means
        infection_P = zeros(num_agents, 1);
        infection_P(1:group2_size) = exprnd(10, group2_size, 1);
        infection_P(group2_size+1:end) = exprnd(4, group1_size, 1);
        mean_infection_P(run) = mean(infection_P);

        % Assign exposure durations (exponential) and compute group-wise means
        exposure_P = zeros(num_agents, 1);
        exposure_P(1:group2_size) = exprnd(4, group2_size, 1);
        exposure_P(group2_size+1:end) = exprnd(1, group1_size, 1);
        mean_exposure_P(run) = mean(exposure_P);

        % Seed one initial exposed case in each group
        S(randi(group2_size)) = 2;
        S(group2_size + randi(group1_size)) = 2;

        % Initialize infection counters and loggers
        I2 = zeros(num_steps+1, 1);  % Infectious in Group 2
        I1 = zeros(num_steps+1, 1);  % Infectious in Group 1
        new2_infections = zeros(num_steps+1, 1);
        new1_infections = zeros(num_steps+1, 1);
        I2(1) = 1; I1(1) = 1;
        new2_infections(1) = 1; new1_infections(1) = 1;

        eft_contacts = zeros(num_agents, 1);  % Placeholder (not used further)

        for t = 1:num_steps
            % Count currently infectious agents per group
            I2(t+1) = sum(S(1:group2_size) == 2);
            I1(t+1) = sum(S(group2_size+1:end) == 2);

            % Update timers for exposed and infectious agents
            exposure_time(S == 2) = exposure_time(S == 2) + 1;
            infection_time(S == 3) = infection_time(S == 3) + 1;

            % Transition from exposed to infectious
            S(exposure_time >= exposure_P & S == 2) = 3;

            % Transition from infectious to recovered
            S(infection_time >= infection_P & S == 3) = 4;

            % Reset positions to original before each movement step
            positions = original_positions;

            % Group 2 (left): move all agents within restricted bounds
            group2_movers = (rand(group2_size,2) - 0.5) * move_range;
            new_positions_group2 = positions(1:group2_size, :) + group2_movers;
            new_positions_group2(:,1) = min(original_positions(1:group2_size,1), new_positions_group2(:,1));
            new_positions_group2(new_positions_group2(:,1) > midline, 1) = original_positions(new_positions_group2(:,1) > midline, 1);
            positions(1:group2_size, :) = new_positions_group2;

            % Group 1 (right): move only non-infectious agents
            group1_moving_mask = (S(group2_size+1:end) ~= 2);
            group1_movers = (rand(sum(group1_moving_mask),2) - 0.5) * move_range;
            new_positions_group1 = positions(group2_size+1:end, :);
            new_positions_group1(group1_moving_mask, :) = new_positions_group1(group1_moving_mask, :) + group1_movers;
            new_positions_group1(:,1) = max(original_positions(group2_size+1:end,1), new_positions_group1(:,1));
            new_positions_group1(new_positions_group1(:,1) < midline, 1) = original_positions(new_positions_group1(:,1) < midline, 1);
            positions(group2_size+1:end, :) = new_positions_group1;

            % Temporarily move a fraction of Group 2 to Group 1 region
            moving_agents = randperm(group2_size, round(movement_fraction * group2_size));
            temp_positions = positions(moving_agents, :);
            positions(moving_agents, 1) = 0.5 + 0.5 * rand(length(moving_agents), 1);

            % Infection transmission step
            new2_infections_t = 0;
            new1_infections_t = 0;

            for i = 1:num_agents
                if S(i) == 3
                    if positions(i,1) < 0.5
                        same_side = positions(:,1) < 0.5;
                    else
                        same_side = positions(:,1) > 0.5;
                    end

                    local_neighbors = find(S == 1 & same_side & vecnorm(positions - positions(i,:), 2, 2) < radius);
                    num_contacts = min(round(C(i) * length(local_neighbors)), length(local_neighbors));
                    contacted = local_neighbors(randperm(length(local_neighbors), num_contacts));

                    for j = 1:length(contacted)
                        if rand < p(contacted(j))
                            S(contacted(j)) = 2;
                            if contacted(j) > group2_size
                                new1_infections_t = new1_infections_t + 1;
                            else
                                new2_infections_t = new2_infections_t + 1;
                            end
                        end
                    end
                end
            end

            % Store new infections for the current time step
            new2_infections(t+1) = new2_infections_t;
            new1_infections(t+1) = new1_infections_t;

            % Return moved agents to original positions
            positions(moving_agents,:) = temp_positions;
        end

        % Store simulation results for each run
        I2_all(run, :) = I2;
        I1_all(run, :) = I1;
        newI2_all(run, :) = new2_infections;
        newI1_all(run, :) = new1_infections;
        run
    end
    % Compute mean
    meanNew_I2 = mean(newI2_all);
    meanNew_I1 = mean(newI1_all);
    meanNew_Tot=mean(newI1_all+newI2_all);
    Data=[meanNew_I1' meanNew_I2' meanNew_Tot'];
    mean_I2m(:,itr1) = mean(I2_all./group2_size)';
    mean_I1m(:,itr1)= mean(I1_all./group1_size)';
    mean_Totm(:,itr1)=mean((I1_all+I2_all)./num_agents)';

    infected_counts=meanNew_Tot;
    valid_indices = find(infected_counts > 0, 30, 'first');
    time_points = 1:length(infected_counts);
    weights = exp(-time_points(valid_indices));
    if length(valid_indices) < 2
        error('Not enough data points with infections to estimate R0.');
    end

    % Perform linear fit on log-transformed infected counts
    fit_params = polyfit(time_points(valid_indices), log(infected_counts(valid_indices)), 1);

    % Extract growth rate r
    r = fit_params(1); % The slope is the estimated growth rate
    mean_infectious_period=mean(mean_infection_P);%10 4
    mean_expose_period=mean(mean_exposure_P);%4 1
    % Compute R0
    % epsilon1=1/mean_expose_period;
    % gamma = 1/mean_infectious_period; % Ensure this is correctly set
    R0_estimated = (1 + r * mean_expose_period)*(1 + r * mean_infectious_period);
    %disp(['Estimated R0: ', num2str(R0_estimated)]);
     EGR(itr1)= R0_estimated;
    itr1
end


%% Plot SCI estimates for different communities

figure;
plot(movement_fractionD, (EGR - 0.96456) ./ EGR, 'k--', 'LineWidth', 2);  % Group 2
hold on;
plot(movement_fractionD, (EGR - 2.5488) ./ EGR, 'k:', 'LineWidth', 2);   % Group 1
xlabel('Group mixing: $m$', 'Interpreter', 'latex');
ylabel('SCI', 'Interpreter', 'latex');
legend('Community 2 (more vulnerable): $\hat{\mathcal{C}}_2$', ...
       'Community 1 (less vulnerable): $\hat{\mathcal{C}}_1$', ...
       'Interpreter', 'latex', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);

%% Infection proportion plots for selected mixing levels

% Time vector
time = 0:num_steps;
movement_fractionD = linspace(0, 0.8, 100);
selected_indices = [1, 20, 41, 61];  % Select movement fractions
movement_fractionD = movement_fractionD(selected_indices);

% Line styles and color sets
line_styles = {'-', '--', '-.', ':'};
colors_I1 = [1 0 0; 0.8 0 0; 0.6 0 0; 0.4 0 0];
colors_I2 = [0 0 1; 0 0 0.8; 0 0 0.6; 0 0 0.4];
colors_Tot = [0 0 0; 0.2 0.2 0.2; 0.4 0.4 0.4; 0.6 0.6 0.6];

plot_data = {mean_I1m, mean_I2m, mean_Totm};
plot_colors = {colors_I1, colors_I2, colors_Tot};

for fig_idx = 1:3
    figure;
    hold on;
    for m_idx = 1:length(movement_fractionD)
        plot(time, plot_data{fig_idx}(:, m_idx), 'LineStyle', line_styles{m_idx}, ...
            'Color', plot_colors{fig_idx}(m_idx, :), 'LineWidth', 2);
    end
    hold off;

    xlabel('Time Step', 'Interpreter', 'latex');
    ylabel('$\frac{I_j}{N_j}$', 'Interpreter', 'latex');
    legend(arrayfun(@(m) sprintf('$m=%.2f$', m), movement_fractionD, 'UniformOutput', false), ...
           'Interpreter', 'latex', 'Location', 'best');
    grid on;
    set(gca, 'FontSize', 12);
end
%% Combined Infection Proportion Plot with Rescaling
close all;
time = 0:num_steps;
movement_fractionD = linspace(0, 0.8, 100);
selected_indices = [1, 7, 88];  % e.g., low, medium, high mixing
movement_fractionD = movement_fractionD(selected_indices);

line_styles = {'-', '--', ':'};
plot_colors = {[0.6, 0.3, 0.2], [0.8, 0.45, 0.3], [0, 0, 0]};  % I1, I2, Total

figure;
hold on;
for m_idx = 1:length(movement_fractionD)
    plot(time, mean_I2m(:, m_idx) * 1000, 'LineStyle', line_styles{m_idx}, ...
         'Color', plot_colors{2}, 'LineWidth', 2);
    plot(time, mean_I1m(:, m_idx) * 1000, 'LineStyle', line_styles{m_idx}, ...
         'Color', plot_colors{1}, 'LineWidth', 2);
    plot(time, mean_Totm(:, m_idx) * 2000, 'LineStyle', line_styles{m_idx}, ...
         'Color', plot_colors{3}, 'LineWidth', 2);
end
hold off;

xlabel('Time Step', 'Interpreter', 'latex');
ylabel('Infection Proportion ($\frac{I}{N}$)', 'Interpreter', 'latex');
legend_entries = {};
for m_idx = 1:length(movement_fractionD)
    legend_entries{end+1} = sprintf('$I_2$, $m=%.2f$', movement_fractionD(m_idx));
    legend_entries{end+1} = sprintf('$I_1$, $m=%.2f$', movement_fractionD(m_idx));
    legend_entries{end+1} = sprintf('$I_{Tot}$, $m=%.2f$', movement_fractionD(m_idx));
end
legend(legend_entries, 'Interpreter', 'latex', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12);
