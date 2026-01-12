% Script for doing analysis on the packaged behavioural data


%% Run this to undock the figures
set(groot, "defaultFigureWindowStyle", "normal");

%% Load variables and define colours

load("C:\Users\havgana\Desktop\DPhil\packaged_data\combined_sig_day_all_protocols_big_stim_and_post_static_filtered_03_12_25.mat") % sig days
load("C:\Users\havgana\Desktop\DPhil\packaged_data\behaviour_structure_all_animals_all_protocols_04_12_25.mat") % behaviour

% set up two custom colors
cLearner    = [0 0.6 0];   % dark green
cNonLearner = [0.6 0 0.6]; % purple

lick_detection_time_window = [-5, 5]; % Time window in seconds
bin_size = 0.15; % Size of time bins in seconds

% Define time bins for PSTH
time_bins = lick_detection_time_window(1):bin_size:lick_detection_time_window(2);
bin_centers = time_bins(1:end-1) + bin_size / 2;

% Define the surround time based on the packed data
surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

%% Optional: Filter out big stim protocol from the behavioural structure

% Provides an option to set whether to filter out post-big regular sized
% static or to keep them just removing big stim days

% ----------------- CONFIG -----------------
opts.protocol_list = { ...
    'visual_operant_lick_two_stim_right_move', ...
    'visual_operant_lick_two_stim_static', ...
    'visual_operant_lick_two_stim_right_move_big_stim', ...
    'visual_operant_lick_two_stim_static_big_stim', ...
    'stim_wheel*' ...
    };

opts.big_names = { ...
    'visual_operant_lick_two_stim_right_move_big_stim', ...
    'visual_operant_lick_two_stim_static_big_stim' ...
    };

% protocols to drop AFTER the first big_stim, if drop_following=true
opts.following_names = { 'visual_operant_lick_two_stim_static' };
% opts.following_names = { 'stim_wheel*' }; % for HA014 for now

opts.drop_following   = true;   % <--- true if to drop following ; false otherwise
opts.error_on_unknown = true;   % error if a workflow isn't in protocol_list

% ------------------------------------------

animal_list = {behaviour_data.animal_id};

for animal_idx = 1:numel(animal_list)
    rd = behaviour_data(animal_idx).recording_day;

    [keptLocalIdx] = ha.helper_func.select_sessions_by_protocol(rd, opts);

    if isempty(keptLocalIdx)
        fprintf('Animal %s: nothing kept (had %d sessions)\n', ...
            animal_list{animal_idx}, numel(rd));
        continue
    end

    % apply to both structures (keep them in sync)
    behaviour_data(animal_idx).recording_day = rd(keptLocalIdx);

    fprintf('\nAnimal %s: kept %d/%d sessions %s', ...
        animal_list{animal_idx}, numel(keptLocalIdx), numel(rd));
end

%% Create indexing variable for workflows

% make protocol index (n all days x workflow number ordered)
workflow_animal = cellfun(@(x) {x.workflow},{behaviour_data.recording_day},'uni',false);
workflow_cat = grp2idx(horzcat(workflow_animal{:}));

% Create a logical learning index variable (n all days x [0,1])
learning_index_animal = vertcat(combined_sig_day_all_protocols{:});

%% Plot a single animal raster plot and PSTH aligned to stim onset/stim in center

% Define parameters
animal_idx = 2; % Choose your animal index here
lick_detection_time_window = [-5 5];

% Modern color scheme
colors = struct();
colors.rewarded = [0.2, 0.5, 0.8];      % Blue for CS+ (rewarded)
colors.nonrewarded = [0.9, 0.3, 0.3];   % Coral red for CS- (non-rewarded)
colors.stim_marker = [0.2, 0.7, 0.3];   % Green for stimulus onset marker

trial_types = {'right', 'left'}; % Right = rewarded, Left = non-rewarded
trial_labels = {'CS+ (Rewarded)', 'CS- (Non-rewarded)'};

% Extract data for selected animal
animal = animal_list{animal_idx};
recordings = behaviour_data(animal_idx).recording_day;

% Filter recordings by first workflow type
all_workflows = {recordings.workflow};
first_workflow = all_workflows{1};
workflow_mask = strcmp(all_workflows, first_workflow);
workflow_indices = find(workflow_mask);

% Check if we have enough days in this workflow
if length(workflow_indices) < 2
    error('Not enough recording days in workflow "%s" for animal %d. Only %d day(s) found.', ...
        first_workflow, animal_idx, length(workflow_indices));
end

% Select early and late days from the first workflow
early_day_idx = workflow_indices(end-1);
late_day_idx = workflow_indices(end);

% Create figure with 2x2 layout per alignment type
figure('Position', [100, 100, 1400, 1000], 'Color', 'w');

for day_position = 1:2
    if day_position == 1
        day_idx = early_day_idx;
        day_label = 'Early Training';
    else
        day_idx = late_day_idx;
        day_label = 'Late Training';
    end

    rec_day = recordings(day_idx).day;
    behaviour = recordings(day_idx);

    for type_idx = 1:length(trial_types)
        trial_type = trial_types{type_idx};
        trial_label = trial_labels{type_idx};

        if strcmp(trial_type, 'right')
            color = colors.rewarded;
            stim_on_times = behaviour.([trial_type, '_stim_on_times']);
            stim_center_times = behaviour.([trial_type, '_stim_off_times']);
            lick_times = behaviour.lick_event_times;

            reward_times = behaviour.reward_times;
            stim_center_times = stim_center_times(1:min(length(reward_times), length(stim_center_times)));
            [~, sorted_idx] = sort(reward_times - stim_center_times, 'descend');
            stim_on_times = stim_on_times(sorted_idx);
            stim_center_times = stim_center_times(sorted_idx);

        else
            color = colors.nonrewarded;
            stim_on_times = behaviour.([trial_type, '_stim_on_times']);
            stim_center_times= behaviour.([trial_type, '_stim_center_times']);

            % Get the times of the first lick from CS- off
            reward_times = behaviour.all_reward_times_including_first_lick_post_cs_minus (behaviour.cs_labels ==0); % get only CS- trials
            lick_times = behaviour.lick_event_times;

            stim_center_times = stim_center_times(1:min(length(reward_times), length(stim_center_times)));
            [~, sorted_idx] = sort(reward_times - stim_center_times, 'descend');
            stim_on_times = stim_on_times(sorted_idx);
            stim_center_times = stim_center_times(sorted_idx);
        end

% --- STIMULUS ONSET ALIGNED ---
subplot(2, 4, (type_idx-1)*2 + day_position);
hold on;

% Calculate stimulus duration range (min and max across trials)
stim_durations = stim_center_times - stim_on_times;
min_stim_duration = min(stim_durations);
max_stim_duration = max(stim_durations);

% Draw grey shaded region for stimulus duration
y_limits = [0, length(stim_on_times)+1];

% Draw stimulus onset marker
xline(0, '--', 'Color', colors.stim_marker, 'LineWidth', 1.5, ...
      'HandleVisibility', 'off');

% Plot licks
for trial = 1:length(stim_on_times)
    relative_licks = lick_times - stim_on_times(trial);
    licks_in_window = relative_licks(relative_licks >= lick_detection_time_window(1) & ...
                                     relative_licks <= lick_detection_time_window(2));
    scatter(licks_in_window, trial * ones(size(licks_in_window)), 15, color, 'filled', ...
           'MarkerFaceAlpha', 0.7);
end

% Formatting
title([trial_label, ' - ', day_label], 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time from Stim Onset (s)', 'FontSize', 10);
ylabel('Trial #', 'FontSize', 10);
xlim(lick_detection_time_window);
ylim(y_limits);
box on;
set(gca, 'FontSize', 9);
hold off;

% --- STIMULUS OFFSET ALIGNED ---
subplot(2, 4, 4 + (type_idx-1)*2 + day_position);
hold on;

% Calculate stimulus duration range (relative to offset)
relative_min_duration = -max_stim_duration;  % Most negative (earliest onset)
relative_max_duration = -min_stim_duration;  % Least negative (latest onset)

% Draw grey shaded region for stimulus duration (backwards from offset at 0)
y_limits = [0, length(stim_center_times)+1];
patch([relative_min_duration, relative_max_duration, relative_max_duration, relative_min_duration], ...
    [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
    [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
    'HandleVisibility', 'off');

% Draw stimulus offset marker (reference point)
xline(0, '--', 'Color', colors.stim_marker, 'LineWidth', 1.5, ...
    'HandleVisibility', 'off');

% Add movement time indicator (relative to stimulus offset)
if strcmp(trial_type, 'right')  % CS+ trials
    xline(-1, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 2, ...
        'Label', 'Movement', 'LabelVerticalAlignment', 'middle', ...
        'FontSize', 8, 'HandleVisibility', 'off');
else  % CS- trials
    
    xline(-1, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, ...
        'Label', 'Putative', 'LabelVerticalAlignment', 'middle', ...
        'FontSize', 8, 'HandleVisibility', 'off');
end

% Plot licks
for trial = 1:length(stim_center_times)
    relative_licks = lick_times - stim_center_times(trial);
    licks_in_window = relative_licks(relative_licks >= lick_detection_time_window(1) & ...
        relative_licks <= lick_detection_time_window(2));
    scatter(licks_in_window, trial * ones(size(licks_in_window)), 15, color, 'filled', ...
        'MarkerFaceAlpha', 0.7);

    % Stimulus onset marker (where stim started relative to offset)
    relative_stim_on = stim_on_times(trial) - stim_center_times(trial);
    if relative_stim_on >= lick_detection_time_window(1) && relative_stim_on <= lick_detection_time_window(2)
        scatter(relative_stim_on, trial, 20, colors.stim_marker, 'filled', ...
            'MarkerFaceAlpha', 0.5);
    end
end

% Formatting
title([trial_label, ' - ', day_label], 'FontSize', 11, 'FontWeight', 'bold');
xlabel('Time from Stim Offset (s)', 'FontSize', 10);
ylabel('Trial #', 'FontSize', 10);
xlim(lick_detection_time_window);
ylim(y_limits);
box on;
set(gca, 'FontSize', 9);
hold off;
    end
end


% Overall title with workflow information
sgtitle(['Animal ', num2str(animal), ' - Workflow: ', first_workflow, ...
    ' (Day ', recordings(early_day_idx).day, ' vs Day ', recordings(late_day_idx).day, ')'], ...
    'FontSize', 14, 'FontWeight', 'bold');

%% Plot Average PSTHs split based on whether anticipatory lick occured or not (to match the modulation plot)

% Define parameters
lick_detection_time_window = [-2 2];
bin_width = 0.05;
time_bins = lick_detection_time_window(1):bin_width:lick_detection_time_window(2);
bin_centers = time_bins(1:end-1) + bin_width/2;

% Color scheme
colors = struct();
colors.cs_plus_with = [0.2, 0.5, 0.8];
colors.cs_plus_no = [0.6, 0.75, 0.9];
colors.cs_minus_with = [0.9, 0.3, 0.3];
colors.cs_minus_no = [0.95, 0.65, 0.65];
colors.stim_marker = [0.2, 0.7, 0.3];

target_workflow = 'visual_operant_lick_two_stim_right_move';

% Initialize storage: [condition, learning_phase, animal]
% Conditions: 1=CS+ with lick, 2=CS+ no lick, 3=CS- with lick, 4=CS- no lick
% Learning phase: 1=pre, 2=post
all_psth_onset = cell(4, 2, length(animal_list));
all_psth_offset = cell(4, 2, length(animal_list));

% Loop through animals
for animal_idx = 1:length(animal_list)
    recordings = behaviour_data(animal_idx).recording_day;
    all_workflows = {recordings.workflow};
    workflow_mask = strcmp(all_workflows, target_workflow);
    workflow_indices = find(workflow_mask);

    if ~ismember(animal_idx,non_learner_indices)
        continue
    end
    
    if isempty(workflow_indices)
        warning('Animal %d: no recordings found for workflow, skipping', animal_idx);
        continue;
    end
    
    % Get learning day
    sig_days_all = combined_sig_day_all_protocols{animal_idx};
    sig_days = sig_days_all(workflow_mask);
    
    ld = find(sig_days, 1, 'first');
    if isempty(ld)
        warning('Animal %d: no learning day found, skipping', animal_idx);
        continue;
    end
    
    % Separate pre and post learning days
    pre_day_indices = workflow_indices(1:ld-1);
    post_day_indices = workflow_indices(ld:end);
    
    if isempty(pre_day_indices) || isempty(post_day_indices)
        warning('Animal %d: insufficient pre or post learning days, skipping', animal_idx);
        continue;
    end
    
    % Process pre and post separately
    for learning_phase = 1:2
        if learning_phase == 1
            day_indices = pre_day_indices;
        else
            day_indices = post_day_indices;
        end
        
        % Initialize accumulators for this phase
        phase_psth_onset = cell(4, 1);
        phase_psth_offset = cell(4, 1);
        
        % Loop through all days in this phase
        for day_idx = day_indices
            behaviour = recordings(day_idx);
            
            cs_plus_mask = behaviour.cs_labels == 1;
            cs_minus_mask = behaviour.cs_labels == 0;
            anticipatory_licks = behaviour.anticipatory_licks;
            
            % Construct aligned timing arrays
            all_stim_on_times = nan(size(cs_plus_mask));
            all_stim_on_times(cs_plus_mask) = behaviour.right_stim_on_times;
            all_stim_on_times(~cs_plus_mask) = behaviour.left_stim_on_times;
            
            all_stim_off_times = nan(size(cs_plus_mask));
            all_stim_off_times(cs_plus_mask) = behaviour.right_stim_off_times;
            all_stim_off_times(~cs_plus_mask) = behaviour.left_stim_off_times;
            
            lick_times = behaviour.lick_event_times;
            
            n_trials = length(behaviour.cs_labels);
            all_stim_on_times = all_stim_on_times(1:n_trials);
            all_stim_off_times = all_stim_off_times(1:n_trials);
            
            % CS+ trials
            cs_plus_indices = find(cs_plus_mask);
            cs_plus_antc_licks = anticipatory_licks(cs_plus_mask);
            cs_plus_with_lick = cs_plus_indices(cs_plus_antc_licks ~= 0);
            cs_plus_no_lick = cs_plus_indices(cs_plus_antc_licks == 0);
            
            % CS- trials
            cs_minus_indices = find(cs_minus_mask);
            cs_minus_antc_licks = anticipatory_licks(cs_minus_mask);
            cs_minus_with_lick = cs_minus_indices(cs_minus_antc_licks ~= 0);
            cs_minus_no_lick = cs_minus_indices(cs_minus_antc_licks == 0);
            
            conditions = {
                cs_plus_with_lick, 1;
                cs_plus_no_lick, 2;
                cs_minus_with_lick, 3;
                cs_minus_no_lick, 4
            };
            
            % Process each condition for this day
            for cond_idx = 1:4
                trial_subset = conditions{cond_idx, 1};
                cond_num = conditions{cond_idx, 2};
                
                if isempty(trial_subset), continue; end
                
                % ONSET ALIGNED
                psth_onset = zeros(length(trial_subset), length(time_bins)-1);
                for i = 1:length(trial_subset)
                    trial = trial_subset(i);
                    relative_licks = lick_times - all_stim_on_times(trial);
                    psth_onset(i, :) = histcounts(relative_licks, time_bins);
                end
                phase_psth_onset{cond_num} = [phase_psth_onset{cond_num}; psth_onset];
                
                % OFFSET ALIGNED
                psth_offset = zeros(length(trial_subset), length(time_bins)-1);
                for i = 1:length(trial_subset)
                    trial = trial_subset(i);
                    relative_licks = lick_times - all_stim_off_times(trial);
                    psth_offset(i, :) = histcounts(relative_licks, time_bins);
                end
                phase_psth_offset{cond_num} = [phase_psth_offset{cond_num}; psth_offset];
            end
        end
        
        % Store accumulated trials for this animal and learning phase
        for cond_num = 1:4
            all_psth_onset{cond_num, learning_phase, animal_idx} = phase_psth_onset{cond_num};
            all_psth_offset{cond_num, learning_phase, animal_idx} = phase_psth_offset{cond_num};
        end
    end
end


% Plot average PSTHs - ONSET ALIGNED ONLY
figure('Position', [100, 100, 1400, 800], 'Color', 'w');

condition_colors = [colors.cs_plus_with; colors.cs_plus_no; 
                    colors.cs_minus_with; colors.cs_minus_no];
condition_labels = {'With anticipatory lick', 'No anticipatory lick'};

day_labels = ["Early Training", "Late Training"];
max_y_lim= 3.9;

% Layout: 2 rows (CS+, CS-) x 2 columns (Early, Late)
for day_position = 1:2
    
    % ROW 1: CS+ 
    subplot(2, 2, day_position);
    hold on;
    
    for cond_idx = 1:2  % CS+ conditions
        pooled = cat(1, all_psth_onset{cond_idx, day_position, :});
        if isempty(pooled), continue; end
        
        mean_psth = mean(pooled, 1) / bin_width;
        sem_psth = std(pooled, 0, 1) / sqrt(size(pooled, 1)) / bin_width;

        fill([bin_centers fliplr(bin_centers)], ...
            [mean_psth+sem_psth fliplr(mean_psth-sem_psth)], ...
            condition_colors(cond_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        patch([1, 2, 2, 1], ...
            [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');

        plot(bin_centers, mean_psth, 'Color', condition_colors(cond_idx,:), ...
             'LineWidth', 2.5, 'DisplayName', condition_labels{cond_idx});
    end
    
    xline(0, '--', 'Color', colors.stim_marker, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Time from Stim Onset (s)', 'FontSize', 12);
    ylabel('Lick Rate (Hz)', 'FontSize', 12);
    title(['CS+ - ', day_labels(day_position)], 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    ylim([0,max_y_lim]);
    xlim(lick_detection_time_window);
    set(gca, 'LineWidth', 1.5);
    
    % ROW 2: CS-
    subplot(2, 2, 2 + day_position);
    hold on;
    
    for cond_idx = 3:4  % CS- conditions
        pooled = cat(1, all_psth_onset{cond_idx, day_position, :});
        if isempty(pooled), continue; end
        
        mean_psth = mean(pooled, 1) / bin_width;
        sem_psth = std(pooled, 0, 1) / sqrt(size(pooled, 1)) / bin_width;
        
        fill([bin_centers fliplr(bin_centers)], ...
             [mean_psth+sem_psth fliplr(mean_psth-sem_psth)], ...
             condition_colors(cond_idx,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
          patch([1, 2, 2, 1], ...
            [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
            [0.7 0.7 0.7], 'FaceAlpha', 0.2, 'EdgeColor', 'none', ...
            'HandleVisibility', 'off');
        plot(bin_centers, mean_psth, 'Color', condition_colors(cond_idx,:), ...
             'LineWidth', 2.5, 'DisplayName', condition_labels{cond_idx - 2});
    end
    
    xline(0, '--', 'Color', colors.stim_marker, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    xlabel('Time from Stim Onset (s)', 'FontSize', 12);
    ylabel('Lick Rate (Hz)', 'FontSize', 12);
    ylim([0,max_y_lim]);
    title(['CS- - ', day_labels(day_position)], 'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', 10);
    xlim(lick_detection_time_window);
    set(gca, 'LineWidth', 1.5);
end

sgtitle('Lick Behavior Aligned to Stimulus Onset', 'FontSize', 16, 'FontWeight', 'bold');


%% Plots indivdual animal PSTH aligned to stim onset/stim center - CS+ and CS- overlaid

% Define parameters
animal_idx = 2; % Choose your animal index here
lick_detection_time_window = [-5 5];

% Color scheme matching your style
colors = struct();
colors.rewarded = [0.2, 0.5, 0.8];      % Blue for CS+ (rewarded)
colors.nonrewarded = [0.9, 0.3, 0.3];   % Coral/Red for CS- (non-rewarded)

% Extract data for selected animal
animal = animal_list{animal_idx};
recordings = behaviour_data(animal_idx).recording_day;

% Filter recordings by first workflow type
all_workflows = {recordings.workflow};
first_workflow = all_workflows{1};
workflow_mask = strcmp(all_workflows, first_workflow);
workflow_indices = find(workflow_mask);

% Check if we have enough days in this workflow
if length(workflow_indices) < 2
    error('Not enough recording days in workflow "%s" for animal %d. Only %d day(s) found.', ...
        first_workflow, animal_idx, length(workflow_indices));
end

% Select early and late days from the first workflow
early_day_idx = workflow_indices(end-1);
late_day_idx = workflow_indices(end);
valid_rec_days = recordings([early_day_idx, late_day_idx]);

fprintf('Animal %d - Workflow: %s\n', animal_idx, first_workflow);
fprintf('PSTH Comparison: Day %s vs Day %s\n', ...
    valid_rec_days(1).day, valid_rec_days(2).day);

% Define alignment types and corresponding field names
alignment_configs = struct();
alignment_configs(1).event_name = 'Stimulus Onset';
alignment_configs(1).rew_field = 'avg_psth_rewarded_stim_on';
alignment_configs(1).non_rew_field = 'avg_psth_non_rewarded_stim_on';

alignment_configs(2).event_name = 'Stimulus Offset';
alignment_configs(2).rew_field = 'avg_psth_rewarded_stim_final_position';
alignment_configs(2).non_rew_field = 'avg_psth_non_rewarded_stim_final_position';

% Create figure - 2 rows (alignments) × 2 columns (early/late days)
num_alignments = length(alignment_configs);
num_days = length(valid_rec_days);
figure('Position', [100, 100, 1200, 800], 'Color', 'w');

for align_idx = 1:num_alignments
    event_name = alignment_configs(align_idx).event_name;
    rew_field = alignment_configs(align_idx).rew_field;
    non_rew_field = alignment_configs(align_idx).non_rew_field;

    for day_idx = 1:num_days
        subplot(num_alignments, num_days, (align_idx-1)*num_days + day_idx);

        rec_day_struct = valid_rec_days(day_idx);

        % Extract PSTH data
        if isfield(rec_day_struct, rew_field) && isfield(rec_day_struct, non_rew_field)
            rewarded_psth_data = rec_day_struct.(rew_field);
            non_rewarded_psth_data = rec_day_struct.(non_rew_field);

            % Handle different data formats
            if isstruct(rewarded_psth_data)
                bin_centers = rewarded_psth_data.time;
                rewarded_psth = rewarded_psth_data.rate;
                non_rewarded_psth = non_rewarded_psth_data.rate;
            else
                % If just vectors, create bin centers
                rewarded_psth = rewarded_psth_data;
                non_rewarded_psth = non_rewarded_psth_data;
                bin_centers = linspace(lick_detection_time_window(1), ...
                    lick_detection_time_window(2), ...
                    length(rewarded_psth));
            end

            % Plot
            plot(bin_centers, rewarded_psth, 'Color', colors.rewarded, 'LineWidth', 2);
            hold on;
            plot(bin_centers, non_rewarded_psth, 'Color', colors.nonrewarded, 'LineWidth', 2);

            % Formatting
            title([event_name, ' - Day ', rec_day_struct.day], 'Interpreter', 'none', ...
                'FontSize', 11, 'FontWeight', 'bold');
            xlabel('Time (s)', 'FontSize', 10);
            ylabel('Lick Rate (licks/s)', 'FontSize', 10);
            xline(0, '--k', 'LineWidth', 1.5);
            xlim([lick_detection_time_window(1), lick_detection_time_window(2)]);

            % Dynamic y-axis
            max_y = max([max(rewarded_psth), max(non_rewarded_psth)]);
            ylim([0, max(1, max_y * 1.1)]);

            legend({'CS+', 'CS-'}, 'Location', 'northwest', ...
                'FontSize', 9);
            box on;
            set(gca, 'FontSize', 9);
            hold off;
        else
            % Missing fields
            text(0.5, 0.5, sprintf('Missing fields:\n%s or %s', rew_field, non_rew_field), ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', ...
                'FontSize', 10);
            title([event_name, ' - Day ', rec_day_struct.day, ' (NO DATA)'], ...
                'FontSize', 11, 'FontWeight', 'bold');
        end
    end
end

% Overall title
sgtitle(['PSTH Comparison - Animal ', animal_list{animal_idx}, ...
    ' - Workflow: ', first_workflow], 'FontSize', 14, 'FontWeight', 'bold');


%% Plot individual PSTH for all animals aligned to stim on / stim final position

animal_list= {behaviour_data.animal_id};
target_workflow = {'visual_operant_lick_two_stim_static'};  % define workflow to work with

% Iterate over animals
for animal_idx = 1:length(animal_list)

    % Filter recording days based on selected workflows
    rec_days = behaviour_data(animal_idx).recording_day;
    valid_day_mask = ismember({rec_days.workflow}, target_workflow);
    valid_rec_days = rec_days(valid_day_mask);
    num_days = length(valid_rec_days);

    % Skip if no days match
    if num_days == 0
        fprintf('No matching days for animal %s\n', animal_list{animal_idx});
        continue;
    end

    % Extract fieldnames from first valid day
    day_fields = fieldnames(valid_rec_days(1));
    psth_fields = day_fields(contains(day_fields, 'avg_psth_'));

    % Identify matched rewarded / non-rewarded pairs
    event_pairs = {};
    for i = 1:length(psth_fields)
        this_field = psth_fields{i};
        if contains(this_field, 'rewarded') && contains(this_field, 'avg_psth_rewarded')
            non_rew_field = strrep(this_field, 'rewarded', 'non_rewarded');
            if any(strcmp(psth_fields, non_rew_field))
                event_pairs{end+1, 1} = this_field;
                event_pairs{end, 2} = non_rew_field;
            end
        end
    end

    % Iterate over each matched pair (event)
    for event_idx = 1:size(event_pairs, 1)
        rew_field = event_pairs{event_idx, 1};
        non_rew_field = event_pairs{event_idx, 2};

        event_name = strrep(rew_field, 'avg_psth_rewarded_', '');
        event_name = strrep(event_name, '_', ' ');
        event_name = regexprep(event_name, '(^.)', '${upper($1)}');

        figure('Position', [100, 100, 1200, 800]);
        num_columns = 3;
        num_rows = ceil(num_days / num_columns);

        for day_idx = 1:num_days
            subplot(num_rows, num_columns, day_idx);
            rec_day_struct = valid_rec_days(day_idx);

            rewarded_psth = rec_day_struct.(rew_field);
            non_rewarded_psth = rec_day_struct.(non_rew_field);

            % Plot
            plot(bin_centers, rewarded_psth, 'Color', [0.2, 0.2, 0.8], 'LineWidth', 2);
            hold on;
            plot(bin_centers, non_rewarded_psth, 'Color', [0.8, 0.2, 0.2], 'LineWidth', 2);

            title([event_name, ' - ', rec_day_struct.day], 'Interpreter', 'none');
            xlabel('Time (s)');
            ylabel('Lick Rate (licks/s)');
            xline(0, '--k', 'LabelHorizontalAlignment', 'left');
            xlim([lick_detection_time_window(1), lick_detection_time_window(2)]);

            max_y = max([max(rewarded_psth), max(non_rewarded_psth)]);
            ylim([0, max(1, max_y * 1.1)]);

            legend({'Rewarded', 'Non-Rewarded'}, 'Location', 'northwest');
            hold off;
        end

        sgtitle([event_name, ' PSTH - Animal ', animal_list{animal_idx}], 'FontSize', 16);
        set(gcf, 'Color', 'w');
    end
end

%% Plots the average PSTH overlayed CS+ and CS-, seperated by learners vs non learners for all days

max_days=7; % define max days you want to plot

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA007','HA006','HA009','HA013','HA014','HA015'};

% Make a list of the animal_ids in cell format
all_ids = {behaviour_data.animal_id};

% Get indices of learners and non-learners in all_ids
learner_indices = find(ismember(all_ids, learner_ids));
non_learner_indices = find(ismember(all_ids, non_learner_ids));

% Loop over each event (e.g., stimulus onset, center, etc.)
for event_idx = 1:size(event_pairs, 1)
    rew_field = event_pairs{event_idx, 1};
    non_rew_field = event_pairs{event_idx, 2};

    event_name = strrep(rew_field, 'avg_psth_rewarded_', '');
    event_name = strrep(event_name, '_', ' ');
    event_name = regexprep(event_name, '(^.)', '${upper($1)}');


    learner_csplus_by_day = cell(1, max_days);
    learner_csminus_by_day = cell(1, max_days);
    nonlearner_csplus_by_day = cell(1, max_days);
    nonlearner_csminus_by_day = cell(1, max_days);


    % --- Learners ---
    for i = learner_indices
        days = behaviour_data(i).recording_day;

        % Filter only days that match the target workflow
        valid_days = days(arrayfun(@(x) isfield(x, 'workflow') && strcmp(x.workflow, target_workflow), days));

        for d = 1:min(length(valid_days), max_days)
            if isfield(valid_days(d), rew_field)
                learner_csplus_by_day{d}{end+1} = valid_days(d).(rew_field);
            end

            if isfield(valid_days(d), non_rew_field)
                learner_csminus_by_day{d}{end+1} = valid_days(d).(non_rew_field);
            end
        end
    end

    % --- Non-learners ---
    for i = non_learner_indices
        days = behaviour_data(i).recording_day;

        % Filter only days that match the target workflow
        valid_days = days(arrayfun(@(x) isfield(x, 'workflow') && strcmp(x.workflow, target_workflow), days));

        for d = 1:min(length(valid_days), max_days)
            if isfield(valid_days(d), rew_field)
                nonlearner_csplus_by_day{d}{end+1} = valid_days(d).(rew_field);
            end

            if isfield(valid_days(d), non_rew_field)
                nonlearner_csminus_by_day{d}{end+1} = valid_days(d).(non_rew_field);
            end
        end
    end

    % Plot for Learners
    figure('Position', [100, 100, 1200, 800]);
    num_days = length(learner_csplus_by_day);
    num_columns = 3;
    num_rows = ceil(num_days / num_columns);

    for day_idx = 1:num_days
        subplot(num_rows, num_columns, day_idx); hold on;

        % --- CS+ ---
        if ~isempty(learner_csplus_by_day{day_idx})
            csplus_data = cat(1, learner_csplus_by_day{day_idx}{:});
            avg_plus = mean(csplus_data, 1);
            sem_plus = std(csplus_data, 0, 1) / sqrt(size(csplus_data, 1));

            fill([bin_centers, fliplr(bin_centers)], ...
                [avg_plus + sem_plus, fliplr(avg_plus - sem_plus)], ...
                [0.2 0.2 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(bin_centers, avg_plus, 'Color', [0.2 0.2 1], 'LineWidth', 2);
        end

        % --- CS- ---
        if ~isempty(learner_csminus_by_day{day_idx})
            csminus_data = cat(1, learner_csminus_by_day{day_idx}{:});
            avg_minus = mean(csminus_data, 1);
            sem_minus = std(csminus_data, 0, 1) / sqrt(size(csminus_data, 1));

            fill([bin_centers, fliplr(bin_centers)], ...
                [avg_minus + sem_minus, fliplr(avg_minus - sem_minus)], ...
                [1 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(bin_centers, avg_minus, 'Color', [1 0 0], 'LineWidth', 2);
        end

        title(sprintf('Day %d', day_idx));
        xlabel('Time (s)');
        ylabel('Lick Rate (licks/s)');
        xline(0, '--k');
        xlim([lick_detection_time_window(1), lick_detection_time_window(2)]);
        ylim auto;
        legend({'CS+', '', 'CS-', ''}, 'Location', 'northwest');
    end
    sgtitle([event_name,' Learners: CS+ vs CS− PSTH'], 'FontSize', 16);
    set(gcf, 'Color', 'w');

    % Plot for Non-learners
    figure('Position', [150, 150, 1200, 800]);
    num_days = length(nonlearner_csplus_by_day);
    num_rows = ceil(num_days / num_columns);

    for day_idx = 1:num_days
        subplot(num_rows, num_columns, day_idx); hold on;

        % --- CS+ ---
        if ~isempty(nonlearner_csplus_by_day{day_idx})
            csplus_data = cat(1, nonlearner_csplus_by_day{day_idx}{:});
            avg_plus = mean(csplus_data, 1);
            sem_plus = std(csplus_data, 0, 1) / sqrt(size(csplus_data, 1));

            fill([bin_centers, fliplr(bin_centers)], ...
                [avg_plus + sem_plus, fliplr(avg_plus - sem_plus)], ...
                [0.2 0.2 1], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(bin_centers, avg_plus, 'Color', [0.2 0.2 1], 'LineWidth', 2);
        end

        % --- CS- ---
        if ~isempty(nonlearner_csminus_by_day{day_idx})
            csminus_data = cat(1, nonlearner_csminus_by_day{day_idx}{:});
            avg_minus = mean(csminus_data, 1);
            sem_minus = std(csminus_data, 0, 1) / sqrt(size(csminus_data, 1));

            fill([bin_centers, fliplr(bin_centers)], ...
                [avg_minus + sem_minus, fliplr(avg_minus - sem_minus)], ...
                [1 0.4 0.4], 'FaceAlpha', 0.2, 'EdgeColor', 'none');
            plot(bin_centers, avg_minus, 'Color', [1 0 0], 'LineWidth', 2);
        end

        title(sprintf('Day %d', day_idx));
        xlabel('Time (s)');
        ylabel('Lick Rate (licks/s)');
        xline(0, '--k');
        xlim([lick_detection_time_window(1), lick_detection_time_window(2)]);
        ylim auto;
        legend({'CS+', '', 'CS-', ''}, 'Location', 'northwest');
    end
    sgtitle([event_name,' Non-learners: CS+ vs CS− PSTH'], 'FontSize', 16);
    set(gcf, 'Color', 'w');
end
%% Plot the average pre-post learning PSTH for learners vs non-learners aligned to defined event

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA007','HA006','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};


target_workflow = {'visual_operant_lick_two_stim_static'};  % define workflow to work with


% Pre‐allocate containers for each group
learner_plus_pre  = {}; learner_plus_post  = {};
learner_minus_pre = {}; learner_minus_post = {};
non_plus_pre      = {}; non_plus_post      = {};
non_minus_pre     = {}; non_minus_post     = {};

% Define the fields
% event_pairs;

% Define the event

rew_field=  'avg_psth_rewarded_stim_on';
non_rew_field= 'avg_psth_non_rewarded_stim_on';

% rew_field= 'avg_psth_rewarded_stim_final_position';
% non_rew_field= 'avg_psth_non_rewarded_stim_final_position';


% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    % % skip non-learners for now
    % if isLearner
    %     continue;
    % end

    % … inside your animal loop …
    days_all   = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid    = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow,target_workflow), days_all);
    validDays  = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};

    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays,1,'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld= numel(validDays)+1; % define post-learning outside the scope
    end

    % 4) now split into pre/post exactly as before
    plus_pre  = []; minus_pre  = [];
    plus_post = []; minus_post = [];


    for d = 1:numel(validDays)
        if d < ld
            plus_pre  = [plus_pre;  validDays(d).(rew_field)];
            minus_pre = [minus_pre; validDays(d).(non_rew_field)];
        else
            plus_post  = [plus_post;  validDays(d).(rew_field)];
            minus_post = [minus_post; validDays(d).(non_rew_field)];
        end
    end

    T = numel(validDays(d).(rew_field));  % expected PSTH length

    % --- safe per-animal day-averages (row, NaN-safe) ---
    if isempty(plus_pre),   avg_plus_pre   = nan(1,T); else, avg_plus_pre   = reshape(mean(plus_pre,  1,'omitnan'),1,[]); end
    if isempty(plus_post),  avg_plus_post  = nan(1,T); else, avg_plus_post  = reshape(mean(plus_post, 1,'omitnan'),1,[]); end
    if isempty(minus_pre),  avg_minus_pre  = nan(1,T); else, avg_minus_pre  = reshape(mean(minus_pre, 1,'omitnan'),1,[]); end
    if isempty(minus_post), avg_minus_post = nan(1,T); else, avg_minus_post = reshape(mean(minus_post,1,'omitnan'),1,[]); end

    % --- store into the appropriate group ---
    if isLearner
        learner_plus_pre{end+1,1}   = avg_plus_pre;
        learner_plus_post{end+1,1}  = avg_plus_post;
        learner_minus_pre{end+1,1}  = avg_minus_pre;
        learner_minus_post{end+1,1} = avg_minus_post;
    else
        non_plus_pre{end+1,1}   = avg_plus_pre;
        non_plus_post{end+1,1}  = avg_plus_post;
        non_minus_pre{end+1,1}  = avg_minus_pre;
        non_minus_post{end+1,1} = avg_minus_post;
    end

end

% Convert to matrices (animals × time)
learner_pre_CS_plus = cat(1, learner_plus_pre{:});
learner_pre_CS_minus = cat(1, learner_minus_pre{:});
learner_post_CS_plus= cat(1, learner_plus_post{:});
learner_post_CS_minus= cat(1, learner_minus_post{:});

non_learner_pre_CS_plus = cat(1, non_plus_pre{:});
non_learner_pre_CS_minus = cat(1, non_minus_pre{:});
non_learner_post_CS_plus= cat(1, non_plus_post{:});
non_learner_post_CS_minus= cat(1, non_minus_post{:});


% Compute grand means and SEMs
mean_learner_pre_CS_plus  = nanmean(learner_pre_CS_plus,1);  sem_learner_pre_CS_plus  = nanstd(learner_pre_CS_plus,0,1)/sqrt(size(learner_pre_CS_plus,1));
mean_learner_pre_CS_minus = nanmean(learner_pre_CS_minus,1); sem_learner_pre_CS_minus = nanstd(learner_pre_CS_minus,0,1)/sqrt(size(learner_pre_CS_minus,1));
mean_learner_post_CS_plus = nanmean(learner_post_CS_plus,1);  sem_learner_post_CS_plus  = nanstd(learner_post_CS_plus,0,1)/sqrt(size(learner_post_CS_plus,1));
mean_learner_post_CS_minus = nanmean(learner_post_CS_minus,1); sem_learner_post_CS_minus = nanstd(learner_post_CS_minus,0,1)/sqrt(size(learner_post_CS_minus,1));

mean_non_learner_pre_CS_plus  = nanmean(non_learner_pre_CS_plus,1);  sem_non_learner_pre_CS_plus  = nanstd(non_learner_pre_CS_plus,0,1)/sqrt(size(non_learner_pre_CS_plus,1));
mean_non_learner_pre_CS_minus = nanmean(non_learner_pre_CS_minus,1); sem_non_learner_pre_CS_minus = nanstd(non_learner_pre_CS_minus,0,1)/sqrt(size(non_learner_pre_CS_minus,1));
mean_non_learner_post_CS_plus  = nanmean(non_learner_post_CS_plus,1);  sem_non_learner_post_CS_plus  = nanstd(non_learner_post_CS_plus,0,1)/sqrt(size(non_learner_post_CS_plus,1));
mean_non_learner_post_CS_minus = nanmean(non_learner_post_CS_minus,1); sem_non_learner_post_CS_minus = nanstd(non_learner_post_CS_minus,0,1)/sqrt(size(non_learner_post_CS_minus,1));

% Create a global y-axis
y_global_max= max([mean_learner_pre_CS_plus;mean_learner_pre_CS_minus;mean_learner_post_CS_plus;mean_learner_post_CS_minus;mean_non_learner_pre_CS_plus;mean_non_learner_pre_CS_minus;...
    mean_non_learner_post_CS_plus;mean_non_learner_post_CS_minus],[],'all') *1.1;

% Finally: plot group‐averages for CS+ and CS−, pre vs post, learners vs non
figure('Color','w','Position',[100 100 1000 600]);
t = bin_centers;

% CS+ panel
subplot(2,1,1); hold on;
ax = gca;
hold(ax,'on');

% % --- Learner PRE (draw only if data exist) ---
% if ~isempty(mean_learner_pre_CS_plus)
%     % ribbon (no legend)
%     fill(ax, [t fliplr(t)], ...
%         [mean_learner_pre_CS_plus+sem_learner_pre_CS_plus, ...
%         fliplr(mean_learner_pre_CS_plus-sem_learner_pre_CS_plus)], ...
%         cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
% 
%     % mean line (will appear in legend)
%     plot(ax, t, mean_learner_pre_CS_plus, ...
%         'Color', cLearner, 'LineWidth', 2, ...
%         'DisplayName', 'mPFC+ Pre',LineStyle='--');
% end


% --- Learner POST ---
if ~isempty(mean_learner_post_CS_plus)
% ribbon (no legend)
fill(ax, [t fliplr(t)], ...
    [mean_learner_post_CS_plus+sem_learner_post_CS_plus, ...
    fliplr(mean_learner_post_CS_plus-sem_learner_post_CS_plus)], ...
    cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

% mean line (will appear in legend)
plot(ax, t, mean_learner_post_CS_plus, ...
    'Color', cLearner, 'LineWidth', 2, ...
    'DisplayName', 'mPFC+ Post');
end
% --- Non-learner PRE ---%
if ~isempty(mean_non_learner_post_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_plus+sem_non_learner_pre_CS_plus, ...
        fliplr(mean_non_learner_pre_CS_plus-sem_non_learner_pre_CS_plus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_pre_CS_plus, ...
        'color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC- Pre',LineStyle='-');
end

% % --- Non-learner POST (draw only if data exist) ---
% if ~isempty(mean_non_learner_post_CS_plus)
%     % ribbon (no legend)
%     fill(ax, [t fliplr(t)], ...
%         [mean_non_learner_post_CS_plus+sem_non_learner_post_CS_plus, ...
%         fliplr(mean_non_learner_post_CS_plus-sem_non_learner_post_CS_plus)], ...
%         cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
% 
%     % mean line (will appear in legend)
%     plot(ax, t, mean_non_learner_post_CS_plus, ...
%         'Color', cNonLearner, 'LineWidth', 2, ...
%         'DisplayName', 'mPFC- post');
% end

% vertical zero line
xline(ax,0,'k--','HandleVisibility','off');

% build legend from whatever was actually plotted
legend(ax, 'show', 'Location', 'best');

xlabel('Time (s)');
ylabel('CS+ PSTH');
ylim([0,y_global_max]);
title(['Stim Off (Reward Times)', ' — CS+']);


% CS− panel
subplot(2,1,2); hold on;
ax = gca;
hold(ax,'on');


% % --- Learner PRE (draw only if data exist) ---
% if ~isempty(mean_learner_pre_CS_minus)
%     % ribbon (no legend)
%     fill(ax, [t fliplr(t)], ...
%         [mean_learner_pre_CS_minus+sem_learner_pre_CS_minus, ...
%         fliplr(mean_learner_pre_CS_minus-sem_learner_pre_CS_minus)], ...
%         cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
% 
%     % mean line (will appear in legend)
%     plot(ax, t, mean_learner_pre_CS_minus, ...
%         'Color', cLearner, 'LineWidth', 2, ...
%         'DisplayName', 'mPFC+ Pre',LineStyle='--');
% end


% --- Learner POST ---
if ~isempty(mean_learner_post_CS_minus)
% ribbon (no legend)
fill(ax, [t fliplr(t)], ...
    [mean_learner_post_CS_minus+sem_learner_post_CS_minus, ...
    fliplr(mean_learner_post_CS_minus-sem_learner_post_CS_minus)], ...
    cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

% mean line (will appear in legend)
plot(ax, t, mean_learner_post_CS_minus, ...
    'Color', cLearner, 'LineWidth', 2, ...
    'DisplayName', 'mPFC+ Post');
end
% --- Non-learner PRE ---
if ~isempty(mean_non_learner_post_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_minus+sem_non_learner_pre_CS_minus, ...
        fliplr(mean_non_learner_pre_CS_minus-sem_non_learner_pre_CS_minus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_pre_CS_minus, ...
        'color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC- Pre',LineStyle='-');
end
% % --- Non-learner POST (draw only if data exist) ---
% if ~isempty(mean_non_learner_post_CS_minus)
%     % ribbon (no legend)
%     fill(ax, [t fliplr(t)], ...
%         [mean_non_learner_post_CS_minus+sem_non_learner_post_CS_minus, ...
%         fliplr(mean_non_learner_post_CS_minus-sem_non_learner_post_CS_minus)], ...
%         cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
% 
%     % mean line (will appear in legend)
%     plot(ax, t, mean_non_learner_post_CS_minus, ...
%         'Color', cNonLearner, 'LineWidth', 2, ...
%         'DisplayName', 'mPFC- post');
% end


% vertical zero line
xline(ax,0,'k--','HandleVisibility','off');
% build legend from whatever was actually plotted
% legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('CS- PSTH');
ylim([0,y_global_max]);
title(['First Lick from Stim Off', ' — CS-']);


sgtitle(['Stim On',' — Group pre/post learning']);

%% Plot Pre-Post PSTHs aligned to reward times

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA007','HA006','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};


target_workflow = {'visual_operant_lick_two_stim_right_move'};  % define workflow to work with


% Pre‐allocate containers for each group
learner_plus_pre  = {}; learner_plus_post  = {};
learner_minus_pre = {}; learner_minus_post = {};
non_plus_pre      = {}; non_plus_post      = {};
non_minus_pre     = {}; non_minus_post     = {};

% Define the fields


% Define the event

rew_field= 'avg_psth_reward_times';
non_rew_field= 'avg_psth_rewarded_stim_final_position';

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    % skip non-learners for now
    if ~isLearner
        continue;
    end

    % … inside your animal loop …
    days_all   = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid    = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow,target_workflow), days_all);
    validDays  = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};

    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays,1,'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld= numel(validDays)+1; % define post-learning outside the scope
    end

    % 4) now split into pre/post exactly as before
    plus_pre  = []; minus_pre  = [];
    plus_post = []; minus_post = [];


    for d = numel(validDays)
        if d < ld
            plus_pre  = [plus_pre;  validDays(d).(rew_field)];
            minus_pre = [minus_pre; validDays(d).(non_rew_field)];
        else
            plus_post  = [plus_post;  validDays(d).(rew_field)];
            minus_post = [minus_post; validDays(d).(non_rew_field)];
        end
    end

    T = numel(validDays(d).(rew_field));  % expected PSTH length

    % --- safe per-animal day-averages (row, NaN-safe) ---
    if isempty(plus_pre),   avg_plus_pre   = nan(1,T); else, avg_plus_pre   = reshape(mean(plus_pre,  1,'omitnan'),1,[]); end
    if isempty(plus_post),  avg_plus_post  = nan(1,T); else, avg_plus_post  = reshape(mean(plus_post, 1,'omitnan'),1,[]); end
    if isempty(minus_pre),  avg_minus_pre  = nan(1,T); else, avg_minus_pre  = reshape(mean(minus_pre, 1,'omitnan'),1,[]); end
    if isempty(minus_post), avg_minus_post = nan(1,T); else, avg_minus_post = reshape(mean(minus_post,1,'omitnan'),1,[]); end

    % --- store into the appropriate group ---
    if isLearner
        learner_plus_pre{end+1,1}   = avg_plus_pre;
        learner_plus_post{end+1,1}  = avg_plus_post;
        learner_minus_pre{end+1,1}  = avg_minus_pre;
        learner_minus_post{end+1,1} = avg_minus_post;
    else
        non_plus_pre{end+1,1}   = avg_plus_pre;
        non_plus_post{end+1,1}  = avg_plus_post;
        non_minus_pre{end+1,1}  = avg_minus_pre;
        non_minus_post{end+1,1} = avg_minus_post;
    end

end

% Convert to matrices (animals × time)
learner_pre_CS_plus = cat(1, learner_plus_pre{:});
learner_pre_CS_minus = cat(1, learner_minus_pre{:});
learner_post_CS_plus= cat(1, learner_plus_post{:});
learner_post_CS_minus= cat(1, learner_minus_post{:});

non_learner_pre_CS_plus = cat(1, non_plus_pre{:});
non_learner_pre_CS_minus = cat(1, non_minus_pre{:});
non_learner_post_CS_plus= cat(1, non_plus_post{:});
non_learner_post_CS_minus= cat(1, non_minus_post{:});


% Compute grand means and SEMs
mean_learner_pre_CS_plus  = nanmean(learner_pre_CS_plus,1);  sem_learner_pre_CS_plus  = nanstd(learner_pre_CS_plus,0,1)/sqrt(size(learner_pre_CS_plus,1));
mean_learner_pre_CS_minus = nanmean(learner_pre_CS_minus,1); sem_learner_pre_CS_minus = nanstd(learner_pre_CS_minus,0,1)/sqrt(size(learner_pre_CS_minus,1));
mean_learner_post_CS_plus = nanmean(learner_post_CS_plus,1);  sem_learner_post_CS_plus  = nanstd(learner_post_CS_plus,0,1)/sqrt(size(learner_post_CS_plus,1));
mean_learner_post_CS_minus = nanmean(learner_post_CS_minus,1); sem_learner_post_CS_minus = nanstd(learner_post_CS_minus,0,1)/sqrt(size(learner_post_CS_minus,1));

mean_non_learner_pre_CS_plus  = nanmean(non_learner_pre_CS_plus,1);  sem_non_learner_pre_CS_plus  = nanstd(non_learner_pre_CS_plus,0,1)/sqrt(size(non_learner_pre_CS_plus,1));
mean_non_learner_pre_CS_minus = nanmean(non_learner_pre_CS_minus,1); sem_non_learner_pre_CS_minus = nanstd(non_learner_pre_CS_minus,0,1)/sqrt(size(non_learner_pre_CS_minus,1));
mean_non_learner_post_CS_plus  = nanmean(non_learner_post_CS_plus,1);  sem_non_learner_post_CS_plus  = nanstd(non_learner_post_CS_plus,0,1)/sqrt(size(non_learner_post_CS_plus,1));
mean_non_learner_post_CS_minus = nanmean(non_learner_post_CS_minus,1); sem_non_learner_post_CS_minus = nanstd(non_learner_post_CS_minus,0,1)/sqrt(size(non_learner_post_CS_minus,1));

% Create a global y-axis
y_global_max= max([mean_learner_pre_CS_plus;mean_learner_pre_CS_minus;mean_learner_post_CS_plus;mean_learner_post_CS_minus;mean_non_learner_pre_CS_plus;mean_non_learner_pre_CS_minus;...
    mean_non_learner_post_CS_plus;mean_non_learner_post_CS_minus],[],'all') *1.1;

% Finally: plot group‐averages for CS+ and CS−, pre vs post, learners vs non
figure('Color','w','Position',[100 100 1000 600]);
t = bin_centers;

% CS+ panel
subplot(2,1,1); hold on;
ax = gca;
hold(ax,'on');

% --- Learner PRE (draw only if data exist) ---
if ~isempty(mean_learner_pre_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_pre_CS_plus+sem_learner_pre_CS_plus, ...
        fliplr(mean_learner_pre_CS_plus-sem_learner_pre_CS_plus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_learner_pre_CS_plus, ...
        'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'Learner Pre',LineStyle='--');
end


% --- Learner POST ---

% ribbon (no legend)
fill(ax, [t fliplr(t)], ...
    [mean_learner_post_CS_plus+sem_learner_post_CS_plus, ...
    fliplr(mean_learner_post_CS_plus-sem_learner_post_CS_plus)], ...
    cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

% mean line (will appear in legend)
plot(ax, t, mean_learner_post_CS_plus, ...
    'Color', cLearner, 'LineWidth', 2, ...
    'DisplayName', 'Learners Post');

% --- Non-learner PRE ---%
if ~isempty(mean_non_learner_post_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_plus+sem_non_learner_pre_CS_plus, ...
        fliplr(mean_non_learner_pre_CS_plus-sem_non_learner_pre_CS_plus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_pre_CS_plus, ...
        'color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'Non-learners Pre',LineStyle='--');
end

% --- Non-learner POST (draw only if data exist) ---
if ~isempty(mean_non_learner_post_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_post_CS_plus+sem_non_learner_post_CS_plus, ...
        fliplr(mean_non_learner_post_CS_plus-sem_non_learner_post_CS_plus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_post_CS_plus, ...
        'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'Non-learners post');
end

% vertical zero line
xline(ax,0,'k--');

% build legend from whatever was actually plotted
legend(ax, 'show', 'Location', 'best');

xlabel('Time (s)');
ylabel('CS+ PSTH');
ylim([0,y_global_max]);
title(['Stim Off (Reward Times)', ' — CS+']);


% CS− panel
subplot(2,1,2); hold on;
ax = gca;
hold(ax,'on');


% --- Learner PRE (draw only if data exist) ---
if ~isempty(mean_learner_pre_CS_minus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_pre_CS_minus+sem_learner_pre_CS_minus, ...
        fliplr(mean_learner_pre_CS_minus-sem_learner_pre_CS_minus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_learner_pre_CS_minus, ...
        'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'Learner Pre',LineStyle='--');
end


% --- Learner POST ---

% ribbon (no legend)
fill(ax, [t fliplr(t)], ...
    [mean_learner_post_CS_minus+sem_learner_post_CS_minus, ...
    fliplr(mean_learner_post_CS_minus-sem_learner_post_CS_minus)], ...
    cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

% mean line (will appear in legend)
plot(ax, t, mean_learner_post_CS_minus, ...
    'Color', cLearner, 'LineWidth', 2, ...
    'DisplayName', 'Learners Post');

% --- Non-learner PRE ---
if ~isempty(mean_non_learner_post_CS_plus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_minus+sem_non_learner_pre_CS_minus, ...
        fliplr(mean_non_learner_pre_CS_minus-sem_non_learner_pre_CS_minus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_pre_CS_minus, ...
        'color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'Non-learners Pre',LineStyle='--');
end
% --- Non-learner POST (draw only if data exist) ---
if ~isempty(mean_non_learner_post_CS_minus)
    % ribbon (no legend)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_post_CS_minus+sem_non_learner_post_CS_minus, ...
        fliplr(mean_non_learner_post_CS_minus-sem_non_learner_post_CS_minus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');

    % mean line (will appear in legend)
    plot(ax, t, mean_non_learner_post_CS_minus, ...
        'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'Non-learners post');
end


% vertical zero line
xline(ax,0,'k--');
% build legend from whatever was actually plotted
% legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('CS- PSTH');
ylim([0,y_global_max]);
title(['First Lick from Stim Off', ' — CS-']);


sgtitle(['Stim in Final Position',' — Group pre/post learning']);




%% Plot mean PSTHs split by mPFC+ and mPFC- seperate panels for Pre and Post learning

% FIGURE 1: PRE-LEARNING
figure('Color','w','Position',[100 100 1400 600]);
t = bin_centers;

% Panel 1: CS+
subplot(1,2,1); hold on;
ax = gca;

% Learners CS+
if ~isempty(mean_learner_pre_CS_plus)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_pre_CS_plus+sem_learner_pre_CS_plus, ...
        fliplr(mean_learner_pre_CS_plus-sem_learner_pre_CS_plus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_learner_pre_CS_plus, 'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'Learners');
    % xline(-1,'k-','Stim Move',HandleVisibility='off') % if aligned to center times
end

% Non-learners CS+
if ~isempty(mean_non_learner_pre_CS_plus)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_plus+sem_non_learner_pre_CS_plus, ...
        fliplr(mean_non_learner_pre_CS_plus-sem_non_learner_pre_CS_plus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_non_learner_pre_CS_plus, 'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'Non-learners');
    % xline(-1,'k-','Stim Move',HandleVisibility='off') % if aligned to center times
end

xline(ax, 0, 'k--', 'HandleVisibility', 'off');
legend(ax, 'show', 'Location', 'best', 'FontSize', 11);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Lick Rate (Hz)', 'FontSize', 12);
ylim([0, y_global_max]);
title('CS+', 'FontSize', 14, 'FontWeight', 'bold');

% Panel 2: CS-
subplot(1,2,2); hold on;
ax = gca;

% Learners CS-
if ~isempty(mean_learner_pre_CS_minus)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_pre_CS_minus+sem_learner_pre_CS_minus, ...
        fliplr(mean_learner_pre_CS_minus-sem_learner_pre_CS_minus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_learner_pre_CS_minus, 'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC+');
end

% Non-learners CS-
if ~isempty(mean_non_learner_pre_CS_minus)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_pre_CS_minus+sem_non_learner_pre_CS_minus, ...
        fliplr(mean_non_learner_pre_CS_minus-sem_non_learner_pre_CS_minus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_non_learner_pre_CS_minus, 'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC-');
end

xline(ax, 0, 'k--', 'HandleVisibility', 'off');
% legend(ax, 'show', 'Location', 'best', 'FontSize', 11);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Lick Rate (Hz)', 'FontSize', 12);
ylim([0, y_global_max/2]);
title('CS-', 'FontSize', 14, 'FontWeight', 'bold');

sgtitle('Mean PSTHs Pre-Learning', 'FontSize', 16, 'FontWeight', 'bold');

% FIGURE 2: POST-LEARNING
figure('Color','w','Position',[100 100 1400 600]);

% Panel 1: CS+
subplot(1,2,1); hold on;
ax = gca;

% Learners CS+
if ~isempty(mean_learner_post_CS_plus)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_post_CS_plus+sem_learner_post_CS_plus, ...
        fliplr(mean_learner_post_CS_plus-sem_learner_post_CS_plus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_learner_post_CS_plus, 'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC+');
        % xline(-1,'k-','Stim Move',HandleVisibility='off') % if aligned to center times

end

% Non-learners CS+
if ~isempty(mean_non_learner_post_CS_plus)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_post_CS_plus+sem_non_learner_post_CS_plus, ...
        fliplr(mean_non_learner_post_CS_plus-sem_non_learner_post_CS_plus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_non_learner_post_CS_plus, 'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC-');
        % xline(-1,'k-','Stim Move',HandleVisibility='off') % if aligned to center times

end

xline(ax, 0, 'k--', 'HandleVisibility', 'off');
% legend(ax, 'show', 'Location', 'best', 'FontSize', 11);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Lick Rate (Hz)', 'FontSize', 12);
ylim([0, y_global_max]);
title('CS+', 'FontSize', 14, 'FontWeight', 'bold');

% Panel 2: CS-
subplot(1,2,2); hold on;
ax = gca;

% Learners CS-
if ~isempty(mean_learner_post_CS_minus)
    fill(ax, [t fliplr(t)], ...
        [mean_learner_post_CS_minus+sem_learner_post_CS_minus, ...
        fliplr(mean_learner_post_CS_minus-sem_learner_post_CS_minus)], ...
        cLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_learner_post_CS_minus, 'Color', cLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC+');

end

% Non-learners CS-
if ~isempty(mean_non_learner_post_CS_minus)
    fill(ax, [t fliplr(t)], ...
        [mean_non_learner_post_CS_minus+sem_non_learner_post_CS_minus, ...
        fliplr(mean_non_learner_post_CS_minus-sem_non_learner_post_CS_minus)], ...
        cNonLearner, 'FaceAlpha', .1, 'EdgeColor', 'none', 'HandleVisibility','off');
    plot(ax, t, mean_non_learner_post_CS_minus, 'Color', cNonLearner, 'LineWidth', 2, ...
        'DisplayName', 'mPFC-');

end

xline(ax, 0, 'k--', 'HandleVisibility', 'off');
% legend(ax, 'show', 'Location', 'best', 'FontSize', 11);
xlabel('Time (s)', 'FontSize', 12);
ylabel('Lick Rate (Hz)', 'FontSize', 12);
ylim([0, y_global_max/2]);
title('CS-', 'FontSize', 14, 'FontWeight', 'bold');

sgtitle('Mean PSTHs Post-Learning', 'FontSize', 16, 'FontWeight', 'bold');


%% Plots seperate PSTHs with individual data points for CS+ and CS- : for each stimulus two panels, first pre and post

target_workflow = {'visual_operant_lick_two_stim_right_move'};  % define workflow to work with

% Pre-allocate containers for all animals
all_plus_pre  = {}; all_plus_post  = {};
all_minus_pre = {}; all_minus_post = {};

% Define the fields
rew_field = 'avg_psth_rewarded_stim_on';
non_rew_field = 'avg_psth_non_rewarded_stim_on';

% Color definitions
cCS_plus = [0.2 0.2 0.8];   % Blue for CS+
cCS_minus = [0.8 0.2 0.2];  % Red for CS-

% Loop over all animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};

    days_all = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow,target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays,1,'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays)+1;
    end

    % 4) now split into pre/post
    plus_pre  = []; minus_pre  = [];
    plus_post = []; minus_post = [];

    for d = 1:numel(validDays)
        if d < ld
            plus_pre  = [plus_pre;  validDays(d).(rew_field)];
            minus_pre = [minus_pre; validDays(d).(non_rew_field)];
        else
            plus_post  = [plus_post;  validDays(d).(rew_field)];
            minus_post = [minus_post; validDays(d).(non_rew_field)];
        end
    end

    T = numel(validDays(d).(rew_field));

    % --- safe per-animal day-averages (row, NaN-safe) ---
    if isempty(plus_pre),   avg_plus_pre   = nan(1,T); else, avg_plus_pre   = reshape(mean(plus_pre,  1,'omitnan'),1,[]); end
    if isempty(plus_post),  avg_plus_post  = nan(1,T); else, avg_plus_post  = reshape(mean(plus_post, 1,'omitnan'),1,[]); end
    if isempty(minus_pre),  avg_minus_pre  = nan(1,T); else, avg_minus_pre  = reshape(mean(minus_pre, 1,'omitnan'),1,[]); end
    if isempty(minus_post), avg_minus_post = nan(1,T); else, avg_minus_post = reshape(mean(minus_post,1,'omitnan'),1,[]); end

    % --- store for all animals ---
    all_plus_pre{end+1,1}   = avg_plus_pre;
    all_plus_post{end+1,1}  = avg_plus_post;
    all_minus_pre{end+1,1}  = avg_minus_pre;
    all_minus_post{end+1,1} = avg_minus_post;
end

% Convert to matrices (animals × time)
all_pre_CS_plus = cat(1, all_plus_pre{:});
all_pre_CS_minus = cat(1, all_minus_pre{:});
all_post_CS_plus = cat(1, all_plus_post{:});
all_post_CS_minus = cat(1, all_minus_post{:});

% Compute grand means and SEMs across all animals
mean_all_pre_CS_plus  = nanmean(all_pre_CS_plus,1);
sem_all_pre_CS_plus  = nanstd(all_pre_CS_plus,0,1)/sqrt(size(all_pre_CS_plus,1));
mean_all_pre_CS_minus = nanmean(all_pre_CS_minus,1);
sem_all_pre_CS_minus = nanstd(all_pre_CS_minus,0,1)/sqrt(size(all_pre_CS_minus,1));
mean_all_post_CS_plus = nanmean(all_post_CS_plus,1);
sem_all_post_CS_plus  = nanstd(all_post_CS_plus,0,1)/sqrt(size(all_post_CS_plus,1));
mean_all_post_CS_minus = nanmean(all_post_CS_minus,1);
sem_all_post_CS_minus = nanstd(all_post_CS_minus,0,1)/sqrt(size(all_post_CS_minus,1));

t = bin_centers;
n_animals = size(all_pre_CS_plus, 1);

% ========== FIGURE 1: CS+ Trials ==========
figure('Color','w','Position',[100 100 1000 800]);

% Top panel: Pre-learning
subplot(2,1,1); hold on;
ax = gca;

% Plot individual animal trajectories (thin lines)
for i = 1:n_animals
    plot(ax, t, all_pre_CS_plus(i,:), '-', ...
        'Color', [cCS_plus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot group average (thick line)
plot(ax, t, mean_all_pre_CS_plus, '-', ...
    'Color', cCS_plus, 'LineWidth', 3, 'DisplayName', 'Average');

% Add SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_pre_CS_plus+sem_all_pre_CS_plus, ...
    fliplr(mean_all_pre_CS_plus-sem_all_pre_CS_plus)], ...
    cCS_plus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5);
legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Lick Rate (Hz)');
title('Pre-Learning');
grid on;

% Bottom panel: Post-learning
subplot(2,1,2); hold on;
ax = gca;

% Plot individual animal trajectories (thin lines)
for i = 1:n_animals
    plot(ax, t, all_post_CS_plus(i,:), '-', ...
        'Color', [cCS_plus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot group average (thick line)
plot(ax, t, mean_all_post_CS_plus, '-', ...
    'Color', cCS_plus, 'LineWidth', 3, 'DisplayName', 'Average');

% Add SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_post_CS_plus+sem_all_post_CS_plus, ...
    fliplr(mean_all_post_CS_plus-sem_all_post_CS_plus)], ...
    cCS_plus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5);
legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Lick Rate (Hz)');
title('Post-Learning');
grid on;

sgtitle(sprintf('CS+ Trials - All Animals (n=%d)', n_animals), 'FontSize', 14, 'FontWeight', 'bold');

% ========== FIGURE 2: CS- Trials ==========
figure('Color','w','Position',[150 150 1000 800]);

% Top panel: Pre-learning
subplot(2,1,1); hold on;
ax = gca;

% Plot individual animal trajectories (thin lines)
for i = 1:n_animals
    plot(ax, t, all_pre_CS_minus(i,:), '-', ...
        'Color', [cCS_minus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot group average (thick line)
plot(ax, t, mean_all_pre_CS_minus, '-', ...
    'Color', cCS_minus, 'LineWidth', 3, 'DisplayName', 'Average');

% Add SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_pre_CS_minus+sem_all_pre_CS_minus, ...
    fliplr(mean_all_pre_CS_minus-sem_all_pre_CS_minus)], ...
    cCS_minus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5);
legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Lick Rate (Hz)');
title('Pre-Learning');
grid on;

% Bottom panel: Post-learning
subplot(2,1,2); hold on;
ax = gca;

% Plot individual animal trajectories (thin lines)
for i = 1:n_animals
    plot(ax, t, all_post_CS_minus(i,:), '-', ...
        'Color', [cCS_minus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot group average (thick line)
plot(ax, t, mean_all_post_CS_minus, '-', ...
    'Color', cCS_minus, 'LineWidth', 3, 'DisplayName', 'Average');

% Add SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_post_CS_minus+sem_all_post_CS_minus, ...
    fliplr(mean_all_post_CS_minus-sem_all_post_CS_minus)], ...
    cCS_minus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5);
% legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Lick Rate (Hz)');
title('Post-Learning');
grid on;

sgtitle(sprintf('CS- Trials - All Animals (n=%d)', n_animals), 'FontSize', 14, 'FontWeight', 'bold');

%% Plots Pre-post PSTHs of CS+ and CS- overlaid with individual data points

target_workflow = {'visual_operant_lick_two_stim_static'};  % define workflow to work with

% Make a list of the animal_ids in cell format
all_ids = {behaviour_data.animal_id};


% Pre-allocate containers for all animals
all_plus_pre  = {}; all_plus_post  = {};
all_minus_pre = {}; all_minus_post = {};

% Define the fields
rew_field = 'avg_psth_rewarded_stim_final_position';
non_rew_field = 'avg_psth_non_rewarded_stim_final_position';

% Color definitions
cCS_plus = [0.2, 0.5, 0.8];   % Blue for CS+
cCS_minus = [0.9, 0.3, 0.3];  % Red for CS-

% Loop over all animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};

    days_all = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow,target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays,1,'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays)+1;
    end

    % 4) now split into pre/post
    plus_pre  = []; minus_pre  = [];
    plus_post = []; minus_post = [];

    for d = 1:numel(validDays)
        if d < ld
            plus_pre  = [plus_pre;  validDays(d).(rew_field)];
            minus_pre = [minus_pre; validDays(d).(non_rew_field)];
        else
            plus_post  = [plus_post;  validDays(d).(rew_field)];
            minus_post = [minus_post; validDays(d).(non_rew_field)];
        end
    end

    T = numel(validDays(d).(rew_field));

    % --- safe per-animal day-averages (row, NaN-safe) ---
    if isempty(plus_pre),   avg_plus_pre   = nan(1,T); else, avg_plus_pre   = reshape(mean(plus_pre,  1,'omitnan'),1,[]); end
    if isempty(plus_post),  avg_plus_post  = nan(1,T); else, avg_plus_post  = reshape(mean(plus_post, 1,'omitnan'),1,[]); end
    if isempty(minus_pre),  avg_minus_pre  = nan(1,T); else, avg_minus_pre  = reshape(mean(minus_pre, 1,'omitnan'),1,[]); end
    if isempty(minus_post), avg_minus_post = nan(1,T); else, avg_minus_post = reshape(mean(minus_post,1,'omitnan'),1,[]); end

    % --- store for all animals ---
    all_plus_pre{end+1,1}   = avg_plus_pre;
    all_plus_post{end+1,1}  = avg_plus_post;
    all_minus_pre{end+1,1}  = avg_minus_pre;
    all_minus_post{end+1,1} = avg_minus_post;
end

% Convert to matrices (animals × time)
all_pre_CS_plus = cat(1, all_plus_pre{:});
all_pre_CS_minus = cat(1, all_minus_pre{:});
all_post_CS_plus = cat(1, all_plus_post{:});
all_post_CS_minus = cat(1, all_minus_post{:});

% Compute grand means and SEMs across all animals
mean_all_pre_CS_plus  = nanmean(all_pre_CS_plus,1);
sem_all_pre_CS_plus  = nanstd(all_pre_CS_plus,0,1)/sqrt(size(all_pre_CS_plus,1));
mean_all_pre_CS_minus = nanmean(all_pre_CS_minus,1);
sem_all_pre_CS_minus = nanstd(all_pre_CS_minus,0,1)/sqrt(size(all_pre_CS_minus,1));
mean_all_post_CS_plus = nanmean(all_post_CS_plus,1);
sem_all_post_CS_plus  = nanstd(all_post_CS_plus,0,1)/sqrt(size(all_post_CS_plus,1));
mean_all_post_CS_minus = nanmean(all_post_CS_minus,1);
sem_all_post_CS_minus = nanstd(all_post_CS_minus,0,1)/sqrt(size(all_post_CS_minus,1));

t = bin_centers;
n_animals = size(all_pre_CS_plus, 1);

% Set ylim based on global max and 0
max_y_global=max([all_pre_CS_plus;all_pre_CS_minus;all_post_CS_plus;all_post_CS_minus],[],'all');
min_y_global= 0;

% ========== Single Figure: CS+ and CS- Overlaid ==========
figure('Color','w','Position',[100 100 1000 800]);

% ===== Top panel: Pre-learning (CS+ and CS- overlaid) =====
subplot(2,1,1); hold on;
ax = gca;

% Plot individual CS+ trajectories (thin blue lines)
for i = 1:n_animals
    plot(ax, t, all_pre_CS_plus(i,:), '-', ...
        'Color', [cCS_plus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot individual CS- trajectories (thin red lines)
for i = 1:n_animals
    plot(ax, t, all_pre_CS_minus(i,:), '-', ...
        'Color', [cCS_minus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot CS+ group average (thick blue line)
plot(ax, t, mean_all_pre_CS_plus, '-', ...
    'Color', cCS_plus, 'LineWidth', 3, 'DisplayName', 'CS+');

% Add CS+ SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_pre_CS_plus+sem_all_pre_CS_plus, ...
    fliplr(mean_all_pre_CS_plus-sem_all_pre_CS_plus)], ...
    cCS_plus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

% Plot CS- group average (thick red line)
plot(ax, t, mean_all_pre_CS_minus, '-', ...
    'Color', cCS_minus, 'LineWidth', 3, 'DisplayName', 'CS-');

% Add CS- SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_pre_CS_minus+sem_all_pre_CS_minus, ...
    fliplr(mean_all_pre_CS_minus-sem_all_pre_CS_minus)], ...
    cCS_minus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5,'HandleVisibility','off');
legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylim([0,max_y_global]);
ylabel('Lick Rate (Hz)');
title('Pre-Learning');
% box on;

% ===== Bottom panel: Post-learning (CS+ and CS- overlaid) =====
subplot(2,1,2); hold on;
ax = gca;

% Plot individual CS+ trajectories (thin blue lines)
for i = 1:n_animals
    plot(ax, t, all_post_CS_plus(i,:), '-', ...
        'Color', [cCS_plus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot individual CS- trajectories (thin red lines)
for i = 1:n_animals
    plot(ax, t, all_post_CS_minus(i,:), '-', ...
        'Color', [cCS_minus 0.3], 'LineWidth', 0.5, 'HandleVisibility', 'off');
end

% Plot CS+ group average (thick blue line)
plot(ax, t, mean_all_post_CS_plus, '-', ...
    'Color', cCS_plus, 'LineWidth', 3, 'DisplayName', 'CS+ ');

% Add CS+ SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_post_CS_plus+sem_all_post_CS_plus, ...
    fliplr(mean_all_post_CS_plus-sem_all_post_CS_plus)], ...
    cCS_plus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

% Plot CS- group average (thick red line)
plot(ax, t, mean_all_post_CS_minus, '-', ...
    'Color', cCS_minus, 'LineWidth', 3, 'DisplayName', 'CS-');

% Add CS- SEM shading
fill(ax, [t fliplr(t)], ...
    [mean_all_post_CS_minus+sem_all_post_CS_minus, ...
    fliplr(mean_all_post_CS_minus-sem_all_post_CS_minus)], ...
    cCS_minus, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility','off');

xline(ax, 0, 'k--', 'LineWidth', 1.5);
legend(ax, 'show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Lick Rate (Hz)');
ylim([0,max_y_global]);
title('Post-Learning');
% box on;


sgtitle(sprintf('CS+ and CS- Trials - Pre vs Post Learning (n=%d animals)', n_animals), 'FontSize', 14, 'FontWeight', 'bold');


%% Plots a heatmap of anticipaotry licks aligned to event for a single animal


% Plot heatmap of anticipatory licks for a single animal across days
% Set parameters
target_animal = 'HA012'; % Define the animal to plot
target_workflow = 'visual_operant_lick_two_stim_right_move'; % Define workflow

% Find the animal index
all_ids = {behaviour_data.animal_id};
animal_idx = find(strcmp(all_ids, target_animal));

if isempty(animal_idx)
    error('Animal %s not found in behaviour_data', target_animal);
end

% Get recording days for this animal
days_all = behaviour_data(animal_idx).recording_day;

% Select only the days matching the workflow
isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
validDays = days_all(isValid);

if isempty(validDays)
    error('Animal %s: no recording days found for workflow %s', target_animal, target_workflow);
end

% Extract PSTH data for each day
n_days = numel(validDays);
n_bins = numel(bin_centers);
psth_matrix = nan(n_days, n_bins);

for d = 1:n_days
    day_data = validDays(d);

    % Get the PSTH for rewarded stim (CS+)
    if isfield(day_data, 'avg_psth_rewarded_stim_final_position')
        psth_data = day_data.avg_psth_non_rewarded_stim_final_position;

        % Ensure correct dimensions
        if length(psth_data) == n_bins
            psth_matrix(d, :) = psth_data;
        else
            warning('Day %d: PSTH length mismatch (%d vs %d expected)', d, length(psth_data), n_bins);
        end
    else
        warning('Day %d: avg_psth_rewarded_stim_final_position field not found', d);
    end
end

% Normalize PSTH data
% Option 1: Normalize each day by its own maximum (better for seeing within-day patterns)
psth_matrix_norm = nan(size(psth_matrix));
for d = 1:n_days
    day_max = max(psth_matrix(d, :), [], 'omitnan');
    if day_max > 0
        psth_matrix_norm(d, :) = psth_matrix(d, :) / day_max;
    else
        psth_matrix_norm(d, :) = psth_matrix(d, :);
    end
end

% Option 2: Z-score normalization per day (uncomment to use)
% psth_matrix_norm = nan(size(psth_matrix));
% for d = 1:n_days
%     day_mean = mean(psth_matrix(d, :), 'omitnan');
%     day_std = std(psth_matrix(d, :), 'omitnan');
%     if day_std > 0
%         psth_matrix_norm(d, :) = (psth_matrix(d, :) - day_mean) / day_std;
%     else
%         psth_matrix_norm(d, :) = psth_matrix(d, :) - day_mean;
%     end
% end

% Option 3: Normalize by global maximum (previous approach, uncomment to use)
% global_max = max(psth_matrix(:), [], 'omitnan');
% if global_max > 0
%     psth_matrix_norm = psth_matrix / global_max;
% else
%     psth_matrix_norm = psth_matrix;
% end

% Create figure
figure('Color', 'w', 'Name', sprintf('%s - PSTH Heatmap - %s', target_animal, target_workflow));

% Plot heatmap with adjusted color limits
imagesc(bin_centers, 1:n_days, psth_matrix_norm);
colorbar;
colormap('parula'); % You can change colormap (e.g., 'jet', 'parula', 'hot', 'turbo')

% Adjust color scale based on anticipatory period only (before stim center)
% Find bins before time = 0 (anticipatory period)
anticipatory_bins = bin_centers < 0;
anticipatory_values = psth_matrix_norm(:, anticipatory_bins);
anticipatory_values = anticipatory_values(~isnan(anticipatory_values));

% Set color limits based on anticipatory period statistics
if ~isempty(anticipatory_values)
    anticipatory_max = max(anticipatory_values);
    anticipatory_95th = prctile(anticipatory_values, 95); % 95th percentile

    % Use 95th percentile or max of anticipatory period (whichever gives better contrast)
    clim([0 anticipatory_95th * 3.1]); % Scale up slightly to avoid saturation at peak
else
    clim([0 0.5]); % Fallback if no anticipatory data
end

% Labels and title
xlabel('Time Relative to Stim Center (s)');
ylabel('Recording Day');
title(sprintf('%s - Normalized Anticipatory Licks Across Days\n%s', target_animal, strrep(target_workflow, '_', ' ')), ...
    'Interpreter', 'none');

% Add vertical line at stim center (time = 0)
hold on;
plot([0 0], [0.5 n_days+0.5], 'w--', 'LineWidth', 2);
hold off;

% Set y-axis to show day numbers (latest days closer to x-axis)
yticks(1:n_days);
yticklabels(1:n_days);

% Improve appearance
set(gca, 'YDir', 'reverse'); % Latest day at bottom, Day 1 at top
set(gca, 'FontSize', 10);
axis tight;


%% Plots a defined behavioural metric seperated by n components for within session analysis- returns 2 figures : CS+ vs CS- overlaid and normalized plot underneath

% Analyze anticipatory licks aligned to learning day for specific protocol - CS+/CS- separated
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move'; % define workflow to work with
n_components = 4; % Number of components to split each session into

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA007','HA006','HA009','HA013','HA014','HA015'};

% Make a list of the animal_ids in cell format
all_ids = {behaviour_data.animal_id};

% Initialize storage - including CS type
component_data = struct('Animal', {}, 'Day', {}, 'ComponentMeans', {}, 'RelativeLearningDay', {}, 'CSType', {});

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    % Skip non-learners
    if ~isLearner
        continue;
    end

    % Get recording days for this animal
    days_all = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays) + 1; % define post-learning outside the scope
    end

    % Process each valid day for this animal
    for d = 1:numel(validDays)
        day_data = validDays(d);

        % Get trial labels (logical)
        cs_labels = day_data.cs_labels;

        % Define the metric you want to plot
        data_metric = day_data.anticipatory_licks;

        if length(cs_labels) ~= length(data_metric)
            warning('CS labels length mismatch with data metric for animal %s day %d', animal_id, d);
            continue;
        end

        % Create trial type masks (logical arrays)
        cs_plus_mask = cs_labels;      % true for CS+, false for CS-
        cs_minus_mask = ~cs_labels;    % false for CS+, true for CS-

        % Process CS+ and CS- separately
        cs_types = [1, 0]; % 1 for CS+, 0 for CS-
        for cs_idx = 1:length(cs_types)
            cs_type_val = cs_types(cs_idx);

            if cs_type_val == 1
                trial_mask = cs_plus_mask;
                cs_type_name = 'CS+';
            else
                trial_mask = cs_minus_mask;
                cs_type_name = 'CS-';
            end

            if sum(trial_mask) < n_components
                continue % Skip if too few trials for this CS type
            end

            % Get data for this CS type
            cs_data = data_metric(trial_mask);

            % Split trials into n components by index
            nT = length(cs_data);
            edges = round(linspace(1, nT+1, n_components+1));
            component_means = nan(1, n_components);

            for c = 1:n_components
                idx = edges(c):edges(c+1)-1;
                component_means(c) = median(cs_data(idx), 'omitnan');
            end

            % Calculate relative learning day
            relative_learning_day = d - ld;

            component_data(end+1) = struct( ...
                'Animal', string(animal_id), ...
                'Day', d, ...
                'ComponentMeans', component_means, ...
                'RelativeLearningDay', relative_learning_day, ...
                'CSType', string(cs_type_name));
        end
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Compute global y-limits across both CS types
all_means = vertcat(component_data.ComponentMeans);
all_means = all_means(~isnan(all_means));
if isempty(all_means)
    return;
end
ymin = max(0, min(all_means));
ymax = max(all_means);

% Get unique CS types
cs_types = unique(string({component_data.CSType})', 'stable');

% Prepare data for both CS types
cs_plus_data = struct();
cs_minus_data = struct();

for cs_i = 1:numel(cs_types)
    cs_type = cs_types(cs_i);

    % Select data for this CS type
    cs_mask = strcmp(string({component_data.CSType})', cs_type);
    relative_days = [component_data(cs_mask).RelativeLearningDay]';
    means = vertcat(component_data(cs_mask).ComponentMeans);

    % Count how many animals per relative learning day and filter days with <3 animals
    unique_rel_days = unique(relative_days);
    day_counts = arrayfun(@(d) sum(relative_days == d), unique_rel_days);

    % Keep only days with 3+ animals
    valid_day_mask = day_counts >= 3;
    valid_rel_days = unique_rel_days(valid_day_mask);

    if isempty(valid_rel_days)
        warning('No relative learning days with >=3 animals for %s', cs_type);
        continue;
    end

    % Average across animals for valid relative learning days only
    avg_means = nan(numel(valid_rel_days), n_components);
    for di = 1:numel(valid_rel_days)
        day_mask = relative_days == valid_rel_days(di);
        avg_means(di, :) = mean(means(day_mask, :), 1, 'omitnan');
    end

    % Sort by relative learning day
    [sorted_rel_days, sort_idx] = sort(valid_rel_days);
    sorted_means = avg_means(sort_idx, :);

    % Store data
    if strcmp(cs_type, 'CS+')
        cs_plus_data.days = sorted_rel_days;
        cs_plus_data.means = sorted_means;
    else
        cs_minus_data.days = sorted_rel_days;
        cs_minus_data.means = sorted_means;
    end
end

% Create figure with two panels
figure('Color','w','Name',sprintf('Learning-Aligned Anticipatory Licks (n=%d) - %s', n_components, target_workflow));
tiledlayout(2, 1, 'TileSpacing','compact');

% Add main title for the entire figure
sgtitle(sprintf('Anticipatory Licks - %s Protocol', strrep(target_workflow, '_', ' ')), 'FontSize', 14, 'FontWeight', 'bold');

% Panel 1: Overlayed CS+ and CS- lines
nexttile;
hold on;
component_spacing = 0.15; % Horizontal spacing between components within a day

% Plot CS+
if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_plus = plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.2 0.8]);

        % Only add to legend for first iteration
        if di == 1
            h_plus.DisplayName = 'CS+';
        else
            h_plus.HandleVisibility = 'off';
        end
    end
end

% Plot CS-
if isfield(cs_minus_data, 'days')
    for di = 1:numel(cs_minus_data.days)
        day_x = cs_minus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_minus = plot(component_x_positions, cs_minus_data.means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.8 0.2 0.2]);

        % Only add to legend for first iteration
        if di == 1
            h_minus.DisplayName = 'CS-';
        else
            h_minus.HandleVisibility = 'off';
        end
    end
end

hold off;
xlabel('Days Relative to Learning');
ylabel('Mean Anticipatory Licks');
ylim([ymin ymax]);
title('CS+ and CS- Trials Overlay');
legend('Location', 'best');
xline(0, '--k', 'Learning Day', 'LineWidth', 2);
grid on;

% Panel 2: Discrimination index (CS+ - CS-) / (CS+ + CS-)
nexttile;

% Find common days between CS+ and CS-
if isfield(cs_plus_data, 'days') && isfield(cs_minus_data, 'days')
    common_days = intersect(cs_plus_data.days, cs_minus_data.days);

    if ~isempty(common_days)
        discrimination_index = nan(numel(common_days), n_components);

        for di = 1:numel(common_days)
            day = common_days(di);

            % Find indices in CS+ and CS- data
            idx_plus = find(cs_plus_data.days == day);
            idx_minus = find(cs_minus_data.days == day);

            cs_plus_vals = cs_plus_data.means(idx_plus, :);
            cs_minus_vals = cs_minus_data.means(idx_minus, :);

            % Calculate discrimination index: (CS+ - CS-) / (CS+ + CS-)
            discrimination_index(di, :) = (cs_plus_vals - cs_minus_vals) ./ (cs_plus_vals + cs_minus_vals);
        end

        % Plot discrimination index
        hold on;
        for di = 1:numel(common_days)
            day_x = common_days(di);
            component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
                component_spacing*(n_components-1)/2, n_components);

            plot(component_x_positions, discrimination_index(di, :), '-o', ...
                'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.4 0.4 0.4]);
        end
        hold off;

        xlabel('Days Relative to Learning');
        ylabel('Discrimination Index');
        title('(CS+ - CS-) / (CS+ + CS-)');
        xline(0, '--k', 'Learning Day', 'LineWidth', 2);
        yline(0, '--', 'No Discrimination', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
        grid on;
    else
        text(0.5, 0.5, 'No common days between CS+ and CS-', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
else
    text(0.5, 0.5, 'Insufficient data for discrimination index', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
end

%% Plots a defined behavioural metric seperated by n components for within session analysis- returns 3 figures : CS+ vs CS- overlaid , log scale for CS+ and normalized plot underneath

% Analyze anticipatory licks aligned to learning day for specific protocol - CS+/CS- separated
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move'; % define workflow to work with
n_components = 4; % Number of components to split each session into

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

% Make a list of the animal_ids in cell format
all_ids = {behaviour_data.animal_id};

% Initialize storage - including CS type
component_data = struct('Animal', {}, 'Day', {}, 'ComponentMeans', {}, 'RelativeLearningDay', {}, 'CSType', {});

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};

    isLearner = ismember(animal_id, learner_ids);

    % Skip non-learners
    % if ~isLearner
    %     continue;
    % end

    % Get recording days for this animal
    days_all = behaviour_data(ai).recording_day;

    % 1) select only the days matching the workflow
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    % get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % 3) find the first significant day
    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays) + 1; % define post-learning outside the scope
    end

    % Process each valid day for this animal
    for d = 1:numel(validDays)
        day_data = validDays(d);

        % Get trial labels (logical)
        cs_labels = day_data.cs_labels;

        % Define the metric you want to plot
        data_metric = day_data.anticipatory_licks;

        if length(cs_labels) ~= length(data_metric)
            warning('CS labels length mismatch with data metric for animal %s day %d', animal_id, d);
            continue;
        end

        % Create trial type masks (logical arrays)
        cs_plus_mask = cs_labels;      % true for CS+, false for CS-
        cs_minus_mask = ~cs_labels;    % false for CS+, true for CS-

        % Process CS+ and CS- separately
        cs_types = [1, 0]; % 1 for CS+, 0 for CS-
        for cs_idx = 1:length(cs_types)
            cs_type_val = cs_types(cs_idx);

            if cs_type_val == 1
                trial_mask = cs_plus_mask;
                cs_type_name = 'CS+';
            else
                trial_mask = cs_minus_mask;
                cs_type_name = 'CS-';
            end

            if sum(trial_mask) < n_components
                continue % Skip if too few trials for this CS type
            end

            % Get data for this CS type
            cs_data = data_metric(trial_mask);

            % Split trials into n components by index
            nT = length(cs_data);
            edges = round(linspace(1, nT+1, n_components+1));
            component_means = nan(1, n_components);

            for c = 1:n_components
                idx = edges(c):edges(c+1)-1;
                component_means(c) = median(cs_data(idx), 'omitnan');
            end

            % Calculate relative learning day
            relative_learning_day = d - ld;

            component_data(end+1) = struct( ...
                'Animal', string(animal_id), ...
                'Day', d, ...
                'ComponentMeans', component_means, ...
                'RelativeLearningDay', relative_learning_day, ...
                'CSType', string(cs_type_name));
        end
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Compute global y-limits across both CS types
all_means = vertcat(component_data.ComponentMeans);
all_means = all_means(~isnan(all_means));
if isempty(all_means)
    return;
end
ymin = max(0, min(all_means));
ymax = max(all_means);

% Get unique CS types
cs_types = unique(string({component_data.CSType})', 'stable');

% Prepare data for both CS types
cs_plus_data = struct();
cs_minus_data = struct();

for cs_i = 1:numel(cs_types)
    cs_type = cs_types(cs_i);

    % Select data for this CS type
    cs_mask = strcmp(string({component_data.CSType})', cs_type);
    relative_days = [component_data(cs_mask).RelativeLearningDay]';
    means = vertcat(component_data(cs_mask).ComponentMeans);

    % Count how many animals per relative learning day and filter days with <3 animals
    unique_rel_days = unique(relative_days);
    day_counts = arrayfun(@(d) sum(relative_days == d), unique_rel_days);

    % Keep only days with 3+ animals
    valid_day_mask = day_counts >= 3;
    valid_rel_days = unique_rel_days(valid_day_mask);

    if isempty(valid_rel_days)
        warning('No relative learning days with >=3 animals for %s', cs_type);
        continue;
    end

    % Average across animals for valid relative learning days only
    avg_means = nan(numel(valid_rel_days), n_components);
    for di = 1:numel(valid_rel_days)
        day_mask = relative_days == valid_rel_days(di);
        avg_means(di, :) = mean(means(day_mask, :), 1, 'omitnan');
    end

    % Sort by relative learning day
    [sorted_rel_days, sort_idx] = sort(valid_rel_days);
    sorted_means = avg_means(sort_idx, :);

    % Store data
    if strcmp(cs_type, 'CS+')
        cs_plus_data.days = sorted_rel_days;
        cs_plus_data.means = sorted_means;
    else
        cs_minus_data.days = sorted_rel_days;
        cs_minus_data.means = sorted_means;
    end
end

% Create figure with three panels
figure('Color','w','Name',sprintf('Learning-Aligned Reaction Time (n=%d) - %s', n_components, target_workflow));
tiledlayout(3, 1, 'TileSpacing','compact');

% Add main title for the entire figure
% sgtitle(sprintf('Reaction Time (Diff from Optimal) - %s Protocol', strrep(target_workflow, '_', ' ')), 'FontSize', 14, 'FontWeight', 'bold');

% Panel 1: Overlayed CS+ and CS- lines
nexttile;
hold on;
component_spacing = 0.15; % Horizontal spacing between components within a day

% Plot CS+
if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_plus = plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'color', cCS_plus);

        % Only add to legend for first iteration
        if di == 1
            h_plus.DisplayName = 'CS+';
        else
            h_plus.HandleVisibility = 'off';
        end
    end
end

% Plot CS-
if isfield(cs_minus_data, 'days')
    for di = 1:numel(cs_minus_data.days)
        day_x = cs_minus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_minus = plot(component_x_positions, cs_minus_data.means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'color', cCS_minus);

        % Only add to legend for first iteration
        if di == 1
            h_minus.DisplayName = 'CS-';
        else
            h_minus.HandleVisibility = 'off';
        end
    end
end

hold off;
xlabel('Days Relative to Learning');
ylabel('Mean Anticipatory Licks');
ylim([ymin ymax]);
xlim([-3,3]);
title('CS+ and CS- Trials Overlay');
legend('Location', 'best');
xline(0, '--k', 'LineWidth', 2);

% Panel 2: CS+ only with log scale for reaction time clarity
nexttile;
hold on;

if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'color', cCS_plus);
    end
end
hold off;

set(gca, 'YScale', 'log'); % Log scale for better visualization
xlabel('Days Relative to Learning');
ylabel('Median Reaction Time (s, log scale)');
title('CS+ Trials Only (Log Scale)');
xline(0, '--k', 'Learning Day', 'LineWidth', 2);

% Panel 3: Discrimination index (CS+ - CS-) / (CS+ + CS-)
nexttile;

% Find common days between CS+ and CS-
if isfield(cs_plus_data, 'days') && isfield(cs_minus_data, 'days')
    common_days = intersect(cs_plus_data.days, cs_minus_data.days);

    if ~isempty(common_days)
        discrimination_index = nan(numel(common_days), n_components);

        for di = 1:numel(common_days)
            day = common_days(di);

            % Find indices in CS+ and CS- data
            idx_plus = find(cs_plus_data.days == day);
            idx_minus = find(cs_minus_data.days == day);

            cs_plus_vals = cs_plus_data.means(idx_plus, :);
            cs_minus_vals = cs_minus_data.means(idx_minus, :);

            % Calculate discrimination index: (CS+ - CS-) / (CS+ + CS-)
            discrimination_index(di, :) = (cs_plus_vals - cs_minus_vals) ./ (cs_plus_vals + cs_minus_vals);
        end

        % Plot discrimination index
        hold on;
        for di = 1:numel(common_days)
            day_x = common_days(di);
            component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
                component_spacing*(n_components-1)/2, n_components);

            plot(component_x_positions, discrimination_index(di, :), '-o', ...
                'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.4 0.4 0.4]);
        end
        hold off;

        xlabel('Days Relative to Learning');
        ylabel('Discrimination Index');
        title('(CS+ - CS-) / (CS+ + CS-)');
        xline(0, '--k', 'Learning Day', 'LineWidth', 2);
        yline(0, '--', 'No Discrimination', 'LineWidth', 1, 'Color', [0.5 0.5 0.5]);
    else
        text(0.5, 0.5, 'No common days between CS+ and CS-', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
else
    text(0.5, 0.5, 'Insufficient data for discrimination index', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
end

%% Seperate behaviour metric by n and plot seperated plots for CS+ and CS- (Good for RT)

% Analyze anticipatory licks aligned to learning day for specific protocol - CS+/CS- separated
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move';
n_components = 4;

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};

% Initialize storage
component_data = struct('Animal', {}, 'Day', {}, 'ComponentMeans', {}, 'RelativeLearningDay', {}, 'CSType', {});

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    if isLearner
        continue;
    end

    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays) + 1;
    end

    for d = 1:numel(validDays)
        day_data = validDays(d);
        cs_labels = day_data.cs_labels;
        data_metric = day_data.all_stim_diff_from_optimal_reward;

        if length(cs_labels) ~= length(data_metric)
            warning('CS labels length mismatch with data metric for animal %s day %d', animal_id, d);
            continue;
        end

        cs_plus_mask = cs_labels;
        cs_minus_mask = ~cs_labels;

        cs_types = [1, 0];
        for cs_idx = 1:length(cs_types)
            cs_type_val = cs_types(cs_idx);

            if cs_type_val == 1
                trial_mask = cs_plus_mask;
                cs_type_name = 'CS+';
            else
                trial_mask = cs_minus_mask;
                cs_type_name = 'CS-';
            end

            if sum(trial_mask) < n_components
                continue;
            end

            cs_data = data_metric(trial_mask);
            nT = length(cs_data);
            edges = round(linspace(1, nT+1, n_components+1));
            component_means = nan(1, n_components);

            for c = 1:n_components
                idx = edges(c):edges(c+1)-1;
                component_means(c) = median(cs_data(idx), 'omitnan');
            end

            relative_learning_day = d - ld;

            component_data(end+1) = struct( ...
                'Animal', string(animal_id), ...
                'Day', d, ...
                'ComponentMeans', component_means, ...
                'RelativeLearningDay', relative_learning_day, ...
                'CSType', string(cs_type_name));
        end
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Get unique CS types
cs_types = unique(string({component_data.CSType})', 'stable');

% Prepare data for both CS types
cs_plus_data = struct();
cs_minus_data = struct();

for cs_i = 1:numel(cs_types)
    cs_type = cs_types(cs_i);

    cs_mask = strcmp(string({component_data.CSType})', cs_type);
    relative_days = [component_data(cs_mask).RelativeLearningDay]';
    means = vertcat(component_data(cs_mask).ComponentMeans);

    unique_rel_days = unique(relative_days);
    day_counts = arrayfun(@(d) sum(relative_days == d), unique_rel_days);

    valid_day_mask = day_counts >= 3;
    valid_rel_days = unique_rel_days(valid_day_mask);

    if isempty(valid_rel_days)
        warning('No relative learning days with >=3 animals for %s', cs_type);
        continue;
    end

    avg_means = nan(numel(valid_rel_days), n_components);
    for di = 1:numel(valid_rel_days)
        day_mask = relative_days == valid_rel_days(di);
        avg_means(di, :) = mean(means(day_mask, :), 1, 'omitnan');
    end

    [sorted_rel_days, sort_idx] = sort(valid_rel_days);
    sorted_means = avg_means(sort_idx, :);

    if strcmp(cs_type, 'CS+')
        cs_plus_data.days = sorted_rel_days;
        cs_plus_data.means = sorted_means;
    else
        cs_minus_data.days = sorted_rel_days;
        cs_minus_data.means = sorted_means;
    end
end

% Create figure with four panels in 2x2 layout
figure('Color','w','Position', [100, 100, 1400, 900], ...
    'Name',sprintf('Learning-Aligned Reaction Time (n=%d) - %s', n_components, target_workflow));
component_spacing = 0.15;

% Overall title
sgtitle(sprintf('Reaction Time (Diff from Optimal) - %s Protocol', strrep(target_workflow, '_', ' ')), ...
    'FontSize', 14, 'FontWeight', 'bold');

% --- PANEL 1: CS+ Linear Scale (Top Left) ---
subplot(2, 2, 1);
hold on;

if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', cCS_plus, 'MarkerFaceColor', cCS_plus);
    end
end
hold off;

xlabel('Days Relative to Learning', 'FontSize', 10);
ylabel('Reaction Time (s)', 'FontSize', 10, 'Color', cCS_plus);
title('CS+ Trials (Linear Scale)', 'FontSize', 11, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
ax1 = gca;
ax1.YColor = cCS_plus;
box on; grid on;

% --- PANEL 2: CS- Linear Scale (Top Right) ---
subplot(2, 2, 2);
hold on;

if isfield(cs_minus_data, 'days')
    for di = 1:numel(cs_minus_data.days)
        day_x = cs_minus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_minus_data.means(di, :), '-o', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', cCS_minus, 'MarkerFaceColor', cCS_minus);
    end
end
hold off;

xlabel('Days Relative to Learning', 'FontSize', 10);
ylabel('Reaction Time (s)', 'FontSize', 10, 'Color', cCS_minus);
title('CS- Trials (Linear Scale)', 'FontSize', 11, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
ax2 = gca;
ax2.YColor = cCS_minus;
box on; grid on;

% --- PANEL 3: CS+ Log Scale (Bottom Left) ---
subplot(2, 2, 3);
hold on;

if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', cCS_plus, 'MarkerFaceColor', cCS_plus);
    end
end
hold off;

set(gca, 'YScale', 'log');
xlabel('Days Relative to Learning', 'FontSize', 10);
ylabel('Reaction Time (s, log scale)', 'FontSize', 10, 'Color', cCS_plus);
title('CS+ Trials (Log Scale)', 'FontSize', 11, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
ax3 = gca;
ax3.YColor = cCS_plus;
box on; grid on;

% --- PANEL 4: Discrimination Index (Bottom Right) ---
subplot(2, 2, 4);

if isfield(cs_plus_data, 'days') && isfield(cs_minus_data, 'days')
    common_days = intersect(cs_plus_data.days, cs_minus_data.days);

    if ~isempty(common_days)
        discrimination_index = nan(numel(common_days), n_components);

        for di = 1:numel(common_days)
            day = common_days(di);

            idx_plus = find(cs_plus_data.days == day);
            idx_minus = find(cs_minus_data.days == day);

            cs_plus_vals = cs_plus_data.means(idx_plus, :);
            cs_minus_vals = cs_minus_data.means(idx_minus, :);

            % Discrimination: (CS- - CS+) / (CS- + CS+)
            % Positive values = CS- slower (longer RT) than CS+, indicating learning
            discrimination_index(di, :) = (cs_minus_vals - cs_plus_vals) ./ (cs_minus_vals + cs_plus_vals);
        end

        hold on;
        for di = 1:numel(common_days)
            day_x = common_days(di);
            component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
                component_spacing*(n_components-1)/2, n_components);

            plot(component_x_positions, discrimination_index(di, :), '-o', ...
                'LineWidth', 2, 'MarkerSize', 7, 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.5 0.5 0.5]);
        end
        hold off;

        xlabel('Days Relative to Learning', 'FontSize', 10);
        ylabel('Discrimination Index', 'FontSize', 10);
        title('(CS- - CS+) / (CS- + CS+)', 'FontSize', 11, 'FontWeight', 'bold');
        xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
        yline(0, '--', 'No Discrimination', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5]);
        box on; grid on;
    else
        text(0.5, 0.5, 'No common days between CS+ and CS-', ...
            'HorizontalAlignment', 'center', 'Units', 'normalized');
    end
else
    text(0.5, 0.5, 'Insufficient data for discrimination index', ...
        'HorizontalAlignment', 'center', 'Units', 'normalized');
end

%% Plots CS+ and CS- RTs on a log scale overlaid
% Analyze anticipatory licks aligned to learning day for specific protocol - CS+/CS- separated
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move';
n_components = 4;

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};

% Initialize storage
component_data = struct('Animal', {}, 'Day', {}, 'ComponentMeans', {}, 'RelativeLearningDay', {}, 'CSType', {});

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    if isLearner
        continue;
    end

    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', behaviour_data(ai).animal_id);
        ld = numel(validDays) + 1;
    end

    for d = 1:numel(validDays)
        day_data = validDays(d);
        cs_labels = day_data.cs_labels;
        data_metric = day_data.all_stim_diff_from_optimal_reward;

        if length(cs_labels) ~= length(data_metric)
            warning('CS labels length mismatch with data metric for animal %s day %d', animal_id, d);
            continue;
        end

        cs_plus_mask = cs_labels;
        cs_minus_mask = ~cs_labels;

        cs_types = [1, 0];
        for cs_idx = 1:length(cs_types)
            cs_type_val = cs_types(cs_idx);

            if cs_type_val == 1
                trial_mask = cs_plus_mask;
                cs_type_name = 'CS+';
            else
                trial_mask = cs_minus_mask;
                cs_type_name = 'CS-';
            end

            if sum(trial_mask) < n_components
                continue;
            end

            cs_data = data_metric(trial_mask);
            nT = length(cs_data);
            edges = round(linspace(1, nT+1, n_components+1));
            component_means = nan(1, n_components);

            for c = 1:n_components
                idx = edges(c):edges(c+1)-1;
                component_means(c) = median(cs_data(idx), 'omitnan');
            end

            relative_learning_day = d - ld;

            component_data(end+1) = struct( ...
                'Animal', string(animal_id), ...
                'Day', d, ...
                'ComponentMeans', component_means, ...
                'RelativeLearningDay', relative_learning_day, ...
                'CSType', string(cs_type_name));
        end
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Get unique CS types
cs_types = unique(string({component_data.CSType})', 'stable');

% Prepare data for both CS types
cs_plus_data = struct();
cs_minus_data = struct();

for cs_i = 1:numel(cs_types)
    cs_type = cs_types(cs_i);

    cs_mask = strcmp(string({component_data.CSType})', cs_type);
    relative_days = [component_data(cs_mask).RelativeLearningDay]';
    means = vertcat(component_data(cs_mask).ComponentMeans);

    unique_rel_days = unique(relative_days);
    day_counts = arrayfun(@(d) sum(relative_days == d), unique_rel_days);

    valid_day_mask = day_counts >= 3;
    valid_rel_days = unique_rel_days(valid_day_mask);

    if isempty(valid_rel_days)
        warning('No relative learning days with >=3 animals for %s', cs_type);
        continue;
    end

    avg_means = nan(numel(valid_rel_days), n_components);
    for di = 1:numel(valid_rel_days)
        day_mask = relative_days == valid_rel_days(di);
        avg_means(di, :) = mean(means(day_mask, :), 1, 'omitnan');
    end

    [sorted_rel_days, sort_idx] = sort(valid_rel_days);
    sorted_means = avg_means(sort_idx, :);

    if strcmp(cs_type, 'CS+')
        cs_plus_data.days = sorted_rel_days;
        cs_plus_data.means = sorted_means;
    else
        cs_minus_data.days = sorted_rel_days;
        cs_minus_data.means = sorted_means;
    end
end

% Create single overlaid plot on log scale
figure('Color','w','Position', [100, 100, 800, 600]);
component_spacing = 0.15;
hold on;

% Plot CS+
if isfield(cs_plus_data, 'days')
    for di = 1:numel(cs_plus_data.days)
        day_x = cs_plus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_plus_data.means(di, :), '-o', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', cCS_plus, 'MarkerFaceColor', cCS_plus);
    end
end

% Plot CS-
if isfield(cs_minus_data, 'days')
    for di = 1:numel(cs_minus_data.days)
        day_x = cs_minus_data.days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, cs_minus_data.means(di, :), '-o', ...
            'LineWidth', 2, 'MarkerSize', 7, 'Color', cCS_minus, 'MarkerFaceColor', cCS_minus);
    end
end

hold off;

set(gca, 'YScale', 'log');
xlabel('Days Relative to Learning', 'FontSize', 12);
ylabel('Reaction Time (s, log scale)', 'FontSize', 12);
title(sprintf('CS+ vs CS- Reaction Time - %s', strrep(target_workflow, '_', ' ')), ...
    'FontSize', 13, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend({'CS+', 'CS-'}, 'Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11);

%% Plots learners vs non-learners overlaid of a behaviour metric split by n components

% Analyze anticipatory licks aligned to learning day - CS+ separated by learner status
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move';
n_components = 5;

learner_ids = {'HA005','HA008','HA010','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};

% Initialize storage - now includes IsLearner
component_data = struct('Animal', {}, 'Day', {}, 'ComponentMeans', {}, ...
    'RelativeLearningDay', {}, 'CSType', {}, 'IsLearner', {});

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording was found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, skipping', behaviour_data(ai).animal_id);
        continue;
    end

    for d = 1:numel(validDays)
        day_data = validDays(d);
        cs_labels = day_data.cs_labels;
        data_metric = day_data.ITI_lick_counts;

        if length(cs_labels) ~= length(data_metric)
            warning('CS labels length mismatch with data metric for animal %s day %d', animal_id, d);
            continue;
        end

        cs_plus_mask = cs_labels;
        cs_minus_mask = ~cs_labels;

        cs_types = [1, 0];
        for cs_idx = 1:length(cs_types)
            cs_type_val = cs_types(cs_idx);

            if cs_type_val == 1
                trial_mask = cs_plus_mask;
                cs_type_name = 'CS+';
            else
                trial_mask = cs_minus_mask;
                cs_type_name = 'CS-';
            end

            if sum(trial_mask) < n_components
                continue;
            end

            cs_data = data_metric(trial_mask);
            nT = length(cs_data);
            edges = round(linspace(1, nT+1, n_components+1));
            component_means = nan(1, n_components);

            for c = 1:n_components
                idx = edges(c):edges(c+1)-1;
                component_means(c) = median(cs_data(idx), 'omitnan');
            end

            relative_learning_day = d - ld;

            component_data(end+1) = struct( ...
                'Animal', string(animal_id), ...
                'Day', d, ...
                'ComponentMeans', component_means, ...
                'RelativeLearningDay', relative_learning_day, ...
                'CSType', string(cs_type_name), ...
                'IsLearner', isLearner);
        end
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Filter for CS+ trials only
cs_plus_mask = strcmp(string({component_data.CSType})', 'CS+');
cs_plus_component_data = component_data(cs_plus_mask);

% Separate learners and non-learners
learner_mask = [cs_plus_component_data.IsLearner]';
learner_component_data = cs_plus_component_data(learner_mask);
non_learner_component_data = cs_plus_component_data(~learner_mask);

% Aggregate data for learners
learner_rel_days = [learner_component_data.RelativeLearningDay]';
learner_means = vertcat(learner_component_data.ComponentMeans);

unique_learner_days = unique(learner_rel_days);
day_counts_learner = arrayfun(@(d) sum(learner_rel_days == d), unique_learner_days);

valid_day_mask_learner = day_counts_learner >= 3;
valid_learner_days = unique_learner_days(valid_day_mask_learner);

learner_avg_means = nan(numel(valid_learner_days), n_components);
for di = 1:numel(valid_learner_days)
    day_mask = learner_rel_days == valid_learner_days(di);
    learner_avg_means(di, :) = mean(learner_means(day_mask, :), 1, 'omitnan');
end

[sorted_learner_days, sort_idx] = sort(valid_learner_days);
sorted_learner_means = learner_avg_means(sort_idx, :);

% Aggregate data for non-learners
non_learner_rel_days = [non_learner_component_data.RelativeLearningDay]';
non_learner_means = vertcat(non_learner_component_data.ComponentMeans);

unique_non_learner_days = unique(non_learner_rel_days);
day_counts_non_learner = arrayfun(@(d) sum(non_learner_rel_days == d), unique_non_learner_days);

valid_day_mask_non_learner = day_counts_non_learner >= 3;
valid_non_learner_days = unique_non_learner_days(valid_day_mask_non_learner);

non_learner_avg_means = nan(numel(valid_non_learner_days), n_components);
for di = 1:numel(valid_non_learner_days)
    day_mask = non_learner_rel_days == valid_non_learner_days(di);
    non_learner_avg_means(di, :) = mean(non_learner_means(day_mask, :), 1, 'omitnan');
end

[sorted_non_learner_days, sort_idx] = sort(valid_non_learner_days);
sorted_non_learner_means = non_learner_avg_means(sort_idx, :);

% Create figure with 3 panels: Learners, Non-learners, Comparison
figure('Color','w','Position', [100, 100, 1800, 600], ...
    'Name',sprintf('CS+ Reaction Time by Group (n=%d) - %s', n_components, target_workflow));
component_spacing = 0.15;

% Overall title
sgtitle('CS+ Reaction Time: mPFC+ vs mPFC- Animals', ...
    'FontSize', 15, 'FontWeight', 'bold');

% --- PANEL 1: Learners (mPFC+) ---
subplot(1, 3, 1);
hold on;

if ~isempty(sorted_learner_days)
    for di = 1:numel(sorted_learner_days)
        day_x = sorted_learner_days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, sorted_learner_means(di, :), '-o', ...
            'LineWidth', 2.5, 'MarkerSize', 8, 'Color', cLearner, 'MarkerFaceColor', cLearner);
    end
end
hold off;

xlabel('Days Relative to Learning', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reaction Time (s)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', cLearner);
title(sprintf('mPFC+ (n=%d)', numel(learner_ids)), 'FontSize', 12, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
ax1 = gca;
ax1.YColor = cLearner;
box on; grid on;
set(gca, 'FontSize', 10);

% --- PANEL 2: Non-learners (mPFC-) ---
subplot(1, 3, 2);
hold on;

if ~isempty(sorted_non_learner_days)
    for di = 1:numel(sorted_non_learner_days)
        day_x = sorted_non_learner_days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        plot(component_x_positions, sorted_non_learner_means(di, :), '-o', ...
            'LineWidth', 2.5, 'MarkerSize', 8, 'Color', cNonLearner, 'MarkerFaceColor', cNonLearner);
    end
end
hold off;

xlabel('Days Relative to Learning', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reaction Time (s)', 'FontSize', 11, 'FontWeight', 'bold', 'Color', cNonLearner);
title(sprintf('mPFC- (n=%d)', numel(non_learner_ids)), 'FontSize', 12, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
ax2 = gca;
ax2.YColor = cNonLearner;
box on; grid on;
set(gca, 'FontSize', 10);

% --- PANEL 3: Overlaid Comparison ---
subplot(1, 3, 3);
hold on;

% Plot learners
if ~isempty(sorted_learner_days)
    for di = 1:numel(sorted_learner_days)
        day_x = sorted_learner_days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_learner = plot(component_x_positions, sorted_learner_means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 5, 'Color', cLearner, 'MarkerFaceColor', cLearner);

        if di == 1
            h_learner.DisplayName = 'mPFC+';
        else
            h_learner.HandleVisibility = 'off';
        end
    end
end

% Plot non-learners
if ~isempty(sorted_non_learner_days)
    for di = 1:numel(sorted_non_learner_days)
        day_x = sorted_non_learner_days(di);
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        h_non_learner = plot(component_x_positions, sorted_non_learner_means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 5, 'Color', cNonLearner, 'MarkerFaceColor', cNonLearner);

        if di == 1
            h_non_learner.DisplayName = 'mPFC-';
        else
            h_non_learner.HandleVisibility = 'off';
        end
    end
end
hold off;

xlabel('Days Relative to Learning', 'FontSize', 11, 'FontWeight', 'bold');
ylabel('Reaction Time (s)', 'FontSize', 11, 'FontWeight', 'bold');
title('Direct Comparison', 'FontSize', 12, 'FontWeight', 'bold');
xline(0, '--k', 'Learning', 'LineWidth', 2, 'LabelVerticalAlignment', 'bottom');
legend('Location', 'best', 'FontSize', 10);
set(gca, 'FontSize', 10);


%% Plot learners and non-learners overlaid for a behavioural metric seperated by n components for joint protocols

n_components = 4;

% Initialize storage
component_data = struct('Animal', {}, 'Day', {}, 'Workflow', {}, 'ComponentMeans', {}, 'IsLearner', {});

for a = 1:numel(behaviour_data)
    rd = behaviour_data(a).recording_day;

    % Get animal ID (handle both string and char consistently)
    if isfield(behaviour_data(a), 'animal_id')
        curr_animal_id = string(behaviour_data(a).animal_id);
    else
        error('animal_id field not found in behaviour_data(%d)', a);
    end

    % Find all stage 1 conditioning days
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    if ~any(isClassical), continue; end

    classical_days = find(isClassical);

    % Get learning day from combined_sig_day_all_protocols
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: no learning day found, skipping\n', curr_animal_id);
        continue;
    end

    % Determine if learner or non-learner
    isLearner = ismember(curr_animal_id, string(learner_ids));

    if ~isLearner
        continue;
    end

    for d = 1:numel(rd)
        wf = string(rd(d).workflow);
        stim = rd(d).right_stim_on_times;
        licks = rd(d).lick_event_times;

        if isempty(stim) || isempty(licks)
            continue
        end

        % Get CS+ and CS- trials
        cs_plus_mask = rd(d).cs_labels;

        % In case of ITI lick rates
        % data_metric = rd(d).ITI_lick_counts ./ rd(d).ITI_actual_duration;
        data_metric= rd(d).anticipatory_licks(cs_plus_mask);

        % Split trials into n components by index
        nT = numel(stim); % if not discrimination

        if nT < n_components
            continue  % Skip sessions with too few trials
        end

        edges = round(linspace(1, nT+1, n_components+1));
        component_means = nan(1, n_components);

        for c = 1:n_components
            idx = edges(c):edges(c+1)-1;
            component_means(c) = mean(data_metric(idx), 'omitnan');
        end

        component_data(end+1) = struct( ...
            'Animal', curr_animal_id, ...
            'Day', d, ...
            'Workflow', wf, ...
            'ComponentMeans', component_means, ...
            'IsLearner', isLearner);
    end
end

% Plot results
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Get unique protocols
all_protocols = unique(string({component_data.Workflow})', 'stable');

% Create figure
figure('Color','w','Position', [100, 100, 1400, 600], ...
    'Name',sprintf('Behavioral Components (n=%d) - Learners vs Non-Learners', n_components));
tiledlayout(1, numel(all_protocols), 'TileSpacing','compact');

% Compute global y-limits
all_means = vertcat(component_data.ComponentMeans);
all_means = all_means(~isnan(all_means));
if isempty(all_means)
    return;
end
ymin = min(0, min(all_means));
ymax = max(all_means);

for wi = 1:numel(all_protocols)
    nexttile;
    wf = all_protocols(wi);

    % Select data for this protocol
    protocol_mask = strcmp(string({component_data.Workflow})', wf);
    protocol_data = component_data(protocol_mask);

    % Separate learners and non-learners
    learner_mask = [protocol_data.IsLearner]';
    learner_data = protocol_data(learner_mask);
    non_learner_data = protocol_data(~learner_mask);

    % Process learners
    if ~isempty(learner_data)
        unique_animals_L = unique(string({learner_data.Animal})');
        workflow_relative_days_L = nan(size(learner_data));

        for ai = 1:numel(unique_animals_L)
            animal_name = unique_animals_L(ai);
            animal_mask = strcmp(string({learner_data.Animal})', animal_name);
            animal_indices = find(animal_mask);
            animal_days = [learner_data(animal_indices).Day];
            [~, sort_order] = sort(animal_days);
            sorted_indices = animal_indices(sort_order);

            for i = 1:numel(sorted_indices)
                workflow_relative_days_L(sorted_indices(i)) = i;
            end
        end

        days_L = workflow_relative_days_L';
        means_L = vertcat(learner_data.ComponentMeans);

        % For second protocol, limit to first 7 days
        if wi == 2
            valid_mask_L = days_L <= 7;
            days_L = days_L(valid_mask_L);
            means_L = means_L(valid_mask_L, :);
        end

        unique_days_L = unique(days_L);
        day_counts_L = arrayfun(@(d) sum(days_L == d), unique_days_L);
        valid_day_mask_L = day_counts_L >= 3;
        valid_days_L = unique_days_L(valid_day_mask_L);

        if ~isempty(valid_days_L)
            avg_means_L = nan(numel(valid_days_L), n_components);
            for di = 1:numel(valid_days_L)
                day_mask = days_L == valid_days_L(di);
                avg_means_L(di, :) = mean(means_L(day_mask, :), 1, 'omitnan');
            end
            [sorted_days_L, sort_idx] = sort(valid_days_L);
            sorted_means_L = avg_means_L(sort_idx, :);
        else
            sorted_days_L = [];
            sorted_means_L = [];
        end
    else
        sorted_days_L = [];
        sorted_means_L = [];
    end

    % Process non-learners
    if ~isempty(non_learner_data)
        unique_animals_NL = unique(string({non_learner_data.Animal})');
        workflow_relative_days_NL = nan(size(non_learner_data));

        for ai = 1:numel(unique_animals_NL)
            animal_name = unique_animals_NL(ai);
            animal_mask = strcmp(string({non_learner_data.Animal})', animal_name);
            animal_indices = find(animal_mask);
            animal_days = [non_learner_data(animal_indices).Day];
            [~, sort_order] = sort(animal_days);
            sorted_indices = animal_indices(sort_order);

            for i = 1:numel(sorted_indices)
                workflow_relative_days_NL(sorted_indices(i)) = i;
            end
        end

        days_NL = workflow_relative_days_NL';
        means_NL = vertcat(non_learner_data.ComponentMeans);

        % For second protocol, limit to first 7 days
        if wi == 2
            valid_mask_NL = days_NL <= 7;
            days_NL = days_NL(valid_mask_NL);
            means_NL = means_NL(valid_mask_NL, :);
        end

        unique_days_NL = unique(days_NL);
        day_counts_NL = arrayfun(@(d) sum(days_NL == d), unique_days_NL);
        valid_day_mask_NL = day_counts_NL >= 3;
        valid_days_NL = unique_days_NL(valid_day_mask_NL);

        if ~isempty(valid_days_NL)
            avg_means_NL = nan(numel(valid_days_NL), n_components);
            for di = 1:numel(valid_days_NL)
                day_mask = days_NL == valid_days_NL(di);
                avg_means_NL(di, :) = mean(means_NL(day_mask, :), 1, 'omitnan');
            end
            [sorted_days_NL, sort_idx] = sort(valid_days_NL);
            sorted_means_NL = avg_means_NL(sort_idx, :);
        else
            sorted_days_NL = [];
            sorted_means_NL = [];
        end
    else
        sorted_days_NL = [];
        sorted_means_NL = [];
    end

    % Plot components overlaid
    hold on;
    component_spacing = 0.15;

    % Plot non-learners
    if ~isempty(sorted_days_NL)
        for di = 1:numel(sorted_days_NL)
            day_x = sorted_days_NL(di);
            component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
                component_spacing*(n_components-1)/2, n_components);

            h_NL = plot(component_x_positions, sorted_means_NL(di, :), '-o', ...
                'LineWidth', 2, 'MarkerSize', 7, 'Color', cNonLearner, 'MarkerFaceColor', cNonLearner);

            if di == 1
                h_NL.DisplayName = sprintf('mPFC- (n=%d)', numel(unique_animals_NL));
            else
                h_NL.HandleVisibility = 'off';
            end
        end
    end

    % Plot learners
    if ~isempty(sorted_days_L)
        for di = 1:numel(sorted_days_L)
            day_x = sorted_days_L(di);
            component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
                component_spacing*(n_components-1)/2, n_components);

            h_L = plot(component_x_positions, sorted_means_L(di, :), '-o', ...
                'LineWidth', 2, 'MarkerSize', 7, 'Color', cLearner, 'MarkerFaceColor', cLearner);

            if di == 1
                h_L.DisplayName = sprintf('mPFC+ (n=%d)', numel(unique_animals_L));
            else
                h_L.HandleVisibility = 'off';
            end
        end
    end

    hold off;

    xlabel('Day (Workflow-Relative)', 'FontSize', 11, 'FontWeight', 'bold');
    ylabel('ITI Lick Rate (licks/s)', 'FontSize', 11, 'FontWeight', 'bold');
    ylim([ymin ymax]);

    % Add protocol-specific title info
    if wi == 2
        title(sprintf('%s (First 7 Days)', strrep(char(wf), '_', ' ')), ...
            'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
    else
        title(sprintf('%s', strrep(char(wf), '_', ' ')), ...
            'Interpreter', 'none', 'FontSize', 12, 'FontWeight', 'bold');
    end

    legend('Location', 'best', 'FontSize', 10);
    box on;
    set(gca, 'FontSize', 10, 'LineWidth', 1.2);
end

sgtitle('ITI Lick Rate Across Training: mPFC+ vs mPFC-', ...
    'FontSize', 14, 'FontWeight', 'bold');

%%

% Analyze behavioral metric aligned to learning day

% Set parameters
learners_group_ID = {'HA005','HA008','HA010','HA011','HA012'}; % Define learners
n_components = 4; % Number of components to split each session into

% Initialize storage - now including relative learning day
component_data = struct('Animal', {}, 'Day', {}, 'Workflow', {}, 'ComponentMeans', {}, 'RelativeLearningDay', {});

for a = 1:numel(behaviour_data)
    rd = behaviour_data(a).recording_day;

    % Get animal ID (handle both string and char consistently)
    if isfield(behaviour_data(a), 'animal_id')
        curr_animal_id = string(behaviour_data(a).animal_id);
    else
        error('animal_id field not found in behaviour_data(%d)', a);
    end

    % Skip non-learners
    if ~ismember(curr_animal_id, string(learners_group_ID))
        continue
    end

    % Find animal index for learning day extraction

    % Get the significance days vector for this animal & protocol
    sig_days_all = combined_sig_day_all_protocols{a};
    sigDays = sig_days_all;

    % Find the first significant day (learning day)
    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, saving just pre-learning days', curr_animal_id);
        ld = numel(validDays) + 1; % define post-learning outside the scope
    end

    for d = 1:numel(rd)
        wf = string(rd(d).workflow);
        stim = rd(d).right_stim_on_times;
        licks = rd(d).lick_event_times;

        if isempty(stim) || isempty(licks)
            continue
        end

        % Count licks in [0,1) s after each stim
        data_metric = arrayfun(@(t) sum(licks >= t & licks < t+1), stim);

        % Split trials into n components by index
        nT = numel(stim);
        if nT < n_components
            continue  % Skip sessions with too few trials
        end

        edges = round(linspace(1, nT+1, n_components+1));
        component_means = nan(1, n_components);

        for c = 1:n_components
            idx = edges(c):edges(c+1)-1;
            component_means(c) = mean(data_metric(idx), 'omitnan');
        end

        % Calculate relative learning day
        relative_learning_day = d - ld;

        component_data(end+1) = struct( ...
            'Animal', curr_animal_id, ...
            'Day', d, ...
            'Workflow', wf, ...
            'ComponentMeans', component_means, ...
            'RelativeLearningDay', relative_learning_day);
    end
end

% Plot results aligned to learning day
if isempty(component_data)
    warning('No data to plot');
    return;
end

% Get unique protocols
all_protocols = unique(string({component_data.Workflow})', 'stable');

% Create figure
figure('Color','w','Name',sprintf('Learning-Aligned Components (n=%d) - Average Across Animals', n_components));
tiledlayout(1, numel(all_protocols), 'TileSpacing','compact');

% Compute global y-limits
all_means = vertcat(component_data.ComponentMeans);
all_means = all_means(~isnan(all_means));
if isempty(all_means)
    return;
end
ymin = max(0, min(all_means));
ymax = max(all_means);

for wi = 1:numel(all_protocols)
    nexttile;
    wf = all_protocols(wi);

    % Select data for this protocol
    protocol_mask = strcmp(string({component_data.Workflow})', wf);
    relative_days = [component_data(protocol_mask).RelativeLearningDay]';
    means = vertcat(component_data(protocol_mask).ComponentMeans);

    % Count how many animals per relative learning day and filter days with <2 animals
    unique_rel_days = unique(relative_days);
    day_counts = arrayfun(@(d) sum(relative_days == d), unique_rel_days);

    % Keep only days with 2+ animals
    valid_day_mask = day_counts >= 2;
    valid_rel_days = unique_rel_days(valid_day_mask);

    if isempty(valid_rel_days)
        warning('No relative learning days with >=2 animals for protocol %s', wf);
        continue;
    end

    % Average across animals for valid relative learning days only
    avg_means = nan(numel(valid_rel_days), n_components);
    for di = 1:numel(valid_rel_days)
        day_mask = relative_days == valid_rel_days(di);
        avg_means(di, :) = mean(means(day_mask, :), 1, 'omitnan');
    end

    % Sort by relative learning day
    [sorted_rel_days, sort_idx] = sort(valid_rel_days);
    sorted_means = avg_means(sort_idx, :);

    % Plot components for each relative day as separate points connected by lines
    hold on;
    component_spacing = 0.15; % Horizontal spacing between components within a day

    for di = 1:numel(sorted_rel_days)
        day_x = sorted_rel_days(di);
        % Create x-positions for components around each day
        component_x_positions = day_x + linspace(-component_spacing*(n_components-1)/2, ...
            component_spacing*(n_components-1)/2, n_components);

        % Plot components for this day connected by a line (same color for all days)
        plot(component_x_positions, sorted_means(di, :), '-o', ...
            'LineWidth', 1.5, 'MarkerSize', 6, 'Color', [0.2 0.4 0.8]);
    end
    hold off;

    xlabel('Days Relative to Learning');
    ylabel('Mean Licks (0-1s)');
    ylim([ymin ymax]);
    title(sprintf('Protocol: %s', wf), 'Interpreter', 'none');
    grid on;

    % Add vertical line at learning day (day 0)
    xline(0, '--r', 'Learning Day', 'LineWidth', 2);
end

%% Wheel velocity analysis
%% Creates a wheel velocity structure and plots the absolute average wheel velocity aligned to stim onset across all days

% create grouping animal index
behaviour_animal_idx = grp2idx(cell2mat(cellfun(@(animal,wf) repmat(animal,length(wf),1), ...
    {behaviour_data.animal_id},{behaviour_data.recording_day},'uni',false)'));

% animals that have not learned the static
group_animals = {'DS017','HA006','HA007','HA009','HA011','HA013','HA014','HA015'};

% Create a list of animal IDs repeated per recording day
animal_ids_all_days = cellfun(@(animal, wf) repmat({animal}, length(wf), 1), ...
    {behaviour_data.animal_id}, {behaviour_data.recording_day}, 'UniformOutput', false);

% Convert to column vector
animal_ids_all_days_stacked = vertcat(animal_ids_all_days{:});

% Create a logical index for whether each row belongs to the group
is_group_animal = ismember(animal_ids_all_days_stacked, group_animals);  % (n_days_all x 1


% Create a logical learning index variable (n all days x [0,1])
learning_index_animal = vertcat(combined_sig_day_all_protocols{:});

% Get all wheel velocity, workflow, and CS+ logical mask across animals x days
stim_on_aligned_wheel_vel_all_trials= cellfun(@(x) {x.stim_on_aligned_wheel_vel},{behaviour_data.recording_day},'uni',false);
stim_on_aligned_wheel_vel_all_flat = [stim_on_aligned_wheel_vel_all_trials{:}]';

cs_labels_all= cellfun(@(x) {x.cs_labels},{behaviour_data.recording_day},'uni',false);
cs_labels_all_flat= [cs_labels_all{:}]';

workflow_animal = cellfun(@(x) {x.workflow},{behaviour_data.recording_day},'uni',false);
workflow_cat = grp2idx(horzcat(workflow_animal{:}));

selected_workflow= workflow_cat==2; % select workflow

% run through and seperate based on CS+ or CS-
N = numel(stim_on_aligned_wheel_vel_all_flat);
cs_plus_flat  = cell(N,1);
cs_minus_flat = cell(N,1);

for d = 1:N
    V = abs(stim_on_aligned_wheel_vel_all_flat{d});   % nTrials×T
    L = cs_labels_all_flat{d};                   % nTrials×1 of "CS+" / "CS-"

    cs_plus_flat{d}  = V(L,:);   % only the CS+ trials
    cs_minus_flat{d} = V(~L,:);   % only the CS- trials
end


% Filter out
behaviour_animal_idx= behaviour_animal_idx(selected_workflow);
is_group_animal= is_group_animal(selected_workflow);
learning_index_animal= learning_index_animal(selected_workflow);
cs_plus_flat= cs_plus_flat(selected_workflow);
cs_minus_flat= cs_minus_flat(selected_workflow);
t= surround_time_points;


animals = unique(behaviour_animal_idx);

for ai = animals(:)'
    % which days belong to this animal?
    d_idx = find(behaviour_animal_idx==ai);
    nDays = numel(d_idx);
    if nDays==0, continue; end

    % pick a colormap of length nDays
    cmap = parula(nDays);

    % build a new figure for this animal
    animalID = behaviour_data(ai).animal_id;
    figure('Name',sprintf('Animal %s — wheel‐vel',animalID),'Color','w',...
        'Position',[200 200 800 400]);
    tlo = tiledlayout(2,1,'TileSpacing','compact');

    % ---- CS+ panel ----
    nexttile;
    hold on;
    for j = 1:nDays
        d = d_idx(j);
        V = cs_plus_flat{d};    % [nTrials×T]
        % 1×T

        m   = mean( V, 1, 'omitnan' );
        sem =  std( V, 0, 1, 'omitnan' ) ./ sqrt(size(V,1));

        % shaded error
        fill([t, fliplr(t)], [m+sem, fliplr(m-sem)], ...
            cmap(j,:), 'FaceAlpha',0.2, 'EdgeColor','none');
        % mean line
        plot(t, m, 'Color',cmap(j,:), 'LineWidth',1.5);
    end
    xline(0,'k--');
    ylabel('CS+ vel');
    title('CS+ across days');
    hold off;

    % ---- CS– panel ----
    nexttile;
    hold on;
    for j = 1:nDays
        d = d_idx(j);
        V = cs_minus_flat{d};

        m   = mean( V, 1, 'omitnan' );
        sem =  std( V, 0, 1, 'omitnan' ) ./ sqrt(size(V,1));

        fill([t, fliplr(t)], [m+sem, fliplr(m-sem)], ...
            cmap(j,:), 'FaceAlpha',0.2, 'EdgeColor','none');
        plot(t, m, 'Color',cmap(j,:), 'LineWidth',1.5);
    end
    xline(0,'k--');
    xlabel('Time (s)');
    ylabel('CS- vel');
    title('CS- across days');
    hold off;

    % shared legend: day 1→nDays
    cb = colorbar('Ticks',linspace(0,1,nDays), ...
        'TickLabels',arrayfun(@num2str,1:nDays,'uni',false));
    cb.Label.String = 'Day #';
    sgtitle(sprintf('Animal %s — wheel velocity by day', animalID));
end


%% Pre-post learning comparison of the wheel velocity

for ai = animals(:)'
    % ——— find which “days” belong to this animal ———
    d_idx = find(behaviour_animal_idx==ai);
    if isempty(d_idx), continue; end

    % ——— split into pre vs post according to learning_index_animal ———
    pre_idx  = d_idx( learning_index_animal(d_idx)==0 );
    post_idx = d_idx( learning_index_animal(d_idx)==1 );

    % require at least one post‐learning day
    if isempty(post_idx)
        warning('Animal %d has no post‐learning days, skipping.', ai);
        continue;
    end


    % ——— gather all trials from those days ———
    Vp_plus  = cellfun(@(V) mean(V,1,'omitnan'), cs_plus_flat(pre_idx),  'UniformOutput',false);
    Vp_plus_mat= vertcat(Vp_plus{:});
    mP_pre= mean(Vp_plus_mat,  1,'omitnan');
    semP_pre= std( Vp_plus_mat,  0, 1, 'omitnan') ./ sqrt(size(Vp_plus_mat,1));

    Vp_minus= cellfun(@(V) mean(V,1,'omitnan'), cs_minus_flat(pre_idx),  'UniformOutput',false);
    Vp_minus_mat= vertcat(Vp_minus{:});
    mM_pre  = mean(Vp_minus_mat, 1, 'omitnan');
    semM_pre= std( Vp_minus_mat, 0, 1, 'omitnan') ./ sqrt(size(Vp_minus_mat,1));


    Vt_plus  = cellfun(@(V) mean(V,1,'omitnan'), cs_plus_flat(post_idx),  'UniformOutput',false);
    Vt_plus_mat= vertcat(Vt_plus{:});
    mP_post  = mean(Vt_plus_mat,  1, 'omitnan');
    semP_post= std( Vt_plus_mat,  0, 1, 'omitnan') ./ sqrt(size(Vt_plus_mat,1));

    Vt_minus = cellfun(@(V) mean(V,1,'omitnan'), cs_minus_flat(post_idx),  'UniformOutput',false);
    Vt_minus_mat= vertcat(Vt_minus{:});
    mM_post  = mean(Vt_minus_mat, 1, 'omitnan');
    semM_post= std( Vt_minus_mat, 0, 1, 'omitnan') ./ sqrt(size(Vt_minus_mat,1));


    % ——— Plot ———
    animalID = behaviour_data(ai).animal_id;
    figure('Name',sprintf('Animal %s — pre vs post',animalID),'Color','w',...
        'Position',[200 200 800 400]);
    tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

    % CS+ panel
    ax1 = nexttile;
    hold(ax1,'on');
    if ~isempty(mP_pre) % in case no pre-learning days
        fill([t,fliplr(t)],[mP_pre+semP_pre,fliplr(mP_pre-semP_pre)], ...
            'b','FaceAlpha',0.2,'EdgeColor','none');
        plot(t,mP_pre,'b-','LineWidth',2);
    end
    fill([t,fliplr(t)],[mP_post+semP_post,fliplr(mP_post-semP_post)], ...
        'b','FaceAlpha',0.1,'EdgeColor','none');
    plot(t,mP_post,'b--','LineWidth',2);
    xline(0,'k--','LineWidth',1);
    ylabel('CS+ wheel vel');
    legend({'Pre ±SEM','Pre mean','Post ±SEM','Post mean'},'Location','best');
    title('CS+ pre vs post');
    hold(ax1,'off');

    % CS– panel
    ax2 = nexttile;
    hold(ax2,'on');
    if ~isempty(mM_pre)
        fill([t,fliplr(t)],[mM_pre+semM_pre,fliplr(mM_pre-semM_pre)], ...
            'r','FaceAlpha',0.2,'EdgeColor','none');
        plot(t,mM_pre,'r-','LineWidth',2);
    end
    fill([t,fliplr(t)],[mM_post+semM_post,fliplr(mM_post-semM_post)], ...
        'r','FaceAlpha',0.1,'EdgeColor','none');
    plot(t,mM_post,'r--','LineWidth',2);
    xline(0,'k--','LineWidth',1);
    xlabel('Time (s)');
    ylabel('CS– wheel vel');
    legend({'Pre ±SEM','Pre mean','Post ±SEM','Post mean'},'Location','best');
    title('CS– pre vs post');
    hold(ax2,'off');

    sgtitle(sprintf('Animal %s — CS+ and CS– pre vs post learning', animalID));
end




%% Pre-post learning wheel velocity seperated by averaging learners vs non-learners

% Filter out based on learners vs non-learners
behaviour_animal_idx= behaviour_animal_idx(is_group_animal);
learning_index_animal= learning_index_animal(is_group_animal);
cs_plus_flat= cs_plus_flat(is_group_animal);
cs_minus_flat= cs_minus_flat(is_group_animal);


% 1) split into pre vs post learning
pre_idx  = find(learning_index_animal==0 );
post_idx = find(learning_index_animal==1 );
assert(~isempty(pre_idx) && ~isempty(post_idx), 'need both pre & post days!');

% 2) build day×time matrices
%    for CS+

M_pre_mat  = cell2mat(cellfun(@(V) mean(V,1,'omitnan'), cs_plus_flat(pre_idx),  'uni',false));
M_post_mat = cell2mat( cellfun(@(V) mean(V,1,'omitnan'), cs_plus_flat(post_idx), 'uni',false));

%    for CS–
N_pre_mat  = cell2mat( cellfun(@(V) mean(V,1,'omitnan'), cs_minus_flat(pre_idx),  'uni',false));
N_post_mat = cell2mat( cellfun(@(V) mean(V,1,'omitnan'), cs_minus_flat(post_idx), 'uni',false));

% 3) now get mean±SEM across *days* (rows)
mu_M_pre   = mean(M_pre_mat,  1, 'omitnan');
sem_M_pre  = std( M_pre_mat,  0, 1, 'omitnan') ./ sqrt(size(M_pre_mat,1));

mu_M_post  = mean(M_post_mat, 1, 'omitnan');
sem_M_post = std( M_post_mat, 0, 1, 'omitnan') ./ sqrt(size(M_post_mat,1));

mu_N_pre   = mean(N_pre_mat,  1, 'omitnan');
sem_N_pre  = std( N_pre_mat,  0, 1, 'omitnan') ./ sqrt(size(N_pre_mat,1));

mu_N_post  = mean(N_post_mat, 1, 'omitnan');
sem_N_post = std( N_post_mat, 0, 1, 'omitnan') ./ sqrt(size(N_post_mat,1));

% 4) plot them all in one figure
figure('Color','w','Position',[200 200 900 500]);
t = surround_time_points;    % your time‐axis vector

% CS+ subplot
subplot(2,1,1); hold on;
fill([t,fliplr(t)], [mu_M_pre+sem_M_pre, fliplr(mu_M_pre-sem_M_pre)], 'b','FaceAlpha',0.2,'EdgeColor','none');
plot(t,mu_M_pre, '-b','LineWidth',2);
fill([t,fliplr(t)], [mu_M_post+sem_M_post, fliplr(mu_M_post-sem_M_post)], 'b','FaceAlpha',0.1,'EdgeColor','none');
plot(t,mu_M_post,'--b','LineWidth',2);
xline(0,'k--');
title('CS+ pre vs post (group‐average)');
ylabel('Mean wheel‐vel');
legend({'Pre ±SEM','Pre','Post ±SEM','Post'},'Location','best');
hold off;

% CS– subplot
subplot(2,1,2); hold on;
fill([t,fliplr(t)], [mu_N_pre+sem_N_pre, fliplr(mu_N_pre-sem_N_pre)], 'r','FaceAlpha',0.2,'EdgeColor','none');
plot(t,mu_N_pre, '-r','LineWidth',2);
fill([t,fliplr(t)], [mu_N_post+sem_N_post, fliplr(mu_N_post-sem_N_post)], 'r','FaceAlpha',0.1,'EdgeColor','none');
plot(t,mu_N_post,'--r','LineWidth',2);
xline(0,'k--');
title('CS– pre vs post (group‐average)');
xlabel('Time (s)');
ylabel('Mean wheel‐vel');
legend({'Pre ±SEM','Pre','Post ±SEM','Post'},'Location','best');
hold off;

sgtitle('Group‐average Non-Learners CS+ and CS– wheel‐velocity pre vs post learning');

%% Plot ITI lick Rate Aligned to Learning Day Split by Learners vs Non-Learners

% Analyze ITI lick rate aligned to learning day
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move';

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};
% non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};

% Initialize storage arrays
all_animal_ids = {};
all_iti_lick_rates = [];
all_relative_days = [];
all_is_learner = [];

% Loop over animals to collect data
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording found, skipping', animal_id);
        continue;
    end

    % Get significance days for learning day determination
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % Find first significant day (learning day)
    ld = find(sigDays, 1, 'first');

    if isempty(ld)
        warning('Animal %s: no learning day found, skipping', animal_id);
        continue;
    end

    % Calculate relative days for this animal
    n = numel(validDays);
    rel_day = (1:n) - ld; % 0 = day of learning

    % Process each valid day
    for d = 1:numel(validDays)
        day_data = validDays(d);

        % Check if required fields exist
        if ~isfield(day_data, 'ITI_lick_counts') || ~isfield(day_data, 'ITI_actual_duration')
            warning('Animal %s day %d: missing ITI fields, skipping', animal_id, d);
            continue;
        end

        iti_licks = day_data.ITI_lick_counts;
        iti_duration = day_data.ITI_actual_duration;

        % %Calculate ITI lick rate (licks per second)
        mean_iti_lick_rate = nanmean(iti_licks ./ iti_duration);

        % % Test RTs for CS+
        % cs_plus_mask= day_data.cs_labels;
        % mean_iti_lick_rate= median(day_data.all_stim_diff_from_optimal_reward(cs_plus_mask));

        % Store data
        all_animal_ids{end+1} = animal_id;
        all_iti_lick_rates(end+1) = mean_iti_lick_rate;
        all_relative_days(end+1) = rel_day(d);
        all_is_learner(end+1) = isLearner;
    end
end

% Check if we have data
if isempty(all_iti_lick_rates)
    warning('No ITI data to plot');
    return;
end

% Convert to column vectors
all_iti_lick_rates = all_iti_lick_rates(:);
all_relative_days = all_relative_days(:);
all_is_learner = all_is_learner(:);

% Separate learners and non-learners
learner_mask = all_is_learner;
learner_rel_days = all_relative_days(learner_mask==1);
learner_rates = all_iti_lick_rates(learner_mask==1);
non_learner_rel_days = all_relative_days(~learner_mask);
non_learner_rates = all_iti_lick_rates(~learner_mask);

% Aggregate learners by relative day (require at least 3 animals)
min_animals = 3;
unique_learner_days = unique(learner_rel_days);
learner_day_counts = arrayfun(@(d) sum(learner_rel_days == d), unique_learner_days);
valid_learner_days = unique_learner_days(learner_day_counts >= min_animals);

learner_mean = nan(size(valid_learner_days));
learner_sem = nan(size(valid_learner_days));

for di = 1:numel(valid_learner_days)
    day_mask = learner_rel_days == valid_learner_days(di);
    day_rates = learner_rates(day_mask);

    learner_mean(di) = mean(day_rates, 'omitnan');
    learner_sem(di) = std(day_rates, 'omitnan') / sqrt(sum(~isnan(day_rates)));
end

% Sort learner data
[learner_days_sorted, sort_idx] = sort(valid_learner_days);
learner_mean = learner_mean(sort_idx);
learner_sem = learner_sem(sort_idx);

% Aggregate non-learners by relative day
unique_non_learner_days = unique(non_learner_rel_days);
non_learner_day_counts = arrayfun(@(d) sum(non_learner_rel_days == d), unique_non_learner_days);
valid_non_learner_days = unique_non_learner_days(non_learner_day_counts >= min_animals);

non_learner_mean = nan(size(valid_non_learner_days));
non_learner_sem = nan(size(valid_non_learner_days));

for di = 1:numel(valid_non_learner_days)
    day_mask = non_learner_rel_days == valid_non_learner_days(di);
    day_rates = non_learner_rates(day_mask);

    non_learner_mean(di) = mean(day_rates, 'omitnan');
    non_learner_sem(di) = std(day_rates, 'omitnan') / sqrt(sum(~isnan(day_rates)));
end

% Sort non-learner data
[non_learner_days_sorted, sort_idx] = sort(valid_non_learner_days);
non_learner_mean = non_learner_mean(sort_idx);
non_learner_sem = non_learner_sem(sort_idx);

% Create figure
figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

% Plot learners
if ~isempty(learner_days_sorted)
    % Error band (SEM)
    fill([learner_days_sorted; flipud(learner_days_sorted)], ...
        [learner_mean + learner_sem; flipud(learner_mean - learner_sem)], ...
        cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Mean line
    plot(learner_days_sorted, learner_mean, '-o', 'Color', cLearner, ...
        'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', cLearner, ...
        'DisplayName', sprintf('mPFC+ (n=%d)', numel(learner_ids)));
end

% Plot non-learners
if ~isempty(non_learner_days_sorted)
    % Error band (SEM)
    fill([non_learner_days_sorted; flipud(non_learner_days_sorted)], ...
        [non_learner_mean + non_learner_sem; flipud(non_learner_mean - non_learner_sem)], ...
        cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Mean line
    plot(non_learner_days_sorted, non_learner_mean, '-o', 'Color', cNonLearner, ...
        'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', cNonLearner, ...
        'DisplayName', sprintf('mPFC- (n=%d)', numel(non_learner_ids)));
end

% Mark learning day
xline(0, '--k', 'Learning Day', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold','HandleVisibility', 'off');

% Formatting
xlabel('Days Relative to Learning', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('ITI Lick Counts (Licks/s)', 'FontSize', 13, 'FontWeight', 'bold');
title('ITI Licks Aligned to Learning Day', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

% Clean axes
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

% Print summary statistics
fprintf('\n=== ITI Lick Rate Summary ===\n');
if ~isempty(learner_days_sorted)
    fprintf('Learners (mPFC+): n=%d animals, %d days\n', numel(learner_ids), numel(learner_days_sorted));
    pre_idx = learner_days_sorted < 0;
    post_idx = learner_days_sorted >= 0;
    if any(pre_idx)
        fprintf('  Pre-learning mean: %.3f licks/s\n', mean(learner_mean(pre_idx)));
    end
    if any(post_idx)
        fprintf('  Post-learning mean: %.3f licks/s\n', mean(learner_mean(post_idx)));
    end
end
if ~isempty(non_learner_days_sorted)
    fprintf('Non-learners (mPFC-): n=%d animals, %d days\n', numel(non_learner_ids), numel(non_learner_days_sorted));
    pre_idx = non_learner_days_sorted < 0;
    post_idx = non_learner_days_sorted >= 0;
    if any(pre_idx)
        fprintf('  Pre-learning mean: %.3f licks/s\n', mean(non_learner_mean(pre_idx)));
    end
    if any(post_idx)
        fprintf('  Post-learning mean: %.3f licks/s\n', mean(non_learner_mean(post_idx)));
    end
end

%% Plot the reward rate split between mPFC+ and mPFC-

% Analyze reward rate aligned to learning day
% Set parameters
target_workflow = 'visual_operant_lick_two_stim_right_move';

learner_ids = {'HA005','HA008','HA010','HA011','HA012'};
non_learner_ids = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};

all_ids = {behaviour_data.animal_id};

% Initialize storage arrays
all_animal_ids = {};
all_reward_rates = [];
all_relative_days = [];
all_is_learner = [];

% Loop over animals to collect data
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, learner_ids);

    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);

    if isempty(validDays)
        warning('Animal %s: no recording found, skipping', animal_id);
        continue;
    end

    % Get significance days for learning day determination
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);

    % Find first significant day (learning day)
    ld = find(sigDays, 1, 'first');

    if isempty(ld)
        warning('Animal %s: no learning day found, skipping', animal_id);
        continue;
    end

    % Calculate relative days for this animal
    n = numel(validDays);
    rel_day = (1:n) - ld; % 0 = day of learning

    % Process each valid day
    for d = 1:numel(validDays)
        day_data = validDays(d);

        % Check if reward_rate field exists
        if ~isfield(day_data, 'reward_rate')
            warning('Animal %s day %d: missing reward_rate field, skipping', animal_id, d);
            continue;
        end

        reward_rate = day_data.reward_rate;

        % Skip if reward rate is invalid
        if isnan(reward_rate) || isinf(reward_rate)
            warning('Animal %s day %d: invalid reward rate, skipping', animal_id, d);
            continue;
        end

        % Store data
        all_animal_ids{end+1} = animal_id;
        all_reward_rates(end+1) = reward_rate;
        all_relative_days(end+1) = rel_day(d);
        all_is_learner(end+1) = isLearner;
    end
end

% Check if we have data
if isempty(all_reward_rates)
    warning('No reward rate data to plot');
    return;
end

% Convert to column vectors
all_reward_rates = all_reward_rates(:);
all_relative_days = all_relative_days(:);
all_is_learner = all_is_learner(:);

% Separate learners and non-learners
learner_mask = all_is_learner;
learner_rel_days = all_relative_days(learner_mask==1);
learner_rates = all_reward_rates(learner_mask==1);
non_learner_rel_days = all_relative_days(~learner_mask);
non_learner_rates = all_reward_rates(~learner_mask);

% Aggregate learners by relative day (require at least 3 animals)
min_animals = 2;
unique_learner_days = unique(learner_rel_days);
learner_day_counts = arrayfun(@(d) sum(learner_rel_days == d), unique_learner_days);
valid_learner_days = unique_learner_days(learner_day_counts >= min_animals);

learner_mean = nan(size(valid_learner_days));
learner_sem = nan(size(valid_learner_days));
learner_n_animals = nan(size(valid_learner_days));

for di = 1:numel(valid_learner_days)
    day_mask = learner_rel_days == valid_learner_days(di);
    day_rates = learner_rates(day_mask);

    learner_mean(di) = mean(day_rates, 'omitnan');
    learner_sem(di) = std(day_rates, 'omitnan') / sqrt(sum(~isnan(day_rates)));
    learner_n_animals(di) = sum(~isnan(day_rates));
end

% Sort learner data
[learner_days_sorted, sort_idx] = sort(valid_learner_days);
learner_mean = learner_mean(sort_idx);
learner_sem = learner_sem(sort_idx);
learner_n_animals = learner_n_animals(sort_idx);

% Aggregate non-learners by relative day
unique_non_learner_days = unique(non_learner_rel_days);
non_learner_day_counts = arrayfun(@(d) sum(non_learner_rel_days == d), unique_non_learner_days);
valid_non_learner_days = unique_non_learner_days(non_learner_day_counts >= min_animals);

non_learner_mean = nan(size(valid_non_learner_days));
non_learner_sem = nan(size(valid_non_learner_days));
non_learner_n_animals = nan(size(valid_non_learner_days));

for di = 1:numel(valid_non_learner_days)
    day_mask = non_learner_rel_days == valid_non_learner_days(di);
    day_rates = non_learner_rates(day_mask);

    non_learner_mean(di) = mean(day_rates, 'omitnan');
    non_learner_sem(di) = std(day_rates, 'omitnan') / sqrt(sum(~isnan(day_rates)));
    non_learner_n_animals(di) = sum(~isnan(day_rates));
end

% Sort non-learner data
[non_learner_days_sorted, sort_idx] = sort(valid_non_learner_days);
non_learner_mean = non_learner_mean(sort_idx);
non_learner_sem = non_learner_sem(sort_idx);
non_learner_n_animals = non_learner_n_animals(sort_idx);

% Create figure
figure('Color', 'w', 'Position', [100, 100, 1000, 650]);
hold on;

% Plot learners
if ~isempty(learner_days_sorted)
    % Error band (SEM)
    fill([learner_days_sorted; flipud(learner_days_sorted)], ...
        [learner_mean + learner_sem; flipud(learner_mean - learner_sem)], ...
        cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Mean line
    plot(learner_days_sorted, learner_mean, '-o', 'Color', cLearner, ...
        'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', cLearner, ...
        'DisplayName', sprintf('mPFC+ (n=%d)', numel(learner_ids)));
end

% Plot non-learners
if ~isempty(non_learner_days_sorted)
    % Error band (SEM)
    fill([non_learner_days_sorted; flipud(non_learner_days_sorted)], ...
        [non_learner_mean + non_learner_sem; flipud(non_learner_mean - non_learner_sem)], ...
        cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % Mean line
    plot(non_learner_days_sorted, non_learner_mean, '-o', 'Color', cNonLearner, ...
        'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', cNonLearner, ...
        'DisplayName', sprintf('mPFC- (n=%d)', numel(non_learner_ids)));
end

% Mark learning day
xline(0, '--k', 'Learning Day', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');

% Formatting
xlabel('Days Relative to Learning', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Reward Rate', 'FontSize', 13, 'FontWeight', 'bold');
title('Reward Rate Aligned to Learning Day', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

% Set y-axis limits (0 to 1 if reward rate is a proportion)
ylim([0 max([non_learner_mean;learner_mean],[],1)*1.2]);  

% Clean axes
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
box on;
hold off;

%% Plot bars that compares lick counts : ITI vs Static time vs Anticipatory licks

% Colors


metrics     = {'static_time_lick_counts','ITI_lick_counts','anticipatory_licks'};
metricNames = {'Static','ITI','Anticipatory'};

for a = 1:numel(behaviour_data)
    days_all = behaviour_data(a).recording_day;
    if isempty(days_all), continue; end

    % protocols present for this animal
    wf_all = string({days_all.workflow})';
    uW = unique(wf_all, 'stable');

    for wi = 1:numel(uW)
        target_workflow = uW(wi);

        if contains(target_workflow,'stim_wheel') % skip wheel protocol for now
            continue;
        end

        % --- filter days by protocol (your pattern) ---
        isValid   = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
        validDays = days_all(isValid);
        if isempty(validDays), continue; end

        % chronological order (if day field exists)
        dayOrder = 1:numel(validDays);
        dayLabels = string(dayOrder);
        if isfield(validDays,'day') && ~isempty([validDays.day])
            try
                if isdatetime([validDays.day])
                    dt = [validDays.day]';
                else
                    dt = datetime(vertcat({validDays.day})', 'InputFormat','yyyy-MM-dd');
                end
                [~, dayOrder] = sort(dt);
                dayLabels = string(cellstr(datestr(dt(dayOrder),'yyyy-mm-dd'))); %#ok<DATST>
            catch
                % keep default order if parsing fails
            end
        end

        K = numel(validDays);
        figure('Color','w','Name',sprintf('Animal %s — %s', behaviour_data(a).animal_id, target_workflow));
        tiledlayout(numel(metrics), 1, 'TileSpacing','compact', 'Padding','compact');

        for m = 1:numel(metrics)
            meanPlus  = nan(K,1); semPlus  = nan(K,1);
            meanMinus = nan(K,1); semMinus = nan(K,1);

            % per-day stats for this metric
            for ii = 1:K
                d = validDays(dayOrder(ii));

                if ~isfield(d,'cs_labels') || isempty(d.cs_labels), continue; end
                static_count = d.cs_labels(:);

                if ~isfield(d, metrics{m}) || isempty(d.(metrics{m})), continue; end
                vals = d.(metrics{m})(:);

                L = min(numel(vals), numel(static_count));   % align lengths safely
                vals = vals(1:L);  static_count = static_count(1:L);

                % CS+ stats
                vp = vals(static_count);
                if ~isempty(vp)
                    meanPlus(ii) = mean(vp, 'omitnan');
                    np = sum(~isnan(vp));
                    semPlus(ii)  = (np>1) * (std(vp,'omitnan')/sqrt(max(1,np)));
                    if np<=1, semPlus(ii) = NaN; end
                end

                % CS- stats
                vm = vals(~static_count);
                if ~isempty(vm)
                    meanMinus(ii) = mean(vm, 'omitnan');
                    nm = sum(~isnan(vm));
                    semMinus(ii)  = (nm>1) * (std(vm,'omitnan')/sqrt(max(1,nm)));
                    if nm<=1, semMinus(ii) = NaN; end
                end
            end

            % y-lims for this metric (within this protocol)
            yvals = [meanPlus+semPlus; meanPlus-semPlus; meanMinus+semMinus; meanMinus-semMinus];
            yvals = yvals(isfinite(yvals));
            if isempty(yvals), yvals = 0; end
            ymin = max(0, min(yvals)); ymax = max(yvals);
            if ymin==ymax, ymax = ymin + 1; end

            % ---- plot grouped bars with SEM ----
            nexttile; hold on;
            M   = [meanPlus,  meanMinus];   % columns: CS+, CS-
            SEM = [semPlus,   semMinus];

            b = bar(M, 'grouped');
            b(1).FaceColor = cCSplus;
            b(2).FaceColor = cCSminus;   % CS- is red

            % error bars on bar endpoints
            for k = 1:2
                xk = b(k).XEndPoints;
                errorbar(xk, M(:,k), SEM(:,k), 'k', 'linestyle','none', 'linewidth', 1);
            end

            xlim([0.5, K+0.5]);
            ylim([ymin, ymax]);
            xticks(1:K); xticklabels(dayLabels); xtickangle(45);
            ylabel('Licks (mean \pm SEM)');
            title(sprintf('%s — %s', strrep(target_workflow,'_',' '), metricNames{m}), 'Interpreter','none');
            grid on; box off;
            if m==1, legend({'CS+','CS-'}, 'Location','best'); end
        end

        % xlabel(tiledlayout(gcf), 'Day (within protocol, chronological)');
        sgtitle(sprintf('Animal %s — %s', behaviour_data(a).animal_id, strrep(target_workflow,'_',' ')), 'Interpreter','none');
    end
end

%% Plot the ITI lick rate aligned to learning day split by learners vs non-learners


%% Plot the average licks for : ITI , Static, Anticipatory seperated by learners vs non learners for CS+ and CS- seperated

% ==== CONFIG ====


metrics     = {'static_time_lick_counts','ITI_lick_counts','anticipatory_licks'};
metricNames = {'Static','ITI','Anticipatory'};

% Build protocol list across dataset (skip wheel protocols)
all_wf = string(horzcat(workflow_animal{:}))';
workflow_list = unique(all_wf(~contains(all_wf,'stim_wheel')), 'stable');

% ---- aggregator: for each protocol & metric & CS(+/-) & group → cell by dayIdx
% proto(wi).(metric).plus.L / plus.NL / minus.L / minus.NL = cell array {dayIdx} -> [values across animals]
proto = struct();
for wi = 1:numel(workflow_list)
    for m = 1:numel(metrics)
        proto(wi).(metrics{m}).plus.L  = {};
        proto(wi).(metrics{m}).plus.NL = {};
        proto(wi).(metrics{m}).minus.L  = {};
        proto(wi).(metrics{m}).minus.NL = {};
    end
end

% ==== GATHER PER-ANIMAL → APPEND TO GROUP CELLS ====
for a = 1:numel(behaviour_data)
    days_all = behaviour_data(a).recording_day;
    if isempty(days_all), continue; end
    aid = string(behaviour_data(a).animal_id);
    isLearner = ismember(aid, learner_ids);

    if isLearner % assign group tag for structure
        grpTag= 'L';
    else
        grpTag= 'NL';
    end


    for wi = 1:numel(workflow_list)
        target_workflow = workflow_list(wi);

        % Filter this animal's days by protocol
        isValid   = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
        validDays = days_all(isValid);
        if isempty(validDays), continue; end

        % Chronological order (if 'day' exists)
        order = 1:numel(validDays);
        if isfield(validDays,'day') && ~isempty([validDays.day])
            try
                if isdatetime([validDays.day]), dt = [validDays.day]';
                else, dt = datetime(vertcat({validDays.day})','InputFormat','yyyy-MM-dd');
                end
                [~, order] = sort(dt);
            catch
            end
        end

        K = numel(validDays);

        for m = 1:numel(metrics)
            met = metrics{m};

            for ii = 1:K
                d = validDays(order(ii));

                % CS labels
                if isfield(d,'cs_labels') && ~isempty(d.cs_labels)
                    static_count = d.cs_labels;
                else
                    continue
                end

                % Metric values
                if ~isfield(d, met) || isempty(d.(met)), continue; end
                vals = d.(met)(:);

                % Align lengths safely
                Lmin = min(numel(vals), numel(static_count));
                vals = vals(1:Lmin);
                static_count   = static_count(1:Lmin);

                % Per-day means (CS+ and CS−)
                vPlus  = vals(static_count);
                vMinus = vals(~static_count);

                mPlus  = mean(vPlus,  'omitnan');
                mMinus = mean(vMinus, 'omitnan');

                % Append to group cells (skip NaNs) — no helper function
                if ~isnan(mPlus)
                    C = proto(wi).(met).plus.(grpTag);
                    if numel(C) < ii, C{ii} = []; end        % pad up to ii
                    C{ii} = [C{ii}; mPlus];                  % append value
                    proto(wi).(met).plus.(grpTag) = C;
                end

                if ~isnan(mMinus)
                    C = proto(wi).(met).minus.(grpTag);
                    if numel(C) < ii, C{ii} = []; end        % pad up to ii
                    C{ii} = [C{ii}; mMinus];                 % append value
                    proto(wi).(met).minus.(grpTag) = C;
                end


            end
        end
    end
end


% ==== REDUCE (cell → mean & SEM per dayIdx) + PLOT ====
for wi = 1:numel(workflow_list)
    figure('Color','w','Name',sprintf('Group averages — %s', workflow_list(wi)));
    t = tiledlayout(numel(metrics), 2, 'TileSpacing','compact','Padding','compact');

    for m = 1:numel(metrics)
        met = metrics{m};
        % ---------- CS+ (Learners vs Non-learners) ----------
        C_L  = proto(wi).(met).plus.L;    % {dayIdx} -> vector of animal means
        C_NL = proto(wi).(met).plus.NL;

        % mean ± SEM per dayIdx
        if isempty(C_L),  L_mu_plus = [];  L_se_plus = [];
        else
            L_mu_plus = cellfun(@(x) mean(x,'omitnan'), C_L)';                         % 1×K
            L_se_plus = cellfun(@(x) std(x,'omitnan') ./ max(1,sqrt(sum(~isnan(x)))), C_L)';
        end
        if isempty(C_NL), NL_mu_plus = []; NL_se_plus = [];
        else
            NL_mu_plus = cellfun(@(x) mean(x,'omitnan'), C_NL)';
            NL_se_plus = cellfun(@(x) std(x,'omitnan') ./ max(1,sqrt(sum(~isnan(x)))), C_NL)';
        end

        Kp = max([numel(L_mu_plus), numel(NL_mu_plus), 1]);

        % y-limits for CS+ (clamp at 0)
        yvalsP = [L_mu_plus+L_se_plus; L_mu_plus-L_se_plus; NL_mu_plus+NL_se_plus; NL_mu_plus-NL_se_plus];
        yvalsP = yvalsP(isfinite(yvalsP)); if isempty(yvalsP), yvalsP = 0; end
        yminP = min(0, min(yvalsP)); ymaxP = max(yvalsP); if yminP==ymaxP, ymaxP=yminP+1; end

        % pad mean/sem to Kp
        L_mu_plus_pad  = nan(Kp,1);  L_mu_plus_pad (1:numel(L_mu_plus))  = L_mu_plus(:);
        NL_mu_plus_pad = nan(Kp,1);  NL_mu_plus_pad(1:numel(NL_mu_plus)) = NL_mu_plus(:);
        L_se_plus_pad  = nan(Kp,1);  L_se_plus_pad (1:numel(L_se_plus))  = L_se_plus(:);
        NL_se_plus_pad = nan(Kp,1);  NL_se_plus_pad(1:numel(NL_se_plus)) = NL_se_plus(:);

        nexttile; hold on;
        off = 0.18; dotSize = 28; meanSize = 60; lw = 1.8;

        % scatter each animal value per dayIdx (slight horizontal offset)
        for d = 1:Kp
            % learners
            if d <= numel(C_L) && ~isempty(C_L{d})
                v = C_L{d}; v = v(~isnan(v));
                if ~isempty(v)
                    scatter(d - off, v, dotSize, 'MarkerFaceColor', cLearner, ...
                        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
                end
            end

            % non-learners
            if d <= numel(C_NL) && ~isempty(C_NL{d})
                v = C_NL{d}; v = v(~isnan(v));
                if ~isempty(v)
                    scatter(d + off, v, dotSize, 'MarkerFaceColor', cNonLearner, ...
                        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
                end
            end
        end

        % overlay group mean ± SEM as bold dots with error bars
        errorbar((1:Kp)-off, L_mu_plus_pad,  L_se_plus_pad,  'o', ...
            'Color', cLearner, 'MarkerFaceColor', cLearner, 'MarkerSize', meanSize/10, 'LineWidth', lw, ...
            'DisplayName','Learners (mean±SEM)');
        errorbar((1:Kp)+off, NL_mu_plus_pad, NL_se_plus_pad, 'o', ...
            'Color', cNonLearner, 'MarkerFaceColor', cNonLearner, 'MarkerSize', meanSize/10, 'LineWidth', lw, ...
            'DisplayName','Non-learners (mean±SEM)');

        xlim([0.5, Kp+0.5]); ylim([yminP, ymaxP]);
        xticks(1:Kp); xticklabels(string(1:Kp)); % day index within protocol
        ylabel([metricNames{m} ' (CS+)']);
        % title(strrep(workflow_list(wi),'_',' '), 'Interpreter','none');
        grid on; box off;
        if m==1, legend('Location','best'); end

        % ---------- CS− (Learners vs Non-learners) ----------
        C_Lm  = proto(wi).(met).minus.L;
        C_NLm = proto(wi).(met).minus.NL;

        if isempty(C_Lm),  L_mu_minus = [];  L_se_minus = [];
        else
            L_mu_minus = cellfun(@(x) mean(x,'omitnan'), C_Lm)';
            L_se_minus = cellfun(@(x) std(x,'omitnan') ./ max(1,sqrt(sum(~isnan(x)))), C_Lm)';
        end
        if isempty(C_NLm), NL_mu_minus = []; NL_se_minus = [];
        else
            NL_mu_minus = cellfun(@(x) mean(x,'omitnan'), C_NLm)';
            NL_se_minus = cellfun(@(x) std(x,'omitnan') ./ max(1,sqrt(sum(~isnan(x)))), C_NLm)';
        end

        Kn = max([numel(L_mu_minus), numel(NL_mu_minus), 1]);

        yvalsN = [L_mu_minus+L_se_minus; L_mu_minus-L_se_minus; NL_mu_minus+NL_se_minus; NL_mu_minus-NL_se_minus];
        yvalsN = yvalsN(isfinite(yvalsN)); if isempty(yvalsN), yvalsN = 0; end
        yminN = min(0, min(yvalsN)); ymaxN = max(yvalsN); if yminN==ymaxN, ymaxN=yminN+1; end

        L_mu_minus_pad  = nan(Kn,1);  L_mu_minus_pad (1:numel(L_mu_minus))  = L_mu_minus(:);
        NL_mu_minus_pad = nan(Kn,1);  NL_mu_minus_pad(1:numel(NL_mu_minus)) = NL_mu_minus(:);
        L_se_minus_pad  = nan(Kn,1);  L_se_minus_pad (1:numel(L_se_minus))  = L_se_minus(:);
        NL_se_minus_pad = nan(Kn,1);  NL_se_minus_pad(1:numel(NL_se_minus)) = NL_se_minus(:);

        nexttile; hold on;

        for d = 1:Kn
            if d <= numel(C_Lm) && ~isempty(C_Lm{d})
                v = C_Lm{d}; v = v(~isnan(v));
                if ~isempty(v)
                    scatter(d - off, v, dotSize, 'MarkerFaceColor', cLearner, ...
                        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
                end
            end
            if d <= numel(C_NLm) && ~isempty(C_NLm{d})
                v = C_NLm{d}; v = v(~isnan(v));
                if ~isempty(v)
                    scatter(d + off, v, dotSize, 'MarkerFaceColor', cNonLearner, ...
                        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
                end
            end
        end

        errorbar((1:Kn)-off, L_mu_minus_pad,  L_se_minus_pad,  'o', ...
            'Color', cLearner, 'MarkerFaceColor', cLearner, 'MarkerSize', meanSize/10, 'LineWidth', lw, ...
            'DisplayName','Learners (mean±SEM)');
        errorbar((1:Kn)+off, NL_mu_minus_pad, NL_se_minus_pad, 'o', ...
            'Color', cNonLearner, 'MarkerFaceColor', cNonLearner, 'MarkerSize', meanSize/10, 'LineWidth', lw, ...
            'DisplayName','Non-learners (mean±SEM)');

        xlim([0.5, Kn+0.5]); ylim([yminN, ymaxN]);
        xticks(1:Kn); xticklabels(string(1:Kn));
        ylabel([metricNames{m} ' (CS−)']);
        grid on; box off;
    end

    xlabel(t, 'Day index within protocol (chronological)');
    sgtitle(sprintf('Group averages — %s', strrep(workflow_list(wi),'_',' ')), 'Interpreter','none');
end


%% Plot only non-learners static protocol for combined CS+ and CS-

% ==== REDUCE (cell → mean & SEM per dayIdx) + PLOT — STATIC NON-LEARNERS ONLY ====
for wi = 1:numel(workflow_list)
    % only keep static (non-big_stim) protocols
    isStaticProto = contains(workflow_list(wi),'static') && ~contains(workflow_list(wi),'big_stim');
    if ~isStaticProto, continue; end

    figure('Color','w','Name',sprintf('Static (NL, CS± combined) — %s', workflow_list(wi)));
    t = tiledlayout(numel(metrics), 1, 'TileSpacing','compact','Padding','compact');

    for m = 1:numel(metrics)
        met = metrics{m};

        % Non-learner cells: CS+ and CS− per day index
        C_NLp = proto(wi).(met).plus.NL;   % {dayIdx} -> values across animals (CS+)
        C_NLm = proto(wi).(met).minus.NL;  % {dayIdx} -> values across animals (CS−)

        % Combine CS± per day index
        Kc = max(numel(C_NLp), numel(C_NLm));
        C_comb = cell(1, Kc);
        for d = 1:Kc
            vals = [];
            if d <= numel(C_NLp) && ~isempty(C_NLp{d}), vals = [vals; C_NLp{d}(:)]; end
            if d <= numel(C_NLm) && ~isempty(C_NLm{d}), vals = [vals; C_NLm{d}(:)]; end
            C_comb{d} = vals;
        end

        % Mean ± SEM per day index (NaN-safe)
        if isempty(C_comb) || all(cellfun(@isempty,C_comb))
            warning('No non-learner values for %s in %s; skipping.', met, workflow_list(wi));
            nexttile; axis off; title([metricNames{m} ' — no data']); continue;
        end
        mu = cellfun(@(x) mean(x,'omitnan'), C_comb)';                    % 1×K
        se = cellfun(@(x) std(x,'omitnan') ./ max(1,sqrt(sum(~isnan(x)))), C_comb)';

        K  = max(numel(mu),1);
        mu_pad = nan(K,1); mu_pad(1:numel(mu)) = mu;
        se_pad = nan(K,1); se_pad(1:numel(se)) = se;

        % y-limits (counts → clamp at 0)
        yvals = [mu+se; mu-se]; yvals = yvals(isfinite(yvals));
        if isempty(yvals), yvals = 0; end
        ymin = min(0, min(yvals)); ymax = max(yvals); if ymin==ymax, ymax = ymin+1; end

        % ---- plot dots (animals) + bold mean±SEM ----
        nexttile; hold on;
        dotSize = 28; meanSize = 60; lw = 1.8;
        jitter = 0.06;

        % individual animal dots per day
        for d = 1:K
            if d <= numel(C_comb) && ~isempty(C_comb{d})
                v = C_comb{d}; v = v(~isnan(v));
                if ~isempty(v)
                    xj = d + (rand(size(v)) - 0.5) * jitter;   % light jitter
                    scatter(xj, v, dotSize, 'MarkerFaceColor', cNonLearner, ...
                        'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
                end
            end
        end

        % group mean ± SEM (bold)
        errorbar(1:K, mu_pad, se_pad, 'o', ...
            'Color', cNonLearner, 'MarkerFaceColor', cNonLearner, ...
            'MarkerSize', meanSize/10, 'LineWidth', lw, ...
            'DisplayName','Non-learners (CS± combined)');

        xlim([0.5, K+0.5]); ylim([ymin, ymax]);
        xticks(1:K); xticklabels(string(1:K));  % day index within protocol
        ylabel(metricNames{m});
        title(sprintf('%s — Non-learners (CS± combined)', strrep(workflow_list(wi),'_',' ')), 'Interpreter','none');
        grid on; box off;
        if m==1, legend('Location','best'); end
    end

    xlabel(t, 'Day index within protocol (chronological)');
    sgtitle(sprintf('Static protocol — Non-learners (CS± combined): %s', strrep(workflow_list(wi),'_',' ')), ...
        'Interpreter','none');
end


%% Plot a modulation index based on two licking rates (e.g. static rate vs ITI rates) split by learners vs non-learners


% ---------- CONFIG ----------
learners_group_ID = ["HA005","HA008","HA010","HA012"];   % learners

dotSize  = 28; meanSize = 60; lw = 1.8; jitter = 0.06;

% ---------- COLLECT MI PER DAY (STATIC, non-big_stim), by group ----------
MI_days_L  = {};   % {dayIdx} -> [MI across learner animals]
MI_days_NL = {};   % {dayIdx} -> [MI across non-learner animals]

for a = 1:numel(behaviour_data)
    aid  = string(behaviour_data(a).animal_id);
    days = behaviour_data(a).recording_day;
    if isempty(days), continue; end

    % keep static and not big_stim
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), days);
    dlist = days(is_stage_2);
    if isempty(dlist), continue; end

    % take order of dayts
    ord = 1:numel(dlist);

    isLearner = ismember(aid, learners_group_ID);

    % compute per-day MI for this animal
    for k = 1:numel(ord)

        D = dlist(ord(k));

        % counts: STATIC vs ITI
        cs_plus_mask= D.cs_labels; % logical mask

        cs = D.static_time_lick_counts(~cs_plus_mask);      % static counts per trial
        ci = D.ITI_lick_counts(~cs_plus_mask);        % ITI counts per trial

        if isempty(cs) || isempty(ci) || nansum(cs)==0 || nansum(ci)==0 % add the sum clause given HA008 has no licks on day 5 static
            continue;
        end

        % --- durations per trial (fallbacks) ---
        if isfield(D,'all_static_times') && ~isempty(D.all_static_times)
            ts = D.all_static_times;

        end

        if isfield(D,'ITI_actual_duration') && ~isempty(D.ITI_actual_duration)
            ti = D.ITI_actual_duration;
        end

        % --- rates and MI ---
        rS = nansum(cs) / max(eps, nansum(ts));   % licks/s during static
        rI = nansum(ci) / max(eps, nansum(ti));   % licks/s during ITI
        MI = (rS - rI) / max(eps, (rS + rI));     % in [-1,1]

        % append to the appropriate group/day cell
        if isLearner
            if numel(MI_days_L)  < k, MI_days_L{k}  = []; end
            MI_days_L{k}  = [MI_days_L{k};  MI];
        else
            if numel(MI_days_NL) < k, MI_days_NL{k} = []; end
            MI_days_NL{k} = [MI_days_NL{k}; MI];
        end
    end
end

% ---------- REDUCE: mean ± SEM per day index, for both groups ----------
K = max(numel(MI_days_L), numel(MI_days_NL)); if K==0, warning('No static data.'); return; end
muL = nan(K,1); seL = nan(K,1); muNL = nan(K,1); seNL = nan(K,1);

for d = 1:K
    % learners
    if d <= numel(MI_days_L) && ~isempty(MI_days_L{d})
        v = MI_days_L{d}; v = v(~isnan(v));
        muL(d) = mean(v);  n = numel(v); seL(d) = (n>1)*std(v)/sqrt(max(1,n)); if n<=1, seL(d)=NaN; end
    end
    % non-learners
    if d <= numel(MI_days_NL) && ~isempty(MI_days_NL{d})
        v = MI_days_NL{d}; v = v(~isnan(v));
        muNL(d) = mean(v); n = numel(v); seNL(d) = (n>1)*std(v)/sqrt(max(1,n)); if n<=1, seNL(d)=NaN; end
    end
end

% ---------- PLOT: both groups on same axes ----------
figure('Color','w','Name','Static protocol (CS± pooled): MI vs day — Learners vs Non-learners');
hold on; off = 0.12;

% per-animal dots
for d = 1:K
    % learners
    if d <= numel(MI_days_L) && ~isempty(MI_days_L{d})
        v = MI_days_L{d}; v = v(~isnan(v));
        if ~isempty(v)
            xj = (d - off) + (rand(size(v)) - 0.5)*jitter;
            scatter(xj, v, dotSize, 'MarkerFaceColor', cLearner, ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
        end
    end
    % non-learners
    if d <= numel(MI_days_NL) && ~isempty(MI_days_NL{d})
        v = MI_days_NL{d}; v = v(~isnan(v));
        if ~isempty(v)
            xj = (d + off) + (rand(size(v)) - 0.5)*jitter;
            scatter(xj, v, dotSize, 'MarkerFaceColor', cNonLearner, ...
                'MarkerEdgeColor','none', 'MarkerFaceAlpha', 0.55, 'HandleVisibility','off');
        end
    end
end

% group means ± SEM
errorbar((1:K)-off, muL,  seL,  'o-', 'Color', cLearner,    'MarkerFaceColor', cLearner,    ...
    'MarkerSize', meanSize/10, 'LineWidth', lw, 'DisplayName','Learners (mean±SEM)');
errorbar((1:K)+off, muNL, seNL, 'o-', 'Color', cNonLearner, 'MarkerFaceColor', cNonLearner, ...
    'MarkerSize', meanSize/10, 'LineWidth', lw, 'DisplayName','Non-learners (mean±SEM)');

xlim([0.5, K+0.5]); ylim([-1, 1]);  % MI bounded
xticks(1:K); xticklabels(string(1:K));
yline(0,'k--','HandleVisibility','off');

xlabel('Day index within static protocol');
ylabel('Modulation index  ( (r_{static}-r_{ITI}) / (r_{static}+r_{ITI}) )');
title('Static CS- — Learners vs Non-learners');
grid on; box off; legend('Location','best');



%% Q: Do the learners learn from the first static trials? Figure 1 (Individual animals)

% Plots individual animal anticipatory licks overlayed with RT for learners for the first 20 trials in first static day

% ---------- CONFIG ----------
learners_group_ID = ["HA005","HA008","HA010","HA012"];  % your learner IDs
Nshow = 20;  % first N trials to visualize

% colors
cAnt = [0.90 0.55 0.10];   % anticipatory (orange)
cSta = [0.20 0.55 0.90];   % static (blue)
cRT  = [0.20 0.20 0.20];   % reaction time (dark)

for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    if ~ismember(aid, learners_group_ID), continue; end

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % ---- pick FIRST "static" (and not big_stim) day ----
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static'), rd);
    if ~any(is_stage_2)
        fprintf('Animal %s: no static day found; skipping.\n', aid);
        continue;
    end

    cand = find(is_stage_2);

    firstIdx = cand(end);

    D = rd(firstIdx);

    % ---- pull series (anticipatory, static) ----
    if ~isfield(D,'anticipatory_licks') || isempty(D.anticipatory_licks),  continue; end
    if ~isfield(D,'static_time_lick_counts') || isempty(D.static_time_lick_counts), continue; end

    cs_plus_mask= D.cs_labels;

    % grab the relevant metrics as lick rates
    anticip = D.anticipatory_licks(cs_plus_mask)/1; % given its a 1s window
    staticC = D.static_time_lick_counts(cs_plus_mask)./ D.all_static_times (cs_plus_mask);
    rt= D.all_stim_diff_from_optimal_reward(cs_plus_mask);


    % ---- trim to first N trials with matched length ----
    nAvail = min([numel(anticip), numel(staticC), max(1, numel(rt))]);
    n = min(Nshow, nAvail);
    anticip = anticip(1:n);
    staticC = staticC(1:n);
    if ~isempty(rt), rt = rt(1:n); end
    trials = 1:n;

    % ---- plot ----
    figure('Color','w','Name',sprintf('First static day — %s', aid));
    hold on;
    yyaxis left
    plot(trials, anticip, '-o', 'Color', cAnt, 'MarkerFaceColor', cAnt, 'LineWidth', 1.8);
    plot(trials, staticC, '-s', 'Color', cSta, 'MarkerFaceColor', cSta, 'LineWidth', 1.8);
    ylabel('Licks (count)');

    yyaxis right
    if ~isempty(rt)
        plot(trials, rt, '-^', 'Color', cRT, 'MarkerFaceColor', cRT, 'LineWidth', 1.8);
        ylabel('Reaction time (s)');
    else
        ylabel('Reaction time (s)');  % empty, but keeps axis labeled
    end

    xlabel('Trial # (first 20)');
    title(sprintf('Animal %s — First static day (CS± pooled)', aid), 'Interpreter','none');
    grid on; box off;

    % legend (dynamic)
    L = {'Anticipatory', 'Static'};
    if ~isempty(rt), L{end+1} = 'Reaction time'; end
    legend(L, 'Location','best');
end


%% Plot Q: Do the learners learn from the first static trials? Figure 2 (Average)

% Plot the mean anticipatory and static lick rate for the first N trials
% for both learners and non-learners

%% Collect data for BOTH learners and non-learners in one loop
Nshow = 50; % number of trials to show

% Color definitions
cAnt = [0.8500 0.3250 0.0980]; % Orange for anticipatory
cSta = [0 0.4470 0.7410];      % Blue for static


% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learners_group_ID));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate cell arrays
A_list_NL = cell(1, n_non_learners); % Non-learners anticipatory
S_list_NL = cell(1, n_non_learners); % Non-learners static
A_list_L = cell(1, n_learners);      % Learners anticipatory
S_list_L = cell(1, n_learners);      % Learners static

% Initialize counters for indexing
idx_L = 0;
idx_NL = 0;

% ===== Single loop to collect both groups =====
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % first "static" day (skip big_stim)
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);
    if ~any(is_stage_2), continue; end

    cand = find(is_stage_2);
    firstIdx = cand(0);
    D = rd(firstIdx);

    % required fields
    if ~isfield(D,'anticipatory_licks') || isempty(D.anticipatory_licks), continue; end
    if ~isfield(D,'static_time_lick_counts') || isempty(D.static_time_lick_counts), continue; end
    if ~isfield(D,'cs_labels') || isempty(D.cs_labels), continue; end

    cs_plus_mask = D.cs_labels;
    anticip = D.anticipatory_licks(cs_plus_mask) / 1;
    staticC = D.static_time_lick_counts(cs_plus_mask) ./ D.all_static_times(cs_plus_mask);

    % trim to first Nshow trials
    n = min(Nshow, numel(anticip));
    if n < 1, continue; end

    anticip = anticip(1:n);
    staticC = staticC(1:n);

    % Conditional: assign to learner or non-learner list
    if ismember(aid, learners_group_ID)
        % Learner
        idx_L = idx_L + 1;
        A_list_L{idx_L} = anticip;
        S_list_L{idx_L} = staticC;
    else
        % Non-learner
        idx_NL = idx_NL + 1;
        A_list_NL{idx_NL} = anticip;
        S_list_NL{idx_NL} = staticC;
    end
end

% ===== Stack into matrices (NON-LEARNERS) =====
nNL = numel(A_list_NL);
A_mat_NL_s1 = nan(Nshow, nNL);
S_mat_NL_s1 = nan(Nshow, nNL);

for i = 1:nNL
    va = A_list_NL{i};
    vs = S_list_NL{i};
    A_mat_NL_s1(1:numel(va), i) = va;
    S_mat_NL_s1(1:numel(vs), i) = vs;
end

% ===== Stack into matrices (LEARNERS) =====
nL = numel(A_list_L);
A_mat_L_s1 = nan(Nshow, nL);
S_mat_L_s1 = nan(Nshow, nL);

for i = 1:nL
    va = A_list_L{i};
    vs = S_list_L{i};
    A_mat_L_s1(1:numel(va), i) = va;
    S_mat_L_s1(1:numel(vs), i) = vs;
end

% ===== Compute mean ± SEM =====
t = (1:Nshow)';

% Non-learners
mA_NL = mean(A_mat_NL_s1, 2, 'omitnan');
nA_NL = sum(~isnan(A_mat_NL_s1), 2);
sA_NL = std(A_mat_NL_s1, 0, 2, 'omitnan') ./ max(1, sqrt(nA_NL));

mS_NL = mean(S_mat_NL_s1, 2, 'omitnan');
nS_NL = sum(~isnan(S_mat_NL_s1), 2);
sS_NL = std(S_mat_NL_s1, 0, 2, 'omitnan') ./ max(1, sqrt(nS_NL));

% Learners
mA_L = mean(A_mat_L_s1, 2, 'omitnan');
nA_L = sum(~isnan(A_mat_L_s1), 2);
sA_L = std(A_mat_L_s1, 0, 2, 'omitnan') ./ max(1, sqrt(nA_L));

mS_L = mean(S_mat_L_s1, 2, 'omitnan');
nS_L = sum(~isnan(S_mat_L_s1), 2);
sS_L = std(S_mat_L_s1, 0, 2, 'omitnan') ./ max(1, sqrt(nS_L));

% ===== Determine shared y-axis limits =====
all_data = [A_mat_NL_s1(:); S_mat_NL_s1(:); A_mat_L_s2(:); S_mat_L_s1(:)];
y_min = min(all_data, [], 'omitnan');
y_max = max(all_data, [], 'omitnan');
y_margin = (y_max - y_min) * 0.1; % 10% margin
ylim_shared = [y_min - y_margin, y_max + y_margin];

% ===== Plot side-by-side =====
figure('Color', 'w', 'Position', [100 100 1200 500], ...
    'Name', 'Lick Rates: Learners vs Non-Learners — First Static Day');

% --- Non-Learners Panel ---
subplot(1, 2, 1);
hold on;

xf = [t; flipud(t)];

% Anticipatory
yfA_NL = [mA_NL + sA_NL; flipud(mA_NL - sA_NL)];
fill(xf, yfA_NL, cAnt, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mA_NL, '-', 'Color', cAnt, 'LineWidth', 2, 'DisplayName', 'Anticipatory');

% Static
yfS_NL = [mS_NL + sS_NL; flipud(mS_NL - sS_NL)];
fill(xf, yfS_NL, cSta, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mS_NL, '-', 'Color', cSta, 'LineWidth', 2, 'DisplayName', 'Static');

ylabel('Lick rate (licks/s)');
xlabel('Trial # (first 20)');
title(sprintf('Non-Learners (n=%d)', nNL));
ylim(ylim_shared);
grid on;
box off;
legend('Location', 'best');

% --- Learners Panel ---
subplot(1, 2, 2);
hold on;

% Anticipatory
yfA_L = [mA_L + sA_L; flipud(mA_L - sA_L)];
fill(xf, yfA_L, cAnt, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mA_L, '-', 'Color', cAnt, 'LineWidth', 2, 'DisplayName', 'Anticipatory');

% Static
yfS_L = [mS_L + sS_L; flipud(mS_L - sS_L)];
fill(xf, yfS_L, cSta, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mS_L, '-', 'Color', cSta, 'LineWidth', 2, 'DisplayName', 'Static');

ylabel('Lick rate (licks/s)');
xlabel('Trial # (first 20)');
title(sprintf('Learners (n=%d)', nL));
ylim(ylim_shared);
grid on;
box off;
legend('Location', 'best');

sgtitle('Last day of stage 1 (CS+): Anticipatory vs Static Lick Rates');

fprintf('\n===== Summary Statistics =====\n');
fprintf('Non-Learners (n=%d):\n', nNL);
fprintf('  Anticipatory: %.2f ± %.2f licks/s\n', mean(mA_NL, 'omitnan'), mean(sA_NL, 'omitnan'));
fprintf('  Static: %.2f ± %.2f licks/s\n', mean(mS_NL, 'omitnan'), mean(sS_NL, 'omitnan'));
fprintf('\nLearners (n=%d):\n', nL);
fprintf('  Anticipatory: %.2f ± %.2f licks/s\n', mean(mA_L, 'omitnan'), mean(sA_L, 'omitnan'));
fprintf('  Static: %.2f ± %.2f licks/s\n', mean(mS_L, 'omitnan'), mean(sS_L, 'omitnan'));

%% Saves both last stage 1 and first stage 2 simultanously - relevant for next plotting for strategy shifts

% Collect data ensuring matched animals across stages

Nshow = 20;

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learners_group_ID));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate with animal IDs to track matching
stage1_data_L = struct('animal_id', {}, 'anticip', {}, 'static', {});
stage2_data_L = struct('animal_id', {}, 'anticip', {}, 'static', {});
stage1_data_NL = struct('animal_id', {}, 'anticip', {}, 'static', {});
stage2_data_NL = struct('animal_id', {}, 'anticip', {}, 'static', {});

idx_L = 0;
idx_NL = 0;

% Single loop collecting BOTH stages per animal
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % ===== Find Stage 1 (last classical day) =====
    is_stage_1 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);

    stage1_anticip = [];
    stage1_static = [];
    if any(is_stage_1)
        cand_classical = find(is_stage_1);
        lastClassicalIdx = cand_classical(end); % LAST classical day
        D1 = rd(lastClassicalIdx);

        % Check required fields
        if isfield(D1,'anticipatory_licks') && ~isempty(D1.anticipatory_licks) && ...
                isfield(D1,'static_time_lick_counts') && ~isempty(D1.static_time_lick_counts) && ...
                isfield(D1,'cs_labels') && ~isempty(D1.cs_labels)

            cs_plus_mask1 = D1.cs_labels;
            anticip1 = D1.anticipatory_licks(cs_plus_mask1) / 1;
            staticC1 = D1.static_time_lick_counts(cs_plus_mask1) ./ D1.all_static_times(cs_plus_mask1);

            % Trim to first Nshow trials
            n1 = min(Nshow, numel(anticip1));
            if n1 >= 1
                stage1_anticip = anticip1(1:n1);
                stage1_static = staticC1(1:n1);
            end
        end
    end

    % ===== Find Stage 2 (first static day) =====
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);

    stage2_anticip = [];
    stage2_static = [];
    if any(is_stage_2)
        cand_static = find(is_stage_2);
        firstStaticIdx = cand_static(end); % FIRST static day
        D2 = rd(firstStaticIdx);

        % Check required fields
        if isfield(D2,'anticipatory_licks') && ~isempty(D2.anticipatory_licks) && ...
                isfield(D2,'static_time_lick_counts') && ~isempty(D2.static_time_lick_counts) && ...
                isfield(D2,'cs_labels') && ~isempty(D2.cs_labels)

            cs_plus_mask2 = D2.cs_labels;
            anticip2 = D2.anticipatory_licks(cs_plus_mask2) / 1;
            staticC2 = D2.static_time_lick_counts(cs_plus_mask2) ./ D2.all_static_times(cs_plus_mask2);

            % Trim to first Nshow trials
            n2 = min(Nshow, numel(anticip2));
            if n2 >= 1
                stage2_anticip = anticip2(1:n2);
                stage2_static = staticC2(1:n2);
            end
        end
    end

    % ===== Only include animals with BOTH stages =====
    if isempty(stage1_anticip) || isempty(stage2_anticip)
        fprintf('Skipping animal %s: missing stage 1 or stage 2 data\n', aid);
        continue;
    end

    % ===== Assign to learner or non-learner =====
    if ismember(aid, learners_group_ID)
        idx_L = idx_L + 1;
        stage1_data_L(idx_L).animal_id = aid;
        stage1_data_L(idx_L).anticip = stage1_anticip;
        stage1_data_L(idx_L).static = stage1_static;

        stage2_data_L(idx_L).animal_id = aid;
        stage2_data_L(idx_L).anticip = stage2_anticip;
        stage2_data_L(idx_L).static = stage2_static;
    else
        idx_NL = idx_NL + 1;
        stage1_data_NL(idx_NL).animal_id = aid;
        stage1_data_NL(idx_NL).anticip = stage1_anticip;
        stage1_data_NL(idx_NL).static = stage1_static;

        stage2_data_NL(idx_NL).animal_id = aid;
        stage2_data_NL(idx_NL).anticip = stage2_anticip;
        stage2_data_NL(idx_NL).static = stage2_static;
    end
end

fprintf('\n===== Data Collection Summary =====\n');
fprintf('Learners with both stages: %d\n', idx_L);
fprintf('Non-learners with both stages: %d\n', idx_NL);

% ===== Stack into matrices =====

% Stage 1 - Non-learners
nNL = length(stage1_data_NL);
A_mat_NL_s1 = nan(Nshow, nNL);
S_mat_NL_s1 = nan(Nshow, nNL);
for i = 1:nNL
    va = stage1_data_NL(i).anticip;
    vs = stage1_data_NL(i).static;
    A_mat_NL_s1(1:numel(va), i) = va;
    S_mat_NL_s1(1:numel(vs), i) = vs;
end

% Stage 2 - Non-learners
A_mat_NL = nan(Nshow, nNL);
S_mat_NL = nan(Nshow, nNL);
for i = 1:nNL
    va = stage2_data_NL(i).anticip;
    vs = stage2_data_NL(i).static;
    A_mat_NL(1:numel(va), i) = va;
    S_mat_NL(1:numel(vs), i) = vs;
end

% Stage 1 - Learners
nL = length(stage1_data_L);
A_mat_L_s1 = nan(Nshow, nL);
S_mat_L_s1 = nan(Nshow, nL);
for i = 1:nL
    va = stage1_data_L(i).anticip;
    vs = stage1_data_L(i).static;
    A_mat_L_s1(1:numel(va), i) = va;
    S_mat_L_s1(1:numel(vs), i) = vs;
end

% Stage 2 - Learners
A_mat_L = nan(Nshow, nL);
S_mat_L = nan(Nshow, nL);
for i = 1:nL
    va = stage2_data_L(i).anticip;
    vs = stage2_data_L(i).static;
    A_mat_L(1:numel(va), i) = va;
    S_mat_L(1:numel(vs), i) = vs;
end

% Verify matching
fprintf('\n===== Verification =====\n');
fprintf('Stage 1 matrices: NL=%dx%d, L=%dx%d\n', size(A_mat_NL_s1,1), size(A_mat_NL_s1,2), size(A_mat_L_s1,1), size(A_mat_L_s1,2));
fprintf('Stage 2 matrices: NL=%dx%d, L=%dx%d\n', size(A_mat_NL,1), size(A_mat_NL,2), size(A_mat_L,1), size(A_mat_L,2));

% Print which animals are included
fprintf('\nLearners included: ');
for i = 1:nL
    fprintf('%s ', stage1_data_L(i).animal_id);
end
fprintf('\n');

fprintf('Non-learners included: ');
for i = 1:nNL
    fprintf('%s ', stage1_data_NL(i).animal_id);
end
fprintf('\n');


%% Plot the average 0-1s from stim oneset for the last 2 days of stage 1 and the individual trials for first day stage 2 

% Parameters
Nshow = 30;  % Number of trials to show from stage 2
lick_window = [0, 1];  % Window from stim onset (seconds)

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learners_group_ID));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate structures for last 2 days of stage 1 + first day of stage 2
stage1_data_L = struct('animal_id', {}, 'day_minus2', {}, 'day_minus1', {});
stage2_data_L = struct('animal_id', {}, 'first_day_trials', {});
stage1_data_NL = struct('animal_id', {}, 'day_minus2', {}, 'day_minus1', {});
stage2_data_NL = struct('animal_id', {}, 'first_day_trials', {});

idx_L = 0;
idx_NL = 0;

% Helper function to calculate lick count in window
calc_lick_count = @(lick_times, ref_time, window) ...
    sum(lick_times >= (ref_time + window(1)) & lick_times < (ref_time + window(2)));

% Single loop collecting BOTH stages per animal
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % ===== Find Stage 1 (last 2 classical days) =====
    is_stage_1 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    
    stage1_day_minus2 = [];
    stage1_day_minus1 = [];
    
    if sum(is_stage_1) >= 2  % Need at least 2 days
        cand_classical = find(is_stage_1);
        lastIdx = cand_classical(end);
        secondLastIdx = cand_classical(end-1);
        
        % Process day -2 (second to last)
        D_minus2 = rd(secondLastIdx);
        if isfield(D_minus2, 'cs_labels') && isfield(D_minus2, 'lick_event_times') && ...
           isfield(D_minus2, 'right_stim_on_times')
            
            
            cs_plus_mask = D_minus2.cs_labels == 1;
            stim_on_times = D_minus2.right_stim_on_times;
            lick_times = D_minus2.lick_event_times;
            
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), ...
                                   stim_on_times);
            stage1_day_minus2 = mean(lick_counts, 'omitnan');
        end
        
        % Process day -1 (last day)
        D_minus1 = rd(lastIdx);
        if isfield(D_minus1, 'cs_labels') && isfield(D_minus1, 'lick_event_times') && ...
           isfield(D_minus1, 'right_stim_on_times')
            
            cs_plus_mask = D_minus1.cs_labels == 1;
            stim_on_times = D_minus1.right_stim_on_times;
            lick_times = D_minus1.lick_event_times;
            
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), ...
                                   stim_on_times);
            stage1_day_minus1 = mean(lick_counts, 'omitnan');
        end
    end

    % ===== Find Stage 2 (first static day, first N trials) =====
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ...
                          ~contains(d.workflow,'big_stim'), rd);
    
    stage2_first_trials = [];
    
    if any(is_stage_2)
        cand_static = find(is_stage_2);
        firstStaticIdx = cand_static(1);  % FIRST static day
        D2 = rd(firstStaticIdx);
        
        if isfield(D2, 'cs_labels') && isfield(D2, 'lick_event_times') && ...
           isfield(D2, 'right_stim_on_times')
            
            cs_plus_mask = D2.cs_labels == 1;
            stim_on_times = D2.right_stim_on_times;
            lick_times = D2.lick_event_times;
            
            % Get first Nshow trials
            n_trials = min(Nshow, length(stim_on_times));
            if n_trials >= 1
                stage2_first_trials = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), ...
                                               stim_on_times(1:n_trials));
            end
        end
    end

    % ===== Only include animals with ALL data =====
    if isempty(stage1_day_minus2) || isempty(stage1_day_minus1) || isempty(stage2_first_trials)
        fprintf('Skipping animal %s: missing data\n', aid);
        continue;
    end

    % ===== Assign to learner or non-learner =====
    if ismember(aid, learners_group_ID)
        idx_L = idx_L + 1;
        stage1_data_L(idx_L).animal_id = aid;
        stage1_data_L(idx_L).day_minus2 = stage1_day_minus2;
        stage1_data_L(idx_L).day_minus1 = stage1_day_minus1;
        
        stage2_data_L(idx_L).animal_id = aid;
        stage2_data_L(idx_L).first_day_trials = stage2_first_trials;
    else
        idx_NL = idx_NL + 1;
        stage1_data_NL(idx_NL).animal_id = aid;
        stage1_data_NL(idx_NL).day_minus2 = stage1_day_minus2;
        stage1_data_NL(idx_NL).day_minus1 = stage1_day_minus1;
        
        stage2_data_NL(idx_NL).animal_id = aid;
        stage2_data_NL(idx_NL).first_day_trials = stage2_first_trials;
    end
end

fprintf('\n===== Data Collection Summary =====\n');
fprintf('Learners with complete data: %d\n', idx_L);
fprintf('Non-learners with complete data: %d\n', idx_NL);

% ===== Prepare plotting data =====
nL = length(stage1_data_L);
nNL = length(stage1_data_NL);

% Collect day -2, -1 averages
L_day_minus2 = arrayfun(@(x) x.day_minus2, stage1_data_L);
L_day_minus1 = arrayfun(@(x) x.day_minus1, stage1_data_L);
NL_day_minus2 = arrayfun(@(x) x.day_minus2, stage1_data_NL);
NL_day_minus1 = arrayfun(@(x) x.day_minus1, stage1_data_NL);

% Collect stage 2 trials into matrix
L_stage2_matrix = nan(Nshow, nL);
NL_stage2_matrix = nan(Nshow, nNL);

for i = 1:nL
    trials = stage2_data_L(i).first_day_trials;
    L_stage2_matrix(1:length(trials), i) = trials;
end

for i = 1:nNL
    trials = stage2_data_NL(i).first_day_trials;
    NL_stage2_matrix(1:length(trials), i) = trials;
end

% Compute trial means across animals
L_trial_means = mean(L_stage2_matrix, 2, 'omitnan');
NL_trial_means = mean(NL_stage2_matrix, 2, 'omitnan');

% Calculate global y-limits
all_data = [L_day_minus2(:); L_day_minus1(:); L_stage2_matrix(:); ...
            NL_day_minus2(:); NL_day_minus1(:); NL_stage2_matrix(:)];
y_min = min(all_data(~isnan(all_data))) - 0.5;
y_max = max(all_data(~isnan(all_data))) + 0.5;

jitter_amount = 0.3;

% ===== FIGURE 1: LEARNERS =====
figure('Position', [100, 100, 900, 700], 'Color', 'w');
hold on;

% IMPROVEMENT 1: Move stage 1 days further left for clarity
day_minus2_pos = -8;
day_minus1_pos = -2;

% Plot day -2 (individual points + mean)
x_minus2 = day_minus2_pos + (rand(nL, 1) - 0.5) * jitter_amount;
scatter(x_minus2, L_day_minus2, 180, cLearner, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cLearner*0.3, 'LineWidth', 1.5);

% Saturated color for mean
plot(day_minus2_pos, mean(L_day_minus2, 'omitnan'), 'o', 'MarkerSize', 18, ...
     'MarkerFaceColor', cLearner*0.8, 'MarkerEdgeColor', cLearner*0.5, 'LineWidth', 2.5);

% Plot day -1 (individual points + mean)
x_minus1 = day_minus1_pos + (rand(nL, 1) - 0.5) * jitter_amount;
scatter(x_minus1, L_day_minus1, 180, cLearner, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cLearner*0.3, 'LineWidth', 1.5);

plot(day_minus1_pos, mean(L_day_minus1, 'omitnan'), 'o', 'MarkerSize', 18, ...
     'MarkerFaceColor', cLearner*0.8, 'MarkerEdgeColor', cLearner*0.5, 'LineWidth', 2.5);

% Connect day -2 to day -1 for each animal
for i = 1:nL
    plot([day_minus2_pos, day_minus1_pos], [L_day_minus2(i), L_day_minus1(i)], '-', ...
         'Color', [cLearner, 0.25], 'LineWidth', 2);
end

% Plot stage 2 trials (jittered scatter) - lighter
for i = 1:nL
    trials = stage2_data_L(i).first_day_trials;
    n_trials = length(trials);
    x_trials = (1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25);
end

% HIGHLIGHTED: Plot group mean line for stage 2 with saturated color
trial_indices = 1:Nshow;
valid_trials = ~isnan(L_trial_means);
smooth_mean = movmean(L_trial_means(valid_trials), 5);  % 5-trial smoothing
plot(trial_indices(valid_trials), smooth_mean, '-', ...
     'Color', cLearner*0.8, 'LineWidth', 4.5);

% Transition line
xline(0, 'k--', 'LineWidth', 3, 'Alpha', 0.7);
text(0, y_max * 1.04, 'Transition', 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', ...
     'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5, ...
     'Margin', 3);

% Set x-ticks for both stage 1 days and stage 2 trials
trial_ticks = 1:5:Nshow;  % Show every 5th trial
xticks([day_minus2_pos, day_minus1_pos, trial_ticks]);
xticklabels([{'Day -2', 'Day -1'}, arrayfun(@(x) sprintf('T%d', x), trial_ticks, 'UniformOutput', false)]);

% Labels
xlabel('Day / Trial Number', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(sprintf('Lick Count [%.1f-%.1f]s', lick_window(1), lick_window(2)), ...
       'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Learners (n=%d)', nL), ...
      'FontSize', 18, 'FontWeight', 'bold', 'Color', cLearner*0.8);

xlim([day_minus2_pos - 1, Nshow + 0.5]);
ylim([y_min, y_max]);

% Clean axes without grid
set(gca, 'FontSize', 13, 'LineWidth', 2, 'Box', 'off');

% Stage labels
text((day_minus2_pos + day_minus1_pos)/2, y_max * 0.92, 'Stage 1', ...
     'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
text(Nshow/2, y_max * 0.92, 'Stage 2', 'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);

% ===== FIGURE 2: NON-LEARNERS =====
figure('Position', [920, 100, 900, 700], 'Color', 'w');
hold on;

% Move stage 1 days further left
day_minus2_pos = -8;
day_minus1_pos = -2;

% Plot day -2 (individual points + mean)
x_minus2 = day_minus2_pos + (rand(nNL, 1) - 0.5) * jitter_amount;
scatter(x_minus2, NL_day_minus2, 180, cNonLearner, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cNonLearner*0.6, 'LineWidth', 1.5);

plot(day_minus2_pos, mean(NL_day_minus2, 'omitnan'), 'o', 'MarkerSize', 18, ...
     'MarkerFaceColor', cNonLearner*0.8, 'MarkerEdgeColor', cNonLearner*0.5, 'LineWidth', 2.5);

% Plot day -1 (individual points + mean)
x_minus1 = day_minus1_pos + (rand(nNL, 1) - 0.5) * jitter_amount;
scatter(x_minus1, NL_day_minus1, 180, cNonLearner, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', cNonLearner*0.6, 'LineWidth', 1.5);

plot(day_minus1_pos, mean(NL_day_minus1, 'omitnan'), 'o', 'MarkerSize', 18, ...
     'MarkerFaceColor', cNonLearner*0.8, 'MarkerEdgeColor', cNonLearner*0.5, 'LineWidth', 2.5);

% Connect day -2 to day -1 for each animal
for i = 1:nNL
    plot([day_minus2_pos, day_minus1_pos], [NL_day_minus2(i), NL_day_minus1(i)], '-', ...
         'Color', [cNonLearner, 0.25], 'LineWidth', 2);
end

% Plot stage 2 trials (jittered scatter)
for i = 1:nNL
    trials = stage2_data_NL(i).first_day_trials;
    n_trials = length(trials);
    x_trials = (1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cNonLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25);
end

% HIGHLIGHTED: Plot group mean line for stage 2 with saturated color
trial_indices = 1:Nshow;
valid_trials = ~isnan(NL_trial_means);
smooth_mean = movmean(NL_trial_means(valid_trials), 5);
plot(trial_indices(valid_trials), smooth_mean, '-', ...
     'Color', cNonLearner*0.8, 'LineWidth', 4.5);

% Transition line
xline(0, 'k--', 'LineWidth', 3, 'Alpha', 0.7);
text(0, y_max * 1.04, 'Transition', 'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'top', 'FontSize', 14, 'FontWeight', 'bold', ...
     'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5, ...
     'Margin', 3);

% Set x-ticks for both stage 1 days and stage 2 trials
trial_ticks = 1:5:Nshow;  % Show every 5th trial
xticks([day_minus2_pos, day_minus1_pos, trial_ticks]);
xticklabels([{'Day -2', 'Day -1'}, arrayfun(@(x) sprintf('T%d', x), trial_ticks, 'UniformOutput', false)]);

% Labels
xlabel('Day / Trial Number', 'FontSize', 16, 'FontWeight', 'bold');
ylabel(sprintf('Lick Count [%.1f-%.1f]s', lick_window(1), lick_window(2)), ...
       'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('Non-learners (n=%d)', nNL), ...
      'FontSize', 18, 'FontWeight', 'bold', 'Color', cNonLearner*0.8);

xlim([day_minus2_pos - 1, Nshow + 0.5]);
ylim([y_min, y_max]);

% Clean axes without grid
set(gca, 'FontSize', 13, 'LineWidth', 2, 'Box', 'off');

% Stage labels
text((day_minus2_pos + day_minus1_pos)/2, y_max * 0.92, 'Stage 1', ...
     'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
text(Nshow/2, y_max * 0.92, 'Stage 2', 'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);

%% Plots a bar plot that compares static licks of last 2 days stage 1 vs first 2 days of stage 2

% Parameters
lick_window = [-3, -2];  % Window from stim onset (seconds)

% Pre-allocate structures for last 2 days of stage 1 + first 2 days of stage 2
stage1_data_L = struct('animal_id', {}, 'day_minus2', {}, 'day_minus1', {});
stage2_data_L = struct('animal_id', {}, 'day1', {}, 'day2', {});
stage1_data_NL = struct('animal_id', {}, 'day_minus2', {}, 'day_minus1', {});
stage2_data_NL = struct('animal_id', {}, 'day1', {}, 'day2', {});

idx_L = 0;
idx_NL = 0;

% Helper function to calculate lick count in window
calc_lick_count = @(lick_times, ref_time, window) ...
    sum(lick_times >= (ref_time + window(1)) & lick_times < (ref_time + window(2)));

% Single loop collecting BOTH stages per animal
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % manually remove HA013
    if aid == "HA013"
        continue;
    end
    % ===== Find Stage 1 (last 2 classical days) =====
    is_stage_1 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    
    stage1_day_minus2 = [];
    stage1_day_minus1 = [];
    
    if sum(is_stage_1) >= 2
        cand_classical = find(is_stage_1);
        lastIdx = cand_classical(end);
        secondLastIdx = cand_classical(end-2);
        
        % Process day -2
        D_minus2 = rd(secondLastIdx);
        if isfield(D_minus2, 'cs_labels') && isfield(D_minus2, 'lick_event_times') && ...
           isfield(D_minus2, 'right_stim_on_times')
            stim_on_times = D_minus2.right_stim_on_times;
            lick_times = D_minus2.lick_event_times;
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), stim_on_times);
            stage1_day_minus2 = mean(lick_counts, 'omitnan');
        end
        
        % Process day -1
        D_minus1 = rd(lastIdx);
        if isfield(D_minus1, 'cs_labels') && isfield(D_minus1, 'lick_event_times') && ...
           isfield(D_minus1, 'right_stim_on_times')
            stim_on_times = D_minus1.right_stim_on_times;
            lick_times = D_minus1.lick_event_times;
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), stim_on_times);
            stage1_day_minus1 = mean(lick_counts, 'omitnan');
        end
    end

    % ===== Find Stage 2 (first 2 static days) =====
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ...
                          ~contains(d.workflow,'big_stim'), rd);
    
    stage2_day1 = [];
    stage2_day2 = [];
    
    if sum(is_stage_2) >= 2
        cand_static = find(is_stage_2);
        firstStaticIdx = cand_static(1);
        secondStaticIdx = cand_static(2);
        
        % Process first static day
        D2_1 = rd(firstStaticIdx);
        if isfield(D2_1, 'cs_labels') && isfield(D2_1, 'lick_event_times') && ...
           isfield(D2_1, 'right_stim_on_times')
            stim_on_times = D2_1.right_stim_on_times;
            lick_times = D2_1.lick_event_times;
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), stim_on_times);
            stage2_day1 = mean(lick_counts, 'omitnan');
        end
        
        % Process second static day
        D2_2 = rd(secondStaticIdx);
        if isfield(D2_2, 'cs_labels') && isfield(D2_2, 'lick_event_times') && ...
           isfield(D2_2, 'right_stim_on_times')
            stim_on_times = D2_2.right_stim_on_times;
            lick_times = D2_2.lick_event_times;
            lick_counts = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), stim_on_times);
            stage2_day2 = mean(lick_counts, 'omitnan');
        end
    end

    % ===== Only include animals with ALL data =====
    if isempty(stage1_day_minus2) || isempty(stage1_day_minus1) || ...
       isempty(stage2_day1) || isempty(stage2_day2)
        fprintf('Skipping animal %s: missing data\n', aid);
        continue;
    end

    % ===== Assign to learner or non-learner =====
    if ismember(aid, learners_group_ID)
        idx_L = idx_L + 1;
        stage1_data_L(idx_L).animal_id = aid;
        stage1_data_L(idx_L).day_minus2 = stage1_day_minus2;
        stage1_data_L(idx_L).day_minus1 = stage1_day_minus1;
        stage2_data_L(idx_L).animal_id = aid;
        stage2_data_L(idx_L).day1 = stage2_day1;
        stage2_data_L(idx_L).day2 = stage2_day2;
    else
        idx_NL = idx_NL + 1;
        stage1_data_NL(idx_NL).animal_id = aid;
        stage1_data_NL(idx_NL).day_minus2 = stage1_day_minus2;
        stage1_data_NL(idx_NL).day_minus1 = stage1_day_minus1;
        stage2_data_NL(idx_NL).animal_id = aid;
        stage2_data_NL(idx_NL).day1 = stage2_day1;
        stage2_data_NL(idx_NL).day2 = stage2_day2;
    end
end

fprintf('\n===== Data Collection Summary =====\n');
fprintf('Learners with complete data: %d\n', idx_L);
fprintf('Non-learners with complete data: %d\n', idx_NL);

% ===== Compute averages per animal =====
nL = length(stage1_data_L);
nNL = length(stage1_data_NL);

% Learners: average last 2 days of stage 1 and first 2 days of stage 2
L_stage1_avg = nan(nL, 1);
L_stage2_avg = nan(nL, 1);
for i = 1:nL
    L_stage1_avg(i) = mean([stage1_data_L(i).day_minus2, stage1_data_L(i).day_minus1]);
    L_stage2_avg(i) = mean([stage2_data_L(i).day1, stage2_data_L(i).day2]);
end

% Non-learners: average last 2 days of stage 1 and first 2 days of stage 2
NL_stage1_avg = nan(nNL, 1);
NL_stage2_avg = nan(nNL, 1);
for i = 1:nNL
    NL_stage1_avg(i) = mean([stage1_data_NL(i).day_minus2, stage1_data_NL(i).day_minus1]);
    NL_stage2_avg(i) = mean([stage2_data_NL(i).day1, stage2_data_NL(i).day2]);
end

% ===== Statistical tests (Wilcoxon signed-rank for paired data) =====
p_L = signrank(L_stage1_avg, L_stage2_avg);
p_NL = signrank(NL_stage1_avg, NL_stage2_avg);

fprintf('\n===== Statistical Tests =====\n');
fprintf('Learners: Stage 1 vs Stage 2, p = %.4f\n', p_L);
fprintf('Non-learners: Stage 1 vs Stage 2, p = %.4f\n', p_NL);

% ===== Plotting =====
figure('Position', [100, 100, 1200, 500], 'Color', 'w');

% Learners panel
subplot(1, 2, 1);
bar_data_L = [median(L_stage1_avg, 'omitnan'), median(L_stage2_avg, 'omitnan')];
b = bar(bar_data_L, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
b.CData = [cLearner*0.6; cLearner*0.9];
hold on;

% Plot individual animal connections
for i = 1:nL
    plot([1 2], [L_stage1_avg(i), L_stage2_avg(i)], ...
        'o-', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'none', 'LineWidth', 1);
end

ylabel(sprintf('Lick Count [%.1f-%.1f]s', lick_window(1), lick_window(2)), 'FontSize', 12);
xticklabels({'Stage 1 (last 2 days)', 'Stage 2 (first 2 days)'});
title(sprintf('Learners (n=%d)\np = %.4f', nL, p_L), 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
ylim([0, max([L_stage1_avg; L_stage2_avg; NL_stage1_avg; NL_stage2_avg])*1.1]);

% Non-learners panel
subplot(1, 2, 2);
bar_data_NL = [median(NL_stage1_avg, 'omitnan'), median(NL_stage2_avg, 'omitnan')];
b = bar(bar_data_NL, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
b.CData = [cNonLearner*0.6; cNonLearner*0.9];
hold on;

% Plot individual animal connections
for i = 1:nNL
    plot([1 2], [NL_stage1_avg(i), NL_stage2_avg(i)], ...
        'o-', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'none', 'LineWidth', 1);
end

ylabel(sprintf('Lick Count [%.1f-%.1f]s', lick_window(1), lick_window(2)), 'FontSize', 12);
xticklabels({'Stage 1 (last 2 days)', 'Stage 2 (first 2 days)'});
title(sprintf('Non-learners (n=%d)\np = %.4f', nNL, p_NL), 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'FontSize', 12);
ylim([0, max([L_stage1_avg; L_stage2_avg; NL_stage1_avg; NL_stage2_avg])*1.1]);

%%

% Parameters
Nshow = 50;  % Number of trials to show
lick_window = [0, 1];  % Window from stim onset (seconds)
iti_window = [-3, -2];  % ITI baseline window (1s before stim)

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learners_group_ID));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate structures
stage1_data_L = struct('animal_id', {}, 'last_day_trials', {}, 'last_day_iti', {});
stage2_data_L = struct('animal_id', {}, 'first_day_trials', {}, 'first_day_iti', {});
stage1_data_NL = struct('animal_id', {}, 'last_day_trials', {}, 'last_day_iti', {});
stage2_data_NL = struct('animal_id', {}, 'first_day_trials', {}, 'first_day_iti', {});

idx_L = 0;
idx_NL = 0;

% Helper function to calculate lick count in window
calc_lick_count = @(lick_times, ref_time, window) ...
    sum(lick_times >= (ref_time + window(1)) & lick_times < (ref_time + window(2)));

% Single loop collecting BOTH stages per animal
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end
    
    % Find classical conditioning days (Stage 1)
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);

    if ~any(isClassical)
        fprintf('Animal %s: No Stage 1 days, skipping\n', aid);
        continue;
    end

    % Check if learned Stage 1
    classical_days = find(isClassical);
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: Never learned Stage 1,skipping\n', aid);
        continue;
    end

    % ===== Find Stage 1 (last day, last N trials) =====
    is_stage_1 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    
    stage1_last_trials = [];
    stage1_last_iti = [];
    
    if any(is_stage_1)
        cand_classical = find(is_stage_1);
        lastIdx = cand_classical(end);
        
        % Process last day
        D_last = rd(lastIdx);
        if isfield(D_last, 'cs_labels') && isfield(D_last, 'lick_event_times') && ...
           isfield(D_last, 'right_stim_on_times')
            
            stim_on_times = D_last.right_stim_on_times;
            lick_times = D_last.lick_event_times;
            
            % Get last Nshow trials
            n_trials = min(Nshow, length(stim_on_times));
            if n_trials >= 1
                % last_n_indices = 1:n_trials;
                
                stage1_last_trials = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), ...
                                              stim_on_times(1:n_trials));
                stage1_last_iti = arrayfun(@(t) calc_lick_count(lick_times, t, iti_window), ...
                                           stim_on_times(1:n_trials));
            end
        end
    end

    % ===== Find Stage 2 (first day, first N trials) =====
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ...
                          ~contains(d.workflow,'big_stim'), rd);
    
    stage2_first_trials = [];
    stage2_first_iti = [];
    
    if any(is_stage_2)
        cand_static = find(is_stage_2);
        firstStaticIdx = cand_static(1);
        
        % First day
        D2_first = rd(firstStaticIdx);
        if isfield(D2_first, 'cs_labels') && isfield(D2_first, 'lick_event_times') && ...
           isfield(D2_first, 'right_stim_on_times')
            
            stim_on_times = D2_first.right_stim_on_times;
            lick_times = D2_first.lick_event_times;
            
            n_trials = min(Nshow, length(stim_on_times));
            if n_trials >= 1
                stage2_first_trials = arrayfun(@(t) calc_lick_count(lick_times, t, lick_window), ...
                                               stim_on_times(1:n_trials));
                stage2_first_iti = arrayfun(@(t) calc_lick_count(lick_times, t, iti_window), ...
                                            stim_on_times(1:n_trials));
            end
        end
    end

    % ===== Only include animals with ALL data =====
    if isempty(stage1_last_trials) || isempty(stage2_first_trials)
        fprintf('Skipping animal %s: missing data\n', aid);
        continue;
    end

    % ===== Assign to learner or non-learner =====
    if ismember(aid, learners_group_ID)
        idx_L = idx_L + 1;
        stage1_data_L(idx_L).animal_id = aid;
        stage1_data_L(idx_L).last_day_trials = stage1_last_trials;
        stage1_data_L(idx_L).last_day_iti = stage1_last_iti;
        
        stage2_data_L(idx_L).animal_id = aid;
        stage2_data_L(idx_L).first_day_trials = stage2_first_trials;
        stage2_data_L(idx_L).first_day_iti = stage2_first_iti;
    else
        idx_NL = idx_NL + 1;
        stage1_data_NL(idx_NL).animal_id = aid;
        stage1_data_NL(idx_NL).last_day_trials = stage1_last_trials;
        stage1_data_NL(idx_NL).last_day_iti = stage1_last_iti;
        
        stage2_data_NL(idx_NL).animal_id = aid;
        stage2_data_NL(idx_NL).first_day_trials = stage2_first_trials;
        stage2_data_NL(idx_NL).first_day_iti = stage2_first_iti;
    end
end

fprintf('\n===== Data Collection Summary =====\n');
fprintf('Learners with complete data: %d\n', idx_L);
fprintf('Non-learners with complete data: %d\n', idx_NL);

% ===== Prepare plotting data =====
nL = length(stage1_data_L);
nNL = length(stage1_data_NL);

% Collect stage 1 last trials into matrix
L_stage1_last = nan(Nshow, nL);
L_stage1_last_iti = nan(Nshow, nL);
NL_stage1_last = nan(Nshow, nNL);
NL_stage1_last_iti = nan(Nshow, nNL);

for i = 1:nL
    trials = stage1_data_L(i).last_day_trials;
    L_stage1_last(1:length(trials), i) = trials;
    iti_trials = stage1_data_L(i).last_day_iti;
    L_stage1_last_iti(1:length(iti_trials), i) = iti_trials;
end

for i = 1:nNL
    trials = stage1_data_NL(i).last_day_trials;
    NL_stage1_last(1:length(trials), i) = trials;
    iti_trials = stage1_data_NL(i).last_day_iti;
    NL_stage1_last_iti(1:length(iti_trials), i) = iti_trials;
end

% Collect stage 2 first trials into matrix
L_stage2_first = nan(Nshow, nL);
L_stage2_first_iti = nan(Nshow, nL);
NL_stage2_first = nan(Nshow, nNL);
NL_stage2_first_iti = nan(Nshow, nNL);

for i = 1:nL
    trials = stage2_data_L(i).first_day_trials;
    L_stage2_first(1:length(trials), i) = trials;
    iti_trials = stage2_data_L(i).first_day_iti;
    L_stage2_first_iti(1:length(iti_trials), i) = iti_trials;
end

for i = 1:nNL
    trials = stage2_data_NL(i).first_day_trials;
    NL_stage2_first(1:length(trials), i) = trials;
    iti_trials = stage2_data_NL(i).first_day_iti;
    NL_stage2_first_iti(1:length(iti_trials), i) = iti_trials;
end

% Compute trial means
L_stage1_means = mean(L_stage1_last, 2, 'omitnan');
L_stage2_means = mean(L_stage2_first, 2, 'omitnan');
L_stage1_iti_means = mean(L_stage1_last_iti, 2, 'omitnan');
L_stage2_iti_means = mean(L_stage2_first_iti, 2, 'omitnan');

NL_stage1_means = mean(NL_stage1_last, 2, 'omitnan');
NL_stage2_means = mean(NL_stage2_first, 2, 'omitnan');
NL_stage1_iti_means = mean(NL_stage1_last_iti, 2, 'omitnan');
NL_stage2_iti_means = mean(NL_stage2_first_iti, 2, 'omitnan');

% Calculate global y-limits
all_data = [L_stage1_last(:); L_stage2_first(:); L_stage1_last_iti(:); L_stage2_first_iti(:); ...
            NL_stage1_last(:); NL_stage2_first(:); NL_stage1_last_iti(:); NL_stage2_first_iti(:)];
y_min = min(all_data(~isnan(all_data))) - 0.5;
y_max = max(all_data(~isnan(all_data))) + 0.5;

jitter_amount = 0.3;

% ===== FIGURE 1: LEARNERS =====
figure('Position', [100, 100, 1200, 700], 'Color', 'w');
hold on;

% X-axis: Stage 1 trials are negative (-50 to -1), Stage 2 trials are positive (1 to 50)
stage1_x = (-Nshow:-1);
stage2_x = (1:Nshow);

% Plot Stage 1 last 50 trials (scatter - NO legend entry)
for i = 1:nL
    trials = stage1_data_L(i).last_day_trials;
    n_trials = length(trials);
    x_trials = stage1_x(1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25, ...
            'HandleVisibility', 'off');
end

valid_s1 = ~isnan(L_stage1_means);
smooth_s1 = movmean(L_stage1_means(valid_s1), 5);
plot(stage1_x(valid_s1), smooth_s1, '-', ...
     'Color', cLearner*0.8, 'LineWidth', 4.5, 'DisplayName', sprintf('Stim [%.1f-%.1f]s', lick_window));

% Plot Stage 2 first 50 trials (scatter - NO legend entry)
for i = 1:nL
    trials = stage2_data_L(i).first_day_trials;
    n_trials = length(trials);
    x_trials = stage2_x(1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25, ...
            'HandleVisibility', 'off');
end

valid_s2 = ~isnan(L_stage2_means);
smooth_s2 = movmean(L_stage2_means(valid_s2), 5);
plot(stage2_x(valid_s2), smooth_s2, '-', ...
     'Color', cLearner*0.7, 'LineWidth', 4.5, 'HandleVisibility', 'off');

% === ITI BASELINE OVERLAY (dashed lines) ===
% Stage 1 ITI smoothed mean
valid_s1_iti = ~isnan(L_stage1_iti_means);
smooth_s1_iti = movmean(L_stage1_iti_means(valid_s1_iti), 5);
plot(stage1_x(valid_s1_iti), smooth_s1_iti, '-', ...
     'Color', [0.7 0.7 0.7], 'LineWidth', 3.5, 'DisplayName', sprintf('ITI [%.1f-%.1f]s', iti_window));

% Stage 2 ITI smoothed mean
valid_s2_iti = ~isnan(L_stage2_iti_means);
smooth_s2_iti = movmean(L_stage2_iti_means(valid_s2_iti), 5);
plot(stage2_x(valid_s2_iti), smooth_s2_iti, '-', ...
     'Color', [0.7 0.7 0.7], 'LineWidth', 3.5, 'HandleVisibility', 'off');

% Transition line at x=0
xline(0, 'k-', 'LineWidth', 3, 'Alpha', 0.9,'HandleVisibility', 'off');

% X-axis labels - showing trial numbers within each stage
trial_ticks_s1 = [-40, -30, -20, -10];  % Corresponds to trials 10, 20, 30, 40 in Stage 1
trial_ticks_s2 = [10, 20, 30, 40];      % Trials 10, 20, 30, 40 in Stage 2
xticks([trial_ticks_s1, 0, trial_ticks_s2]);
xticklabels({'10', '20', '30', '40', '0', '10', '20', '30', '40', '50'});
xlabel('Trial Number', 'FontSize', 16, 'FontWeight', 'bold');

% xlabel('Trial Number (relative to stage transition)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Lick Count', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('mPFC+ (n=%d)', nL), 'FontSize', 18, 'FontWeight', 'bold');

xlim([-Nshow-2, Nshow+2]);
ylim([y_min, y_max]);
set(gca, 'FontSize', 13, 'LineWidth', 2, 'Box', 'off');

% Stage labels
text(-Nshow/2, y_max * 0.92, 'Stage 1', ...
     'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
text(Nshow/2, y_max * 0.92, 'Stage 2', 'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);

legend('Location', 'best', 'FontSize', 11);

% ===== FIGURE 2: NON-LEARNERS =====
figure('Position', [920, 100, 1200, 700], 'Color', 'w');
hold on;

% Plot Stage 1 last 50 trials (scatter - NO legend entry)
for i = 1:nNL
    trials = stage1_data_NL(i).last_day_trials;
    n_trials = length(trials);
    x_trials = stage1_x(1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cNonLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25, ...
            'HandleVisibility', 'off');
end

valid_s1 = ~isnan(NL_stage1_means);
smooth_s1 = movmean(NL_stage1_means(valid_s1), 5);
plot(stage1_x(valid_s1), smooth_s1, '-', ...
     'Color', cNonLearner*0.8, 'LineWidth', 4.5, 'DisplayName', sprintf('Stim [%.1f-%.1f]s', lick_window));

% Plot Stage 2 first 50 trials (scatter - NO legend entry)
for i = 1:nNL
    trials = stage2_data_NL(i).first_day_trials;
    n_trials = length(trials);
    x_trials = stage2_x(1:n_trials) + (rand(1, n_trials) - 0.5) * jitter_amount * 0.8;
    scatter(x_trials, trials, 40, cNonLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.25, ...
            'HandleVisibility', 'off');
end

valid_s2 = ~isnan(NL_stage2_means);
smooth_s2 = movmean(NL_stage2_means(valid_s2), 5);
plot(stage2_x(valid_s2), smooth_s2, '-', ...
     'Color', cNonLearner*0.8, 'LineWidth', 4.5, 'HandleVisibility', 'off');

% === ITI BASELINE OVERLAY ===
valid_s1_iti = ~isnan(NL_stage1_iti_means);
smooth_s1_iti = movmean(NL_stage1_iti_means(valid_s1_iti), 5);
plot(stage1_x(valid_s1_iti), smooth_s1_iti, '-', ...
     'Color', [0.7 0.7 0.7], 'LineWidth', 3.5, 'DisplayName', sprintf('ITI [%.1f-%.1f]s', iti_window));

valid_s2_iti = ~isnan(NL_stage2_iti_means);
smooth_s2_iti = movmean(NL_stage2_iti_means(valid_s2_iti), 5);
plot(stage2_x(valid_s2_iti), smooth_s2_iti, '-', ...
     'Color', [0.7 0.7 0.7], 'LineWidth', 3.5, 'HandleVisibility', 'off');

% Transition line
xline(0, 'k-', 'LineWidth', 3, 'Alpha', 0.9,'HandleVisibility', 'off');
% text(0, y_max * 1.04, 'Transition', 'HorizontalAlignment', 'center', ...
%      'VerticalAlignment', 'top', 'FontSize', 12, 'FontWeight', 'bold', ...
%      'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'LineWidth', 1.5, ...
%      'Margin', 3);

% X-axis labels
trial_ticks_s1 = [-40, -30, -20, -10];  % Corresponds to trials 10, 20, 30, 40 in Stage 1
trial_ticks_s2 = [10, 20, 30, 40];      % Trials 10, 20, 30, 40 in Stage 2
xticks([trial_ticks_s1, 0, trial_ticks_s2]);
xticklabels({'10', '20', '30', '40', '0', '10', '20', '30', '40', '50'});
xlabel('Trial Number', 'FontSize', 16, 'FontWeight', 'bold');

ylabel('Lick Count', 'FontSize', 16, 'FontWeight', 'bold');
title(sprintf('mPFC- (n=%d)', nNL), 'FontSize', 18, 'FontWeight', 'bold');

xlim([-Nshow-2, Nshow+2]);
ylim([y_min, y_max]);
set(gca, 'FontSize', 13, 'LineWidth', 2, 'Box', 'off');

% Stage labels
text(-Nshow/2, y_max * 0.92, 'Stage 1', ...
     'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);
text(Nshow/2, y_max * 0.92, 'Stage 2', 'FontSize', 13, 'FontWeight', 'bold', ...
     'HorizontalAlignment', 'center', 'Color', [0.3 0.3 0.3]);

legend('Location', 'best', 'FontSize', 11);


%% Plots a bar plot that compares the startegy shift between las
% t day stage 1 and first day stage 2 split between learners vs-non-learners

% ===== STEP 1: Calculate normalized ratios for each animal =====

% Normalized ratio: (A - S) / (A + S)
% This ranges from -1 (only static) to +1 (only anticipatory)

% Stage 1 - average across trials for each animal
stage1_anticip_NL = mean(A_mat_NL_s1, 1, 'omitnan')'; % per animal
stage1_static_NL = mean(S_mat_NL_s1, 1, 'omitnan')';
stage1_ratio_NL = (stage1_anticip_NL - stage1_static_NL) ./ (stage1_anticip_NL + stage1_static_NL + eps);

stage1_anticip_L = mean(A_mat_L_s1, 1, 'omitnan')';
stage1_static_L = mean(S_mat_L_s1, 1, 'omitnan')';
stage1_ratio_L = (stage1_anticip_L - stage1_static_L) ./ (stage1_anticip_L + stage1_static_L + eps);

% Stage 2 - average across trials for each animal
stage2_anticip_NL = mean(A_mat_NL, 1, 'omitnan')';
stage2_static_NL = mean(S_mat_NL, 1, 'omitnan')';
stage2_ratio_NL = (stage2_anticip_NL - stage2_static_NL) ./ (stage2_anticip_NL + stage2_static_NL + eps);

stage2_anticip_L = mean(A_mat_L, 1, 'omitnan')';
stage2_static_L = mean(S_mat_L, 1, 'omitnan')';
stage2_ratio_L = (stage2_anticip_L - stage2_static_L) ./ (stage2_anticip_L + stage2_static_L + eps);

% ===== STEP 2: Calculate strategy shift (delta) =====

% Change in normalized ratio from Stage 1 to Stage 2
delta_ratio_NL = stage2_ratio_NL - stage1_ratio_NL;
delta_ratio_L = stage2_ratio_L - stage1_ratio_L;

% ===== STEP 3: Statistical comparison =====

fprintf('\n===== Strategy Shift Analysis (Normalized Ratio) =====\n');
fprintf('Normalized ratio = (Anticip - Static) / (Anticip + Static)\n');
fprintf('  +1 = only anticipatory, 0 = equal, -1 = only static\n\n');

% Test 1: Compare delta between groups
[h1, p1] = ttest2(delta_ratio_NL, delta_ratio_L);
fprintf('T-test: Strategy shift differs between groups\n');
fprintf('  p-value = %.4f\n', p1);
fprintf('  Non-learners Δratio: %.3f ± %.3f\n', mean(delta_ratio_NL), std(delta_ratio_NL));
fprintf('  Learners Δratio: %.3f ± %.3f\n', mean(delta_ratio_L), std(delta_ratio_L));

% Test 2: Within-group changes
[h2, p2] = ttest(stage1_ratio_NL, stage2_ratio_NL);
fprintf('\nNon-learners: Stage 1 vs Stage 2 ratio\n');
fprintf('  p-value = %.4f\n', p2);
fprintf('  Stage 1: %.3f ± %.3f\n', mean(stage1_ratio_NL), std(stage1_ratio_NL));
fprintf('  Stage 2: %.3f ± %.3f\n', mean(stage2_ratio_NL), std(stage2_ratio_NL));

[h3, p3] = ttest(stage1_ratio_L, stage2_ratio_L);
fprintf('\nLearners: Stage 1 vs Stage 2 ratio\n');
fprintf('  p-value = %.4f\n', p3);
fprintf('  Stage 1: %.3f ± %.3f\n', mean(stage1_ratio_L), std(stage1_ratio_L));
fprintf('  Stage 2: %.3f ± %.3f\n', mean(stage2_ratio_L), std(stage2_ratio_L));

% Effect size (Cohen's d)
cohens_d = (mean(delta_ratio_NL) - mean(delta_ratio_L)) / ...
    sqrt((std(delta_ratio_NL)^2 + std(delta_ratio_L)^2) / 2);
fprintf('\nEffect size (Cohen''s d): %.3f\n', cohens_d);

fprintf('\n===== Direction of shift =====\n');
fprintf('Non-learners with decreased ratio (more static): %d/%d (%.0f%%)\n', ...
    sum(delta_ratio_NL < 0), length(delta_ratio_NL), ...
    100*sum(delta_ratio_NL < 0)/length(delta_ratio_NL));
fprintf('Learners with decreased ratio (more static): %d/%d (%.0f%%)\n', ...
    sum(delta_ratio_L < 0), length(delta_ratio_L), ...
    100*sum(delta_ratio_L < 0)/length(delta_ratio_L));

% ===== STEP 4: Visualization =====

figure('Color', 'w', 'Position', [100 100 600 500]);
hold on;

% Prepare data for grouped bar plot
group_means = [mean(stage1_ratio_L), mean(stage1_ratio_NL); ...
    mean(stage2_ratio_L), mean(stage2_ratio_NL)];
group_sem = [std(stage1_ratio_L)/sqrt(length(stage1_ratio_L)), std(stage1_ratio_NL)/sqrt(length(stage1_ratio_NL)); ...
    std(stage2_ratio_L)/sqrt(length(stage2_ratio_L)), std(stage2_ratio_NL)/sqrt(length(stage2_ratio_NL))];

% Create bar plot
b = bar(group_means);
b(1).FaceColor = [0.8500 0.3250 0.0980]; % Learners
b(2).FaceColor = [0 0.4470 0.7410];      % Non-learners
b(1).FaceAlpha = [0.5]; % Stage 1 lighter, Stage 2 darker
b(2).FaceAlpha = [0.5];

% Add error bars
x_positions = [1 2]; % Stage 1, Stage 2
offset = 0.15; % Offset for grouped bars
for i = 1:2 % Stages
    errorbar(x_positions(i) - offset, group_means(i,1), group_sem(i,1), ...
        'k', 'LineStyle', 'none', 'LineWidth', 2);
    errorbar(x_positions(i) + offset, group_means(i,2), group_sem(i,2), ...
        'k', 'LineStyle', 'none', 'LineWidth', 2);
end

% Add individual data points
jitter = 0.05;
% Stage 1 Learners
scatter(ones(size(stage1_ratio_L)) * (1 - offset) + jitter*randn(size(stage1_ratio_L)), ...
    stage1_ratio_L, 50, [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha', 0.4);
% Stage 1 Non-learners
scatter(ones(size(stage1_ratio_NL)) * (1 + offset) + jitter*randn(size(stage1_ratio_NL)), ...
    stage1_ratio_NL, 50, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.4);
% Stage 2 Learners
scatter(ones(size(stage2_ratio_L)) * (2 - offset) + jitter*randn(size(stage2_ratio_L)), ...
    stage2_ratio_L, 50, [0.8500 0.3250 0.0980], 'filled', 'MarkerFaceAlpha', 0.4);
% Stage 2 Non-learners
scatter(ones(size(stage2_ratio_NL)) * (2 + offset) + jitter*randn(size(stage2_ratio_NL)), ...
    stage2_ratio_NL, 50, [0 0.4470 0.7410], 'filled', 'MarkerFaceAlpha', 0.4);

yline(0, 'k--', 'LineWidth', 1.5);
xticks([1 2]);
xticklabels({'Stage 1 (Last Day)', 'Stage 2 (First Day)'});
ylabel('Normalized Ratio: (A-S)/(A+S)');
title('Licking Strategy: Stage 1 → Stage 2');
legend(b, {'Learners', 'Non-learners'}, 'Location', 'best');
grid on;
box off;


%% ===== CS+ vs CS- Discrimination Analysis =====

% For LAST day of Stage 1 (when both groups have "learned")

% Initialize storage
discrimination_metrics = struct();
idx = 0;  % Counter for valid animals only

for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % Find last classical conditioning day
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    if ~any(isClassical), continue; end

    % Create another conditional to skip animals that haven't learned the first stage
    classical_days = find(isClassical);
    sigDays = combined_sig_day_all_protocols{a}(classical_days); % get the sig days for relevant protocol
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: no learning day found, skipping\n', behaviour_data(a).animal_id);
        continue
    end

    cand = find(isClassical);
    lastIdx = cand(end);
    D = rd(lastIdx);

    % Required fields
    if ~isfield(D,'anticipatory_licks') || ~isfield(D,'cs_labels'), continue; end

    % Create a labels mask
    cs_plus_mask = D.cs_labels;

    % === Metric 1: Anticipatory lick rate ===
    anticip_cs_plus = mean(D.anticipatory_licks(cs_plus_mask), 'omitnan');
    anticip_cs_minus = mean(D.anticipatory_licks(~cs_plus_mask), 'omitnan');

    % === Metric 2: Static period lick rate ===
    static_cs_plus = mean(D.static_time_lick_counts(cs_plus_mask) ./ ...
        D.all_static_times(cs_plus_mask), 'omitnan');
    static_cs_minus = mean(D.static_time_lick_counts(~cs_plus_mask) ./ ...
        D.all_static_times(~cs_plus_mask), 'omitnan');

    % === Metric 3: Total licks per trial ===
    total_licks_cs_plus = mean(D.anticipatory_licks(cs_plus_mask) + ...
        D.static_time_lick_counts(cs_plus_mask), 'omitnan');
    total_licks_cs_minus = mean(D.anticipatory_licks(~cs_plus_mask) + ...
        D.static_time_lick_counts(~cs_plus_mask), 'omitnan');

    % === Metric 4: Discrimination index ===
    % (CS+ - CS-) / (CS+ + CS-) for anticipatory licking
    discrimination_idx = (anticip_cs_plus - anticip_cs_minus) / ...
        (anticip_cs_plus + anticip_cs_minus + eps);

    % Store metrics - INCREMENT idx only for valid animals
    idx = idx + 1;  % This is the key fix!

    discrimination_metrics(idx).animal_id = aid;
    discrimination_metrics(idx).is_learner = ismember(aid, learner_ids);
    discrimination_metrics(idx).anticip_cs_plus = anticip_cs_plus;
    discrimination_metrics(idx).anticip_cs_minus = anticip_cs_minus;
    discrimination_metrics(idx).static_cs_plus = static_cs_plus;
    discrimination_metrics(idx).static_cs_minus = static_cs_minus;
    discrimination_metrics(idx).total_cs_plus = total_licks_cs_plus;
    discrimination_metrics(idx).total_cs_minus = total_licks_cs_minus;
    discrimination_metrics(idx).discrimination_idx = discrimination_idx;
end

% Verification
fprintf('\n===== Data Collection Summary =====\n');
fprintf('Total animals processed: %d\n', idx);
fprintf('Learners: %d\n', sum([discrimination_metrics.is_learner]));
fprintf('Non-learners: %d\n', sum(~[discrimination_metrics.is_learner]));

% ===== Statistical Comparison =====

learners_mask = [discrimination_metrics.is_learner];


% Test 1: Anticipatory licking to CS+
anticip_plus_L = [discrimination_metrics(learners_mask).anticip_cs_plus];
anticip_plus_NL = [discrimination_metrics(~learners_mask).anticip_cs_plus];
[~, p1] = ttest2(anticip_plus_L, anticip_plus_NL);

% Test 2: Anticipatory licking to CS-
anticip_minus_L = [discrimination_metrics(learners_mask).anticip_cs_minus];
anticip_minus_NL = [discrimination_metrics(~learners_mask).anticip_cs_minus];
[~, p2] = ttest2(anticip_minus_L, anticip_minus_NL);

% Test 3: Discrimination index
discrim_idx_L = [discrimination_metrics(learners_mask).discrimination_idx];
discrim_idx_NL = [discrimination_metrics(~learners_mask).discrimination_idx];
[~, p3] = ttest2(discrim_idx_L, discrim_idx_NL);

% Test 4: Within-group discrimination (CS+ > CS-)
[~, p4_L] = ttest(anticip_plus_L, anticip_minus_L);
[~, p4_NL] = ttest(anticip_plus_NL, anticip_minus_NL);


% ===== Visualization =====

figure('Color', 'w', 'Position', [100 100 1200 400]);

% Panel 1: CS+ vs CS- anticipatory licking
subplot(1, 3, 1);
hold on;

% Learners
bar(1, mean(anticip_plus_L), 'FaceColor', cCS_plus, 'FaceAlpha', 0.5,'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);
bar(2, mean(anticip_minus_L), 'FaceColor', cCS_minus, 'FaceAlpha', 0.5,'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);

% Non-learners
bar(4, mean(anticip_plus_NL), 'FaceColor', cCS_plus, 'FaceAlpha', 0.5,'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);
bar(5, mean(anticip_minus_NL), 'FaceColor', cCS_minus, 'FaceAlpha', 0.7,'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);

% Error bars
errorbar([1 2 4 5], ...
    [mean(anticip_plus_L), mean(anticip_minus_L), mean(anticip_plus_NL), mean(anticip_minus_NL)], ...
    [std(anticip_plus_L)/sqrt(length(anticip_plus_L)), ...
    std(anticip_minus_L)/sqrt(length(anticip_minus_L)), ...
    std(anticip_plus_NL)/sqrt(length(anticip_plus_NL)), ...
    std(anticip_minus_NL)/sqrt(length(anticip_minus_NL))], ...
    'k', 'LineStyle', 'none', 'LineWidth', 2);

% Individual points
scatter(ones(size(anticip_plus_L)), anticip_plus_L, 50, cCS_plus, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(2*ones(size(anticip_minus_L)), anticip_minus_L, 50, cCS_minus, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(4*ones(size(anticip_plus_NL)), anticip_plus_NL, 50, cCS_plus, 'filled', 'MarkerFaceAlpha', 0.5);
scatter(5*ones(size(anticip_minus_NL)), anticip_minus_NL, 50, cCS_minus, 'filled', 'MarkerFaceAlpha', 0.5);

xticks([1.5 4.5]);
xticklabels({'mPFC+', 'mPFC-'});
ylabel('Anticipatory lick rate (licks/s)');
title('CS+ vs CS- Discrimination');
legend({'CS+', 'CS-'}, 'Location', 'best');

% Panel 2: Discrimination index
subplot(1, 3, 2);
hold on;

data_for_bar = [mean(discrim_idx_L), mean(discrim_idx_NL)];
b = bar([1 2], data_for_bar, 'FaceColor', 'flat','FaceAlpha',0.5,'EdgeColor', [0.2 0.2 0.2], ...
    'LineWidth', 1.8);
b.CData(1,:) = cLearner;
b.CData(2,:) = cNonLearner;

errorbar([1 2], data_for_bar, ...
    [std(discrim_idx_L)/sqrt(length(discrim_idx_L)), ...
    std(discrim_idx_NL)/sqrt(length(discrim_idx_NL))], ...
    'k', 'LineStyle', 'none', 'LineWidth', 2);

scatter(ones(size(discrim_idx_L)), discrim_idx_L, 80, cLearner, 'filled', 'MarkerFaceAlpha', 0.6);
scatter(2*ones(size(discrim_idx_NL)), discrim_idx_NL, 80, cNonLearner, 'filled', 'MarkerFaceAlpha', 0.6);

yline(0, 'k--', 'LineWidth', 1.5);
xticks([1 2]);
xticklabels({'mPFC+', 'mPFC-'});
ylabel('Discrimination Index');
ylim([0,1.2]);
title(sprintf('(CS+-CS-)/(CS++CS-) [p=%.3f]', p3));
box off;

% Panel 3: Total licks per trial
subplot(1, 3, 3);
hold on;

total_plus_L = [discrimination_metrics(learners_mask).total_cs_plus];
total_plus_NL = [discrimination_metrics(~learners_mask).total_cs_plus];
total_minus_L = [discrimination_metrics(learners_mask).total_cs_minus];
total_minus_NL = [discrimination_metrics(~learners_mask).total_cs_minus];

bar(1, mean(total_plus_L), 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.7);
bar(2, mean(total_minus_L), 'FaceColor', [0.8500 0.3250 0.0980], 'FaceAlpha', 0.3);
bar(4, mean(total_plus_NL), 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.7);
bar(5, mean(total_minus_NL), 'FaceColor', [0 0.4470 0.7410], 'FaceAlpha', 0.3);

errorbar([1 2 4 5], ...
    [mean(total_plus_L), mean(total_minus_L), mean(total_plus_NL), mean(total_minus_NL)], ...
    [std(total_plus_L)/sqrt(length(total_plus_L)), ...
    std(total_minus_L)/sqrt(length(total_minus_L)), ...
    std(total_plus_NL)/sqrt(length(total_plus_NL)), ...
    std(total_minus_NL)/sqrt(length(total_minus_NL))], ...
    'k', 'LineStyle', 'none', 'LineWidth', 2);

xticks([1.5 4.5]);
xticklabels({'Learners', 'Non-learners'});
ylabel('Total licks per trial');
title('Overall Licking');
grid on;
box off;

sgtitle('Stage 1 Behavioral Equivalence (Last Day)');



%% ===== Learning Speed Analysis =====

learning_speed = struct();
idx = 0;

for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % Find all classical conditioning days
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    if ~any(isClassical), continue; end

    classical_days = find(isClassical);

    % Get learning day from combined_sig_day_all_protocols
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: no learning day found, skipping\n', aid);
        continue;
    end

    % ===== Calculate cumulative trials until learning =====
    cumulative_trials = 0;
    cumulative_cs_plus_trials = 0;

    for day_idx = 1:learned_day
        D = rd(classical_days(day_idx));

        % Check if cs_labels field exists
        if isfield(D, 'cs_labels') && ~isempty(D.cs_labels)
            cs_plus_mask = D.cs_labels;

            % Count trials for this day
            n_total_trials_this_day = length(cs_plus_mask);
            n_cs_plus_trials_this_day = sum(cs_plus_mask);

            cumulative_trials = cumulative_trials + n_total_trials_this_day;
            cumulative_cs_plus_trials = cumulative_cs_plus_trials + n_cs_plus_trials_this_day;
        end
    end

    % Store
    idx = idx + 1;
    learning_speed(idx).animal_id = aid;
    learning_speed(idx).is_learner = ismember(aid, learner_ids);
    learning_speed(idx).days_to_criterion = learned_day;
    learning_speed(idx).total_classical_days = length(classical_days);
    learning_speed(idx).cumulative_trials = cumulative_trials;
    learning_speed(idx).cumulative_cs_plus_trials = cumulative_cs_plus_trials;
end

% Statistical comparison
learners_mask = [learning_speed.is_learner];

% Days to criterion
days_L = [learning_speed(learners_mask).days_to_criterion];
days_NL = [learning_speed(~learners_mask).days_to_criterion];

% Cumulative trials
trials_L = [learning_speed(learners_mask).cumulative_trials];
trials_NL = [learning_speed(~learners_mask).cumulative_trials];

% Cumulative CS+ trials
cs_trials_L = [learning_speed(learners_mask).cumulative_cs_plus_trials];
cs_trials_NL = [learning_speed(~learners_mask).cumulative_cs_plus_trials];

% Statistical tests - Days
[~, p_days] = ttest2(days_L, days_NL);
p_days_ranksum = ranksum(days_L, days_NL);

% Statistical tests - Total trials
[~, p_trials] = ttest2(trials_L, trials_NL);
p_trials_ranksum = ranksum(trials_L, trials_NL);

% Statistical tests - CS+ trials
[~, p_cs_trials] = ttest2(cs_trials_L, cs_trials_NL);
p_cs_trials_ranksum = ranksum(cs_trials_L, cs_trials_NL);

% ===== Professional Visualization: Boxplots with Overlaid Points =====

figure('Color', 'w', 'Position', [100 100 1500 600]);

% Define positions and styling
positions = [1 2];
data_sets = {
    {days_L, days_NL},
    {trials_L, trials_NL},
    {cs_trials_L, cs_trials_NL}
    };
y_labels = {'Days to Criterion', 'Total Trials to Criterion', 'CS+ Trials to Criterion'};
p_values = [p_days, p_trials, p_cs_trials];
titles = {'Days to Criterion', 'Total Trials', 'CS+ Trials'};

% Color scheme - more saturated, professional
colors = {cLearner, cNonLearner};

for panel = 1:3
    subplot(1, 3, panel);
    hold on;

    data_to_plot = data_sets{panel};

    % Plot each group
    for i = 1:2
        data = data_to_plot{i};
        pos = positions(i);
        col = colors{i};

        % Quartiles
        q25 = prctile(data, 25);
        q50 = median(data);
        q75 = prctile(data, 75);
        iqr_val = q75 - q25;

        % Whiskers (1.5*IQR, standard Tukey)
        whisker_low = max(min(data), q25 - 1.5*iqr_val);
        whisker_high = min(max(data), q75 + 1.5*iqr_val);

        % Draw box
        rectangle('Position', [pos-0.15, q25, 0.3, q75-q25], ...
            'FaceColor', col, 'EdgeColor', [0.2 0.2 0.2], ...
            'LineWidth', 1.8, 'FaceAlpha', 0.5);

        % Median line (thicker, white for contrast)
        plot([pos-0.15, pos+0.15], [q50, q50], 'w-', 'LineWidth', 3);
        plot([pos-0.15, pos+0.15], [q50, q50], 'k-', 'LineWidth', 2);

        % Whiskers
        plot([pos, pos], [whisker_low, q25], 'k-', 'LineWidth', 1.5);
        plot([pos, pos], [q75, whisker_high], 'k-', 'LineWidth', 1.5);
        plot([pos-0.08, pos+0.08], [whisker_low, whisker_low], 'k-', 'LineWidth', 1.5);
        plot([pos-0.08, pos+0.08], [whisker_high, whisker_high], 'k-', 'LineWidth', 1.5);

        % Individual points with jitter - smaller, more transparent
        jitter = 0.06 * randn(size(data));
        scatter(pos + jitter, data, 45, col, 'filled', ...
            'MarkerFaceAlpha', 0.4, 'MarkerEdgeColor', 'none');
    end

    % Add significance annotation
    p_val = p_values(panel);
    y_max = max([data_to_plot{1}; data_to_plot{2}],[],'all');
    y_range = range([data_to_plot{1}; data_to_plot{2}]);

    % Determine significance level
    if p_val < 0.001
        sig_text = '***';
    elseif p_val < 0.01
        sig_text = '**';
    elseif p_val < 0.05
        sig_text = '*';
    else
        sig_text = 'n.s.';
    end


    % Formatting
    xlim([0.4 2.6]);
    xticks([1 2]);
    xticklabels({'mPFC+', 'mPFC-'});
    ylabel(y_labels{panel}, 'FontSize', 12, 'FontWeight', 'bold');
    title(titles{panel}, 'FontSize', 13, 'FontWeight', 'bold');

    % Set y-axis limits with padding
    if panel == 1  % Days - force integers
        y_min = floor(min([data_to_plot{1}; data_to_plot{2}],[],'all'));
        y_max_data = ceil(max([data_to_plot{1}; data_to_plot{2}],[],'all'));
        ylim([min(y_min-1), y_max_data+1]);
        yticks(y_min:1:y_max_data);
    else
        ylim([min([data_to_plot{1}; data_to_plot{2}],[],'all')*0.7, y_max*1.2]);
    end


    % Clean axes
    set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
    box off;

    hold off;
end

% Overall title - more subtle
sgtitle('Learning Speed: mPFC+ vs mPFC- Animals', 'FontSize', 15, 'FontWeight', 'bold');

% Tighten layout
set(gcf, 'Units', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperPositionMode', 'Auto', 'PaperUnits', 'inches', 'PaperSize', [pos(3), pos(4)]);



%% Collect RT data for BOTH learners and non-learners

Nshow = 20; % number of trials to show

% Color definitions
cRT = [0.4660 0.6740 0.1880]; % Green for RT

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learner_ids));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate cell arrays
R_list_NL = cell(1, n_non_learners); % Non-learners RT
R_list_L = cell(1, n_learners);      % Learners RT

% Initialize counters for indexing
idx_L = 0;
idx_NL = 0;

% ===== Single loop to collect both groups =====
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % first "static" day (skip big_stim)
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);
    if ~any(is_stage_2), continue; end

    cand = find(is_stage_2);
    firstIdx = cand(3);
    D = rd(firstIdx);

    % required fields
    if ~isfield(D,'all_stim_diff_from_optimal_reward') || isempty(D.all_stim_diff_from_optimal_reward), continue; end
    if ~isfield(D,'cs_labels') || isempty(D.cs_labels), continue; end

    cs_plus_mask = D.cs_labels;
    rt = D.all_stim_diff_from_optimal_reward(cs_plus_mask);

    % trim to first Nshow trials
    n = min(Nshow, numel(rt));
    if n < 1, continue; end

    rt = rt(1:n);

    % Conditional: assign to learner or non-learner list
    if ismember(aid, learner_ids)
        idx_L = idx_L + 1;
        R_list_L{idx_L} = rt;
    else
        idx_NL = idx_NL + 1;
        R_list_NL{idx_NL} = rt;
    end
end

% Trim unused cells
% R_list_L = R_list_L(1:idx_L);
% R_list_NL = R_list_NL(1:idx_NL);

% ===== Stack into matrices (NON-LEARNERS) =====
nNL = numel(R_list_NL);
R_mat_NL = nan(Nshow, nNL);

for i = 1:nNL
    vr = R_list_NL{i};
    R_mat_NL(1:numel(vr), i) = vr;
end

% ===== Stack into matrices (LEARNERS) =====
nL = numel(R_list_L);
R_mat_L = nan(Nshow, nL);

for i = 1:nL
    vr = R_list_L{i};
    R_mat_L(1:numel(vr), i) = vr;
end

% ===== Compute mean ± SEM =====
t = (1:Nshow)';

% Non-learners
mR_NL = median(R_mat_NL, 2, 'omitnan');
nR_NL = sum(~isnan(R_mat_NL), 2);
sR_NL = std(R_mat_NL, 0, 2, 'omitnan') ./ max(1, sqrt(nR_NL));

% Learners
mR_L = median(R_mat_L, 2, 'omitnan');
nR_L = sum(~isnan(R_mat_L), 2);
sR_L = std(R_mat_L, 0, 2, 'omitnan') ./ max(1, sqrt(nR_L));

% ===== Determine shared y-axis limits =====
all_data = [R_mat_NL(:); R_mat_L(:)];
y_min = min(all_data, [], 'omitnan');
y_max = max(all_data, [], 'omitnan');
y_margin = (y_max - y_min) * 0.1; % 10% margin
ylim_shared = [y_min - y_margin, y_max + y_margin];

% ===== Plot side-by-side =====
figure('Color', 'w', 'Position', [100 100 1200 500], ...
    'Name', 'Reaction Time: Learners vs Non-Learners — First Static Day');

% --- Non-Learners Panel ---
subplot(1, 2, 1);
hold on;

xf = [t; flipud(t)];

% RT
yfR_NL = [mR_NL + sR_NL; flipud(mR_NL - sR_NL)];
fill(xf, yfR_NL, cRT, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mR_NL, '-', 'Color', cRT, 'LineWidth', 2, 'DisplayName', 'Reaction Time');

% Add zero line
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

ylabel('RT / diff from optimal (s)');
xlabel('Trial # (first 20)');
title(sprintf('Non-Learners (n=%d)', nNL));
ylim(ylim_shared);
grid on;
box off;
legend('Location', 'best');

% --- Learners Panel ---
subplot(1, 2, 2);
hold on;

% RT
yfR_L = [mR_L + sR_L; flipud(mR_L - sR_L)];
fill(xf, yfR_L, cRT, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(t, mR_L, '-', 'Color', cRT, 'LineWidth', 2, 'DisplayName', 'Reaction Time');

% Add zero line
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

ylabel('RT / diff from optimal (s)');
xlabel('Trial # (first 20)');
title(sprintf('Learners (n=%d)', nL));
ylim(ylim_shared);
grid on;
box off;
legend('Location', 'best');

sgtitle('First Static Day (CS+): Reaction Time');

fprintf('\n===== Summary Statistics =====\n');
fprintf('Non-Learners (n=%d):\n', nNL);
fprintf('  RT: %.3f ± %.3f s\n', mean(mR_NL, 'omitnan'), mean(sR_NL, 'omitnan'));
fprintf('\nLearners (n=%d):\n', nL);
fprintf('  RT: %.3f ± %.3f s\n', mean(mR_L, 'omitnan'), mean(sR_L, 'omitnan'));


%% Collect RT data for BOTH learners and non-learners

Nshow = 40; % number of trials to show

% Color definitions
cRT_L = cLearner;  % Orange for learners
cRT_NL = cNonLearner;      % Blue for non-learners

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learner_ids));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate cell arrays
R_list_NL = cell(1, n_non_learners); % Non-learners RT
R_list_L = cell(1, n_learners);      % Learners RT
animal_ids_NL = cell(1, n_non_learners);
animal_ids_L = cell(1, n_learners);

% Initialize counters for indexing
idx_L = 0;
idx_NL = 0;

% ===== Single loop to collect both groups =====
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % first "static" day (skip big_stim)
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);
    if ~any(is_stage_2), continue; end

    cand = find(is_stage_2);
    firstIdx = cand(1);
    D = rd(firstIdx);

    % required fields
    if ~isfield(D,'all_stim_diff_from_optimal_reward') || isempty(D.all_stim_diff_from_optimal_reward), continue; end
    if ~isfield(D,'cs_labels') || isempty(D.cs_labels), continue; end

    cs_plus_mask = D.cs_labels;
    rt = D.all_stim_diff_from_optimal_reward(cs_plus_mask);

    % trim to first Nshow trials
    n = min(Nshow, numel(rt));
    if n < 1, continue; end

    rt = rt(1:n);

    % Conditional: assign to learner or non-learner list
    if ismember(aid, learner_ids)
        idx_L = idx_L + 1;
        R_list_L{idx_L} = rt;
        animal_ids_L{idx_L} = aid;
    else
        idx_NL = idx_NL + 1;
        R_list_NL{idx_NL} = rt;
        animal_ids_NL{idx_NL} = aid;
    end
end

% Trim unused cells
R_list_L = R_list_L(1:idx_L);
R_list_NL = R_list_NL(1:idx_NL);
animal_ids_L = animal_ids_L(1:idx_L);
animal_ids_NL = animal_ids_NL(1:idx_NL);

% ===== Stack into matrices =====
nNL = numel(R_list_NL);
R_mat_NL = nan(Nshow, nNL);

for i = 1:nNL
    vr = R_list_NL{i};
    R_mat_NL(1:numel(vr), i) = vr;
end

nL = numel(R_list_L);
R_mat_L = nan(Nshow, nL);

for i = 1:nL
    vr = R_list_L{i};
    R_mat_L(1:numel(vr), i) = vr;
end

% ===== Compute mean for overlay =====
t = (1:Nshow)';
mR_NL = mean(R_mat_NL, 2, 'omitnan');
mR_L = mean(R_mat_L, 2, 'omitnan');

% ===== Plot with scatter for individual animals =====
figure('Color', 'w', 'Position', [100 100 1400 500], ...
    'Name', 'Reaction Time: Learners vs Non-Learners — First Static Day');

% --- Non-Learners Panel ---
subplot(1, 2, 1);
hold on;

% Plot individual animals as scatter
for i = 1:nNL
    scatter(t, R_mat_NL(:, i), 40, cRT_NL, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'HandleVisibility', 'off');
end

% Plot mean as thick line
plot(t, mR_NL, '-', 'Color', cRT_NL, 'LineWidth', 3, 'DisplayName', 'Mean RT');

% Add zero line
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Optimal timing');

ylabel('RT / diff from optimal (s)');
xlabel('Trial # (first 20)');
title(sprintf('Non-Learners (n=%d)', nNL));
grid on;
box off;
legend('Location', 'best');
xlim([0 Nshow+1]);

% --- Learners Panel ---
subplot(1, 2, 2);
hold on;

% Plot individual animals as scatter
for i = 1:nL
    scatter(t, R_mat_L(:, i), 40, cRT_L, 'filled', 'MarkerFaceAlpha', 0.3, ...
        'HandleVisibility', 'off');
end

% Plot mean as thick line
plot(t, mR_L, '-', 'Color', cRT_L, 'LineWidth', 3, 'DisplayName', 'Mean RT');

% Add zero line
yline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Optimal timing');

ylabel('RT / diff from optimal (s)');
xlabel('Trial # (first 20)');
title(sprintf('Learners (n=%d)', nL));
grid on;
box off;
legend('Location', 'best');
xlim([0 Nshow+1]);

% Link y-axes for easier comparison
linkaxes([subplot(1,2,1), subplot(1,2,2)], 'y');

sgtitle('First Static Day (CS+): Reaction Time (Individual Animals + Mean)');

fprintf('\n===== Summary Statistics =====\n');
fprintf('Non-Learners (n=%d):\n', nNL);
fprintf('  Trial 1 RT: %.2f ± %.2f s\n', mR_NL(1), std(R_mat_NL(1,:), 'omitnan'));
fprintf('  Trial 20 RT: %.2f ± %.2f s\n', mR_NL(end), std(R_mat_NL(end,:), 'omitnan'));
fprintf('  Overall mean RT: %.2f ± %.2f s\n', mean(mR_NL, 'omitnan'), std(R_mat_NL(:), 'omitnan'));

fprintf('\nLearners (n=%d):\n', nL);
fprintf('  Trial 1 RT: %.2f ± %.2f s\n', mR_L(1), std(R_mat_L(1,:), 'omitnan'));
fprintf('  Trial 20 RT: %.2f ± %.2f s\n', mR_L(end), std(R_mat_L(end,:), 'omitnan'));
fprintf('  Overall mean RT: %.2f ± %.2f s\n', mean(mR_L, 'omitnan'), std(R_mat_L(:), 'omitnan'));



%% Plot individual trial trajectory and bar plot of behavioural component for the first n trials of first static day

% to play around
learner_ids= {'HA005','HA008','HA010','HA012'};

Nshow = 30; % number of trials to show

% Color definitions
cRT_L = cLearner;  % Learners
cRT_NL = cNonLearner; % Non-learners

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learner_ids));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate cell arrays
R_list_NL = cell(1, n_non_learners);
R_list_L = cell(1, n_learners);
animal_ids_NL = cell(1, n_non_learners);
animal_ids_L = cell(1, n_learners);

% Initialize counters
idx_L = 0;
idx_NL = 0;

% ===== Collect data for both groups =====
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % Find all classical conditioning days
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    if ~any(isClassical), continue; end

    classical_days = find(isClassical);

    % Get learning day from combined_sig_day_all_protocols
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: no learning day found, skipping\n', aid);
        continue;
    end


    % First "static" day (skip big_stim)
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);
    if ~any(is_stage_2), continue; end

    cand = find(is_stage_2);
    firstIdx = cand(1);
    D = rd(firstIdx);

    % Required fields
    if ~isfield(D,'all_stim_diff_from_optimal_reward') || isempty(D.all_stim_diff_from_optimal_reward), continue; end
    if ~isfield(D,'cs_labels') || isempty(D.cs_labels), continue; end

    cs_plus_mask = D.cs_labels;

    rt = D.all_stim_diff_from_optimal_reward(cs_plus_mask);


    % Trim to first Nshow trials
    n = min(Nshow, numel(rt));
    if n < 1, continue; end

    rt = rt(1:n);

    % Assign to learner or non-learner list
    if ismember(aid, learner_ids)
        idx_L = idx_L + 1;
        R_list_L{idx_L} = rt;
        animal_ids_L{idx_L} = aid;
    else
        idx_NL = idx_NL + 1;
        R_list_NL{idx_NL} = rt;
        animal_ids_NL{idx_NL} = aid;
    end
end

% Trim unused cells
R_list_L = R_list_L(1:idx_L);
R_list_NL = R_list_NL(1:idx_NL);
animal_ids_L = animal_ids_L(1:idx_L);
animal_ids_NL = animal_ids_NL(1:idx_NL);

% ===== Stack into matrices =====
nNL = numel(R_list_NL);
R_mat_NL = nan(Nshow, nNL);

for i = 1:nNL
    vr = R_list_NL{i};
    R_mat_NL(1:numel(vr), i) = vr;
end

nL = numel(R_list_L);
R_mat_L = nan(Nshow, nL);

for i = 1:nL
    vr = R_list_L{i};
    R_mat_L(1:numel(vr), i) = vr;
end

% ===== Compute statistics =====
t = (1:Nshow)';

% Non-learners - trial-by-trial median
medR_NL = median(R_mat_NL, 2, 'omitnan');

% Learners - trial-by-trial median
medR_L = median(R_mat_L, 2, 'omitnan');

% For bar plot: calculate each animal's median
animal_avg_NL = mean(R_mat_NL, 1, 'omitnan')';  % [nNL × 1]
animal_avg_L = mean(R_mat_L, 1, 'omitnan')';    % [nL × 1]

% Group statistics for bar plot
group_mean_NL = mean(animal_avg_NL, 'omitnan');
group_sem_NL = std(animal_avg_NL, 'omitnan') / sqrt(sum(~isnan(animal_avg_NL)));

group_mean_L = mean(animal_avg_L, 'omitnan');
group_sem_L = std(animal_avg_L, 'omitnan') / sqrt(sum(~isnan(animal_avg_L)));

% ===== Two-panel figure =====
figure('Color', 'w', 'Position', [100 100 1400, 600]);

% --- PANEL 1: Trial-by-trial trajectory (LEFT) ---
subplot(1, 2, 1);
hold on;

% Plot median lines only (no error bands)
plot(t, medR_NL, '-', 'Color', cRT_NL, 'LineWidth', 3.5, ...
    'DisplayName', sprintf('mPFC- (n=%d)', nNL));
plot(t, medR_L, '-', 'Color', cRT_L, 'LineWidth', 3.5, ...
    'DisplayName', sprintf('mPFC+ (n=%d)', nL));
%
% % Add optimal timing reference
% yline(0, 'k--', 'Optimal', 'LineWidth', 2, ...
%     'LabelVerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');

% Formatting
xlabel('Trial Number', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Reaction Time (s from optimal)', 'FontSize', 13, 'FontWeight', 'bold');
title('Trial-by-Trial Trajectory', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11);
xlim([0 Nshow+1]);
% ylim([0 4]);
grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

% --- PANEL 2: Summary bar plot with individual points (RIGHT) ---
subplot(1, 2, 2);
hold on;

positions = [1 2];
bar_width = 0.6;

% Non-learners bar
bar(positions(1), group_mean_NL, bar_width, 'FaceColor', cRT_NL, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, 'FaceAlpha', 0.5);

% Learners bar
bar(positions(2), group_mean_L, bar_width, 'FaceColor', cRT_L, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, 'FaceAlpha', 0.5);

% Error bars (SEM)
errorbar(positions(1), group_mean_NL, group_sem_NL, 'k', 'LineWidth', 2, ...
    'CapSize', 10, 'HandleVisibility', 'off');
errorbar(positions(2), group_mean_L, group_sem_L, 'k', 'LineWidth', 2, ...
    'CapSize', 10, 'HandleVisibility', 'off');

% Individual animal points with jitter
jitter_NL = 0.08 * randn(size(animal_avg_NL));
scatter(positions(1) + jitter_NL, animal_avg_NL, 50, cRT_NL, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

jitter_L = 0.08 * randn(size(animal_avg_L));
scatter(positions(2) + jitter_L, animal_avg_L, 50, cRT_L, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

% Add optimal timing reference
yline(0, 'k--', 'Optimal', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold');

% Formatting
xlim([0.4 2.6]);
xticks(positions);
xticklabels({'mPFC-', 'mPFC+'});
ylabel('Mean RT (s from optimal)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Average RT (First %d Trials)', Nshow), 'FontSize', 13, 'FontWeight', 'bold');
% ylim([0 4]);
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

% Overall title
sgtitle('First Static Day (CS+ Trials): Reaction Time Analysis', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ===== Statistical test =====
[h, p_ttest] = ttest2(animal_avg_NL, animal_avg_L);
[p_ranksum] = ranksum(animal_avg_NL, animal_avg_L);

% Add significance annotation to bar plot
subplot(1, 2, 2);
hold on;
y_max = max([animal_avg_NL; animal_avg_L]);
bracket_y = y_max + 1;

if p_ranksum < 0.001
    sig_text = '***';
elseif p_ranksum < 0.01
    sig_text = '**';
elseif p_ranksum < 0.05
    sig_text = '*';
else
    sig_text = 'n.s.';
end

% Draw significance bracket
plot([1, 1, 2, 2], [bracket_y-0.3, bracket_y, bracket_y, bracket_y-0.3], ...
    'k-', 'LineWidth', 1.5);
text(1.5, bracket_y + 0.2, sig_text, ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
hold off;

% ===== Summary Statistics =====
fprintf('\n===== Summary Statistics =====\n');
fprintf('Non-Learners (mPFC-, n=%d):\n', nNL);
fprintf('  Mean RT: %.2f ± %.2f s (SEM)\n', group_mean_NL, group_sem_NL);
fprintf('  Median RT: %.2f s\n', median(animal_avg_NL, 'omitnan'));

fprintf('\nLearners (mPFC+, n=%d):\n', nL);
fprintf('  Mean RT: %.2f ± %.2f s (SEM)\n', group_mean_L, group_sem_L);
fprintf('  Median RT: %.2f s\n', median(animal_avg_L, 'omitnan'));

fprintf('\n===== Statistical Tests =====\n');
fprintf('Two-sample t-test: p = %.4f %s\n', p_ttest, sig_text);
fprintf('Wilcoxon rank-sum: p = %.4f\n', p_ranksum);
fprintf('Effect size (Cohen''s d): %.2f\n', ...
    (group_mean_L - group_mean_NL) / sqrt((std(animal_avg_L)^2 + std(animal_avg_NL)^2) / 2));

%%  Plot individual trial trajectory and bar plot of discrimination index (CS+ - CS- / CS+ + CS-)

Nshow = 30; % number of trials to show

% Color definitions
cRT_L = cLearner;  % Learners
cRT_NL = cNonLearner; % Non-learners

% Count learners and non-learners
n_learners = sum(ismember(string({behaviour_data.animal_id}), learner_ids));
n_non_learners = numel(behaviour_data) - n_learners;

% Pre-allocate cell arrays
DI_list_NL = cell(1, n_non_learners);
DI_list_L = cell(1, n_learners);
animal_ids_NL = cell(1, n_non_learners);
animal_ids_L = cell(1, n_learners);

% Initialize counters
idx_L = 0;
idx_NL = 0;

% ===== Collect data for both groups =====
for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);

    rd = behaviour_data(a).recording_day;
    if isempty(rd), continue; end

    % Find all classical conditioning days
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    if ~any(isClassical), continue; end

    classical_days = find(isClassical);

    % Get learning day from combined_sig_day_all_protocols
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');

    if isempty(learned_day)
        fprintf('Animal %s: no learning day found, skipping\n', aid);
        continue;
    end


    % First "static" day (skip big_stim)
    is_stage_2 = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static') && ~contains(d.workflow,'big_stim'), rd);
    if ~any(is_stage_2), continue; end

    cand = find(is_stage_2);
    firstIdx = cand(1);
    D = rd(firstIdx);

    % Required fields
    if ~isfield(D,'anticipatory_licks') || isempty(D.anticipatory_licks), continue; end
    if ~isfield(D,'cs_labels') || isempty(D.cs_labels), continue; end

    % Get CS+ and CS- trials
    cs_plus_mask = D.cs_labels;
    cs_plus_licks = D.anticipatory_licks(cs_plus_mask);
    cs_minus_licks = D.anticipatory_licks(~cs_plus_mask);

    % Determine how many paired trials we can make
    n_pairs = min([length(cs_plus_licks), length(cs_minus_licks), Nshow]);

    if n_pairs < 1
        warning('Animal %s: insufficient trials, skipping', aid);
        continue;
    end

    % Take first Nshow trials of each type
    cs_plus_licks = cs_plus_licks(1:n_pairs);
    cs_minus_licks = cs_minus_licks(1:n_pairs);

    % Calculate discrimination index for each paired trial
    % DI = (CS+ - CS-) / (CS+ + CS-)
    discrimination_index = nan(n_pairs, 1);
    for trial = 1:n_pairs
        cs_plus_val = cs_plus_licks(trial);
        cs_minus_val = cs_minus_licks(trial);

        % Calculate DI (avoid division by zero)
        if (cs_plus_val + cs_minus_val) == 0
            discrimination_index(trial) = 0;
        else
            discrimination_index(trial) = (cs_plus_val - cs_minus_val) / (cs_plus_val + cs_minus_val);
        end
    end

    % Assign to learner or non-learner list
    if ismember(aid, learner_ids)
        idx_L = idx_L + 1;
        DI_list_L{idx_L} = discrimination_index;
        animal_ids_L{idx_L} = aid;
    else
        idx_NL = idx_NL + 1;
        DI_list_NL{idx_NL} = discrimination_index;
        animal_ids_NL{idx_NL} = aid;
    end
end

% Trim unused cells
DI_list_L = DI_list_L(1:idx_L);
DI_list_NL = DI_list_NL(1:idx_NL);
animal_ids_L = animal_ids_L(1:idx_L);
animal_ids_NL = animal_ids_NL(1:idx_NL);

% ===== Stack into matrices =====
nNL = numel(DI_list_NL);
DI_mat_NL = nan(Nshow, nNL);

for i = 1:nNL
    vr = DI_list_NL{i};
    DI_mat_NL(1:numel(vr), i) = vr;
end

nL = numel(DI_list_L);
DI_mat_L = nan(Nshow, nL);

for i = 1:nL
    vr = DI_list_L{i};
    DI_mat_L(1:numel(vr), i) = vr;
end

% ===== Compute statistics =====
t = (1:Nshow)';

% Non-learners - trial-by-trial mean
meanDI_NL = mean(DI_mat_NL, 2, 'omitnan');

% Learners - trial-by-trial mean
meanDI_L = mean(DI_mat_L, 2, 'omitnan');

% For bar plot: calculate each animal's average DI across all trials
animal_avg_NL = mean(DI_mat_NL, 1, 'omitnan')';  % [nNL × 1]
animal_avg_L = mean(DI_mat_L, 1, 'omitnan')';    % [nL × 1]

% Group statistics for bar plot
group_mean_NL = mean(animal_avg_NL, 'omitnan');
group_sem_NL = std(animal_avg_NL, 'omitnan') / sqrt(sum(~isnan(animal_avg_NL)));

group_mean_L = mean(animal_avg_L, 'omitnan');
group_sem_L = std(animal_avg_L, 'omitnan') / sqrt(sum(~isnan(animal_avg_L)));

% ===== Two-panel figure =====
figure('Color', 'w', 'Position', [100 100 1400, 600]);

% --- PANEL 1: Trial-by-trial trajectory (LEFT) ---
subplot(1, 2, 1);
hold on;

% Plot mean lines
plot(t, meanDI_NL, '-', 'Color', cRT_NL, 'LineWidth', 3.5, ...
    'DisplayName', sprintf('mPFC- (n=%d)', nNL));
plot(t, meanDI_L, '-', 'Color', cRT_L, 'LineWidth', 3.5, ...
    'DisplayName', sprintf('mPFC+ (n=%d)', nL));

% Add zero line (no discrimination)
yline(0, 'k--', 'No Discrimination', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');

% Formatting
xlabel('Trial Number', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Discrimination Index', 'FontSize', 13, 'FontWeight', 'bold');
title('Trial-by-Trial CS+/CS- Discrimination', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
xlim([0 Nshow+1]);
ylim([-1 1]);
grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

% --- PANEL 2: Summary bar plot with individual points (RIGHT) ---
subplot(1, 2, 2);
hold on;

positions = [1 2];
bar_width = 0.6;

% Non-learners bar
bar(positions(1), group_mean_NL, bar_width, 'FaceColor', cRT_NL, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, 'FaceAlpha', 0.5);

% Learners bar
bar(positions(2), group_mean_L, bar_width, 'FaceColor', cRT_L, ...
    'EdgeColor', [0.2 0.2 0.2], 'LineWidth', 1.5, 'FaceAlpha', 0.5);

% Error bars (SEM)
errorbar(positions(1), group_mean_NL, group_sem_NL, 'k', 'LineWidth', 2, ...
    'CapSize', 10, 'HandleVisibility', 'off');
errorbar(positions(2), group_mean_L, group_sem_L, 'k', 'LineWidth', 2, ...
    'CapSize', 10, 'HandleVisibility', 'off');

% Individual animal points with jitter
jitter_NL = 0.08 * randn(size(animal_avg_NL));
scatter(positions(1) + jitter_NL, animal_avg_NL, 50, cRT_NL, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

jitter_L = 0.08 * randn(size(animal_avg_L));
scatter(positions(2) + jitter_L, animal_avg_L, 50, cRT_L, 'filled', ...
    'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');

% Add zero line
yline(0, 'k--', 'No Discrimination', 'LineWidth', 2, ...
    'LabelVerticalAlignment', 'bottom', 'FontSize', 10, 'FontWeight', 'bold');

% Formatting
xlim([0.4 2.6]);
xticks(positions);
xticklabels({'mPFC-', 'mPFC+'});
ylabel('Mean Discrimination Index', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('Average DI (First %d Trials)', Nshow), 'FontSize', 13, 'FontWeight', 'bold');
ylim([-1 1]);
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

% Overall title
sgtitle('First Static Day: CS+/CS- Discrimination Index (Anticipatory Licks)', ...
    'FontSize', 15, 'FontWeight', 'bold');

% ===== Statistical test =====
[h, p_ttest] = ttest2(animal_avg_NL, animal_avg_L);
[p_ranksum] = ranksum(animal_avg_NL, animal_avg_L);

% Test if each group is significantly different from zero (one-sample t-test)
[~, p_NL_vs_zero] = ttest(animal_avg_NL, 0);
[~, p_L_vs_zero] = ttest(animal_avg_L, 0);

% Add significance annotation to bar plot
subplot(1, 2, 2);
hold on;
y_max = max([animal_avg_NL; animal_avg_L]);
y_min = min([animal_avg_NL; animal_avg_L]);
y_range = y_max - y_min;
bracket_y = y_max + 0.15;

if p_ranksum < 0.001
    sig_text = '***';
elseif p_ranksum < 0.01
    sig_text = '**';
elseif p_ranksum < 0.05
    sig_text = '*';
else
    sig_text = 'n.s.';
end

% Draw significance bracket between groups
plot([1, 1, 2, 2], [bracket_y-0.05, bracket_y, bracket_y, bracket_y-0.05], ...
    'k-', 'LineWidth', 1.5);
text(1.5, bracket_y + 0.05, sig_text, ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
hold off;

% ===== Summary Statistics =====
fprintf('\n===== Summary Statistics =====\n');
fprintf('Discrimination Index = (CS+ - CS-) / (CS+ + CS-)\n');
fprintf('Positive values = more licking to CS+ than CS-\n\n');

fprintf('Non-Learners (mPFC-, n=%d):\n', nNL);
fprintf('  Mean DI: %.3f ± %.3f (SEM)\n', group_mean_NL, group_sem_NL);
fprintf('  Median DI: %.3f\n', median(animal_avg_NL, 'omitnan'));
fprintf('  Different from zero? p = %.4f\n', p_NL_vs_zero);

fprintf('\nLearners (mPFC+, n=%d):\n', nL);
fprintf('  Mean DI: %.3f ± %.3f (SEM)\n', group_mean_L, group_sem_L);
fprintf('  Median DI: %.3f\n', median(animal_avg_L, 'omitnan'));
fprintf('  Different from zero? p = %.4f\n', p_L_vs_zero);

fprintf('\n===== Statistical Tests =====\n');
fprintf('Between groups:\n');
fprintf('  Two-sample t-test: p = %.4f %s\n', p_ttest, sig_text);
fprintf('  Wilcoxon rank-sum: p = %.4f\n', p_ranksum);
fprintf('  Effect size (Cohen''s d): %.2f\n', ...
    (group_mean_L - group_mean_NL) / sqrt((std(animal_avg_L)^2 + std(animal_avg_NL)^2) / 2));