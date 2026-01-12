%% What this script includes:

% Fixes:

% Examine how you can create the events for the PSTHs in a better way
% Modify the code so it can accomodate for 3 groups as well


%% Run this to undock the figures
set(groot, "defaultFigureWindowStyle", "normal");

%% Load variables and define colours

load("C:\Users\havgana\Desktop\DPhil\packaged_data\behaviour_structure_all_animals.mat") % sig days
load("C:\Users\havgana\Desktop\DPhil\packaged_data\combined_sig_day_all_protocols.mat") % behaviour

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

%% Create indexing variable for workflows

% make protocol index (n all days x workflow number ordered)
workflow_animal = cellfun(@(x) {x.workflow},{behaviour_data.recording_day},'uni',false);
workflow_cat = grp2idx(horzcat(workflow_animal{:}));

% Create a logical learning index variable (n all days x [0,1])
learning_index_animal = vertcat(combined_sig_day_all_protocols{:});



%% Plots MEAN pre-post PSTHs seperated by groups for CS+ and CS- and non-learners 


mPFC_plus_grp = {'HA005','HA008','HA010','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA009','HA011','HA014','HA015'};
non_learners_grp = {'AP030','AP031','AP032'};

all_ids = {behaviour_data.animal_id};

target_workflow = {'visual_operant_lick_two_stim_right_move_big_stim'};
 

% Pre-allocate containers for THREE groups
% Store session-level data (no averaging yet)
mPFC_plus_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});
mPFC_minus_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});
non_learner_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});

% Define fields
rew_field = 'avg_psth_rewarded_stim_on';
non_rew_field = 'avg_psth_non_rewarded_stim_on';

% rew_field= 'avg_psth_rewarded_stim_final_position';
% non_rew_field= 'avg_psth_non_rewarded_stim_final_position';




% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    
    % Determine group membership
    is_mPFC_plus = ismember(animal_id, mPFC_plus_grp);
    is_mPFC_minus = ismember(animal_id, mPFC_minus_grp);
    is_non_learner = ismember(animal_id, non_learners_grp);
    
    days_all = behaviour_data(ai).recording_day;
    
    % Select days matching workflow
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow,target_workflow), days_all);
    validDays = days_all(isValid);
    
    if isempty(validDays)
        warning('Animal %s: no recording found, skipping', behaviour_data(ai).animal_id);
        continue;
    end
    
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);
    
    % Find first significant day
    ld = find(sigDays, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found', behaviour_data(ai).animal_id);
        ld = numel(validDays) + 1;
    end
    
    % Collect SESSION-LEVEL data (no averaging yet)
    plus_pre_sessions = [];
    plus_post_sessions = [];
    minus_pre_sessions = [];
    minus_post_sessions = [];
    
    for d = 1:numel(validDays)
        if d < ld
            plus_pre_sessions = [plus_pre_sessions; validDays(d).(rew_field)];
            minus_pre_sessions = [minus_pre_sessions; validDays(d).(non_rew_field)];
        else
            plus_post_sessions = [plus_post_sessions; validDays(d).(rew_field)];
            minus_post_sessions = [minus_post_sessions; validDays(d).(non_rew_field)];
        end
    end
    
    % Store session-level data for this animal
    animal_data.plus_pre = plus_pre_sessions;
    animal_data.plus_post = plus_post_sessions;
    animal_data.minus_pre = minus_pre_sessions;
    animal_data.minus_post = minus_post_sessions;
    
    % Assign to appropriate group
    if is_mPFC_plus
        mPFC_plus_sessions(end+1) = animal_data;
    elseif is_mPFC_minus
        mPFC_minus_sessions(end+1) = animal_data;
    elseif is_non_learner
        non_learner_sessions(end+1) = animal_data;
    end
end


% Compute stats for all groups and conditions
[mean_mPFC_plus_pre_CS_plus, sem_mPFC_plus_pre_CS_plus] = ha.helper_func.compute_group_stats(mPFC_plus_sessions, 'plus_pre');
[mean_mPFC_plus_post_CS_plus, sem_mPFC_plus_post_CS_plus] = ha.helper_func.compute_group_stats(mPFC_plus_sessions, 'plus_post');
[mean_mPFC_plus_pre_CS_minus, sem_mPFC_plus_pre_CS_minus] = ha.helper_func.compute_group_stats(mPFC_plus_sessions, 'minus_pre');
[mean_mPFC_plus_post_CS_minus, sem_mPFC_plus_post_CS_minus] = ha.helper_func.compute_group_stats(mPFC_plus_sessions, 'minus_post');

[mean_mPFC_minus_pre_CS_plus, sem_mPFC_minus_pre_CS_plus] = ha.helper_func.compute_group_stats(mPFC_minus_sessions, 'plus_pre');
[mean_mPFC_minus_post_CS_plus, sem_mPFC_minus_post_CS_plus] = ha.helper_func.compute_group_stats(mPFC_minus_sessions, 'plus_post');
[mean_mPFC_minus_pre_CS_minus, sem_mPFC_minus_pre_CS_minus] = ha.helper_func.compute_group_stats(mPFC_minus_sessions, 'minus_pre');
[mean_mPFC_minus_post_CS_minus, sem_mPFC_minus_post_CS_minus] = ha.helper_func.compute_group_stats(mPFC_minus_sessions, 'minus_post');

[mean_non_learner_pre_CS_plus, sem_non_learner_pre_CS_plus] = ha.helper_func.compute_group_stats(non_learner_sessions, 'plus_pre');
[mean_non_learner_post_CS_plus, sem_non_learner_post_CS_plus] = ha.helper_func.compute_group_stats(non_learner_sessions, 'plus_post');
[mean_non_learner_pre_CS_minus, sem_non_learner_pre_CS_minus] = ha.helper_func.compute_group_stats(non_learner_sessions, 'minus_pre');
[mean_non_learner_post_CS_minus, sem_non_learner_post_CS_minus] = ha.helper_func.compute_group_stats(non_learner_sessions, 'minus_post');

% Plotting
figure('Color','w','Position',[100 100 1200 600]);
t = bin_centers;

% Determine global y-axis
y_global_max_cs_plus = max([mean_mPFC_plus_pre_CS_plus(:); mean_mPFC_plus_post_CS_plus(:); ...
                    mean_mPFC_minus_pre_CS_plus(:); mean_mPFC_minus_post_CS_plus(:); ...
                    mean_non_learner_pre_CS_plus(:); mean_non_learner_post_CS_plus(:)], [], 'omitnan') * 1.1;

y_global_max_cs_minus = max([mean_mPFC_plus_pre_CS_plus(:); mean_mPFC_plus_post_CS_plus(:); ...
                    mean_mPFC_minus_pre_CS_plus(:); mean_mPFC_minus_post_CS_plus(:); ...
                    mean_non_learner_pre_CS_plus(:); mean_non_learner_post_CS_plus(:)], [], 'omitnan') * 1.1;
% CS+ panel
subplot(2,1,1); hold on;

% % mPFC+ Pre
% if ~isempty(mean_mPFC_plus_pre_CS_plus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_plus_pre_CS_plus + sem_mPFC_plus_pre_CS_plus, ...
%           fliplr(mean_mPFC_plus_pre_CS_plus - sem_mPFC_plus_pre_CS_plus)], ...
%          cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_plus_pre_CS_plus, 'Color', cLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC+ Pre',LineStyle='--');
% end

% 
% % mPFC+ Post
% if ~isempty(mean_mPFC_plus_post_CS_plus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_plus_post_CS_plus + sem_mPFC_plus_post_CS_plus, ...
%           fliplr(mean_mPFC_plus_post_CS_plus - sem_mPFC_plus_post_CS_plus)], ...
%          cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_plus_post_CS_plus, 'Color', cLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC+ Post');
% end


% % mPFC- Pre
% if ~isempty(mean_mPFC_minus_pre_CS_plus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_minus_pre_CS_plus + sem_mPFC_minus_pre_CS_plus, ...
%           fliplr(mean_mPFC_minus_pre_CS_plus - sem_mPFC_minus_pre_CS_plus)], ...
%          cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_minus_pre_CS_plus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC- Pre',LineStyle='--');
% end

% % mPFC- Post
% if ~isempty(mean_mPFC_minus_post_CS_plus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_minus_post_CS_plus + sem_mPFC_minus_post_CS_plus, ...
%           fliplr(mean_mPFC_minus_post_CS_plus - sem_mPFC_minus_post_CS_plus)], ...
%          cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_minus_post_CS_plus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC- Post');
% end


% Non-learners Pre
if ~isempty(mean_non_learner_pre_CS_plus) && any(~isnan(mean_non_learner_pre_CS_plus))
    fill([t fliplr(t)], ...
         [mean_non_learner_pre_CS_plus + sem_non_learner_pre_CS_plus, ...
          fliplr(mean_non_learner_pre_CS_plus - sem_non_learner_pre_CS_plus)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_non_learner_pre_CS_plus, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
         'DisplayName', 'Non-learners Pre',LineStyle='--');
end

% Non-learners Post
if ~isempty(mean_non_learner_post_CS_plus)
    fill([t fliplr(t)], ...
         [mean_non_learner_post_CS_plus + sem_non_learner_post_CS_plus, ...
          fliplr(mean_non_learner_post_CS_plus - sem_non_learner_post_CS_plus)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_non_learner_post_CS_plus, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
         'DisplayName', 'Non-learners Post');
end

xline(0, 'k--', 'HandleVisibility', 'off');
legend('show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Mean Lick Rate (Hz)');
ylim([0, y_global_max_cs_plus]);
title('CS+ Responses');

% CS- panel (similar structure)
subplot(2,1,2); hold on;
% 
% % mPFC+ Pre
% if ~isempty(mean_mPFC_plus_pre_CS_minus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_plus_pre_CS_minus + sem_mPFC_plus_pre_CS_minus, ...
%           fliplr(mean_mPFC_plus_pre_CS_minus - sem_mPFC_plus_pre_CS_minus)], ...
%          cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_plus_pre_CS_minus, 'Color', cLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC+ Pre',LineStyle='--');
% end


% mPFC+ Post
if ~isempty(mean_mPFC_plus_post_CS_minus)
    fill([t fliplr(t)], ...
         [mean_mPFC_plus_post_CS_minus + sem_mPFC_plus_post_CS_minus, ...
          fliplr(mean_mPFC_plus_post_CS_minus - sem_mPFC_plus_post_CS_minus)], ...
         cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_mPFC_plus_post_CS_minus, 'Color', cLearner, 'LineWidth', 2.5, ...
         'DisplayName', 'mPFC+ Post');
end

% % mPFC- Pre
% if ~isempty(mean_mPFC_minus_pre_CS_minus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_minus_pre_CS_minus + sem_mPFC_minus_pre_CS_minus, ...
%           fliplr(mean_mPFC_minus_pre_CS_minus - sem_mPFC_minus_pre_CS_minus)], ...
%          cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_minus_pre_CS_minus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC- Pre',LineStyle='--');
% end

% % mPFC- Post
% if ~isempty(mean_mPFC_minus_post_CS_minus)
%     fill([t fliplr(t)], ...
%          [mean_mPFC_minus_post_CS_minus + sem_mPFC_minus_post_CS_minus, ...
%           fliplr(mean_mPFC_minus_post_CS_minus - sem_mPFC_minus_post_CS_minus)], ...
%          cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
%     plot(t, mean_mPFC_minus_post_CS_minus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
%          'DisplayName', 'mPFC- Post');
% end
% %

% Non-learners Pre
if ~isempty(mean_non_learner_pre_CS_minus)
    fill([t fliplr(t)], ...
         [mean_non_learner_pre_CS_minus + sem_non_learner_pre_CS_minus, ...
          fliplr(mean_non_learner_pre_CS_minus - sem_non_learner_pre_CS_minus)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_non_learner_pre_CS_minus, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
         'DisplayName', 'Non-learners Pre',LineStyle='--');
end

% Non-learners Post
if ~isempty(mean_non_learner_post_CS_minus)
    fill([t fliplr(t)], ...
         [mean_non_learner_post_CS_minus + sem_non_learner_post_CS_minus, ...
          fliplr(mean_non_learner_post_CS_minus - sem_non_learner_post_CS_minus)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_non_learner_post_CS_minus, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
         'DisplayName', 'Non-learners Post');
end

xline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Mean Lick Rate (Hz)');
ylim([0, y_global_max_cs_minus]);
title('CS- Responses');

sgtitle('Group Comparisons: Pre vs Post Learning', 'FontSize', 14, 'FontWeight', 'bold');

%% Plot ITI lick Rate Aligned to Learning Day Split by mPFC+ , mPFC- and non-learners

mPFC_plus_grp = {'HA005','HA008','HA010','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA009','HA011','HA014','HA015'};
non_learners_grp = {'HA006','HA013','AP030','AP031','AP032'};

% Set parameters
target_workflow = 'visual_operant_lick_two_stim_static';
all_ids = {behaviour_data.animal_id};

% Step 1: Collect session-level data for each animal
% Data structure: animal_sessions{animal_idx}(session_idx) = [rel_day, value, group]

animal_sessions = cell(numel(all_ids), 1);

for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    
    % Determine group membership
    if ismember(animal_id, mPFC_plus_grp)
        group_label = 1;
    elseif ismember(animal_id, mPFC_minus_grp)
        group_label = 2;
    elseif ismember(animal_id, non_learners_grp)
        group_label = 3;
    else
        continue;
    end
    
    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);
    
    if isempty(validDays)
        continue;
    end
    
    % Get learning day
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);
    ld = find(sigDays, 1, 'first');
    
    if isempty(ld)
        ld = ceil(numel(validDays) / 2);
        fprintf('Animal %s: using middle day (%d) for alignment\n', animal_id, ld);
    end
    
    % Calculate relative days
    n = numel(validDays);
    rel_days = (1:n) - ld;
    
    % Collect session data with mean and variance
    session_data = [];
    
    for d = 1:numel(validDays)
        day_data = validDays(d);
        
        if ~isfield(day_data, 'ITI_lick_counts') || ~isfield(day_data, 'ITI_actual_duration')
            continue;
        end
        
        iti_licks = day_data.ITI_lick_counts;
        iti_duration = day_data.ITI_actual_duration;
        
        cs_plus_mask= day_data.cs_labels==1;

        anticipatory_licks= day_data.anticipatory_licks(cs_plus_mask);
        RTs= day_data.all_stim_diff_from_optimal_reward(cs_plus_mask);

        % Calculate mean and variance for this session
        session_mean = nanmean(iti_licks./iti_duration);
        session_var = nanvar(iti_licks./iti_duration);
        

        % Store: [rel_day, mean, variance, group]
        session_data = [session_data; rel_days(d), session_mean, session_var, group_label];
    end
    
    animal_sessions{ai} = session_data;
end

% Step 2: Aggregate to group-level
% Data structure: group_data{group_idx}(day_idx) = [rel_day, mean, sem, n_animals]

group_data = cell(3, 1);
min_animals = 3;  % Minimum number of animals required

for g = 1:3
    % Collect all relative days for this group
    all_rel_days = [];
    
    for ai = 1:numel(all_ids)
        if isempty(animal_sessions{ai})
            continue;
        end
        
        session_data = animal_sessions{ai};
        if session_data(1, 4) == g  % Check group (4th column)
            all_rel_days = [all_rel_days; session_data(:, 1)];
        end
    end
    
    unique_days = unique(all_rel_days);
    group_summary = [];
    
    for di = 1:length(unique_days)
        day_val = unique_days(di);
        
        % Collect animal means and variances for this relative day
        animal_means = [];
        animal_vars = [];
        
        for ai = 1:numel(all_ids)
            if isempty(animal_sessions{ai})
                continue;
            end
            
            session_data = animal_sessions{ai};
            
            % Check if this animal is in the current group
            if session_data(1, 4) == g

                % Find data for this relative day
                day_idx = find(session_data(:, 1) == day_val);
                
                if ~isempty(day_idx) % only if it has that relative day
                    animal_means = [animal_means; session_data(day_idx, 2)];  % Mean (column 2)
                    animal_vars = [animal_vars; session_data(day_idx, 3)];    % Variance (column 3)
                end
            end
        end
        
        % Only include this day if we have enough animals
        n_animals_this_day = length(animal_means);
        if n_animals_this_day >= min_animals
            % Mean across animals
            group_mean = mean(animal_means, 'omitnan');

            % Between-animal SEM (standard error of the mean across animals)
            group_sem = std(animal_means, 'omitnan') / sqrt(n_animals_this_day);

            % Store: [rel_day, mean, sem, n_animals]
            group_summary = [group_summary; day_val, group_mean, group_sem, n_animals_this_day];
        end
    end
    
     % Sort by relative day (only if we have data)
    if ~isempty(group_summary)
        [~, sort_idx] = sort(group_summary(:, 1));
        group_data{g} = group_summary(sort_idx, :);
    end
    
end

% Plotting

groups = {mPFC_plus_grp, mPFC_minus_grp, non_learners_grp};
group_names = {'mPFC+', 'mPFC-', 'Non-learners'};
group_colors = {cLearner, cNonLearner, [0.5 0.5 0.5]};

figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

for g = 1:3
    if ~isempty(group_data{g})
        days = group_data{g}(:, 1);
        means = group_data{g}(:, 2);
        sems = group_data{g}(:, 3);
        
        % Error band
        fill([days; flipud(days)], ...
             [means + sems; flipud(means - sems)], ...
             group_colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Mean line
        plot(days, means, '-o', 'Color', group_colors{g}, ...
             'LineWidth', 2.5, 'MarkerSize', 8, 'MarkerFaceColor', group_colors{g}, ...
             'DisplayName', sprintf('%s (n=%d)', group_names{g}, numel(groups{g})));
    end
end

xline(0, '--k', 'Learning Day', 'LineWidth', 2, ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 11, 'FontWeight', 'bold', 'HandleVisibility', 'off');

xlabel('Days Relative to Learning', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('ITI Lick Rate (Licks/s)', 'FontSize', 13, 'FontWeight', 'bold');
title('ITI Licking Aligned to Learning Day', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);

set(gca, 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out');
hold off;

