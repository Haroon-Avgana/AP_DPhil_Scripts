%% What this script includes:

% Fixes:

% Examine how you can create the events for the PSTHs in a better way
% Modify the code so it can accomodate for 3 groups as well


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
% bin_centers = time_bins(1:end-1) + bin_size / 2;

% Define the surround time based on the packed data
surround_time = [-5,5];
surround_sample_rate = 100;
surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);



%%

mPFC_plus_grp = {'HA005','HA008','HA010','HA011','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA009','HA014','HA015'};
non_learners_grp = {'HA006','HA013','AP030','AP031','AP032'};

all_ids = {behaviour_data.animal_id};

target_workflow = {'visual_operant_lick_two_stim_static'};

% Pre-allocate containers for THREE groups
% Store session-level data (no averaging yet)
mPFC_plus_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});
mPFC_minus_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});
non_learner_sessions = struct('plus_pre', {}, 'plus_post', {}, 'minus_pre', {}, 'minus_post', {});

% Define fields
rew_field = 'avg_psth_rewarded_stim_on';
non_rew_field = 'avg_psth_non_rewarded_stim_on';

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

%% Hierarchical averaging: Sessions → Pre/Post → Animals → Group

function [group_mean, group_sem] = compute_group_stats(group_sessions, condition)
    % condition: 'plus_pre', 'plus_post', 'minus_pre', 'minus_post'
    
    n_animals = length(group_sessions);
    animal_means = [];
    animal_sems = [];
    
    for a = 1:n_animals
        sessions = group_sessions(a).(condition);
        
        if isempty(sessions)
            continue;
        end
        
        % Step 1: Average across sessions for this animal
        animal_mean = mean(sessions, 1, 'omitnan');
        
        % Step 2: SEM across sessions for this animal
        animal_sem = std(sessions, 0, 1, 'omitnan') / sqrt(size(sessions, 1));
        
        animal_means = [animal_means; animal_mean];
        animal_sems = [animal_sems; animal_sem];
    end
    
    % Step 3: Average across animals
    group_mean = mean(animal_means, 1, 'omitnan');
    
    % Step 4: Propagate SEM across animals
    % Pooled variance: mean of squared SEMs (converts SEM back to variance, averages, converts back)
    pooled_var = mean(animal_sems.^2, 1, 'omitnan');
    group_sem = sqrt(pooled_var) / sqrt(size(animal_means, 1));
end

% Compute stats for all groups and conditions
[mean_mPFC_plus_pre_CS_plus, sem_mPFC_plus_pre_CS_plus] = compute_group_stats(mPFC_plus_sessions, 'plus_pre');
[mean_mPFC_plus_post_CS_plus, sem_mPFC_plus_post_CS_plus] = compute_group_stats(mPFC_plus_sessions, 'plus_post');
[mean_mPFC_plus_pre_CS_minus, sem_mPFC_plus_pre_CS_minus] = compute_group_stats(mPFC_plus_sessions, 'minus_pre');
[mean_mPFC_plus_post_CS_minus, sem_mPFC_plus_post_CS_minus] = compute_group_stats(mPFC_plus_sessions, 'minus_post');

[mean_mPFC_minus_pre_CS_plus, sem_mPFC_minus_pre_CS_plus] = compute_group_stats(mPFC_minus_sessions, 'plus_pre');
[mean_mPFC_minus_post_CS_plus, sem_mPFC_minus_post_CS_plus] = compute_group_stats(mPFC_minus_sessions, 'plus_post');
[mean_mPFC_minus_pre_CS_minus, sem_mPFC_minus_pre_CS_minus] = compute_group_stats(mPFC_minus_sessions, 'minus_pre');
[mean_mPFC_minus_post_CS_minus, sem_mPFC_minus_post_CS_minus] = compute_group_stats(mPFC_minus_sessions, 'minus_post');

[mean_non_learner_pre_CS_plus, sem_non_learner_pre_CS_plus] = compute_group_stats(non_learner_sessions, 'plus_pre');
[mean_non_learner_post_CS_plus, sem_non_learner_post_CS_plus] = compute_group_stats(non_learner_sessions, 'plus_post');
[mean_non_learner_pre_CS_minus, sem_non_learner_pre_CS_minus] = compute_group_stats(non_learner_sessions, 'minus_pre');
[mean_non_learner_post_CS_minus, sem_non_learner_post_CS_minus] = compute_group_stats(non_learner_sessions, 'minus_post');

% Plotting
figure('Color','w','Position',[100 100 1200 600]);
t = bin_centers;

% Determine global y-axis
y_global_max = max([mean_mPFC_plus_pre_CS_plus(:); mean_mPFC_plus_post_CS_plus(:); ...
                    mean_mPFC_minus_pre_CS_plus(:); mean_mPFC_minus_post_CS_plus(:); ...
                    mean_non_learner_pre_CS_plus(:); mean_non_learner_post_CS_plus(:)], [], 'omitnan') * 1.1;

% CS+ panel
subplot(2,1,1); hold on;

% mPFC+ Post
if ~isempty(mean_mPFC_plus_post_CS_plus)
    fill([t fliplr(t)], ...
         [mean_mPFC_plus_post_CS_plus + sem_mPFC_plus_post_CS_plus, ...
          fliplr(mean_mPFC_plus_post_CS_plus - sem_mPFC_plus_post_CS_plus)], ...
         cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_mPFC_plus_post_CS_plus, 'Color', cLearner, 'LineWidth', 2.5, ...
         'DisplayName', 'mPFC+ Post');
end

% mPFC- Pre
if ~isempty(mean_mPFC_minus_pre_CS_plus)
    fill([t fliplr(t)], ...
         [mean_mPFC_minus_pre_CS_plus + sem_mPFC_minus_pre_CS_plus, ...
          fliplr(mean_mPFC_minus_pre_CS_plus - sem_mPFC_minus_pre_CS_plus)], ...
         cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_mPFC_minus_pre_CS_plus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
         'DisplayName', 'mPFC- Pre');
end

% Non-learners (if exists)
if ~isempty(mean_non_learner_pre_CS_plus) && any(~isnan(mean_non_learner_pre_CS_plus))
    fill([t fliplr(t)], ...
         [mean_non_learner_pre_CS_plus + sem_non_learner_pre_CS_plus, ...
          fliplr(mean_non_learner_pre_CS_plus - sem_non_learner_pre_CS_plus)], ...
         [0.5 0.5 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_non_learner_pre_CS_plus, 'Color', [0.5 0.5 0.5], 'LineWidth', 2.5, ...
         'DisplayName', 'Non-learners Pre');
end

xline(0, 'k--', 'HandleVisibility', 'off');
legend('show', 'Location', 'best');
xlabel('Time (s)');
ylabel('Activity (CS+)');
ylim([0, y_global_max]);
title('CS+ Responses');

% CS- panel (similar structure)
subplot(2,1,2); hold on;

% mPFC+ Post
if ~isempty(mean_mPFC_plus_post_CS_minus)
    fill([t fliplr(t)], ...
         [mean_mPFC_plus_post_CS_minus + sem_mPFC_plus_post_CS_minus, ...
          fliplr(mean_mPFC_plus_post_CS_minus - sem_mPFC_plus_post_CS_minus)], ...
         cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_mPFC_plus_post_CS_minus, 'Color', cLearner, 'LineWidth', 2.5, ...
         'DisplayName', 'mPFC+ Post');
end

% mPFC- Pre
if ~isempty(mean_mPFC_minus_pre_CS_minus)
    fill([t fliplr(t)], ...
         [mean_mPFC_minus_pre_CS_minus + sem_mPFC_minus_pre_CS_minus, ...
          fliplr(mean_mPFC_minus_pre_CS_minus - sem_mPFC_minus_pre_CS_minus)], ...
         cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(t, mean_mPFC_minus_pre_CS_minus, 'Color', cNonLearner, 'LineWidth', 2.5, ...
         'DisplayName', 'mPFC- Pre');
end

xline(0, 'k--', 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Activity (CS-)');
ylim([0, y_global_max]);
title('CS- Responses');

sgtitle('Group Comparisons: Pre vs Post Learning', 'FontSize', 14, 'FontWeight', 'bold');

%% for pre-post learning PSTHs with propr SEMs


mPFC_plus_grp = {'HA005','HA008','HA010','HA011','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA006','HA009','HA013','HA014','HA015'};
non_learners_grp= {};

all_ids = {behaviour_data.animal_id};

target_workflow = {'visual_operant_lick_two_stim_static'};

% Pre-allocate containers for each group
% Store both means AND variances for proper error propagation
learner_plus_pre_means  = {}; learner_plus_pre_vars  = {};
learner_plus_post_means = {}; learner_plus_post_vars = {};
learner_minus_pre_means = {}; learner_minus_pre_vars = {};
learner_minus_post_means = {}; learner_minus_post_vars = {};

non_plus_pre_means  = {}; non_plus_pre_vars  = {};
non_plus_post_means = {}; non_plus_post_vars = {};
non_minus_pre_means = {}; non_minus_pre_vars = {};
non_minus_post_means = {}; non_minus_post_vars = {};

% Define fields
rew_field = 'avg_psth_rewarded_stim_on';
non_rew_field = 'avg_psth_non_rewarded_stim_on';

% Loop over animals
for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    isLearner = ismember(animal_id, mPFC_plus_grp);
    
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
    
    % Split into pre/post
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
    
    % --- Compute per-animal mean AND variance (across sessions) ---
    if isempty(plus_pre)
        avg_plus_pre = nan(1, T);
        var_plus_pre = nan(1, T);
    else
        avg_plus_pre = mean(plus_pre, 1, 'omitnan');
        var_plus_pre = var(plus_pre, 0, 1, 'omitnan');
    end
    
    if isempty(plus_post)
        avg_plus_post = nan(1, T);
        var_plus_post = nan(1, T);
    else
        avg_plus_post = mean(plus_post, 1, 'omitnan');
        var_plus_post = var(plus_post, 0, 1, 'omitnan');
    end
    
    if isempty(minus_pre)
        avg_minus_pre = nan(1, T);
        var_minus_pre = nan(1, T);
    else
        avg_minus_pre = mean(minus_pre, 1, 'omitnan');
        var_minus_pre = var(minus_pre, 0, 1, 'omitnan');
    end
    
    if isempty(minus_post)
        avg_minus_post = nan(1, T);
        var_minus_post = nan(1, T);
    else
        avg_minus_post = mean(minus_post, 1, 'omitnan');
        var_minus_post = var(minus_post, 0, 1, 'omitnan');
    end
    
    % Store means and variances
    if isLearner
        learner_plus_pre_means{end+1, 1} = avg_plus_pre;
        learner_plus_pre_vars{end+1, 1} = var_plus_pre;
        learner_plus_post_means{end+1, 1} = avg_plus_post;
        learner_plus_post_vars{end+1, 1} = var_plus_post;
        learner_minus_pre_means{end+1, 1} = avg_minus_pre;
        learner_minus_pre_vars{end+1, 1} = var_minus_pre;
        learner_minus_post_means{end+1, 1} = avg_minus_post;
        learner_minus_post_vars{end+1, 1} = var_minus_post;
    else
        non_plus_pre_means{end+1, 1} = avg_plus_pre;
        non_plus_pre_vars{end+1, 1} = var_plus_pre;
        non_plus_post_means{end+1, 1} = avg_plus_post;
        non_plus_post_vars{end+1, 1} = var_plus_post;
        non_minus_pre_means{end+1, 1} = avg_minus_pre;
        non_minus_pre_vars{end+1, 1} = var_minus_pre;
        non_minus_post_means{end+1, 1} = avg_minus_post;
        non_minus_post_vars{end+1, 1} = var_minus_post;
    end
end

% Convert to matrices
learner_pre_CS_plus = cat(1, learner_plus_pre_means{:});
learner_pre_CS_plus_vars = cat(1, learner_plus_pre_vars{:});
learner_post_CS_plus = cat(1, learner_plus_post_means{:});
learner_post_CS_plus_vars = cat(1, learner_plus_post_vars{:});
learner_pre_CS_minus = cat(1, learner_minus_pre_means{:});
learner_pre_CS_minus_vars = cat(1, learner_minus_pre_vars{:});
learner_post_CS_minus = cat(1, learner_minus_post_means{:});
learner_post_CS_minus_vars = cat(1, learner_minus_post_vars{:});

non_learner_pre_CS_plus = cat(1, non_plus_pre_means{:});
non_learner_pre_CS_plus_vars = cat(1, non_plus_pre_vars{:});
non_learner_post_CS_plus = cat(1, non_plus_post_means{:});
non_learner_post_CS_plus_vars = cat(1, non_plus_post_vars{:});
non_learner_pre_CS_minus = cat(1, non_minus_pre_means{:});
non_learner_pre_CS_minus_vars = cat(1, non_minus_pre_vars{:});
non_learner_post_CS_minus = cat(1, non_minus_post_means{:});
non_learner_post_CS_minus_vars = cat(1, non_minus_post_vars{:});

% Compute grand means and SEMs (pooling variances properly)
n_learner = size(learner_pre_CS_plus, 1);
n_non_learner = size(non_learner_pre_CS_plus, 1);

% Grand means
mean_learner_pre_CS_plus = mean(learner_pre_CS_plus, 1, 'omitnan');
mean_learner_post_CS_plus = mean(learner_post_CS_plus, 1, 'omitnan');
mean_learner_pre_CS_minus = mean(learner_pre_CS_minus, 1, 'omitnan');
mean_learner_post_CS_minus = mean(learner_post_CS_minus, 1, 'omitnan');

mean_non_learner_pre_CS_plus = mean(non_learner_pre_CS_plus, 1, 'omitnan');
mean_non_learner_post_CS_plus = mean(non_learner_post_CS_plus, 1, 'omitnan');
mean_non_learner_pre_CS_minus = mean(non_learner_pre_CS_minus, 1, 'omitnan');
mean_non_learner_post_CS_minus = mean(non_learner_post_CS_minus, 1, 'omitnan');

% Pooled SEMs (averaging variances, then dividing by n_animals)
sem_learner_pre_CS_plus = sqrt(mean(learner_pre_CS_plus_vars, 1, 'omitnan')) / sqrt(n_learner);
sem_learner_post_CS_plus = sqrt(mean(learner_post_CS_plus_vars, 1, 'omitnan')) / sqrt(n_learner);
sem_learner_pre_CS_minus = sqrt(mean(learner_pre_CS_minus_vars, 1, 'omitnan')) / sqrt(n_learner);
sem_learner_post_CS_minus = sqrt(mean(learner_post_CS_minus_vars, 1, 'omitnan')) / sqrt(n_learner);

sem_non_learner_pre_CS_plus = sqrt(mean(non_learner_pre_CS_plus_vars, 1, 'omitnan')) / sqrt(n_non_learner);
sem_non_learner_post_CS_plus = sqrt(mean(non_learner_post_CS_plus_vars, 1, 'omitnan')) / sqrt(n_non_learner);
sem_non_learner_pre_CS_minus = sqrt(mean(non_learner_pre_CS_minus_vars, 1, 'omitnan')) / sqrt(n_non_learner);
sem_non_learner_post_CS_minus = sqrt(mean(non_learner_post_CS_minus_vars, 1, 'omitnan')) / sqrt(n_non_learner);

% [Rest of plotting code remains the same...]