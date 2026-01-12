% Script

%% 11/12/2025 Files that I can load to work on cohort until HA015

% new_packaging_test_task_widefield
% behaviour_structure_all_animals_11_12_25
% the singificance one is incorrect as it set for the 1.5s anticipaotry
% window

%% Load variables

load('C:\Users\havgana\Desktop\DPhil\packaged_data\task_widefield_08_01_26.mat'); % wf task
load('C:\Users\havgana\Desktop\DPhil\packaged_data\behaviour_structure_all_animals_11_12_25.mat') % behaviour
load('C:\Users\havgana\Desktop\DPhil\packaged_data\combined_sig_day_all_protocols_07_30.mat') % sig days

% load master_U
wf_U= plab.wf.load_master_U;

added_time = -1:0.03:3;
added_time_Kernel=  fliplr((-10:100)/30); % kernel

%% Filter our big stim protocols

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

% opts.big_names = { ...
%   'visual_operant_lick_two_stim_right_move_big_stim', ...
%   'visual_operant_lick_two_stim_static_big_stim', ...
%   'stim_wheel*'
% };


opts.big_names = { ...
  'stim_wheel*'
};

% protocols to drop AFTER the first big_stim, if drop_following=true
% opts.following_names = {'visual_operant_lick_two_stim_static'};
opts.following_names = {'stim_wheel*'};  % For HA014 for now

opts.drop_following   = false;   % <--- true if to drop following ; false otherwise
opts.error_on_unknown = true;   % error if a workflow isn't in protocol_list

% ------------------------------------------

animal_list = {behaviour_data.animal_id};

for animal_idx = 1:numel(animal_list)
    rd = behaviour_data(animal_idx).recording_day;
    wf = task_data(animal_idx).widefield;   % optional, if present
    

    [keptLocalIdx] = ha.helper_func.select_sessions_by_protocol(rd, opts);

    if isempty(keptLocalIdx)
        fprintf('Animal %s: nothing kept (had %d sessions)\n', ...
            animal_list{animal_idx}, numel(rd));
        continue
    end

    % apply to both structures (keep them in sync)
    behaviour_data(animal_idx).recording_day = rd(keptLocalIdx);
    if numel(wf) == numel(rd)
        % task_data(animal_idx).widefield = wf(keptLocalIdx);
    else
        warning('Animal %s: widefield length %d != recording_day length %d. Skipping WF filter.', ...
            animal_list{animal_idx}, numel(wf), numel(rd));
    end

    fprintf('\nAnimal %s: kept %d/%d sessions %s', ...
        animal_list{animal_idx}, numel(keptLocalIdx), numel(rd));
end


%% Set up variables for averaged activity 

% only using n components
n_components=200;
kernel_n_components=100;


% make protocol index (n all days x workflow number ordered)
workflow_animal = cellfun(@(x) {x.workflow},{behaviour_data.recording_day},'uni',false);
workflow_cat = grp2idx(horzcat(workflow_animal{:}));

% Create a logical learning index variable (n all days x [0,1])
learning_index_animal = vertcat(combined_sig_day_all_protocols{:});

% Creates an animal index (n all days x animal number ordered)
widefield_animal_idx = grp2idx(cell2mat(cellfun(@(animal,wf) repmat(animal,length(wf),1), ...
    {task_data.animal},{task_data.widefield},'uni',false)'));

% get all the widefield data
widefield_cat= cat(2,task_data.widefield); % concat

% for V - create a V x T x all days variable - rewarded stim
rewarded_stim_v = cellfun(@(x) mean(x, 3), {widefield_cat.rewarded_stim_on_aligned_V}, 'UniformOutput', false); % average across trials
rewarded_stim_v_stacked_data = cat(3, rewarded_stim_v{:});

% for V - non rewarded
non_rewarded_stim_v = cellfun(@(x) mean(x, 3), {widefield_cat.non_rewarded_stim_onset_aligned_V}, 'UniformOutput', false); % average across trials
non_rewarded_stim_v_stacked_data = cat(3, non_rewarded_stim_v{:});

% for V- rewarded stim start move (only relevant to right move protocol)
rewarded_stim_start_move_v= cellfun(@(x) mean(x, 3), {widefield_cat.rewarded_stim_start_to_move_aligned_V}, 'UniformOutput', false); % average across trials
rewarded_stim_start_move_v_stacked= cat(3, rewarded_stim_start_move_v{:});

% for V- reward times
reward_times= cellfun(@(x) mean(x, 3), {widefield_cat.reward_times_aligned_V}, 'UniformOutput', false); % average across trials
reward_times_stacked= cat(3, reward_times{:});


% For kernels
rewarded_stim_kernel= cat(3,widefield_cat.rewarded_stim_on_aligned_kernel);
non_rewarded_stim_kernel= cat(3,widefield_cat.non_rewarded_stim_on_aligned_kernel);
rewarded_stim_start_move_kernel= cat(3,widefield_cat.rewarded_stim_start_to_move_aligned_kernel);
rewarded_stim_final_position_kernel= cat(3,widefield_cat.rewarded_stim_final_position_aligned_kernel);

reward_times_kernel= cat(3,widefield_cat.reward_times_aligned_kernel);


animal_list= {behaviour_data.animal_id}; % get animal list
unique_workflow_labels= unique(horzcat(workflow_animal{:}),"stable"); % unique workflows by order


% ROI masks
mPFC_data = load('C:\Users\havgana\Desktop\DPhil\Coding Scripts\ROIs\mPFC_ROI.mat');
left_mPFC_ROI_mask = mPFC_data.new_pfc_roi_mask;
right_mPFC_ROI_mask = fliplr(left_mPFC_ROI_mask);

ViS_ROI_data = load("C:\Users\havgana\Desktop\DPhil\Coding Scripts\ROIs\right_ViS_mask.mat");
left_ViS_ROI_mask = ViS_ROI_data.visual_cortex_mask;
right_ViS_ROI_mask = fliplr(left_ViS_ROI_mask);

% animals that have not learned the static but learned the right
% group_animals = {'DS017','HA006','HA007','HA009','HA011','HA013','HA014','HA015'};
group_animals = {'DS017','HA006','HA007','HA009','HA013','HA014','HA015'};


% Create a list of animal IDs repeated per recording day
animal_ids_all_days = cellfun(@(animal, wf) repmat({animal}, length(wf), 1), ...
    {task_data.animal}, {task_data.widefield}, 'UniformOutput', false);

% Convert to column vector
animal_ids_all_days_stacked = vertcat(animal_ids_all_days{:});

% Create a logical index for whether each row belongs to the group
is_group_animal = ismember(animal_ids_all_days_stacked, group_animals);  % (n_days_all x 1)


% set up two custom colors for learners and non-learners
cLearner    = [0 0.6 0];   % dark green
cNonLearner = [0.6 0 0.6]; % purple


%% Plot CCF videos of pre-post learning seperated by learners vs non-learners

% a logical mask for protocol 1
wf1_idx = ((workflow_cat==1));% length K


% pre post of a selected workflow
pre = nanmean(rewarded_stim_v_stacked_data(:,:,workflow_cat==1 & learning_index_animal==0 & is_group_animal==0),3);
post = nanmean(rewarded_stim_v_stacked_data(:,:,workflow_cat==1 & learning_index_animal==1 & is_group_animal==0),3);

ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:n_components),pre),plab.wf.svd2px(wf_U(:,:,1:n_components),post)],added_time);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;

figure('Color', 'w');
plot(added_time, roi.trace, 'k-', 'LineWidth', 2.5);

xlabel('Time (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ΔF/F', 'FontSize', 12, 'FontWeight', 'bold');
title('Example ROI from ViS Aligned to CS+', 'FontSize', 14, 'FontWeight', 'bold');

% Add CS+ onset reference line
xline(0, 'r--', 'LineWidth', 2, 'DisplayName', 'CS+ onset');
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


% Define event and baseline time windows
event_times = [0, 0.35];
baseline_time_windows = [-0.1, 0];

% Find the corresponding indices in added_time for the time window
time_window_idx = (added_time >= event_times(1) & added_time <= event_times(2));
baseline_mask = (added_time >= baseline_time_windows(1) & added_time <= baseline_time_windows(2));


pre_time_window=pre(:,time_window_idx);
post_time_window=post(:,time_window_idx);

pre_baseline= mean(pre(:,baseline_mask),2);
post_baseline=mean(post(:,baseline_mask),2);

pre_baseline_sub=mean(pre_time_window- pre_baseline,2);
post_baseline_sub= mean(post_time_window- post_baseline,2);


% Plot learners
figure;
imagesc(plab.wf.svd2px(wf_U(:,:,1:n_components),pre_baseline_sub));
axis image off;
colormap(ap.colormap('PWG'));
title(('Pre-Learning %s V activity — Learners'));
ap.wf_draw('ccf');

figure;
imagesc(plab.wf.svd2px(wf_U(:,:,1:n_components),post_baseline_sub));
axis image off;
colormap(ap.colormap('PWG'));
title(('Post-Learning %s V activity — Learners'));
ap.wf_draw('ccf');

%% Plot CCF videos for individual animals

animal_list;

% pre post of workflow 3 (right move)
pre = nanmean(rewarded_stim_kernel(:,:,workflow_cat==1 & learning_index_animal==0 & widefield_animal_idx==8),3);
post = nanmean(rewarded_stim_kernel(:,:,workflow_cat==1 & learning_index_animal==1 & widefield_animal_idx==8),3);


ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:kernel_n_components),pre),plab.wf.svd2px(wf_U(:,:,1:kernel_n_components),post)],added_time_Kernel);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;

%% Plot indivudal animal plot
%% Plot ROI traces for learners and non-learners of averages including individual lines aligned to different event - also includes peak activity changes

% define the variables to plot
data_to_plot= rewarded_stim_v_stacked_data;
ROI_to_plot= left_mPFC_ROI_mask;
time_to_plot= added_time;
curr_components= n_components;

% Get unique animal IDs for both groups
non_learner_animal_ids = unique(widefield_animal_idx(is_group_animal == 1));
learner_animal_ids = unique(widefield_animal_idx(is_group_animal == 0));

% Colors
cNonLearnerLight = [cNonLearner 0.3];
cLearnerLight = [cLearner 0.3];

% ========== NON-LEARNERS ==========
% Storage for individual animal traces
individual_nl_pre_traces = cell(length(non_learner_animal_ids), 1);
individual_nl_post_traces = cell(length(non_learner_animal_ids), 1);

% Extract traces for each non-learner animal
for i = 1:length(non_learner_animal_ids)
    animal_id = non_learner_animal_ids(i);
    
    % Pre-learning
    pre_data = nanmean(data_to_plot(:,:, ...
        workflow_cat==1 & learning_index_animal==0 & widefield_animal_idx==animal_id), 3);
    
    if ~all(isnan(pre_data(:)))
        individual_nl_pre_traces{i} = permute( ...
            ap.wf_roi(wf_U(:,:,1:curr_components), pre_data,[],[], ROI_to_plot), ...
            [3,2,1] );
    else
        individual_nl_pre_traces{i} = nan(1, size(data_to_plot, 2));
    end
    
    % Post-learning
    post_data = nanmean(data_to_plot(:,:, ...
        workflow_cat==1 & learning_index_animal==1 & widefield_animal_idx==animal_id), 3);
    
    if ~all(isnan(post_data(:)))
        individual_nl_post_traces{i} = permute( ...
            ap.wf_roi(wf_U(:,:,1:curr_components), post_data,[],[], ROI_to_plot), ...
            [3,2,1] );
    else
        individual_nl_post_traces{i} = nan(1, size(data_to_plot, 2));
    end
end

% Calculate group average for non-learners
non_learner_pre = nanmean(data_to_plot(:,:, ...
    workflow_cat==1 & learning_index_animal==0 & is_group_animal==1), 3);
non_learner_post = nanmean(data_to_plot(:,:, ...
    workflow_cat==1 & learning_index_animal==1 & is_group_animal==1), 3);

non_learner_pre_avg = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_pre,[],[], ROI_to_plot), ...
    [3,2,1] );

non_learner_post_avg = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_post,[],[], ROI_to_plot), ...
    [3,2,1] );

% ========== LEARNERS ==========
% Storage for individual animal traces
individual_l_pre_traces = cell(length(learner_animal_ids), 1);
individual_l_post_traces = cell(length(learner_animal_ids), 1);

% Extract traces for each learner animal
for i = 1:length(learner_animal_ids)
    animal_id = learner_animal_ids(i);
    
    % Pre-learning
    pre_data = nanmean(data_to_plot(:,:, ...
        workflow_cat==1 & learning_index_animal==0 & widefield_animal_idx==animal_id), 3);
    
    if ~all(isnan(pre_data(:)))
        individual_l_pre_traces{i} = permute( ...
            ap.wf_roi(wf_U(:,:,1:curr_components), pre_data,[],[], ROI_to_plot), ...
            [3,2,1] );
    else
        individual_l_pre_traces{i} = nan(1, size(data_to_plot, 2));
    end
    
    % Post-learning
    post_data = nanmean(data_to_plot(:,:, ...
        workflow_cat==1 & learning_index_animal==1 & widefield_animal_idx==animal_id), 3);
    
    if ~all(isnan(post_data(:)))
        individual_l_post_traces{i} = permute( ...
            ap.wf_roi(wf_U(:,:,1:curr_components), post_data,[],[], ROI_to_plot), ...
            [3,2,1] );
    else
        individual_l_post_traces{i} = nan(1, size(data_to_plot, 2));
    end
end

% Calculate group average for learners
learner_pre = nanmean(data_to_plot(:,:, ...
    workflow_cat==1 & learning_index_animal==0 & is_group_animal==0), 3);
learner_post = nanmean(data_to_plot(:,:, ...
    workflow_cat==1 & learning_index_animal==1 & is_group_animal==0), 3);

learner_pre_avg = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_pre,[],[], ROI_to_plot), ...
    [3,2,1] );

learner_post_avg = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_post,[],[], ROI_to_plot), ...
    [3,2,1] );

% Find peak times and values for non-learners
[nl_pre_peak_val, nl_pre_peak_idx] = max(non_learner_pre_avg);
nl_pre_peak_time = time_to_plot(nl_pre_peak_idx);

[nl_post_peak_val, nl_post_peak_idx] = max(non_learner_post_avg);
nl_post_peak_time = time_to_plot(nl_post_peak_idx);

nl_peak_shift = nl_post_peak_time - nl_pre_peak_time;

% Find peak times and values for learners
[l_pre_peak_val, l_pre_peak_idx] = max(learner_pre_avg);
l_pre_peak_time = time_to_plot(l_pre_peak_idx);

[l_post_peak_val, l_post_peak_idx] = max(learner_post_avg);
l_post_peak_time = time_to_plot(l_post_peak_idx);

l_peak_shift = l_post_peak_time - l_pre_peak_time;



% Determine shared y-axis limits across all data
all_traces = [vertcat(individual_nl_pre_traces{:}); vertcat(individual_nl_post_traces{:}); ...
              vertcat(individual_l_pre_traces{:}); vertcat(individual_l_post_traces{:}); ...
              non_learner_pre_avg; non_learner_post_avg; learner_pre_avg; learner_post_avg];
y_min = min(all_traces(:), [], 'omitnan');
y_max = max(all_traces(:), [], 'omitnan');
y_limits = [y_min, y_max];

% Create figure with 2 rows, 2 columns
figure('Color','w','Position',[100 100 1200 900]);

% ========== ROW 1: NON-LEARNERS ==========
% Subplot 1: Non-Learners Pre-learning
subplot(2,2,1); hold on;
title('Non-Learners - Pre-Learning');

% Plot individual animals (thin lines)
for i = 1:length(non_learner_animal_ids)
    if ~all(isnan(individual_nl_pre_traces{i}))
        plot(time_to_plot, individual_nl_pre_traces{i}, '-', ...
            'Color', cNonLearnerLight, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Plot group average (thick line)
plot(time_to_plot, non_learner_pre_avg, '-', ...
    'Color', cNonLearner, 'LineWidth', 3, 'DisplayName', 'Group Average');

% Mark peak
plot(nl_pre_peak_time, nl_pre_peak_val, 'v', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Peak @ %.2fs', nl_pre_peak_time));
xline(nl_pre_peak_time, ':', 'Color', 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Mark reward availability time
% xline(0.0, '--', 'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'DisplayName', 'Reward Available');

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Left mPFC Activity (a.u.)');
ylim(y_limits);
legend('Location', 'best');
grid on;
hold off;

% Subplot 2: Non-Learners Post-learning
subplot(2,2,2); hold on;
title('Non-Learners - Post-Learning');

% Plot individual animals (thin lines)
for i = 1:length(non_learner_animal_ids)
    if ~all(isnan(individual_nl_post_traces{i}))
        plot(time_to_plot, individual_nl_post_traces{i}, '-', ...
            'Color', cNonLearnerLight, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Plot group average (thick line)
plot(time_to_plot, non_learner_post_avg, '-', ...
    'Color', cNonLearner, 'LineWidth', 3, 'DisplayName', 'Group Average');

% Mark peak
plot(nl_post_peak_time, nl_post_peak_val, 'v', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Peak @ %.2fs', nl_post_peak_time));
xline(nl_post_peak_time, ':', 'Color', 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
% 
% % Mark reward availability time
% xline(0.0, '--', 'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'DisplayName', 'Reward Available');

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Left mPFC Activity (a.u.)');
ylim(y_limits);
legend('Location', 'best');
grid on;
hold off;

% ========== ROW 2: LEARNERS ==========
% Subplot 3: Learners Pre-learning
subplot(2,2,3); hold on;
title('Learners - Pre-Learning');

% Plot individual animals (thin lines)
for i = 1:length(learner_animal_ids)
    if ~all(isnan(individual_l_pre_traces{i}))
        plot(time_to_plot, individual_l_pre_traces{i}, '-', ...
            'Color', cLearnerLight, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Plot group average (thick line)
plot(time_to_plot, learner_pre_avg, '-', ...
    'Color', cLearner, 'LineWidth', 3, 'DisplayName', 'Group Average');

% Mark peak
plot(l_pre_peak_time, l_pre_peak_val, 'v', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Peak @ %.2fs', l_pre_peak_time));
xline(l_pre_peak_time, ':', 'Color', 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% % Mark reward availability time
% xline(0.0, '--', 'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'DisplayName', 'Reward Available');

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Left mPFC Activity (a.u.)');
ylim(y_limits);
legend('Location', 'best');
grid on;
hold off;

% Subplot 4: Learners Post-learning
subplot(2,2,4); hold on;
title('Learners - Post-Learning');

% Plot individual animals (thin lines)
for i = 1:length(learner_animal_ids)
    if ~all(isnan(individual_l_post_traces{i}))
        plot(time_to_plot, individual_l_post_traces{i}, '-', ...
            'Color', cLearnerLight, 'LineWidth', 0.8, 'HandleVisibility', 'off');
    end
end

% Plot group average (thick line)
plot(time_to_plot, learner_post_avg, '-', ...
    'Color', cLearner, 'LineWidth', 3, 'DisplayName', 'Group Average');

% Mark peak
plot(l_post_peak_time, l_post_peak_val, 'v', 'MarkerSize', 12, 'MarkerFaceColor', 'r', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'DisplayName', sprintf('Peak @ %.2fs', l_post_peak_time));
xline(l_post_peak_time, ':', 'Color', 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% % Mark reward availability time
% xline(0.0, '--', 'Color', [0.2 0.4 0.8], 'LineWidth', 2, 'DisplayName', 'Reward Available');

xline(0, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('mPFC Activity (a.u.)');
ylim(y_limits);
legend('Location', 'best');
grid on;
hold off;

sgtitle(sprintf('Right mPFC Activity (Kernel) - \nNon-Learners (n=%d): Peak shift %.0fms | Learners (n=%d): Peak shift %.0fms', ...
    length(non_learner_animal_ids), -nl_peak_shift*1000, length(learner_animal_ids), -l_peak_shift*1000), ...
    'FontSize', 14, 'FontWeight', 'bold');

% % Print animal IDs for reference
% fprintf('\nNon-learner animals plotted (n=%d):\n', length(non_learner_animal_ids));
% for i = 1:length(non_learner_animal_ids)
%     fprintf('  Animal %d\n', non_learner_animal_ids(i));
% end
% 
% fprintf('\nLearner animals plotted (n=%d):\n', length(learner_animal_ids));
% for i = 1:length(learner_animal_ids)
%     fprintf('  Animal %d\n', learner_animal_ids(i));
% end

%% Plot ROI traces for mPFC and visual cortex - learners vs non-learners overlaid

% The script includes different masks for indexing different variables -
% given that some variables only appear in some of the protocols


% Create a logical masks for protocols 
wf1_idx = ((workflow_cat==1));  % protocol 1
wf1_2_idx= find(workflow_cat ~=3); % protocol 1&2

% Create indices that only apply to days with start_move data
stage_1_workflow = workflow_cat(wf1_idx);
stage_1_learning = learning_index_animal(wf1_idx);
stage_1_group = is_group_animal(wf1_idx);
stage_1_widefield_animal= widefield_animal_idx(wf1_idx);

stages_1_2_workflow= workflow_cat(wf1_2_idx);
stages_1_2_learning= learning_index_animal(wf1_2_idx);
stages_1_2_group= is_group_animal(wf1_2_idx);
stages_1_2_widefield_animal= widefield_animal_idx(wf1_2_idx);


% % If for just stage one (rewarded_stim_start_move)
% learner_pre_idx = (stage_1_workflow==1) & (stage_1_learning==0) & (stage_1_group==0);
% learner_post_idx = (stage_1_workflow==1) & (stage_1_learning==1) & (stage_1_group==0);
% 
% % Non-learners (is_group_animal==1)
% non_learner_pre_idx = (stage_1_workflow==1) & (stage_1_learning==0) & (stage_1_group==1);
% non_learner_post_idx = (stage_1_workflow==1) & (stage_1_learning==1) & (stage_1_group==1);

% If for stages 1&2 if (rewarded_stim_final_position)
% learner_pre_idx = (stages_1_2_workflow==1) & (stages_1_2_learning==0) & (stages_1_2_group==0);
% learner_post_idx = (stages_1_2_workflow==1) & (stages_1_2_learning==1) & (stages_1_2_group==0);
% 
% % Non-learners (is_group_animal==1)
% non_learner_pre_idx = (stages_1_2_workflow==1) & (stages_1_2_learning==0) & (stages_1_2_group==1);
% non_learner_post_idx = (stages_1_2_workflow==1) & (stages_1_2_learning==1) & (stages_1_2_group==1);


% deafult: for all stages together
learner_pre_idx= (workflow_cat==1) & (learning_index_animal==0);
learner_post_idx=(workflow_cat==1) & (learning_index_animal==1); 

non_learner_pre_idx= (workflow_cat==1) & (learning_index_animal==0) & (is_group_animal==1);
non_learner_post_idx= (workflow_cat==1) & (learning_index_animal==1) & (is_group_animal==1);


% Count unique animals per group
unique_learner_animals = unique(animal_ids_all_days_stacked(~is_group_animal));
unique_non_learner_animals = unique(animal_ids_all_days_stacked(is_group_animal));

n_learner = length(unique_learner_animals);
n_non_learner = length(unique_non_learner_animals);

% define the variables to plot
data_to_plot= rewarded_stim_v_stacked_data;
ROI_to_plot= left_mPFC_ROI_mask;
time_to_plot= added_time;
curr_components= n_components;

% Learners (is_group_animal==0)
learner_pre = nanmean(data_to_plot(:,:,learner_pre_idx), 3);
learner_pre_std = nanstd(data_to_plot(:,:,learner_pre_idx), 0, 3);
learner_post = nanmean(data_to_plot(:,:,learner_post_idx), 3);
learner_post_std = nanstd(data_to_plot(:,:,learner_post_idx), 0, 3);

% Non-learners (is_group_animal==1)
non_learner_pre = nanmean(data_to_plot(:,:,non_learner_pre_idx), 3);
non_learner_pre_std = nanstd(data_to_plot(:,:,non_learner_pre_idx), 0, 3);
non_learner_post = nanmean(data_to_plot(:,:,non_learner_post_idx), 3);
non_learner_post_std = nanstd(data_to_plot(:,:,non_learner_post_idx), 0, 3);

% Extract left mPFC ROI traces for learners
left_mPFC_learner_pre = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_pre,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_learner_pre_std = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_pre_std,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_learner_pre_sem = left_mPFC_learner_pre_std / sqrt(n_learner);

left_mPFC_learner_post = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_post,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_learner_post_std = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), learner_post_std,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_learner_post_sem = left_mPFC_learner_post_std / sqrt(n_learner);

% Extract left mPFC ROI traces for non-learners
left_mPFC_non_learner_pre = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_pre,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_non_learner_pre_std = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_pre_std,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_non_learner_pre_sem = left_mPFC_non_learner_pre_std / sqrt(n_non_learner);

left_mPFC_non_learner_post = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_post,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_non_learner_post_std = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), non_learner_post_std,[],[], ROI_to_plot), ...
    [3,2,1] );
left_mPFC_non_learner_post_sem = left_mPFC_non_learner_post_std / sqrt(n_non_learner);

% Create figure with overlaid traces
figure('Color','w','Position',[100 100 1000 600]);
hold on;

% Plot SEM shading for Learner Pre (dashed blue)
fill([time_to_plot, fliplr(time_to_plot)], ...
     [left_mPFC_learner_pre + left_mPFC_learner_pre_sem, fliplr(left_mPFC_learner_pre - left_mPFC_learner_pre_sem)], ...
     cLearner, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Learner pre (dashed blue)
plot(time_to_plot, left_mPFC_learner_pre, '--', ...
    'Color', cLearner, 'LineWidth', 2.5, 'DisplayName', 'mPFC+ Pre');

% Plot SEM shading for Learner Post (solid blue)
fill([time_to_plot, fliplr(time_to_plot)], ...
     [left_mPFC_learner_post + left_mPFC_learner_post_sem, fliplr(left_mPFC_learner_post - left_mPFC_learner_post_sem)], ...
     cLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Learner post (solid blue)
plot(time_to_plot, left_mPFC_learner_post, '-', ...
    'Color', cLearner, 'LineWidth', 2.5, 'DisplayName', 'mPFC+ Post');

% Plot SEM shading for Non-learner Pre (dashed orange)
fill([time_to_plot, fliplr(time_to_plot)], ...
     [left_mPFC_non_learner_pre + left_mPFC_non_learner_pre_sem, fliplr(left_mPFC_non_learner_pre - left_mPFC_non_learner_pre_sem)], ...
     cNonLearner, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Non-learner pre (dashed orange)
plot(time_to_plot, left_mPFC_non_learner_pre, '--', ...
    'Color', cNonLearner, 'LineWidth', 2.5, 'DisplayName', 'mPFC- Pre');

% Plot SEM shading for Non-learner Post (solid orange)
fill([time_to_plot, fliplr(time_to_plot)], ...
     [left_mPFC_non_learner_post + left_mPFC_non_learner_post_sem, fliplr(left_mPFC_non_learner_post - left_mPFC_non_learner_post_sem)], ...
     cNonLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Non-learner post (solid orange)
plot(time_to_plot, left_mPFC_non_learner_post, '-', ...
    'Color', cNonLearner, 'LineWidth', 2.5, 'DisplayName', 'mPFC- Post');

xline(0, 'k--', 'LineWidth', 1.5,'HandleVisibility', 'off');
xlabel('Time (s)');
ylabel('Activity (a.u.)');
title('Left mPFC- Aligned to CS+ Movement');
legend('Location', 'best');
hold off;

%% Creates snippets of CCF aligned to an event (currently to CS+ move)

% pre post of a selected workflow
pre_learners = nanmean(rewarded_stim_start_move_v_stacked(:,:,learner_pre_idx),3);
post_learners = nanmean(rewarded_stim_start_move_v_stacked(:,:,learner_post_idx),3);

pre_non_learners = nanmean(rewarded_stim_start_move_v_stacked(:,:,non_learner_pre_idx),3);
post_non_learners = nanmean(rewarded_stim_start_move_v_stacked(:,:,non_learner_post_idx),3);

% ===== Create CCF Activity Snapshots Across Time =====

% Parameters
time_start = -0.5;  % Start time (s)
time_end = 0.5;     % End time (s)
time_step = 0.1;    % Time step (s)

% Generate timepoints
timepoints = time_start:time_step:time_end;
n_timepoints = length(timepoints);

% Find indices in added_time for each timepoint
time_indices = nan(1, n_timepoints);
for t = 1:n_timepoints
    [~, time_indices(t)] = min(abs(added_time - timepoints(t)));
end

% 
pre_learners_px = plab.wf.svd2px(wf_U(:,:,1:n_components), pre_learners);
post_learners_px = plab.wf.svd2px(wf_U(:,:,1:n_components), post_learners);
pre_non_learners_px = plab.wf.svd2px(wf_U(:,:,1:n_components), pre_non_learners);
post_non_learners_px = plab.wf.svd2px(wf_U(:,:,1:n_components), post_non_learners);

% Prepare data
group_data = {pre_learners, post_learners, pre_non_learners, post_non_learners};
group_names = {'Pre-Learning mPFC(+)', 'Post-Learning mPFC(+)', ...
               'Pre-Learning mPFC(-)', 'Post-Learning mPFC(-)'};
group_colors = {cLearner, cLearner, cNonLearner, cNonLearner};

% Convert to pixel space
group_data_px = cell(1, size(group_data,2));
for g = 1:size(group_data,2)
    group_data_px{g} = plab.wf.svd2px(wf_U(:,:,1:n_components), group_data{g});
end

% Get image dimensions
[img_height, img_width, n_frames] = size(pre_learners_px);

% Determine color limits (use max across all data for consistency)
all_data = cat(3, pre_learners_px(:), post_learners_px(:), ...
                  pre_non_learners_px(:), post_non_learners_px(:));
clim_max = max(abs(all_data(:)));
clim_range = 0.6*[-clim_max, clim_max];

% ===== Alternative: Separate Figure for Each Group (Cleaner) =====

for g = 1:4
    figure('Position', [50 + (g-1)*50, 50 + (g-1)*50, 1600, 800], 'Color', 'w');
    
    % Calculate optimal rows/columns (try to make square-ish)
    n_cols = ceil(sqrt(n_timepoints));
    n_rows = ceil(n_timepoints / n_cols);
    
    t_layout = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for t = 1:n_timepoints
        nexttile;
        
        frame = group_data_px{g}(:, :, time_indices(t));
        imagesc(frame);
        axis image off;
        clim(clim_range);
        colormap(ap.colormap('PWG'));
        ap.wf_draw('ccf');

        
        title(sprintf('t = %.1fs', timepoints(t)), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        if abs(timepoints(t)) < 0.01
            hold on;
            rectangle('Position', [0.5, 0.5, img_width, img_height], ...
                'EdgeColor', 'r', 'LineWidth', 4);
        end
    end
    
    cb = colorbar;
    cb.Layout.Tile = 'east';
    cb.Label.String = 'ΔF/F';
    cb.Label.FontSize = 14;
    cb.Label.FontWeight = 'bold';
    
    title(t_layout, group_names{g}, 'FontSize', 18, 'FontWeight', 'bold', ...
        'Color', group_colors{g});
end

%% Creates and saves a video/gif of a CCF video aligned to an event

% ===== Create GIF and Video of CCF Activity =====

% Parameters
output_dir = 'CCF_videos';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Video parameters
frame_rate = 10; % Frames per second (adjust for speed)
gif_delay = 0.1; % Delay between frames in seconds for GIF

% Prepare data
group_data = {pre_learners, post_learners, pre_non_learners, post_non_learners};
group_names = {'Pre-Learning mPFC(+)', 'Post-Learning mPFC(+)', ...
               'Pre-Learning mPFC(-)', 'Post-Learning mPFC(-)'};
group_colors = {cLearner, cLearner, cNonLearner, cNonLearner};

% Convert to pixel space
group_data_px = cell(1, size(group_data,2));
for g = 1:size(group_data,2)
    group_data_px{g} = plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), group_data{g});
end

% Determine color limits
all_data = cat(3, group_data_px{:});
clim_max = max(abs(all_data(:)));
clim_range = [-clim_max, clim_max];


% ===== Create Video and GIF for Each Group =====

for g = 1:size(group_data,2)
    
    fprintf('\nProcessing group %d: %s\n', g, group_names{g});
    
    % Get data
    data_px = group_data_px{g};
    [img_height, img_width, n_frames] = size(data_px);
    
    % Output filenames
     video_filename = fullfile(output_dir, sprintf('group%d_%s_Kernels.mp4', g, ...
        strrep(group_names{g}, ' ', '_')));
    
    % Create Video Writer
    video_writer = VideoWriter(video_filename, 'MPEG-4');
    video_writer.FrameRate = frame_rate;
    video_writer.Quality = 95;
    open(video_writer);
    
    % Create figure for rendering
    fig = figure('Position', [100 100 900 700], 'Color', 'w', 'Visible', 'off');
    
    for frame_idx = n_frames:-1:1 % goes backwards in case of Kernels
        
        % Clear figure
        clf(fig);
        
        % Display frame
        imagesc(data_px(:, :, frame_idx));
        axis image off;
        clim(clim_range);
        colormap(ap.colormap('PWG'));
        ap.wf_draw('ccf')
        
        % Add colorbar
        cb = colorbar;
        cb.Label.String = 'ΔF/F';
        cb.Label.FontSize = 14;
        cb.Label.FontWeight = 'bold';
        
        % Add time indicator
        current_time = added_time_Kernel(frame_idx); 
        title(sprintf('%s\nTime: %.2f s', group_names{g}, current_time), ...
            'FontSize', 16, 'FontWeight', 'bold', 'Color', group_colors{g});
        
        % Highlight CS+ onset (t = 0)
        if abs(current_time) < 0.05
            hold on;
            rectangle('Position', [0.5, 0.5, img_width, img_height], ...
                'EdgeColor', 'r', 'LineWidth', 5);
            text(img_width/2, 10, 'CS+ ONSET', ...
                'FontSize', 18, 'FontWeight', 'bold', 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.7]);
        end
        
        % Add time progress bar at bottom
        hold on;
        progress = (n_frames - frame_idx + 1) / n_frames; 
        bar_width = img_width * 0.8;
        bar_x = img_width * 0.1;
        bar_y = img_height - 10;
        rectangle('Position', [bar_x, bar_y, bar_width, 5], ...
            'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        rectangle('Position', [bar_x, bar_y, bar_width * progress, 5], ...
            'FaceColor', group_colors{g}, 'EdgeColor', 'none');
        
        drawnow;
        
        % Capture frame for video
        frame_img = getframe(fig);
        writeVideo(video_writer, frame_img);
        
        % % Progress indicator
        % if mod(frame_idx, 10) == 0
        %     fprintf('  Frame %d/%d\n', frame_idx, n_frames);
        % end
    end
    
    % Close video writer
    close(video_writer);
    close(fig);
    
    fprintf('✓ Saved video: %s\n', video_filename);
end



fprintf('✓ Saved comparison video: %s\n', comparison_filename);



%% Plot mean/max mPFC response aligned to learning day seperated by learners vs non-learners

% define the variables to plot
data_to_plot = rewarded_stim_start_move_kernel;
ROI_to_plot = right_mPFC_ROI_mask;
time_to_plot = added_time_Kernel;
curr_components = kernel_n_components;
protocol_idx = 1; % workflow to analyze

% Define event and baseline time windows
event_times = [0, 0.3];
baseline_time_windows = [0, 0];

% Create logical masks for both
event_time_window_mask = (time_to_plot > event_times(1) & time_to_plot < event_times(2));
baseline_mask = (time_to_plot >= baseline_time_windows(1) & time_to_plot <= baseline_time_windows(2));


left_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,event_time_window_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % same dims

left_mPFC_baseline = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,baseline_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % Ndays_total × T_baseline

left_mPFC_ROI_trace = left_mPFC_ROI_trace - mean(left_mPFC_baseline, 2); %substract baseline

% unique animals
animals = unique(stage_1_widefield_animal);

% Initialize to store all relative days
rel_day_all = nan(size(stage_1_widefield_animal));

% Find the relative day for each individual animal
for curr_animal = 1:length(animals)
    sel = (stage_1_widefield_animal == curr_animal) & (stage_1_workflow == protocol_idx);
    days_idx = find(sel);
    
    if isempty(days_idx)
        continue;
    end
    
    ld = find(stage_1_learning(days_idx) == 1, 1, 'first');
    
    if isempty(ld)
        continue;
    end
    
    n = numel(days_idx);
    rel_day = (1:n) - ld; % 0 = day of learning
    rel_day_all(sel) = rel_day;
end


% Find max activity for each day
% left_mPFC_ROI_trace_max = max(left_mPFC_ROI_trace, [], 2); % n×1
left_mPFC_ROI_trace_max = mean(left_mPFC_ROI_trace, 2); % n×1

% Calculate group averages
[group_avg, groups] = ap.groupfun(@mean, left_mPFC_ROI_trace_max, [rel_day_all, stage_1_group]);

% Calculate SEM for each group
sem_func = @(x) std(x, 'omitnan') / sqrt(sum(~isnan(x)));
[group_sem, groups_sem] = ap.groupfun(sem_func, left_mPFC_ROI_trace_max, [rel_day_all, stage_1_group]);

% To count animals per relative day per group
[group_counts, groups_count] = ap.groupfun(@length, ones(size(stage_1_widefield_animal)), [rel_day_all, stage_1_group]);

% Filter to keep only days with >2 animals
min_animals = 2;

% Separate learners and non-learners
learner_mask = groups(:,2) == 0;
non_learner_mask = groups(:,2) == 1;

% Filter learners
learner_valid = group_counts(learner_mask) > min_animals;
learner_groups_filtered = groups(learner_mask, :);
learner_groups_filtered = learner_groups_filtered(learner_valid, :);
learner_avg_filtered = group_avg(learner_mask);
learner_avg_filtered = learner_avg_filtered(learner_valid);
learner_sem_filtered = group_sem(learner_mask);
learner_sem_filtered = learner_sem_filtered(learner_valid);

% Sort learners by relative day
[learner_groups_filtered, learner_sort_idx] = sortrows(learner_groups_filtered, 1);
learner_avg_filtered = learner_avg_filtered(learner_sort_idx);
learner_sem_filtered = learner_sem_filtered(learner_sort_idx);

% Filter non-learners
non_learner_valid = group_counts(non_learner_mask) > min_animals;
non_learner_groups_filtered = groups(non_learner_mask, :);
non_learner_groups_filtered = non_learner_groups_filtered(non_learner_valid, :);
non_learner_avg_filtered = group_avg(non_learner_mask);
non_learner_avg_filtered = non_learner_avg_filtered(non_learner_valid);
non_learner_sem_filtered = group_sem(non_learner_mask);
non_learner_sem_filtered = non_learner_sem_filtered(non_learner_valid);

% Sort non-learners by relative day
[non_learner_groups_filtered, non_learner_sort_idx] = sortrows(non_learner_groups_filtered, 1);
non_learner_avg_filtered = non_learner_avg_filtered(non_learner_sort_idx);
non_learner_sem_filtered = non_learner_sem_filtered(non_learner_sort_idx);

% Plot filtered data with SEM
figure; hold on;

% Plot learner SEM shading
learner_days = learner_groups_filtered(:, 1);
fill([learner_days; flipud(learner_days)], ...
     [learner_avg_filtered + learner_sem_filtered; flipud(learner_avg_filtered - learner_sem_filtered)], ...
     cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot learner mean line
plot(learner_days, learner_avg_filtered, '-o', ...
    'LineWidth', 2, 'DisplayName', 'Learners', 'Color', cLearner, 'MarkerFaceColor', cLearner);

% Plot non-learner SEM shading
non_learner_days = non_learner_groups_filtered(:, 1);
fill([non_learner_days; flipud(non_learner_days)], ...
     [non_learner_avg_filtered + non_learner_sem_filtered; flipud(non_learner_avg_filtered - non_learner_sem_filtered)], ...
     cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot non-learner mean line
plot(non_learner_days, non_learner_avg_filtered, '-o', ...
    'LineWidth', 2, 'DisplayName', 'Non-Learners', 'Color', cNonLearner, 'MarkerFaceColor', cNonLearner);

xline(0, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Learning Day');
xlabel('Days Relative to Learning');
ylabel('Max mPFC Activity (Baseline-subtracted)');
title(sprintf('Right Stimulus Left mPFC Activity (Kernels) (days with >%d animals per group)', min_animals));
legend('Location', 'best');
grid on;
hold off;

%% Summary plot of mPFC aligned to trasnition day across all protocols

% ===== mPFC Activity Across Protocol Transitions =====

% Parameters
ROI_to_extract = right_mPFC_ROI_mask;
curr_components = kernel_n_components;

% Response window
time_window_start = 0;
time_window_end = 0.35;

% Minimum animals per day threshold
min_animals_per_day = 2;

% Define protocol keywords
protocol_1_keyword = 'right_move';  % Classical conditioning
protocol_2_keyword = 'static';      % Static reward
protocol_3_keyword = 'stim_wheel*';  % Wheel task


% ===== Step 1: Extract mPFC response for each animal-day =====

n_animals = numel(behaviour_data);
animal_mPFC_data = struct();

for a = 1:n_animals
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    
    if isempty(rd)
        fprintf('Animal %s: No recording days, skipping\n', aid);
        continue;
    end
    
    % Find matching widefield indices for this animal
    matching_indices = find(widefield_animal_idx == a);
    n_days = length(matching_indices);
    
    if isempty(matching_indices)
        fprintf('Animal %s: No widefield data, skipping\n', aid);
        continue;
    end
    
    % Initialize storage for this animal
    animal_mPFC_data(a).animal_id = aid;
    animal_mPFC_data(a).n_days = n_days;
    animal_mPFC_data(a).mPFC_response = nan(1, n_days);
    animal_mPFC_data(a).protocol = strings(1, n_days);
    animal_mPFC_data(a).day_number = nan(1, n_days);
    
    % Extract mPFC for each day
    for d = 1:n_days
        day_idx = matching_indices(d);
        
        % Get workflow/protocol
        if d <= length(rd) && isfield(rd(d), 'workflow')
            workflow = rd(d).workflow;
            
            % Classify protocol
            if contains(workflow, protocol_1_keyword)
                protocol_name = 'Protocol_1';
            elseif contains(workflow, protocol_2_keyword)
                protocol_name = 'Protocol_2';
            elseif contains(workflow, protocol_3_keyword)
                protocol_name = 'Protocol_3';
            else
                protocol_name = 'Unknown';
            end
        else
            protocol_name = 'Unknown';
        end
        
        animal_mPFC_data(a).protocol(d) = protocol_name;
        animal_mPFC_data(a).day_number(d) = d;
        
        % Extract mean mPFC activity for the session
        V_data = rewarded_stim_kernel(:,:,day_idx);
        
        if isempty(V_data)
            continue;
        end

        % Extract ROI trace
        roi_trace = ap.wf_roi(wf_U(:,:,1:curr_components), V_data, [], [], ROI_to_extract);

        % Calculate response
        response_window_idx = added_time_Kernel >= time_window_start & added_time_Kernel <= time_window_end;
        baseline_window_idx = added_time_Kernel >= -0.1 & added_time_Kernel < 0;

        baseline = mean(roi_trace(baseline_window_idx), 'omitnan');
        response = mean(roi_trace(response_window_idx),2);

        % Save baseline sub for session
        animal_mPFC_data(a).mPFC_response(d) = response- baseline;
    end
    
    fprintf('Animal %s: Extracted %d days\n', aid, n_days);
end


% ===== Create aligned transition data structure (similar to animal_mPFC_data) =====


% Initialize structure array (1 x n_animals)
n_animals = numel(behaviour_data);
aligned_transition_data = struct();

for a = 1:n_animals
 
    
    % Copy basic info
    aligned_transition_data(a).animal_id = animal_mPFC_data(a).animal_id;
    
    % Determine group membership
    if ~ismember(a, group_animals)
        aligned_transition_data(a).group = 'mPFC_plus';
    elseif ismember(a, mPFC_minus_animals)
        aligned_transition_data(a).group = 'mPFC_minus';
    else
        aligned_transition_data(a).group = 'unknown';
    end
    
    protocols = animal_mPFC_data(a).protocol;
    mPFC_response = animal_mPFC_data(a).mPFC_response;
    
    % Find protocol days
    prot1_days = find(protocols == "Protocol_1");
    prot2_days = find(protocols == "Protocol_2");
    prot3_days = find(protocols == "Protocol_3");
    
    % Initialize aligned data arrays
    aligned_transition_data(a).transition_1_to_2.days_relative = [];
    aligned_transition_data(a).transition_1_to_2.mPFC = [];
    aligned_transition_data(a).transition_1_to_2.has_data = false;
    
    aligned_transition_data(a).transition_2_to_3.days_relative = [];
    aligned_transition_data(a).transition_2_to_3.mPFC = [];
    aligned_transition_data(a).transition_2_to_3.has_data = false;
    
    % ===== Transition 1→2 =====
    if ~isempty(prot1_days) && ~isempty(prot2_days)
        if min(prot2_days) == max(prot1_days) + 1
            % Get transition days - last 2 days protocl 1 and +2 days from
            % protocol 2
            transition_days_1 = [prot1_days(max(1, end-2):end), ...
                                prot2_days(1:min(4, length(prot2_days)))];
            
            % Align to transition (day 0 = first day of Protocol 2)
            aligned_days = transition_days_1 - min(prot2_days);
            
            % Store aligned data
            aligned_transition_data(a).transition_1_to_2.days_relative = aligned_days;
            aligned_transition_data(a).transition_1_to_2.mPFC = mPFC_response(transition_days_1);
            aligned_transition_data(a).transition_1_to_2.has_data = true;
            aligned_transition_data(a).transition_1_to_2.n_days = length(aligned_days);
            
            fprintf('Animal %d (%s): Transition 1→2 - %d days\n', ...
                a, aligned_transition_data(a).animal_id, length(aligned_days));
        end
    end
    
    % ===== Transition 2→3 =====
    if ~isempty(prot2_days) && ~isempty(prot3_days)
        if min(prot3_days) == max(prot2_days) + 1
            % Get transition days
            transition_days_2 = [prot2_days(max(1, end-3):end), ...
                                prot3_days(1:min(3, length(prot3_days)))];
            
            % Align to transition (day 0 = first day of Protocol 3)
            aligned_days = transition_days_2 - min(prot3_days);
            
            % Store aligned data
            aligned_transition_data(a).transition_2_to_3.days_relative = aligned_days;
            aligned_transition_data(a).transition_2_to_3.mPFC = mPFC_response(transition_days_2);
            aligned_transition_data(a).transition_2_to_3.has_data = true;
            aligned_transition_data(a).transition_2_to_3.n_days = length(aligned_days);
            
            fprintf('Animal %d (%s): Transition 2→3 - %d days\n', ...
                a, aligned_transition_data(a).animal_id, length(aligned_days));
        end
    end
end

% ===== Aggregate by group and relative day =====

% Separate by group
mPFC_plus_idx = strcmp({aligned_transition_data.group}, 'mPFC_plus');
mPFC_minus_idx = strcmp({aligned_transition_data.group}, 'mPFC_minus');

aggregated_data = struct();

for g = 1:2
    if g == 1
        group_name = 'mPFC_plus';
        group_idx = mPFC_plus_idx;
    else
        group_name = 'mPFC_minus';
        group_idx = mPFC_minus_idx;
    end
    
    % Get animal indices in this group
    group_animal_indices = find(group_idx);

    % Aggregate Transition 1→2
    all_days = [];
    all_mPFC = [];
    all_animals = [];
    
    % Store individual animal data
    individual_animals_1_to_2 = struct('animal_idx', {}, 'animal_id', {}, ...
                                       'days_relative', {}, 'mPFC', {});
    
    animal_counter = 0;
    for a = group_animal_indices
        if aligned_transition_data(a).transition_1_to_2.has_data
            n_days = length(aligned_transition_data(a).transition_1_to_2.days_relative);

            all_days = [all_days, aligned_transition_data(a).transition_1_to_2.days_relative];
            all_mPFC = [all_mPFC, aligned_transition_data(a).transition_1_to_2.mPFC];
            all_animals = [all_animals, repmat(a, 1, n_days)];
            
            % Store individual animal trajectory
            animal_counter = animal_counter + 1;
            individual_animals_1_to_2(animal_counter).animal_idx = a;
            individual_animals_1_to_2(animal_counter).animal_id = aligned_transition_data(a).animal_id;
            individual_animals_1_to_2(animal_counter).days_relative = aligned_transition_data(a).transition_1_to_2.days_relative;
            individual_animals_1_to_2(animal_counter).mPFC = aligned_transition_data(a).transition_1_to_2.mPFC;
        end
    end
    
    % Get unique days and compute mean/sem
    if ~isempty(all_days)
        unique_days = unique(all_days);
        means = nan(size(unique_days));
        sems = nan(size(unique_days));
        ns = nan(size(unique_days));
        
        for d = 1:length(unique_days)
            day_mask = all_days == unique_days(d);
            day_mPFC = all_mPFC(day_mask);
            day_animals = all_animals(day_mask);
            
            % Count unique animals
            n_unique_animals = length(unique(day_animals));
            
            if n_unique_animals >= min_animals_per_day
                means(d) = mean(day_mPFC, 'omitnan');
                sems(d) = std(day_mPFC, 'omitnan') / sqrt(n_unique_animals);
                ns(d) = n_unique_animals;
            end
        end
        
        % Store aggregated data (remove NaN entries)
        valid_idx = ~isnan(means);
        aggregated_data.(group_name).trans_1_to_2.days = unique_days(valid_idx);
        aggregated_data.(group_name).trans_1_to_2.mean = means(valid_idx);
        aggregated_data.(group_name).trans_1_to_2.sem = sems(valid_idx);
        aggregated_data.(group_name).trans_1_to_2.n = ns(valid_idx);
        aggregated_data.(group_name).trans_1_to_2.individual_animals = individual_animals_1_to_2;
        
        fprintf('%s Transition 1→2: %d days with n≥%d, %d animals\n', ...
            group_name, sum(valid_idx), min_animals_per_day, length(individual_animals_1_to_2));
    end
    
    % Aggregate Transition 2→3
    all_days = [];
    all_mPFC = [];
    all_animals = [];
    
    % Store individual animal data
    individual_animals_2_to_3 = struct('animal_idx', {}, 'animal_id', {}, ...
                                       'days_relative', {}, 'mPFC', {});
    
    animal_counter = 0;
    for a = group_animal_indices
        if aligned_transition_data(a).transition_2_to_3.has_data
            n_days = length(aligned_transition_data(a).transition_2_to_3.days_relative);

            all_days = [all_days, aligned_transition_data(a).transition_2_to_3.days_relative];
            all_mPFC = [all_mPFC, aligned_transition_data(a).transition_2_to_3.mPFC];
            all_animals = [all_animals, repmat(a, 1, n_days)];
            
            % Store individual animal trajectory
            animal_counter = animal_counter + 1;
            individual_animals_2_to_3(animal_counter).animal_idx = a;
            individual_animals_2_to_3(animal_counter).animal_id = aligned_transition_data(a).animal_id;
            individual_animals_2_to_3(animal_counter).days_relative = aligned_transition_data(a).transition_2_to_3.days_relative;
            individual_animals_2_to_3(animal_counter).mPFC = aligned_transition_data(a).transition_2_to_3.mPFC;
        end
    end
    
    if ~isempty(all_days)
        unique_days = unique(all_days);
        means = nan(size(unique_days));
        sems = nan(size(unique_days));
        ns = nan(size(unique_days));
        
        for d = 1:length(unique_days)
            day_mask = all_days == unique_days(d);
            day_mPFC = all_mPFC(day_mask);
            day_animals = all_animals(day_mask);
            
            n_unique_animals = length(unique(day_animals));
            
            if n_unique_animals >= min_animals_per_day
                means(d) = mean(day_mPFC, 'omitnan');
                sems(d) = std(day_mPFC, 'omitnan') / sqrt(n_unique_animals);
                ns(d) = n_unique_animals;
            end
        end
        
        valid_idx = ~isnan(means);
        aggregated_data.(group_name).trans_2_to_3.days = unique_days(valid_idx);
        aggregated_data.(group_name).trans_2_to_3.mean = means(valid_idx);
        aggregated_data.(group_name).trans_2_to_3.sem = sems(valid_idx);
        aggregated_data.(group_name).trans_2_to_3.n = ns(valid_idx);
        aggregated_data.(group_name).trans_2_to_3.individual_animals = individual_animals_2_to_3;
        
        fprintf('%s Transition 2→3: %d days with n≥%d, %d animals\n', ...
            group_name, sum(valid_idx), min_animals_per_day, length(individual_animals_2_to_3));
    end
end

% ===== Plot Summary with Individual Animal Traces =====

figure('Position', [100 100 1400 600], 'Color', 'w');

% Colors
color_mPFC_plus = cLearner;
color_mPFC_minus = cNonLearner;
color_mPFC_plus_light = color_mPFC_plus + (1 - color_mPFC_plus) * 0.6;  % Lighter version
color_mPFC_minus_light = color_mPFC_minus + (1 - color_mPFC_minus) * 0.6;

% Subplot 1: Transition 1→2
subplot(1, 2, 1);
hold on;

% Plot individual animals for mPFC+ group
if isfield(aggregated_data.mPFC_plus, 'trans_1_to_2')
    individual_animals = aggregated_data.mPFC_plus.trans_1_to_2.individual_animals;
    
    for i = 1:length(individual_animals)
        plot(individual_animals(i).days_relative, individual_animals(i).mPFC, ...
            'o-', 'Color', color_mPFC_plus_light, 'LineWidth', 1, ...
            'MarkerSize', 4, 'MarkerFaceColor', color_mPFC_plus_light, ...
            'HandleVisibility', 'off');
    end
end

% Plot individual animals for mPFC- group
if isfield(aggregated_data.mPFC_minus, 'trans_1_to_2')
    individual_animals = aggregated_data.mPFC_minus.trans_1_to_2.individual_animals;
    
    for i = 1:length(individual_animals)
        plot(individual_animals(i).days_relative, individual_animals(i).mPFC, ...
            'o-', 'Color', color_mPFC_minus_light, 'LineWidth', 1, ...
            'MarkerSize', 4, 'MarkerFaceColor', color_mPFC_minus_light, ...
            'HandleVisibility', 'off');
    end
end

% Plot group averages with SEM on top
if isfield(aggregated_data.mPFC_plus, 'trans_1_to_2')
    data = aggregated_data.mPFC_plus.trans_1_to_2;
    
    fill([data.days, fliplr(data.days)], ...
         [data.mean + data.sem, fliplr(data.mean - data.sem)], ...
         color_mPFC_plus, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    
    plot(data.days, data.mean, 'o-', 'Color', color_mPFC_plus, ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color_mPFC_plus, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('mPFC+ (n=%d)', length(aggregated_data.mPFC_plus.trans_1_to_2.individual_animals)));
end

if isfield(aggregated_data.mPFC_minus, 'trans_1_to_2')
    data = aggregated_data.mPFC_minus.trans_1_to_2;
    
    fill([data.days, fliplr(data.days)], ...
         [data.mean + data.sem, fliplr(data.mean - data.sem)], ...
         color_mPFC_minus, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    
    plot(data.days, data.mean, 'o-', 'Color', color_mPFC_minus, ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color_mPFC_minus, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('mPFC- (n=%d)', length(aggregated_data.mPFC_minus.trans_1_to_2.individual_animals)));
end

xline(0, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Protocol transition');
yline(0, 'k:', 'LineWidth', 1);

xlabel('Days Relative to Transition', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('Protocol 1 → Protocol 2 Transition', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

% Subplot 2: Transition 2→3
subplot(1, 2, 2);
hold on;

% Plot individual animals for mPFC+ group
if isfield(aggregated_data.mPFC_plus, 'trans_2_to_3')
    individual_animals = aggregated_data.mPFC_plus.trans_2_to_3.individual_animals;
    
    for i = 1:length(individual_animals)
        plot(individual_animals(i).days_relative, individual_animals(i).mPFC, ...
            'o-', 'Color', color_mPFC_plus_light, 'LineWidth', 1, ...
            'MarkerSize', 4, 'MarkerFaceColor', color_mPFC_plus_light, ...
            'HandleVisibility', 'off');
    end
end

% Plot individual animals for mPFC- group
if isfield(aggregated_data.mPFC_minus, 'trans_2_to_3')
    individual_animals = aggregated_data.mPFC_minus.trans_2_to_3.individual_animals;
    
    for i = 1:length(individual_animals)
        plot(individual_animals(i).days_relative, individual_animals(i).mPFC, ...
            'o-', 'Color', color_mPFC_minus_light, 'LineWidth', 1, ...
            'MarkerSize', 4, 'MarkerFaceColor', color_mPFC_minus_light, ...
            'HandleVisibility', 'off');
    end
end

% Plot group averages with SEM on top
if isfield(aggregated_data.mPFC_plus, 'trans_2_to_3')
    data = aggregated_data.mPFC_plus.trans_2_to_3;
    
    fill([data.days, fliplr(data.days)], ...
         [data.mean + data.sem, fliplr(data.mean - data.sem)], ...
         color_mPFC_plus, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    
    plot(data.days, data.mean, 'o-', 'Color', color_mPFC_plus, ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color_mPFC_plus, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('mPFC+ (n=%d)', length(aggregated_data.mPFC_plus.trans_2_to_3.individual_animals)));
end

if isfield(aggregated_data.mPFC_minus, 'trans_2_to_3')
    data = aggregated_data.mPFC_minus.trans_2_to_3;
    
    fill([data.days, fliplr(data.days)], ...
         [data.mean + data.sem, fliplr(data.mean - data.sem)], ...
         color_mPFC_minus, 'FaceAlpha', 0.3, 'EdgeColor', 'none', ...
         'HandleVisibility', 'off');
    
    plot(data.days, data.mean, 'o-', 'Color', color_mPFC_minus, ...
        'LineWidth', 3, 'MarkerSize', 10, 'MarkerFaceColor', color_mPFC_minus, ...
        'MarkerEdgeColor', 'k', 'LineWidth', 2.5, ...
        'DisplayName', sprintf('mPFC- (n=%d)', length(aggregated_data.mPFC_minus.trans_2_to_3.individual_animals)));
end

xline(0, 'k--', 'LineWidth', 2.5, 'DisplayName', 'Protocol transition');
yline(0, 'k:', 'LineWidth', 1);

xlabel('Days Relative to Transition', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('Protocol 2 → Protocol 3 Transition', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 12, 'LineWidth', 1.5);

sgtitle('Mean Right mPFC Activity (Kernel) Across Protocol Transitions', 'FontSize', 18, 'FontWeight', 'bold');


%% Plot mean/max mPFC response aligned to learning day seperated by learners vs non-learners (includes individual animal trajectories)

% Switch color codes to grey
% cLearner =[0.5 0.5 0.5]; 
% cNonLearner= [0.5 0.5 0.5]; 

% define the variables to plot
data_to_plot = rewarded_stim_v_stacked_data;
ROI_to_plot = left_mPFC_ROI_mask;
time_to_plot = added_time;
curr_components = n_components;
protocol_idx = 1; % workflow to analyze

% Define event and baseline time windows
event_times = [0, 0.35];
baseline_time_windows = [-0.1, 0];

% Create logical masks for both
event_time_window_mask = (time_to_plot > event_times(1) & time_to_plot < event_times(2));
baseline_mask = (time_to_plot >= baseline_time_windows(1) & time_to_plot <= baseline_time_windows(2));

left_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,event_time_window_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % same dims

left_mPFC_baseline = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,baseline_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % Ndays_total × T_baseline

left_mPFC_ROI_trace = left_mPFC_ROI_trace - mean(left_mPFC_baseline, 2); %substract baseline

% unique animals
animals = unique(widefield_animal_idx);

% Initialize to store all relative days
rel_day_all = nan(size(widefield_animal_idx));

% Find the relative day for each individual animal
for curr_animal = 1:length(animals)
    sel = (widefield_animal_idx == curr_animal) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    
    if isempty(days_idx)
        continue;
    end
    
    ld = find(learning_index_animal(days_idx) == 1, 1, 'first');
    
    if isempty(ld)
        continue;
    end
    
    n = numel(days_idx);
    rel_day = (1:n) - ld; % 0 = day of learning
    rel_day_all(sel) = rel_day;
end

% Find mean/max activity for each day
left_mPFC_ROI_trace_max = mean(left_mPFC_ROI_trace,2); % n×1

% Store individual animal trajectories
learner_animal_ids = unique(widefield_animal_idx(is_group_animal == 0));
non_learner_animal_ids = unique(widefield_animal_idx(is_group_animal == 1));

learner_trajectories = cell(length(learner_animal_ids), 1);
non_learner_trajectories = cell(length(non_learner_animal_ids), 1);

% Extract learner trajectories
for i = 1:length(learner_animal_ids)
    animal_id = learner_animal_ids(i);
    sel = (widefield_animal_idx == animal_id) & (workflow_cat == protocol_idx) & ~isnan(rel_day_all);
    
    if any(sel)
        learner_trajectories{i} = struct(...
            'days', rel_day_all(sel), ...
            'values', left_mPFC_ROI_trace_max(sel));
    end
end

% Extract non-learner trajectories
for i = 1:length(non_learner_animal_ids)
    animal_id = non_learner_animal_ids(i);
    sel = (widefield_animal_idx == animal_id) & (workflow_cat == protocol_idx) & ~isnan(rel_day_all);
    
    if any(sel)
        non_learner_trajectories{i} = struct(...
            'days', rel_day_all(sel), ...
            'values', left_mPFC_ROI_trace_max(sel));
    end
end

% Calculate group averages

[group_avg, groups] = ap.groupfun(@mean, left_mPFC_ROI_trace_max, [rel_day_all, is_group_animal]);

% Calculate SEM for each group
sem_func = @(x) std(x, 'omitnan') / sqrt(sum(~isnan(x)));
[group_sem, groups_sem] = ap.groupfun(sem_func, left_mPFC_ROI_trace_max, [rel_day_all, is_group_animal]);

% To count animals per relative day per group
[group_counts, groups_count] = ap.groupfun(@length, ones(size(widefield_animal_idx)), [rel_day_all, is_group_animal]);

% Filter to keep only days with >2 animals
min_animals = 2;

% Separate learners and non-learners
learner_mask = groups(:,2) == 0;
non_learner_mask = groups(:,2) == 1;

% Filter learners
learner_valid = group_counts(learner_mask) > min_animals;
learner_groups_filtered = groups(learner_mask, :);
learner_groups_filtered = learner_groups_filtered(learner_valid, :);
learner_avg_filtered = group_avg(learner_mask);
learner_avg_filtered = learner_avg_filtered(learner_valid);
learner_sem_filtered = group_sem(learner_mask);
learner_sem_filtered = learner_sem_filtered(learner_valid);

% Sort learners by relative day
[learner_groups_filtered, learner_sort_idx] = sortrows(learner_groups_filtered, 1);
learner_avg_filtered = learner_avg_filtered(learner_sort_idx);
learner_sem_filtered = learner_sem_filtered(learner_sort_idx);

% Filter non-learners
non_learner_valid = group_counts(non_learner_mask) > min_animals;
non_learner_groups_filtered = groups(non_learner_mask, :);
non_learner_groups_filtered = non_learner_groups_filtered(non_learner_valid, :);
non_learner_avg_filtered = group_avg(non_learner_mask);
non_learner_avg_filtered = non_learner_avg_filtered(non_learner_valid);
non_learner_sem_filtered = group_sem(non_learner_mask);
non_learner_sem_filtered = non_learner_sem_filtered(non_learner_valid);

% Sort non-learners by relative day
[non_learner_groups_filtered, non_learner_sort_idx] = sortrows(non_learner_groups_filtered, 1);
non_learner_avg_filtered = non_learner_avg_filtered(non_learner_sort_idx);
non_learner_sem_filtered = non_learner_sem_filtered(non_learner_sort_idx);

% Plot filtered data with individual trajectories and SEM
figure('Color','w','Position',[100 100 900 600]); hold on;

% Plot individual LEARNER trajectories (thin lines, semi-transparent)
for i = 1:length(learner_trajectories)
    if ~isempty(learner_trajectories{i})
        plot(learner_trajectories{i}.days, learner_trajectories{i}.values, '-', ...
            'Color', [cLearner 0.6], 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
end

% Plot individual NON-LEARNER trajectories (thin lines, semi-transparent)
for i = 1:length(non_learner_trajectories)
    if ~isempty(non_learner_trajectories{i})
        plot(non_learner_trajectories{i}.days, non_learner_trajectories{i}.values, '-', ...
            'Color', [cNonLearner 0.6], 'LineWidth',1.2, 'HandleVisibility', 'off');
    end
end

%Plot learner SEM shading
learner_days = learner_groups_filtered(:, 1);
fill([learner_days; flipud(learner_days)], ...
     [learner_avg_filtered + learner_sem_filtered; flipud(learner_avg_filtered - learner_sem_filtered)], ...
     cLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot learner mean line (BOLD)
plot(learner_days, learner_avg_filtered, '-o', ...
    'LineWidth', 3, 'DisplayName', sprintf('Learners (n=%d)', length(learner_animal_ids)), ...
    'Color', cLearner, 'MarkerFaceColor', cLearner, 'MarkerSize', 8);

% Plot non-learner SEM shading
non_learner_days = non_learner_groups_filtered(:, 1);
fill([non_learner_days; flipud(non_learner_days)], ...
     [non_learner_avg_filtered + non_learner_sem_filtered; flipud(non_learner_avg_filtered - non_learner_sem_filtered)], ...
     cNonLearner, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Plot non-learner mean line (BOLD)
plot(non_learner_days, non_learner_avg_filtered, '-o', ...
    'LineWidth', 3, 'DisplayName', sprintf('Non-Learners (n=%d)', length(non_learner_animal_ids)), ...
    'Color', cNonLearner, 'MarkerFaceColor', cNonLearner, 'MarkerSize', 8);

xline(0, 'k--', 'LineWidth', 2, 'DisplayName', 'Learning Day');
xlabel('Days Relative to Learning', 'FontSize', 12);
ylabel('Mean mPFC Activity (Baseline-subtracted)', 'FontSize', 12);
title(sprintf('Right mPFC Activity (Kernel) Aligned to Learning Day (days with >%d animals per group)', min_animals), ...
    'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
hold off;

fprintf('Plotted trajectories:\n');
fprintf('  Learners: %d animals\n', length(learner_animal_ids));
fprintf('  Non-learners: %d animals\n', length(non_learner_animal_ids));
fprintf('  Days range: %d to +%d\n', min(learner_days(1), non_learner_days(1)), ...
    max(learner_days(end), non_learner_days(end)));

%% Plot max/mean baseline substracted CCF maps for V

% Choose 'max' or 'mean' to plot
statistic = 'mean';

% Define event and baseline time windows
event_times = [0, 0.25];
baseline_time_windows = [-0.1, 0];

% define data to plot
data_to_plot= rewarded_stim_v_stacked_data;

% Find the corresponding indices in added_time for the time window
time_window_idx = (added_time >= event_times(1) & added_time <= event_times(2));
baseline_mask = (added_time >= baseline_time_windows(1) & added_time <= baseline_time_windows(2));

% pre post of right move workflow split by learning days and learners vs
% non-learners
mean_pre_learning_learners = nanmean(data_to_plot(:,:,learner_pre_idx ),3);
mean_post_learning_learners = nanmean(data_to_plot(:,:,learner_post_idx),3);

mean_pre_learning_non_learners = nanmean(data_to_plot(:,:,non_learner_pre_idx),3);
mean_post_learning_non_learners = nanmean(data_to_plot(:,:,non_learner_post_idx),3);
 
% substract baseline
baseline_substracted_pre_learning_LR= mean_pre_learning_learners- nanmean(mean_pre_learning_learners(:,baseline_mask),2);
baseline_substracted_post_learning_LR= mean_post_learning_learners- nanmean(mean_post_learning_learners(:,baseline_mask),2);

baseline_substracted_pre_learning_NLR= mean_pre_learning_non_learners- nanmean(mean_pre_learning_non_learners(:,baseline_mask),2);
baseline_substracted_post_learning_NLR= mean_post_learning_non_learners- nanmean(mean_post_learning_non_learners(:,baseline_mask),2);

% get the specific time window
timewindow_mean_pre_learning_learners= plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_substracted_pre_learning_LR(:,time_window_idx));
timewindow_mean_post_learning_learners= plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_substracted_post_learning_LR(:,time_window_idx));

timewindow_mean_pre_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_substracted_pre_learning_NLR(:,time_window_idx));
timewindow_mean_post_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_substracted_post_learning_NLR(:,time_window_idx));



% Compute the image to show
switch statistic
    case 'max'
        img_pre_learning_learners = max(timewindow_mean_pre_learning_learners, [], 3);
        img_post_learning_learners= max(timewindow_mean_post_learning_learners, [], 3);
        img_pre_learning_non_learners   = max(timewindow_mean_pre_learning_non_learners, [], 3);
        img_post_learning_non_learners   = max(timewindow_mean_post_learning_non_learners, [], 3);

    case 'mean'
        img_pre_learning_learners = mean(timewindow_mean_pre_learning_learners, 3);
        img_post_learning_learners= mean(timewindow_mean_post_learning_learners, 3);
        img_pre_learning_non_learners = mean(timewindow_mean_pre_learning_non_learners, 3);
        img_post_learning_non_learners = mean(timewindow_mean_post_learning_non_learners,3);

    otherwise
        error('statistic must be ''max'' or ''mean''');
end

% Compute a common color‐limit
all_vals = [timewindow_mean_pre_learning_learners(:); timewindow_mean_pre_learning_non_learners(:)];
switch statistic
    case 'max'
        ref_vals = max(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),[],2);
    case 'mean'
        ref_vals = mean(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),2);
end
clim_val = 0.9 * max(abs(ref_vals));

% Plot learners
figure;
imagesc(img_pre_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s V activity — mPFC+', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s V activity — mPFC+', statistic));
ap.wf_draw('ccf');

% Plot non learners
figure;
imagesc(img_pre_learning_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s V activity — mPFC-', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s V activity — mPFC-', statistic));
ap.wf_draw('ccf');


%% Plot max/mean baseline substracted CCF maps for Kernels

% Choose 'max' or 'mean' to plot
statistic = 'max';

% define data to plot
data_to_plot= rewarded_stim_kernel;

added_time_Kernel=  fliplr((-10:52)/30);

% Find the corresponding indices in added_time for the time window
time_window_idx = (added_time_Kernel >= event_times(1) & added_time_Kernel <= event_times(2));
baseline_mask = (added_time_Kernel >= baseline_time_windows(1) & added_time_Kernel <= baseline_time_windows(2));


% pre post of right move workflow split by learning days and learners vs
% non-learners

mean_pre_learning_learners = nanmean(data_to_plot(:,:,learner_pre_idx),3);
mean_post_learning_learners = nanmean(data_to_plot(:,:,learner_post_idx),3);

mean_pre_learning_non_learners = nanmean(data_to_plot(:,:,non_learner_pre_idx),3);
mean_post_learning_non_learners = nanmean(data_to_plot(:,:,non_learner_post_idx),3);

% substract baseline
baseline_substracted_pre_learning_LR= mean_pre_learning_learners- nanmean(mean_pre_learning_learners(:,baseline_mask),2);
baseline_substracted_post_learning_LR= mean_post_learning_learners- nanmean(mean_post_learning_learners(:,baseline_mask),2);

baseline_substracted_pre_learning_NLR= mean_pre_learning_non_learners- nanmean(mean_pre_learning_non_learners(:,baseline_mask),2);
baseline_substracted_post_learning_NLR= mean_post_learning_non_learners- nanmean(mean_post_learning_non_learners(:,baseline_mask),2);

% get the specific time window
timewindow_mean_pre_learning_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), baseline_substracted_pre_learning_LR(:,time_window_idx));
timewindow_mean_post_learning_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), baseline_substracted_post_learning_LR(:,time_window_idx));

timewindow_mean_pre_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), baseline_substracted_pre_learning_NLR(:,time_window_idx));
timewindow_mean_post_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), baseline_substracted_post_learning_NLR(:,time_window_idx));

% Compute the image to show
switch statistic
    case 'max'
        img_pre_learning_learners = max(timewindow_mean_pre_learning_learners, [], 3);
        img_post_learning_learners= max(timewindow_mean_post_learning_learners, [], 3);
        img_pre_learning_non_learners   = max(timewindow_mean_pre_learning_non_learners, [], 3);
        img_post_learning_non_learners   = max(timewindow_mean_post_learning_non_learners, [], 3);

    case 'mean'
        img_pre_learning_learners = mean(timewindow_mean_pre_learning_learners, 3);
        img_post_learning_learners= mean(timewindow_mean_post_learning_learners, 3);
        img_pre_learning_non_learners = mean(timewindow_mean_pre_learning_non_learners, 3);
        img_post_learning_non_learners = mean(timewindow_mean_post_learning_non_learners,3);

    otherwise
        error('statistic must be ''max'' or ''mean''');
end

% Compute a common color‐limit
all_vals = [timewindow_mean_pre_learning_learners(:); timewindow_mean_pre_learning_non_learners(:)];
switch statistic
    case 'max'
        ref_vals = max(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),[],2);
    case 'mean'
        ref_vals = mean(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),2);
end
clim_val = 0.8 * max(abs(ref_vals));
% 
% clim_holder= clim;
% clim_val=clim_holder(2);


% Plot learners
figure;
imagesc(img_pre_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s Kernels — mPFC+', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s Kernels — mPFC+', statistic));
ap.wf_draw('ccf');

% Plot non learners
figure;
imagesc(img_pre_learning_non_learners);
axis image off;
colormap(ap.colormap('GWP'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s Kernels— mPFC-', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_non_learners);
axis image off;
colormap(ap.colormap('GWP'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s Kernels — mPFC-', statistic));
ap.wf_draw('ccf');

%% Plot max/mean baseline substracted CCF maps for Kernels comparing first and last learning days

% Choose 'max' or 'mean' to plot
statistic = 'max';

% Define data to plot
data_to_plot = rewarded_stim_v_stacked_data;

added_time_Kernel = fliplr((-10:52)/30);

% Find the corresponding indices in added_time for the time window
time_window_idx = (added_time >= event_times(1) & added_time <= event_times(2));
baseline_mask = (added_time >= baseline_time_windows(1) & added_time <= baseline_time_windows(2));

% Create protocol masks
protocol_1_mask = (workflow_cat == 1);
protocol_2_mask = (workflow_cat == 2);

% Get unique animals
unique_animals = unique(stage_1_widefield_animal);
n_total_days = size(data_to_plot, 3);

% Initialize masks for first and last learning days
% Protocol 1
first_learning_p1_learners_mask = false(1, n_total_days);
last_learning_p1_learners_mask = false(1, n_total_days);
first_learning_p1_non_learners_mask = false(1, n_total_days);
last_learning_p1_non_learners_mask = false(1, n_total_days);

% Protocol 2
first_learning_p2_learners_mask = false(1, n_total_days);
last_learning_p2_learners_mask = false(1, n_total_days);
first_learning_p2_non_learners_mask = false(1, n_total_days);
last_learning_p2_non_learners_mask = false(1, n_total_days);

% Loop through each animal
for curr_animal = 1:length(unique_animals)
    % Get mask for this animal
    animal_mask = (widefield_animal_idx == curr_animal);
    
    % Determine if this animal is a learner
    isLearner = ~ismember(animal_list{curr_animal}, group_animals);
    
    % === PROTOCOL 1 ===
    % Combine animal mask with protocol 1 mask
    animal_p1_mask = animal_mask & protocol_1_mask;
    
    if any(animal_p1_mask)
        % Get learning days within protocol 1 for this animal
        learning_days_p1 = animal_p1_mask & learning_index_animal;
        
        % Get all indices for this animal in protocol 1
        animal_p1_indices = find(animal_p1_mask);
        
        if any(learning_days_p1)
            % Has learning days in protocol 1
            learning_p1_indices = find(learning_days_p1);
            first_day_p1 = learning_p1_indices(1);
            last_day_p1 = learning_p1_indices(end);
        else
            % No learning days - use first and second day of protocol 1
            if length(animal_p1_indices) >= 2
                first_day_p1 = animal_p1_indices(1);
                last_day_p1 = animal_p1_indices(2);
            elseif length(animal_p1_indices) == 1
                first_day_p1 = animal_p1_indices(1);
                last_day_p1 = animal_p1_indices(1);
            else
                continue;
            end
        end
        
        % Store in appropriate mask
        if isLearner
            first_learning_p1_learners_mask(first_day_p1) = true;
            last_learning_p1_learners_mask(last_day_p1) = true;
        else
            first_learning_p1_non_learners_mask(first_day_p1) = true;
            last_learning_p1_non_learners_mask(last_day_p1) = true;
        end
    end
    
    % === PROTOCOL 2 ===
    % Combine animal mask with protocol 2 mask
    animal_p2_mask = animal_mask & protocol_2_mask;
    
    if any(animal_p2_mask)
        % Get learning days within protocol 2 for this animal
        learning_days_p2 = animal_p2_mask & learning_index_animal;
        
        % Get all indices for this animal in protocol 2
        animal_p2_indices = find(animal_p2_mask);
        
        if any(learning_days_p2)
            % Has learning days in protocol 2
            learning_p2_indices = find(learning_days_p2);
            first_day_p2 = learning_p2_indices(1);
            last_day_p2 = learning_p2_indices(end);
        else
            % No learning days - use first and second day of protocol 2
            if length(animal_p2_indices) >= 2
                first_day_p2 = animal_p2_indices(1);
                last_day_p2 = animal_p2_indices(2);
            elseif length(animal_p2_indices) == 1
                first_day_p2 = animal_p2_indices(1);
                last_day_p2 = animal_p2_indices(1);
            else
                continue;
            end
        end
        
        % Store in appropriate mask
        if isLearner
            first_learning_p2_learners_mask(first_day_p2) = true;
            last_learning_p2_learners_mask(last_day_p2) = true;
        else
            first_learning_p2_non_learners_mask(first_day_p2) = true;
            last_learning_p2_non_learners_mask(last_day_p2) = true;
        end
    end
end

% ===== COMPUTE MEANS FOR PROTOCOL 1 =====
% Learners
mean_first_p1_learners = nanmean(data_to_plot(:,:,first_learning_p1_learners_mask), 3);
mean_last_p1_learners = nanmean(data_to_plot(:,:,last_learning_p1_learners_mask), 3);

% Non-learners
mean_first_p1_non_learners = nanmean(data_to_plot(:,:,first_learning_p1_non_learners_mask), 3);
mean_last_p1_non_learners = nanmean(data_to_plot(:,:,last_learning_p1_non_learners_mask), 3);

% Subtract baseline
baseline_sub_first_p1_learners = mean_first_p1_learners - nanmean(mean_first_p1_learners(:,baseline_mask), 2);
baseline_sub_last_p1_learners = mean_last_p1_learners - nanmean(mean_last_p1_learners(:,baseline_mask), 2);
baseline_sub_first_p1_non_learners = mean_first_p1_non_learners - nanmean(mean_first_p1_non_learners(:,baseline_mask), 2);
baseline_sub_last_p1_non_learners = mean_last_p1_non_learners - nanmean(mean_last_p1_non_learners(:,baseline_mask), 2);

% Get time window
timewindow_first_p1_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_first_p1_learners(:,time_window_idx));
timewindow_last_p1_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_last_p1_learners(:,time_window_idx));
timewindow_first_p1_non_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_first_p1_non_learners(:,time_window_idx));
timewindow_last_p1_non_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_last_p1_non_learners(:,time_window_idx));

% ===== COMPUTE MEANS FOR PROTOCOL 2 =====
% Learners
mean_first_p2_learners = nanmean(data_to_plot(:,:,first_learning_p2_learners_mask), 3);
mean_last_p2_learners = nanmean(data_to_plot(:,:,last_learning_p2_learners_mask), 3);

% Non-learners
mean_first_p2_non_learners = nanmean(data_to_plot(:,:,first_learning_p2_non_learners_mask), 3);
mean_last_p2_non_learners = nanmean(data_to_plot(:,:,last_learning_p2_non_learners_mask), 3);

% Subtract baseline
baseline_sub_first_p2_learners = mean_first_p2_learners - nanmean(mean_first_p2_learners(:,baseline_mask), 2);
baseline_sub_last_p2_learners = mean_last_p2_learners - nanmean(mean_last_p2_learners(:,baseline_mask), 2);
baseline_sub_first_p2_non_learners = mean_first_p2_non_learners - nanmean(mean_first_p2_non_learners(:,baseline_mask), 2);
baseline_sub_last_p2_non_learners = mean_last_p2_non_learners - nanmean(mean_last_p2_non_learners(:,baseline_mask), 2);

% Get time window
timewindow_first_p2_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_first_p2_learners(:,time_window_idx));
timewindow_last_p2_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_last_p2_learners(:,time_window_idx));
timewindow_first_p2_non_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baseline_sub_first_p2_non_learners(:,time_window_idx));
timewindow_last_p2_non_learners = plab.wf.svd2px(wf_U(:,:,1:n_components), baselineb_sub_last_p2_non_learners(:,time_window_idx));

% ===== COMPUTE IMAGES =====
switch statistic
    case 'max'
        % Protocol 1
        img_first_p1_learners = max(timewindow_first_p1_learners, [], 3);
        img_last_p1_learners = max(timewindow_last_p1_learners, [], 3);
        img_first_p1_non_learners = max(timewindow_first_p1_non_learners, [], 3);
        img_last_p1_non_learners = max(timewindow_last_p1_non_learners, [], 3);
        
        % Protocol 2
        img_first_p2_learners = max(timewindow_first_p2_learners, [], 3);
        img_last_p2_learners = max(timewindow_last_p2_learners, [], 3);
        img_first_p2_non_learners = max(timewindow_first_p2_non_learners, [], 3);
        img_last_p2_non_learners = max(timewindow_last_p2_non_learners, [], 3);

    case 'mean'
        % Protocol 1
        img_first_p1_learners = mean(timewindow_first_p1_learners, 3);
        img_last_p1_learners = mean(timewindow_last_p1_learners, 3);
        img_first_p1_non_learners = mean(timewindow_first_p1_non_learners, 3);
        img_last_p1_non_learners = mean(timewindow_last_p1_non_learners, 3);
        
        % Protocol 2
        img_first_p2_learners = mean(timewindow_first_p2_learners, 3);
        img_last_p2_learners = mean(timewindow_last_p2_learners, 3);
        img_first_p2_non_learners = mean(timewindow_first_p2_non_learners, 3);
        img_last_p2_non_learners = mean(timewindow_last_p2_non_learners, 3);

    otherwise
        error('statistic must be ''max'' or ''mean''');
end

% Compute common color limit across ALL conditions
all_vals = [timewindow_first_p1_learners(:); timewindow_last_p1_learners(:); ...
            timewindow_first_p1_non_learners(:); timewindow_last_p1_non_learners(:); ...
            timewindow_first_p2_learners(:); timewindow_last_p2_learners(:); ...
            timewindow_first_p2_non_learners(:); timewindow_last_p2_non_learners(:)];

switch statistic
    case 'max'
        ref_vals = max(reshape(all_vals,[],size(timewindow_first_p1_learners,3)),[],2);
    case 'mean'
        ref_vals = mean(reshape(all_vals,[],size(timewindow_first_p1_learners,3)),2);
end
clim_val = 0.8 * max(abs(ref_vals));

% Print summary
fprintf('\n===== Summary =====\n');
fprintf('PROTOCOL 1:\n');
fprintf('  mPFC+ First: n=%d animals\n', sum(first_learning_p1_learners_mask));
fprintf('  mPFC+ Last: n=%d animals\n', sum(last_learning_p1_learners_mask));
fprintf('  mPFC- First: n=%d animals\n', sum(first_learning_p1_non_learners_mask));
fprintf('  mPFC- Last: n=%d animals\n', sum(last_learning_p1_non_learners_mask));
fprintf('\nPROTOCOL 2:\n');
fprintf('  mPFC+ First: n=%d animals\n', sum(first_learning_p2_learners_mask));
fprintf('  mPFC+ Last: n=%d animals\n', sum(last_learning_p2_learners_mask));
fprintf('  mPFC- First: n=%d animals\n', sum(first_learning_p2_non_learners_mask));
fprintf('  mPFC- Last: n=%d animals\n', sum(last_learning_p2_non_learners_mask));

% ===== PLOTTING =====
% Figure 1: Protocol 1 (2x2 layout)
figure('Position', [100, 100, 1400, 1200], 'Color', 'w', 'Name', 'Protocol 1');

subplot(2, 2, 1);
imagesc(img_first_p1_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 1: mPFC+ First Learning Day\n%s (n=%d)', statistic, sum(first_learning_p1_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 2);
imagesc(img_last_p1_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 1: mPFC+ Last Learning Day\n%s (n=%d)', statistic, sum(last_learning_p1_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 3);
imagesc(img_first_p1_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 1: mPFC- First Learning Day\n%s (n=%d)', statistic, sum(first_learning_p1_non_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 4);
imagesc(img_last_p1_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 1: mPFC- Last Learning Day\n%s (n=%d)', statistic, sum(last_learning_p1_non_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

sgtitle('Protocol 1: Widefield Activity (First vs Last Learning Day)', ...
        'FontSize', 15, 'FontWeight', 'bold');

% Figure 2: Protocol 2 (2x2 layout)
figure('Position', [150, 150, 1400, 1200], 'Color', 'w', 'Name', 'Protocol 2');

subplot(2, 2, 1);
imagesc(img_first_p2_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 2: mPFC+ First Learning Day\n%s (n=%d)', statistic, sum(first_learning_p2_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 2);
imagesc(img_last_p2_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 2: mPFC+ Last Learning Day\n%s (n=%d)', statistic, sum(last_learning_p2_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 3);
imagesc(img_first_p2_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 2: mPFC- First Learning Day\n%s (n=%d)', statistic, sum(first_learning_p2_non_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

subplot(2, 2, 4);
imagesc(img_last_p2_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Protocol 2: mPFC- Last Learning Day\n%s (n=%d)', statistic, sum(last_learning_p2_non_learners_mask)), ...
      'FontSize', 12, 'FontWeight', 'bold');
ap.wf_draw('ccf');

sgtitle('Protocol 2: Widefield Activity (First vs Last Learning Day)', ...
        'FontSize', 15, 'FontWeight', 'bold');



%% Explorartory : Plot mPFC response between learners vs non-learners for right/centre stimulus


% logical mask for stage 1
wf1_idx = ((workflow_cat==1));% length K

% days belonging to protocol 1 & 2 (in case of final position)
wf1_2_idx= find(workflow_cat ~=3);

K = size(rewarded_stim_start_move_kernel,3);

assert(K == numel(wf1_idx), 'Kernel 3rd dim must equal #wf==1 days');

% build masks within the first stage subset
learning_idx_sub   = learning_index_animal(wf1_idx);    % K×1
grp_sub  = is_group_animal(wf1_idx);          % 1 non-learners / 0 learners
workflow_sub= workflow_cat(wf1_idx);

% pre_non_Learners = (workflow_sub==1) & (learning_idx_sub==0) & (grp_sub==0);        % K×1 logical
% post_non_Learners = (workflow_sub==1) &(learning_idx_sub==1) & (grp_sub==0); 
% 
% pre_Learners= (workflow_sub==1) &(learning_idx_sub==0) & (grp_sub==1);  
% post_Learners= (workflow_sub==1) &(learning_idx_sub==1) & (grp_sub==1);  

pre_non_Learners = (learning_idx_sub==0) & (grp_sub==1);        % K×1 logical
post_non_Learners = (learning_idx_sub==1) & (grp_sub==1); 

pre_Learners= (learning_idx_sub==0) & (grp_sub==0);  
post_Learners= (learning_idx_sub==1) & (grp_sub==0);  

% pre post of a selected workflow
pre = nanmean(rewarded_stim_v_stacked_data(:,:,pre_non_Learners),3);
post = nanmean(rewarded_stim_v_stacked_data(:,:,post_non_Learners),3);

ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:n_components),pre),plab.wf.svd2px(wf_U(:,:,1:n_components),post)],added_time);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;




%% Plot CCF average/max mPFC activity (Kernels) comparing learners vs non-learners

% Define the time window in seconds (0 to 200ms)
start_time = 0; % 0s (0ms)
end_time = 0.3; % 0.2s (200ms)

% Choose 'max' or 'mean' to plot
statistic = 'mean';

% define data to plot
data_to_plot= rewarded_stim_start_move_kernel;

% define protocol
protocol_idx=1;

added_time_Kernel=  fliplr((-10:52)/30);

% Find the corresponding indices in added_time for the time window
time_window_idx = find(added_time_Kernel >= start_time & added_time_Kernel <= end_time);


% pre post of right move workflow split by learning days and learners vs
% non-learners

mean_pre_learning_learners = nanmean(data_to_plot(:,:,pre_Learners),3);
mean_post_learning_learners = nanmean(data_to_plot(:,:,post_Learners),3);

mean_pre_learning_non_learners = nanmean(data_to_plot(:,:,pre_non_Learners),3);
mean_post_learning_non_learners = nanmean(data_to_plot(:,:,post_non_Learners),3);

% get the specific time window
timewindow_mean_pre_learning_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), mean_pre_learning_learners(:,time_window_idx));
timewindow_mean_post_learning_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), mean_post_learning_learners(:,time_window_idx));

timewindow_mean_pre_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), mean_pre_learning_non_learners(:,time_window_idx));
timewindow_mean_post_learning_non_learners= plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), mean_post_learning_non_learners(:,time_window_idx));

% Compute the image to show
switch statistic
    case 'max'
        img_pre_learning_learners = max(timewindow_mean_pre_learning_learners, [], 3);
        img_post_learning_learners= max(timewindow_mean_post_learning_learners, [], 3);
        img_pre_learning_non_learners   = max(timewindow_mean_pre_learning_non_learners, [], 3);
        img_post_learning_non_learners   = max(timewindow_mean_post_learning_non_learners, [], 3);

    case 'mean'
        img_pre_learning_learners = mean(timewindow_mean_pre_learning_learners, 3);
        img_post_learning_learners= mean(timewindow_mean_post_learning_learners, 3);
        img_pre_learning_non_learners = mean(timewindow_mean_pre_learning_non_learners, 3);
        img_post_learning_non_learners = mean(timewindow_mean_post_learning_non_learners,3);

    otherwise
        error('statistic must be ''max'' or ''mean''');
end

% Compute a common color‐limit
all_vals = [timewindow_mean_pre_learning_learners(:); timewindow_mean_pre_learning_non_learners(:)];
switch statistic
    case 'max'
        ref_vals = max(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),[],2);
    case 'mean'
        ref_vals = mean(reshape(all_vals,[],size(timewindow_mean_pre_learning_learners,3)),2);
end
clim_val = 0.5 * max(abs(ref_vals));

% Plot learners
figure;
imagesc(img_pre_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s Kernels — Learners', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s Kernels — Learners', statistic));
ap.wf_draw('ccf');

% Plot non learners
figure;
imagesc(img_pre_learning_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Pre-Learning %s Kernels— Non-learners', statistic));
ap.wf_draw('ccf');

figure;
imagesc(img_post_learning_non_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s Kernels — Non-learners', statistic));
ap.wf_draw('ccf');

%% Plot line variables setup basline substracted

event_times= [-0.5,0.75];
baseline_time_windows = [-0.1, 0];

% data to index:  V×T×Ndays_total
aligned_all = rewarded_stim_kernel;

% define the added time 
added_time_V= added_time;
added_time_Kernel;

% Find the corresponding indices in added_time for the time window
time_window_idx = find(added_time_Kernel >= event_times(1) & added_time_Kernel <= event_times(2));
baseline_mask_idx = (added_time_Kernel >= baseline_time_windows(1) & added_time_Kernel <= baseline_time_windows(2));

% % 1×T win vector in seconds
t_for_plot = added_time_Kernel(time_window_idx);

n_components_to_use= kernel_n_components; % kernel_n_components or n_components

right_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,time_window_idx,:),[],[], right_mPFC_ROI_mask), ...
    [3,2,1] );   % gives Ndays_total × T

right_mPFC_ROI_baseline=permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,baseline_mask_idx,:),[],[], right_mPFC_ROI_mask), ...
    [3,2,1] );   % gives Ndays_total × T

right_mPFC_ROI_trace= right_mPFC_ROI_trace- mean(right_mPFC_ROI_baseline,2);

left_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,time_window_idx,:),[],[], left_mPFC_ROI_mask), ...
    [3,2,1] );   % same dims

left_mPFC_ROI_baseline=permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,baseline_mask_idx,:),[],[], left_mPFC_ROI_mask), ...
    [3,2,1] );   % gives Ndays_total × T

left_mPFC_ROI_trace= left_mPFC_ROI_trace- mean(left_mPFC_ROI_baseline,2);

% For visual cortex
right_ViS_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,time_window_idx,:),[],[], right_ViS_ROI_mask), ...
    [3,2,1] );   % gives Ndays_total × T

left_ViS_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:n_components_to_use), aligned_all(:,time_window_idx,:),[],[], left_ViS_ROI_mask), ...
    [3,2,1] );   % gives Ndays_total × T





%% Plot individual animal plots for all days colorcoded baseline substracted


protocol_idx = 1;  % define protocol

% Get unique animals
animals = unique(widefield_animal_idx);
n_animals = length(animals);

% Create figure with subplots
fig = figure('Color', 'w', 'Position', [100, 100, 1800, 1200]);

% Define colormap for progression (early days = cool, late days = warm)
cmap = cool(10);  % Adjust max number of days as needed

subplot_idx = 1;

for ai = animals(:)'
    % ——— 1) Pick out this animal & workflow ———
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    n_days = length(days_idx);
    
    % ——— 2) Locate the learning day ———
    learn_local = learning_index_animal(days_idx) == 1;
    ld = find(learn_local, 1, 'first');
    
    if isempty(ld)
        ld = n_days + 1;  % All pre-learning
    end
    
    % ——— 3) Pull out the day×time traces ———
    TR_pre = right_mPFC_ROI_trace(days_idx, :);   % n_days × T
    TL_pre = left_mPFC_ROI_trace(days_idx, :);    % n_days × T
    TL_baseline= left_mPFC_ROI_baseline(days_idx,:); 
    
    % ——— 4) Create subplot for this animal ———
    subplot(ceil(sqrt(n_animals)), ceil(n_animals/ceil(sqrt(n_animals))), subplot_idx);
    hold on;
    
    % Plot each day with color coding
    for d = 1:n_days
        % Color based on day progression
        day_color = cmap(min(d, size(cmap, 1)), :);
        
        % Line style: dashed before learning, solid after
        if d < ld
            line_style = '--';
            line_width = 1.5;
        else
            line_style = '-';
            line_width = 2.5;
        end
        
        % Plot right mPFC (this day)
        plot(t_for_plot, TR_pre(d, :), ...
             'Color', day_color, ...
             'LineWidth', line_width, ...
             'LineStyle', line_style);
    end
    
    % Mark learning day if exists
    if ld <= n_days
        % Plot learning day prominently
        plot(t_for_plot, TR_pre(ld, :), ...
             'Color', 'g', ...
             'LineWidth', 3, ...
             'LineStyle', '-');
    end
    
    % Mark event time
    xline(0, 'k--', 'LineWidth', 1.5);
    
    xlabel('Time (s)', 'FontSize', 10);
    ylabel('mPFC Activity', 'FontSize', 10);
    title(sprintf('%s (n=%d days)', animal_list{ai}, n_days), ...
          'FontSize', 11, 'FontWeight', 'bold');
    grid on;
    
    % Add text indicating learning day
    if ld <= n_days
        text(0.02, 0.98, sprintf('LD: Day %d', ld), ...
             'Units', 'normalized', ...
             'VerticalAlignment', 'top', ...
             'FontSize', 9, ...
             'BackgroundColor', 'w');
    end
    
    subplot_idx = subplot_idx + 1;
end

% Add colorbar for day progression
cb = colorbar('Position', [0.92, 0.3, 0.02, 0.4]);
colormap(cmap);
clim([1, max(arrayfun(@(a) sum((widefield_animal_idx==a) & (workflow_cat==protocol_idx)), animals))]);
cb.Label.String = 'Day Number';
cb.Label.FontSize = 12;

% After all subplots are created, add:
all_axes = findobj(fig, 'Type', 'axes', '-not', 'Tag', 'Colorbar');
y_limits = [min(arrayfun(@(ax) ax.YLim(1), all_axes)), ...
            max(arrayfun(@(ax) ax.YLim(2), all_axes))];
set(all_axes, 'YLim', y_limits);

% Overall title
sgtitle(sprintf('Right mPFC Activity Across Days (Protocol %d) - Dashed=Pre-Learning, Solid=Post-Learning', protocol_idx), ...
        'FontSize', 14, 'FontWeight', 'bold');



%% Plot pre and post Scatter for Right vs Left mPFC responses for all animals

protocol_idx = 1;
animals = unique(widefield_animal_idx);

% Storage
animal_mPFC_right_pre = []; animal_mPFC_left_pre = [];
animal_mPFC_right_post = []; animal_mPFC_left_post = [];

for ai = animals(:)'
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    if isempty(days_idx), continue; end
    
    learn_local = learning_index_animal(days_idx) == 1;
    ld = find(learn_local, 1, 'first');
    if isempty(ld), continue; end
    
    pre_days = find((1:length(days_idx)) < ld);
    post_days = find((1:length(days_idx)) >= ld);
    if isempty(post_days) || isempty(pre_days), continue; end
    
    response_window = t_for_plot >= 0 & t_for_plot <= 0.35;
    
    % Pre-learning
    TR_pre = right_mPFC_ROI_trace(days_idx(pre_days), :);
    TL_pre = left_mPFC_ROI_trace(days_idx(pre_days), :);

    % average across days
    TR_pre_mean= mean(TR_pre, 1, 'omitnan');
    TL_pre_mean= mean(TL_pre, 1, 'omitnan');

    TR_pre_peak= max(TR_pre_mean(response_window));
    TL_pre_peak= max(TL_pre_mean(response_window));

    animal_mPFC_right_pre = [animal_mPFC_right_pre; TR_pre_peak];
    animal_mPFC_left_pre = [animal_mPFC_left_pre; TL_pre_peak];
    
    % Post-learning
    TR_post = right_mPFC_ROI_trace(days_idx(post_days), :);
    TL_post = left_mPFC_ROI_trace(days_idx(post_days), :);

     % average across days
    TR_post_mean= mean(TR_post, 1, 'omitnan');
    TL_post_mean= mean(TL_post, 1, 'omitnan');

    TR_post_peak= max(TR_post_mean(response_window));
    TL_post_peak= max(TL_post_mean(response_window));

    animal_mPFC_right_post = [animal_mPFC_right_post; TR_post_peak];
    animal_mPFC_left_post = [animal_mPFC_left_post; TL_post_peak];
end


% Plot
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
max_val = max([animal_mPFC_right_pre; animal_mPFC_left_pre; ...
               animal_mPFC_right_post; animal_mPFC_left_post]);

for sp = 1:2
    subplot(1, 2, sp);
    if sp == 1
        data_x = animal_mPFC_right_pre; data_y = animal_mPFC_left_pre;
        tit = 'Pre-Learning';
    else
        data_x = animal_mPFC_right_post; data_y = animal_mPFC_left_post;
        tit = 'Post-Learning';
    end
    
    scatter(data_x, data_y, 150, [0.3 0.5 0.8], 'filled', ...
            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'LineWidth', 1.5);
    hold on;
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 2);
    
    xlabel('Right mPFC Activity', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Left mPFC Activity', 'FontSize', 13, 'FontWeight', 'bold');
    title(tit, 'FontSize', 14, 'FontWeight', 'bold');
    axis equal; xlim([0, max_val*1.05]); ylim([0, max_val*1.05]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
end

%% Scatter plot: Individual animal mPFC responses (post-learning)

protocol_idx = 1;
animals = unique(widefield_animal_idx);

% Storage
animal_mPFC_post = [];
animal_names = {};
animal_group_labels = [];  % 1=mPFC+, 2=mPFC-, 3=Non-learner

% Define groups (adjust these to your actual groupings)
mPFC_plus_grp = {'HA005','HA008','HA010','HA011','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA009','HA014','HA015'};
non_learners_grp = {'HA006','HA013','AP030','AP031','AP032'};

for ai = animals(:)'
    % Get this animal's data
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    % Get learning day
    learn_local = learning_index_animal(days_idx) == 1;
    ld = find(learn_local, 1, 'first');
    
    if isempty(ld)
        % No learning day - skip or treat as all pre-learning
        continue;
    end
    
    % Get POST-learning days only
    post_days = find((1:length(days_idx)) >= ld);
    pre_days=  find((1:length(days_idx)) < ld);
    if isempty(post_days), continue; end
    
    % Get traces
    TR_pre = left_mPFC_ROI_trace(days_idx(pre_days), :);
    
    % Average across post-learning days
    mean_trace = mean(TR_pre, 1, 'omitnan');
    
    % Define response window (0 to 0.5s based on your plots)
    response_window = t_for_plot >= 0 & t_for_plot <= 0.35;
    
    % Compute summary metric (peak or mean)
    peak_response = max(mean_trace(response_window));
    % OR area under curve:
    % peak_response = trapz(t_for_plot(response_window), mean_trace(response_window));
    
    % Store
    animal_mPFC_post = [animal_mPFC_post; peak_response];
    animal_names{end+1} = animal_list{ai};
    
    % Assign group
    if ismember(animal_list{ai}, mPFC_plus_grp)
        animal_group_labels = [animal_group_labels; 1];
    elseif ismember(animal_list{ai}, mPFC_minus_grp)
        animal_group_labels = [animal_group_labels; 2];
    elseif ismember(animal_list{ai}, non_learners_grp)
        animal_group_labels = [animal_group_labels; 3];
    else
        animal_group_labels = [animal_group_labels; 0];  % Unknown
    end
end

% Create scatter plot

group_names = {'mPFC+', 'mPFC-', 'Non-learners', 'Unknown'};
group_colors = {cLearner, cNonLearner, [0.5 0.5 0.5], [0.8 0.8 0.8]};

figure('Color', 'w', 'Position', [100, 100, 1000, 700]);

hold on;

% Plot each group
for g = 1:3
    group_mask = animal_group_labels == g;
    
    if sum(group_mask) == 0, continue; end
    
    % X positions with jitter for visibility
    x_positions = g * ones(sum(group_mask), 1);
    x_jitter = x_positions + 0.15 * randn(sum(group_mask), 1);
    
    % Scatter individual animals
    scatter(x_jitter, animal_mPFC_post(group_mask), 120, ...
            'MarkerFaceColor', group_colors{g}, ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceAlpha', 0.7, ...
            'LineWidth', 1.5);
    
    % Add animal labels
    for i = find(group_mask)'
        text(x_jitter(sum(group_mask(1:i)) == sum(group_mask(1:i))), ...
             animal_mPFC_post(i), ...
             ['  ' animal_names{i}], ...
             'FontSize', 9, ...
             'HorizontalAlignment', 'left');
    end
    
    % Add group mean line
    group_mean = mean(animal_mPFC_post(group_mask));
    plot([g-0.3, g+0.3], [group_mean, group_mean], ...
         'Color', group_colors{g}, 'LineWidth', 4);
    
    % Add SEM error bar
    group_sem = std(animal_mPFC_post(group_mask)) / sqrt(sum(group_mask));
    errorbar(g, group_mean, group_sem, ...
             'Color', group_colors{g}, 'LineWidth', 2.5, ...
             'CapSize', 10, 'HandleVisibility', 'off');
end

% Formatting
xlim([0.5, 3.5]);
xticks(1:3);
xticklabels(group_names(1:3));
ylabel('Peak mPFC Activity (Post-Learning)', 'FontSize', 13, 'FontWeight', 'bold');
title('Individual Variability in Post-Learning mPFC Responses', ...
      'FontSize', 14, 'FontWeight', 'bold');
grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.2);

% Add statistics annotation
fprintf('\n=== Post-Learning mPFC Activity ===\n');
for g = 1:3
    group_mask = animal_group_labels == g;
    if sum(group_mask) > 0
        fprintf('%s (n=%d): Mean=%.5f, SD=%.5f, Range=[%.5f, %.5f]\n', ...
                group_names{g}, sum(group_mask), ...
                mean(animal_mPFC_post(group_mask)), ...
                std(animal_mPFC_post(group_mask)), ...
                min(animal_mPFC_post(group_mask)), ...
                max(animal_mPFC_post(group_mask)));
    end
end

% Statistical test
if sum(animal_group_labels == 1) > 0 && sum(animal_group_labels == 2) > 0
    [p_12, ~] = ranksum(animal_mPFC_post(animal_group_labels == 1), ...
                        animal_mPFC_post(animal_group_labels == 2));
    fprintf('\nmPFC+ vs mPFC-: p = %.4f\n', p_12);
end

if sum(animal_group_labels == 1) > 0 && sum(animal_group_labels == 3) > 0
    [p_13, ~] = ranksum(animal_mPFC_post(animal_group_labels == 1), ...
                        animal_mPFC_post(animal_group_labels == 3));
    fprintf('mPFC+ vs Non-learners: p = %.4f\n', p_13);
end

%% Simple scatter: All animals ranked by response

% Sort by response magnitude
[sorted_responses, sort_idx] = sort(animal_mPFC_post, 'descend');
sorted_names = animal_names(sort_idx);

figure('Color', 'w', 'Position', [100, 100, 1000, 700]);

% Create scatter
scatter(1:length(sorted_responses), sorted_responses, 150, ...
        sorted_responses, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5);

% Colormap
colormap(hot);
% colorbar('Label', 'Peak Activity');

% Add animal labels
for i = 1:length(sorted_names)
    text(i, sorted_responses(i) + 0.0002, sorted_names{i}, ...
         'Rotation', 45, 'FontSize', 10, 'HorizontalAlignment', 'left');
end

% Formatting
xlabel('Animals (ranked by response)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Peak mPFC Activity (Post-Learning)', 'FontSize', 13, 'FontWeight', 'bold');
title('Individual Variability in Post-Learning mPFC Responses', ...
      'FontSize', 14, 'FontWeight', 'bold');

xlim([0, length(sorted_responses) + 1]);
ylim([0,max(sorted_responses)+0.0005]);
set(gca, 'FontSize', 12);

% Add statistics
text(0.7, 0.95, sprintf('Range: %.4f - %.4f\nCV: %.1f%%', ...
                        min(sorted_responses), max(sorted_responses), ...
                        std(sorted_responses)/mean(sorted_responses)*100), ...
     'Units', 'normalized', 'FontSize', 11, 'BackgroundColor', 'w', ...
     'EdgeColor', 'k', 'LineWidth', 1);

%% Plot pre-post average split between learners and non-learners for both hemispheres

protocol_idx=1;  % define protocol
animals= grp2idx(animal_list);

% ————— Prepare containers —————
Rpre_L = [];  Rpost_L = [];
Lpre_L = [];  Lpost_L = [];
Rpre_N = [];  Rpost_N = [];
Lpre_N = [];  Lpost_N = [];

for ai = animals(:)'
    % ——— pick out this animal & protocol days ———
    sel      = (widefield_animal_idx==ai) & (workflow_cat==protocol_idx);
    days_idx = find(sel);
    if isempty(days_idx), continue; end

    % ——— find the learning day within those ———
    ld = find(learning_index_animal(days_idx)==1 , 1,'first' );
    if isempty(ld), ld=length(days_idx)+1; end

    % ——— grab the day×time ROI traces ———
    TR_pre = right_mPFC_ROI_trace(days_idx, :);  % n_days×T
    TL_pre = left_mPFC_ROI_trace(days_idx,  :);  % n_days×T

    % ——— split into pre vs post days ———
    pre_mask  = (1:size(TR_pre,1)) < ld;
    post_mask = ~pre_mask;

    R_pre  = mean(TR_pre(pre_mask, :), 1);
    R_post = mean(TR_pre(post_mask,:), 1);
    L_pre  = mean(TL_pre(pre_mask, :), 1);
    L_post = mean(TL_pre(post_mask,:), 1);

    % ——— decide which group this animal is in ———
    grp = is_group_animal(days_idx(1));   % 0 = learner, 1 = non-learner

    if grp == 0
        Rpre_L(end+1,:)  = R_pre;
        Rpost_L(end+1,:) = R_post;
        Lpre_L(end+1,:)  = L_pre;
        Lpost_L(end+1,:) = L_post;
    else
        Rpre_N(end+1,:)  = R_pre;
        Rpost_N(end+1,:) = R_post;
        Lpre_N(end+1,:)  = L_pre;
        Lpost_N(end+1,:) = L_post;
    end
end

% ————— compute group means & SEMs —————
% learners
nLeaners = size(Rpre_L,1);
grpR_pre_L  = nanmean(Rpre_L,1);    semR_pre_L  = nanstd(Rpre_L,0,1)/sqrt(nLeaners);
grpR_post_L = nanmean(Rpost_L,1);   semR_post_L = nanstd(Rpost_L,0,1)/sqrt(nLeaners);
grpL_pre_L  = nanmean(Lpre_L,1);    semL_pre_L  = nanstd(Lpre_L,0,1)/sqrt(nLeaners);
grpL_post_L = nanmean(Lpost_L,1);   semL_post_L = nanstd(Lpost_L,0,1)/sqrt(nLeaners);

% non‐learners
nNon_Leaners = size(Rpre_N,1);
grpR_pre_N  = nanmean(Rpre_N,1);    semR_pre_N  = nanstd(Rpre_N,0,1)/sqrt(nNon_Leaners);
grpR_post_N = nanmean(Rpost_N,1);   semR_post_N = nanstd(Rpost_N,0,1)/sqrt(nNon_Leaners);
grpL_pre_N  = nanmean(Lpre_N,1);    semL_pre_N  = nanstd(Lpre_N,0,1)/sqrt(nNon_Leaners);
grpL_post_N = nanmean(Lpost_N,1);   semL_post_N = nanstd(Lpost_N,0,1)/sqrt(nNon_Leaners);

% ————— plot with SEM shading —————
figure('Color','w','Position',[200 200 800 600]);

% Right mPFC
ax1 = subplot(2,1,1); hold(ax1,'on');

% Right mPFC
subplot(2,1,1); hold on;
% learner pre
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpR_pre_L+semR_pre_L, fliplr(grpR_pre_L-semR_pre_L)], ...
    cLearner,'FaceAlpha',0.2,'EdgeColor','none');
% learner post
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpR_post_L+semR_post_L, fliplr(grpR_post_L-semR_post_L)], ...
    cLearner,'FaceAlpha',0.1,'EdgeColor','none');
% nonlearner pre
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpR_pre_N+semR_pre_N, fliplr(grpR_pre_N-semR_pre_N)], ...
    cNonLearner,'FaceAlpha',0.2,'EdgeColor','none');
% nonlearner post
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpR_post_N+semR_post_N, fliplr(grpR_post_N-semR_post_N)], ...
    cNonLearner,'FaceAlpha',0.1,'EdgeColor','none');

% now lines on top
plot(t_for_plot, grpR_pre_L,  'color',cLearner,'LineWidth',2,LineStyle='--');
plot(t_for_plot, grpR_post_L, 'color',cLearner,'LineWidth',2,LineStyle='-');
plot(t_for_plot, grpR_pre_N,  'color',cNonLearner,'LineWidth',2,LineStyle='--');
plot(t_for_plot, grpR_post_N, 'color',cNonLearner,'LineWidth',2,LineStyle='-');

xline(ax1,0,'k--','LineWidth',1.5);
title(ax1,'Right mPFC (Kernels)');
ylabel(ax1,'Mean activity');
legend(ax1,{'','','',''...
    'Learners Pre','Learners Post','Non-Learners Pre','Non-Learners Post'},'Location','best');grid(ax1,'on');

% Left mPFC
ax2 = subplot(2,1,2); hold(ax2,'on');

% Left mPFC
subplot(2,1,2); hold on;
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpL_pre_L+semL_pre_L, fliplr(grpL_pre_L-semL_pre_L)], ...
    cLearner,'FaceAlpha',0.2,'EdgeColor','none');
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpL_post_L+semL_post_L, fliplr(grpL_post_L-semL_post_L)], ...
    cLearner,'FaceAlpha',0.1,'EdgeColor','none');
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpL_pre_N+semL_pre_N, fliplr(grpL_pre_N-semL_pre_N)], ...
    cNonLearner,'FaceAlpha',0.2,'EdgeColor','none');
fill([t_for_plot, fliplr(t_for_plot)], ...
    [grpL_post_N+semL_post_N, fliplr(grpL_post_N-semL_post_N)], ...
    cNonLearner,'FaceAlpha',0.1,'EdgeColor','none');

plot(t_for_plot, grpL_pre_L,  'color',cLearner,'LineWidth',2, LineStyle='--');
plot(t_for_plot, grpL_post_L, 'color',cLearner,'LineWidth',2, LineStyle='-');
plot(t_for_plot, grpL_pre_N, 'color',cNonLearner,'LineWidth',2, LineStyle='--');
plot(t_for_plot, grpL_post_N, 'color',cNonLearner,'LineWidth',2, LineStyle='-');

xline(ax2,0,'k--','LineWidth',1.5);
title(ax2,'Left mPFC (Kernels)');
xlabel(ax2,'Time (s)');
ylabel(ax2,'Mean activity');
legend(ax2,{'','','',''...
    'Learners Pre','Learners Post','Non-Learners Pre','Non-Learners Post'},'Location','best');
grid(ax2,'on');

% ——— Make y-limits identical ———
yl1 = ylim(ax1);
yl2 = ylim(ax2);
ymin = min(yl1(1), yl2(1));
ymax = max(yl1(2), yl2(2));
set([ax1, ax2], 'YLim', [ymin ymax]);

sgtitle('Group‐average pre vs post learning with SEM');
 

%% Plot the mPFC activity for the first 20 trials in the second stage-
%% ===== Step 1: Stack ALL trial data and create comprehensive logical masks =====


fprintf('\n===== Creating trial-level logical masks =====\n');

% Pre-allocate arrays for metadata (much faster than structure growth)


animal_ids = nan(n_total_trials, 1);
session_indices = nan(n_total_trials, 1);
trial_indices_in_session = nan(n_total_trials, 1);
learning_stages = nan(n_total_trials, 1);
is_nonlearner = false(n_total_trials, 1);



% Loop once to populate metadata arrays and to ensure no empty trials
trial_counter = 0;
for idx = 1:length(widefield_cat)
    
    day_V_data = widefield_cat(idx).rewarded_stim_on_aligned_V;
    
    if isempty(day_V_data)
        continue;
    end
    
    if workflow_cat(idx)~=1 % skip stage 2 and 3 for stim move analysis
        continue;
    end
    n_trials_this_session = size(day_V_data, 3);
    
    % Store metadata for all trials in this session
    trial_range = trial_counter + (1:n_trials_this_session);
    
    animal_ids(trial_range) = widefield_animal_idx(idx);
    session_indices(trial_range) = idx;
    trial_indices_in_session(trial_range) = 1:n_trials_this_session;
    learning_stages(trial_range) = learning_index_animal(idx);
    is_nonlearner(trial_range) = is_group_animal(idx);
    workflows(trial_range) = workflow_cat(idx);
    
    trial_counter = trial_counter + n_trials_this_session;
end

% Trim to actual size (in case some sessions were empty)
animal_ids = animal_ids(1:trial_counter);
session_indices = session_indices(1:trial_counter);
trial_indices_in_session = trial_indices_in_session(1:trial_counter);
learning_stages = learning_stages(1:trial_counter);
is_nonlearner = is_nonlearner(1:trial_counter);
workflows = workflows(1:trial_counter);

n_total_trials = trial_counter;

fprintf('Actual trials after filtering: %d\n', n_total_trials);

% ===== Step 2: Create boolean masks =====

% Learning group masks
is_learner = ~is_nonlearner;

% Learning stage masks
is_pre_learning = (learning_stages == 0);
is_post_learning = (learning_stages == 1);

% Workflow masks
is_workflow_1 = (workflows == 1); % Stage 1
is_workflow_2 = (workflows == 2); % Stage 2


% ===== Step 3: Create specific mask for learners Stage 2 =====

% Mask: Learners + Stage 2 (workflow 2) + Post-learning
learners_stage2_mask = is_learner & is_workflow_2 & is_post_learning;

% Get unique animals in this group
learner_animal_ids = unique(animal_ids(learners_stage2_mask));
n_learner_animals = length(learner_animal_ids);


% ===== Step 4: For each learner, find first day of Stage 2 and extract first N trials =====

% Parameters
n_trials_to_extract = 50; % Number of trials to extract
ROI_to_extract = left_ViS_ROI_mask; % Which ROI to plot
curr_components = n_components;

% Storage for extracted data
learner_first_day_data = struct();

for a = 1:n_learner_animals
    animal_id = learner_animal_ids(a);
    
    % Mask: This animal + Stage 2 learner trials
    this_animal_stage2_mask = learners_stage2_mask & (animal_ids == animal_id);
    
    % Find all Stage 2 sessions for this animal
    animal_sessions = unique(session_indices(this_animal_stage2_mask));
    
    % Get FIRST Stage 2 session
    first_session_idx = animal_sessions(1);
    
    % Mask: This animal + First Stage 2 day only
    first_day_mask = this_animal_stage2_mask & (session_indices == first_session_idx);
    
    % Get trial indices (global indices into metadata arrays)
    first_day_trial_indices = find(first_day_mask);
    
    % Get either the n trials or the lower
    n_trials_to_use = min(n_trials_to_extract, length(first_day_trial_indices));
    
    % Select first N trials
    trials_to_extract = first_day_trial_indices(1:n_trials_to_use);
    
  
    
    % Store metadata
    learner_first_day_data(a).animal_id = animal_id;
    learner_first_day_data(a).first_session_idx = first_session_idx;
    learner_first_day_data(a).n_trials = n_trials_to_use;
    learner_first_day_data(a).global_trial_indices = trials_to_extract;
    
    % Pre-allocate for ROI traces
    learner_first_day_data(a).roi_traces = nan(n_trials_to_use, length(added_time));
    
    % ===== NOW loop to extract ROI traces for these specific trials =====
    for t = 1:n_trials_to_use
        
        global_trial_idx = trials_to_extract(t);
        
        % Get session and trial-within-session for this global trial
        sess_idx = session_indices(global_trial_idx);
        trial_in_sess = trial_indices_in_session(global_trial_idx);
        
        % Extract V data for this specific trial
        trial_V = widefield_cat(sess_idx).rewarded_stim_on_aligned_V(:, :, trial_in_sess);
        learner_first_day_data(a).trial_V(:,:,t) = trial_V;
        % Store trial V as well

        % Extract ROI activity
        roi_trace = ap.wf_roi(wf_U(:,:,1:curr_components), trial_V, [], [], ROI_to_extract);
        roi_trace = squeeze(roi_trace); % Make 1D
        
        % Store
        learner_first_day_data(a).roi_traces(t, :) = roi_trace;
    end
    
    % Calculate average trace
    learner_first_day_data(a).roi_avg = mean(learner_first_day_data(a).roi_traces, 1, 'omitnan');
   
end


% ===== Step 5: Plot the results =====

figure('Position', [100 100 1400 300*n_learner_animals], 'Color', 'w');

for a = 1:n_learner_animals
    subplot(n_learner_animals, 1, a);
    hold on;
    
    % Plot individual trials
    for trial = 1:learner_first_day_data(a).n_trials
        plot(added_time, learner_first_day_data(a).roi_traces(trial, :), ...
            'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    
    % Plot average
    plot(added_time, learner_first_day_data(a).roi_avg, ...
        'Color', cLearner, 'LineWidth', 3);
    
    % Reference lines
    xline(0, 'k--', 'LineWidth', 2);
    yline(0, 'k:', 'LineWidth', 1);
    
    % Labels
    ylabel('ΔF/F', 'FontSize', 12);
    title(sprintf('Animal %d - First %d trials of Stage 2 Day 1 (Session %d)', ...
        learner_first_day_data(a).animal_id, ...
        learner_first_day_data(a).n_trials, ...
        learner_first_day_data(a).first_session_idx), ...
        'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    
    if a == n_learner_animals
        xlabel('Time from CS+ onset (s)', 'FontSize', 12);
    end
    
    if a == 1
        legend({'Individual trials', 'Average'}, 'Location', 'best', 'FontSize', 10);
    end
    
    set(gca, 'FontSize', 11, 'LineWidth', 1.5);
end

sgtitle('Learners: First N Trials of Stage 2 Day 1 - Left mPFC', ...
    'FontSize', 16, 'FontWeight', 'bold');

% % CCF videos
% animal_1_example= learner_first_day_data(3).trial_V;
% animal_1_example= mean(animal_1_example,3);
% 
% animal_2_example= learner_first_day_data(4).trial_V;
% animal_2_example= mean(animal_2_example,3);
% 
% ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:n_components),animal_1_example),plab.wf.svd2px(wf_U(:,:,1:n_components),animal_2_example)],added_time);
% clim([-max(abs(clim)), max(abs(clim))]);
% colormap(ap.colormap( ...
%     'PWG'));
% axis image;


%% Create Behavioural Metadata

% Get all trials (all protocols)
widefield_V_all_trials = cat(3, widefield_cat(~cellfun(@isempty, {widefield_cat.rewarded_stim_on_aligned_V})).rewarded_stim_on_aligned_V);

% Initialize with NaNs
widefield_V_start_to_move = nan(size(widefield_V_all_trials));

% Get indices where protocol 1 data exists
protocol1_idx = find(~cellfun(@isempty, {widefield_cat.rewarded_stim_start_to_move_aligned_V}));

% Concatenate only protocol 1 data
protocol1_data = cat(3, widefield_cat(protocol1_idx).rewarded_stim_start_to_move_aligned_V);

% Get trial counts for all sessions up to protocol 1 sessions
trial_counts = cellfun(@(x) size(x,3), {widefield_cat.rewarded_stim_on_aligned_V});
trial_counts(cellfun(@isempty, {widefield_cat.rewarded_stim_on_aligned_V})) = 0;
cumsum_trials = cumsum(trial_counts);

% Map protocol 1 data to correct positions
start_idx = [0 cumsum_trials(protocol1_idx(1:end-1))];

for i = 1:length(protocol1_idx)
    n_trials = size(widefield_cat(protocol1_idx(i)).rewarded_stim_start_to_move_aligned_V, 3);
    widefield_V_start_to_move(:, :, start_idx(i)+(1:n_trials)) = ...
        widefield_cat(protocol1_idx(i)).rewarded_stim_start_to_move_aligned_V;
end

% Get all CS- activity 
% 
% % Initialize with NaNs
% widefield_V_cs_minus_on = nan(size(widefield_V_all_trials));
% 
% % Get indices where protocol 1 data exists
% protocol_1_2_idx = find(~cellfun(@isempty, {widefield_cat.non_rewarded_stim_onset_aligned_V}));
% 
% % Concatenate only protocol 1 data
% protocol_1_2_data = cat(3, widefield_cat(protocol_1_2_idx).non_rewarded_stim_onset_aligned_V);
% 
% % Get trial counts for all sessions up to protocol 1 sessions
% trial_counts = cellfun(@(x) size(x,3), {widefield_cat.rewarded_stim_on_aligned_V});
% trial_counts(cellfun(@isempty, {widefield_cat.rewarded_stim_on_aligned_V})) = 0;
% cumsum_trials = cumsum(trial_counts);
% 
% % Map protocol 1 data to correct positions
% start_idx = [0 cumsum_trials(protocol_1_2_idx(1:end-1))];
% 
% for i = 1:length(protocol_1_2_idx)
%     n_trials = size(widefield_cat(protocol_1_2_idx(i)).non_rewarded_stim_onset_aligned_V, 3);
%     widefield_V_cs_minus_on(:, :, start_idx(i)+(1:n_trials)) = ...
%         widefield_cat(protocol_1_2_idx(i)).non_rewarded_stim_onset_aligned_V;
% end

 
% When working with stim move its easier to populate it within 
% CORRECTED: Build session_indices that match widefield_cat indices

n_total_trials = size(widefield_V_all_trials, 3);
session_indices_corrected = nan(n_total_trials, 1);
animal_ids_corrected = nan(n_total_trials, 1);
workflows_all_trials= nan(n_total_trials, 1);
learning_stages = nan(n_total_trials, 1);
is_nonlearner = false(n_total_trials, 1);

% Behavioural masks
anticipatory_licks = nan(n_total_trials, 1);
cs_labels = nan(n_total_trials, 1);  % CS+ vs CS-
diff_from_optimal = nan(n_total_trials, 1);
latency_from_stim_move= nan(n_total_trials,1);

% Also track if behavioral data exists for each trial
has_behavioral_data = false(n_total_trials, 1);

% To skip workflows work with all workflows
all_workflows = nan(sum(arrayfun(@(x) size(x.rewarded_stim_on_aligned_V, 3), widefield_cat)), 1);
trial_counter = 0;

% Fill workflows
for idx = 1:length(widefield_cat)
    
    day_V_data = widefield_cat(idx).rewarded_stim_on_aligned_V;
    
    if isempty(day_V_data)
        continue;
    end

    n_trials = size(day_V_data, 3);
    trial_range = trial_counter + (1:n_trials);

    all_workflows(trial_range) = workflow_cat(idx);
    trial_counter = trial_counter + n_trials;
end


trial_counter = 0;

for idx = 1:length(widefield_cat)
    day_V_data = widefield_cat(idx).rewarded_stim_on_aligned_V;
    
    if isempty(day_V_data)
        continue;
    end

    n_trials = size(day_V_data, 3);
    trial_range = trial_counter + (1:n_trials);
    

    % Store metadata - idx now correctly refers to widefield_cat position
    session_indices_corrected(trial_range) = idx;
    animal_ids_corrected(trial_range) = widefield_animal_idx(idx);
    workflows_all_trials(trial_range)= workflow_cat(idx);
    learning_stages(trial_range) = learning_index_animal(idx);
    is_nonlearner(trial_range) = is_group_animal(idx);


    % Find corresponding animal and day in behaviour_data
    animal_idx = widefield_animal_idx(idx);
    animal_behaviour = behaviour_data(animal_idx);


    widefield_indices_for_animal = find(widefield_animal_idx == animal_idx);
    day_in_animal = find(widefield_indices_for_animal == idx);

    if isempty(day_in_animal) || day_in_animal > length(animal_behaviour.recording_day)
        fprintf('Warning: Cannot find matching day for animal %d, session %d\n', animal_idx, idx);
        trial_counter = trial_counter + n_trials;
        continue;
    end

    day_data = animal_behaviour.recording_day(day_in_animal);

    
% Get CS labels
if isfield(day_data, 'cs_labels') && ~isempty(day_data.cs_labels)
    
    % Get CS+ trials (rewarded stimulus trials)
    cs_plus_trials = find(day_data.cs_labels == 1);
    
    % Make sure we have the right number of trials (match widefield data)
    n_behavioral_trials = min(length(cs_plus_trials), n_trials);
    
    if n_behavioral_trials > 0
        
        % Store CS labels (should all be 1 for CS+ trials)
        cs_labels(trial_range) = ...
            day_data.cs_labels(cs_plus_trials(1:n_behavioral_trials));
        
        % THEN: Extract anticipatory licks for those CS+ trials
        if isfield(day_data, 'anticipatory_licks') && ~isempty(day_data.anticipatory_licks)

            cs_plus_anticipatory_licks = day_data.anticipatory_licks(day_data.cs_labels == 1);

            % Store in the main array
            anticipatory_licks(trial_range)= ...
                cs_plus_anticipatory_licks(1:n_behavioral_trials);
        else
            fprintf('  Warning: No anticipatory_lick field for session %d\n', idx);
        end

        %Extract diff_from_optimal for those CS+ trials
        if isfield(day_data, 'all_stim_diff_from_optimal_reward') && ~isempty(day_data.all_stim_diff_from_optimal_reward)
            cs_plus_diff_from_optimal = day_data.all_stim_diff_from_optimal_reward(day_data.cs_labels == 1);

            % Store in the main array
            diff_from_optimal(trial_range) = ... 
                cs_plus_diff_from_optimal(1:n_behavioral_trials);
        else
            fprintf('  Warning: No diff_from_optimal field for session %d\n', idx);
        end

        %Extract latency to lick from stim move

        if isfield(day_data, 'cs_plus_latency_lick_to_move') && ~isempty(day_data.cs_plus_latency_lick_to_move)
            cs_plus_latency_lick_from_move= day_data.cs_plus_latency_lick_to_move;

            % Store in the main array
            latency_from_stim_move(trial_range) = ...
                cs_plus_latency_lick_from_move(1:n_behavioral_trials);
        else
            fprintf('  Warning: No latency to lick from stim movement field for session %d\n', idx);
        end


        % Mark these trials as having behavioral data
        has_behavioral_data(trial_range(1:n_behavioral_trials)) = true;
        

        % Optional: Report if there's a mismatch
        if length(cs_plus_trials) ~= n_trials_this_session
            fprintf('  Note: %d CS+ trials in behavior, %d trials in widefield\n', ...
                length(cs_plus_trials), n_trials);
        end

    else
        fprintf('  Warning: No CS+ trials found for session %d\n', idx);
    end
    
else
    fprintf('  Warning: No cs_label field for session %d (Animal %d, Day %d)\n', ...
        idx, animal_idx, day_in_animal);

end

    trial_counter = trial_counter + n_trials;
end




%% Test whether the mapping is correct

% Test that the mapping is now correct
fprintf('=== Verifying Corrected Mapping ===\n');

animal_id = 7;
pre_session_list = find(workflow_cat == 1 & learning_index_animal == 1 & widefield_animal_idx == animal_id);

for sess = pre_session_list'
    % Trials in original
    n_expected = size(widefield_cat(sess).rewarded_stim_on_aligned_V,3);
    
    % Trials assigned in corrected mapping
    n_assigned = sum(session_indices_corrected == sess);
    
    fprintf('Session %d: Expected=%d, Assigned=%d, Match=%d\n', ...
            sess, n_expected, n_assigned, n_expected == n_assigned);
    
    % Compare averaged data
    if n_assigned > 0
        % Average from widefield_V_all_trials using corrected indices
        trials_mask = session_indices_corrected == sess;
        V_reconstructed = mean(widefield_V_all_trials(:, :, trials_mask), 3);
        
        % Compare to rewarded_stim_v_stacked_data (which uses direct widefield_cat indexing)
        V_original = rewarded_stim_v_stacked_data(:, :, sess);
        
        diff = max(abs(V_reconstructed(:) - V_original(:)));
        % fprintf('  Max difference: %.6f %s\n', diff, ...
        %         ternary(diff < 1e-6, '✓ MATCH', '✗ MISMATCH'));
    end
end


% Check anticipatory licks (only for trials with behavioral data)
valid_anticipatory = has_behavioral_data & ~isnan(anticipatory_licks);

% Check diff_from_optimal
valid_diff = has_behavioral_data & ~isnan(diff_from_optimal);

valid_lick_stim_move_latency= has_behavioral_data & ~isnan(latency_from_stim_move);

% ===== Create useful masks for downstream analysis =====
% Group masks (can apply to all stages)
is_learner = ~is_nonlearner;


% % Stages 1+2 only masks
% stage_1_2_valid = stage_1_2_mask & has_behavioral_data;

% Get valid trials with diff_from_optimal data
valid_diff_trials = (all_workflows == 1) & ~isnan(diff_from_optimal);

fprintf('Total valid trials with diff_from_optimal: %d\n', sum(valid_diff_trials));

% ===== Median Split for Differece from Optimal =====

diff_values = diff_from_optimal(valid_diff_trials);

% Calculate median
median_diff = median(diff_values);

% Create masks
fast_trials_median = valid_diff_trials & diff_from_optimal <= median_diff;
slow_trials_median = valid_diff_trials & diff_from_optimal > median_diff;


%% Analysis: mPFC+ Post-Learning by Anticipatory Licks (CORRECTED)

% Define ROI
ROI_to_extract= left_mPFC_ROI_mask;

% Use corrected indices
animal_stage_1_idx = animal_ids_corrected;
session_indices = session_indices_corrected;

% Define valid trials
stage_1_valid = ~isnan(session_indices) & has_behavioral_data & (workflows_all_trials == 1);

% Find RTs for stage 1 for pre-post seperatly so the median has
% correspondence to the learning stage 
stage_1_diff_from_optimal_mask= ~isnan(valid_diff) & (workflows_all_trials == 1) & is_post_learning;

stage_1_diff_from_optimal_values= nan(size(diff_from_optimal)); % initalize first

stage_1_diff_from_optimal_values(stage_1_diff_from_optimal_mask)= diff_from_optimal(stage_1_diff_from_optimal_mask);


% stage_1_fast_trials= stage_1_diff_from_optimal_values<= nanmedian(diff_from_optimal(stage_1_diff_from_optimal_mask));
% stage_1_slow_trials= stage_1_diff_from_optimal_values> nanmedian(diff_from_optimal(stage_1_diff_from_optimal_mask));

stage_1_fast_trials= stage_1_diff_from_optimal_values<= 0.2;
stage_1_slow_trials= stage_1_diff_from_optimal_values> 0.2;

% Learning stage masks
is_pre_learning = learning_stages == 0;
is_post_learning = learning_stages == 1;


% Create masks for anticipatory licks
stage_1_valid_plus_post_no_licks_mask = stage_1_valid & is_learner & is_post_learning & stage_1_slow_trials ;
stage_1_valid_plus_post_with_licks_mask = stage_1_valid & is_learner & is_post_learning & stage_1_fast_trials;

% Get unique animals
animals_post = unique(animal_stage_1_idx(stage_1_valid & is_learner & is_post_learning));
n_animals_post = length(animals_post);

fprintf('=== Analysis Setup (CORRECTED) ===\n');
fprintf('Total valid trials: %d\n', sum(stage_1_valid));
fprintf('Post-learning mPFC+ trials:\n');
fprintf('  No licks: %d\n', sum(stage_1_valid_plus_post_no_licks_mask));
fprintf('  With licks: %d\n', sum(stage_1_valid_plus_post_with_licks_mask));
fprintf('  Animals: %d\n', n_animals_post);

% Compute per-animal averages (session-averaged approach)

% Pre-allocate for means and SEMs
animal_means_post_no_licks = nan(size(widefield_V_all_trials,1), length(added_time), n_animals_post);
animal_means_post_with_licks = nan(size(widefield_V_all_trials,1), length(added_time), n_animals_post);
animal_var_post_no_licks = nan(size(widefield_V_all_trials,1), length(added_time), n_animals_post);
animal_var_post_with_licks = nan(size(widefield_V_all_trials,1), length(added_time), n_animals_post);

for a = 1:n_animals_post
    animal_id = animals_post(a);
    
    % Get sessions for this animal
    animal_sessions = find(widefield_animal_idx == animal_id & workflow_cat == 1 & learning_index_animal == 1);
    
    % For each session, average trials with/without licks separately
    session_avg_no_licks = [];
    session_avg_with_licks = [];
    
    for sess = animal_sessions'
        % Trials from this session
        sess_trials_no_licks = session_indices_corrected == sess & stage_1_slow_trials & has_behavioral_data;
        sess_trials_with_licks = session_indices_corrected == sess & stage_1_fast_trials ~= 0 & has_behavioral_data;
        
        % Average trials within session
        if sum(sess_trials_no_licks) > 0
            sess_V_no_licks = mean(widefield_V_all_trials(:, :, sess_trials_no_licks), 3);
            session_avg_no_licks = cat(3, session_avg_no_licks, sess_V_no_licks);
        end
        
        if sum(sess_trials_with_licks) > 0
            sess_V_with_licks = mean(widefield_V_all_trials(:, :, sess_trials_with_licks), 3);
            session_avg_with_licks = cat(3, session_avg_with_licks, sess_V_with_licks);
        end
    end
    
    % Calculate mean AND SEM across sessions for this animal
    if ~isempty(session_avg_no_licks)
        V_no_licks_mean = mean(session_avg_no_licks, 3);
        V_no_licks_var = var(session_avg_no_licks, 0, 3);
        
        animal_means_post_no_licks(:, :, a) = V_no_licks_mean;
        animal_var_post_no_licks(:, :, a) = V_no_licks_var;
    end
    
    if ~isempty(session_avg_with_licks)
        V_with_licks_mean = mean(session_avg_with_licks, 3);
        V_with_licks_var = var(session_avg_with_licks, 0, 3);
        
        animal_means_post_with_licks(:, :, a) = V_with_licks_mean;
        animal_var_post_with_licks(:, :, a) = V_with_licks_var;
    end
    
    fprintf('Animal %d: %d sessions (no-lick=%d, with-lick=%d)\n', ...
            animal_id, length(animal_sessions), size(session_avg_no_licks, 3), size(session_avg_with_licks, 3));
end


% Compute grand means and VARIANCES (not SEMs yet!)
grand_mean_post_no_licks = mean(animal_means_post_no_licks, 3, 'omitnan');
grand_mean_post_with_licks = mean(animal_means_post_with_licks, 3, 'omitnan');

% Pooled VARIANCE = mean(individual variances)
grand_var_post_no_licks = mean(animal_var_post_no_licks, 3, 'omitnan'); 
grand_var_post_with_licks = mean(animal_var_post_with_licks, 3, 'omitnan');

% Transform MEAN to ROI space (normal)
grand_mean_post_no_licks_roi = ap.wf_roi(wf_U(:,:,1:curr_components), grand_mean_post_no_licks, [], [], ROI_to_extract);
grand_mean_post_with_licks_roi = ap.wf_roi(wf_U(:,:,1:curr_components), grand_mean_post_with_licks, [], [], ROI_to_extract);

% Transform VARIANCE to ROI space (use SQUARED weights)
grand_var_post_no_licks_roi = ap.wf_roi(wf_U(:,:,1:curr_components).^2, grand_var_post_no_licks, [], [], ROI_to_extract);
grand_var_post_with_licks_roi = ap.wf_roi(wf_U(:,:,1:curr_components).^2, grand_var_post_with_licks, [], [], ROI_to_extract);

% SEM across animals
grand_sem_post_no_licks = sqrt(grand_var_post_no_licks_roi) / sqrt(n_animals_post);
grand_sem_post_with_licks= sqrt(grand_var_post_with_licks_roi) / sqrt(n_animals_post);

% % Rename for clarity
grand_mean_post_no_licks = grand_mean_post_no_licks_roi;
grand_mean_post_with_licks = grand_mean_post_with_licks_roi;

% Had a bug with fill function so forced everything into a column vector

% Force everything to column vectors
added_time_col = added_time(:);
mean_no = grand_mean_post_no_licks(:);
sem_no = grand_sem_post_no_licks(:);
mean_with = grand_mean_post_with_licks(:);
sem_with = grand_sem_post_with_licks(:);

% Verify same length
assert(length(added_time_col) == length(mean_no), 'Length mismatch!');
assert(length(mean_no) == length(sem_no), 'SEM length mismatch!');

% Construct SEM bands
upper_no = mean_no + sem_no;
lower_no = mean_no - sem_no;
upper_with = mean_with + sem_with;
lower_with = mean_with - sem_with;

% Plot
figure('Position', [100, 100, 900, 600]); hold on;

ylim_max= max([upper_with;lower_with;upper_no;lower_no],[],'all');
ylim_min= min([upper_with;lower_with;upper_no;lower_no],[],'all');

% Calculate ylims

% No licks
fill([added_time_col; flipud(added_time_col)], ...
     [upper_no; flipud(lower_no)], ...
     cNonLearner*0.85, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time_col, mean_no, 'Color', cNonLearner*0.85, 'LineWidth', 3, ...
     'DisplayName', sprintf('Slow RT > 200ms (n trials=%d)', sum(stage_1_valid_plus_post_no_licks_mask)));

% With licks
fill([added_time_col; flipud(added_time_col)], ...
     [upper_with; flipud(lower_with)], ...
     cNonLearner*0.55, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time_col, mean_with, 'Color', cNonLearner*0.55, 'LineWidth', 3, ...
     'DisplayName', sprintf('Fast RT < 200ms (n trials=%d)', sum(stage_1_valid_plus_post_with_licks_mask)));

xline(0, 'k--', 'LineWidth', 4,'HandleVisibility', 'off');
xlabel('Time from Stim Onset (s)', 'FontSize', 25, 'FontWeight', 'bold');
ylabel('Left ViS Activity (%ΔF/F)', 'FontSize', 25, 'FontWeight', 'bold');
title('mPFC+ Pre-Learning: By RTs', 'FontSize', 20, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 14);
xlim([added_time_col(1), added_time_col(end)]);
ylim([ylim_min,ylim_max]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


%% Latency to Peak Analysis and Amplitude Split by fast vs slow trials 

% Define analysis window
baseline_win = added_time >= -0.1 & added_time < 0;
response_win = added_time >= 0 & added_time <= 0.35;
t_response = added_time(response_win);

% Pre-allocate
peak_latency_no_licks = nan(n_animals_post, 1);
peak_latency_with_licks = nan(n_animals_post, 1);
peak_amp_no_licks = nan(n_animals_post, 1);
peak_amp_with_licks = nan(n_animals_post, 1);

for a = 1:n_animals_post

    % Extract ROI-specific traces FIRST
    roi_trace_no = ap.wf_roi(wf_U(:,:,1:curr_components), ...
        animal_means_post_no_licks(:,:,a), [], [], ROI_to_extract);
    roi_trace_with = ap.wf_roi(wf_U(:,:,1:curr_components), ...
        animal_means_post_with_licks(:,:,a), [], [], ROI_to_extract);


    %Baseline-corrected traces
    baseline_no = mean(roi_trace_no(:, baseline_win), 2);
    baseline_with = mean(roi_trace_with(:, baseline_win), 2);

    base_sub_roi_trace_no= roi_trace_no-baseline_no;
    base_sub_roi_trace_with= roi_trace_with- baseline_with;

    % Now find peak in 1D trace
    [peak_amp_no_licks(a), idx_no] = max(base_sub_roi_trace_no(response_win));
    peak_latency_no_licks(a) = t_response(idx_no);

    [peak_amp_with_licks(a), idx_with] = max(base_sub_roi_trace_with(response_win));
    peak_latency_with_licks(a) = t_response(idx_with);

    peak_amp_no_licks(a) = max(base_sub_roi_trace_no(:));
    peak_amp_with_licks(a) = max(base_sub_roi_trace_with(:));
    
    peak_latency_no_licks(a) = t_response(idx_no);
    peak_latency_with_licks(a) = t_response(idx_with);
end

% Non-parametric paired tests (Wilcoxon signed-rank)
p_latency = signrank(peak_latency_no_licks, peak_latency_with_licks);
p_amp = signrank(peak_amp_no_licks, peak_amp_with_licks);

% Bar plot with paired points
figure;
tiledlayout(1, 2);
set(groot, 'DefaultAxesFontSize', 14);

% Latency comparison
nexttile;
bar_data = [median(peak_latency_no_licks, 'omitnan'), median(peak_latency_with_licks, 'omitnan')];
b = bar(bar_data, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
b.CData = [cNonLearner*0.55; cNonLearner*0.85];
hold on;
% Paired lines
for a = 1:n_animals_post
    plot([1 2], [peak_latency_no_licks(a), peak_latency_with_licks(a)], ...
        'o-', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'none', 'LineWidth', 1);
end
ylabel('Peak Latency (s)');
ylim([0, max([peak_latency_no_licks; peak_latency_with_licks])*1.1]);
xticklabels({'Slow RT', 'Fast RT'});
title(sprintf('p = %.3f', p_latency));

% Amplitude comparison
nexttile;
bar_data = [median(peak_amp_no_licks, 'omitnan'), median(peak_amp_with_licks, 'omitnan')];
b = bar(bar_data, 'FaceColor', 'flat', 'FaceAlpha', 0.6);
b.CData = [cNonLearner*0.55; cNonLearner*0.85];
hold on;
for a = 1:n_animals_post
    plot([1 2], [peak_amp_no_licks(a), peak_amp_with_licks(a)], ...
        'o-', 'Color', [0.3 0.3 0.3], 'MarkerFaceColor', [0.2 0.2 0.2], ...
        'MarkerEdgeColor', 'none', 'LineWidth', 1);
end
ylabel('Peak Amplitude (∆F/F)');
xticklabels({'Slow RT', 'Fast RT'});
title(sprintf('p = %.3f', p_amp));


%% Do jittered swarm plot for RTs vs mPFC Split by Pre and Post Learning 

% Define response window
response_win = added_time >= 0 & added_time <= 0.35;

% Create a mask for all RTs  just for mPFC+
% stage_1_valid_plus = stage_1_valid & is_nonlearner & animal_ids_corrected ~=3 & animal_ids_corrected ~=10;
stage_1_valid_plus = stage_1_valid & is_learner;

stage_1_diff_from_optimal_plus_mask= ~isnan(valid_diff) & stage_1_valid_plus; % just for mPFC+
stage_1_diff_from_optimal_plus_values= nan(size(diff_from_optimal)); % initalize first

stage_1_diff_from_optimal_plus_values(stage_1_diff_from_optimal_plus_mask)= diff_from_optimal(stage_1_diff_from_optimal_plus_mask);


mPFC_activity = squeeze(ap.wf_roi(wf_U(:,:,1:curr_components), ...
    widefield_V_all_trials(:,response_win,stage_1_valid_plus), [], [], ROI_to_extract)); % ROI Tx All trials

% Extract valid trials
animal_ids_valid = animal_ids_corrected(stage_1_valid_plus);
session_ids_valid = session_indices_corrected(stage_1_valid_plus);

% Get RTs for valid trials only
RT_valid = stage_1_diff_from_optimal_plus_values(stage_1_valid_plus);

% Concatenate horizontally to create Nx2
combo_matrix = [animal_ids_valid, session_ids_valid];

% Get unique rows
unique_combos = unique(combo_matrix, 'rows');

mPFC_per_session = nan(size(unique_combos, 1), 1);
RT_per_session = nan(size(unique_combos, 1), 1);


figure;
tiledlayout(1, 2);

% Loop through pre and post learning
learning_stages_combined = {is_pre_learning, is_post_learning};
stage_names = {'Pre-Learning', 'Post-Learning'};

for stage_idx = 1:2
    % Create mask for current stage
    stage_1_valid_plus = stage_1_valid & is_learner & learning_stages_combined{stage_idx};
    
    % Extract mPFC activity
    mPFC_activity = squeeze(ap.wf_roi(wf_U(:,:,1:curr_components), ...
        widefield_V_all_trials(:,response_win,stage_1_valid_plus), [], [], ROI_to_extract));
    
    % Extract valid trials
    animal_ids_valid = animal_ids_corrected(stage_1_valid_plus);
    session_ids_valid = session_indices_corrected(stage_1_valid_plus);
    RT_valid = stage_1_diff_from_optimal_plus_values(stage_1_valid_plus);
    
    % Get unique animal-session combinations
    combo_matrix = [animal_ids_valid, session_ids_valid];
    unique_combos = unique(combo_matrix, 'rows');
    
    mPFC_per_session = nan(size(unique_combos, 1), 1);
    RT_per_session = nan(size(unique_combos, 1), 1);
    
    for i = 1:size(unique_combos, 1)
        animal_id = unique_combos(i, 1);
        session_id = unique_combos(i, 2);
        
        session_mask = animal_ids_valid == animal_id & session_ids_valid == session_id;
        session_data = mPFC_activity(:, session_mask);
        mPFC_per_session(i) = mean(mPFC_activity(:, session_mask), 'all', 'omitnan');
        RT_per_session(i) = median(RT_valid(session_mask), 'omitnan');
    end
    
    % Plot
    nexttile;
    swarmchart(RT_per_session, mPFC_per_session, 50, cLearner, 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    
    % Fit and plot regression line
    valid_idx = ~isnan(RT_per_session) & ~isnan(mPFC_per_session);
    RT_clean = RT_per_session(valid_idx);
    mPFC_clean = mPFC_per_session(valid_idx);
    
    p = polyfit(RT_clean, mPFC_clean, 1);
    RT_fit = linspace(min(RT_clean), max(RT_clean), 100);
    mPFC_fit = polyval(p, RT_fit);
    plot(RT_fit, mPFC_fit, 'k-', 'LineWidth', 2, 'Color',cLearner);
    
    % Statistics
    [r, p_val] = corr(RT_clean, mPFC_clean, 'Type', 'Spearman');
    
    xlabel('Response Time (s)');
    ylabel('mPFC Activity (∆F/F)');
    title(sprintf('%s\nSlope = %.3f, r = %.2f, p = %.3f', stage_names{stage_idx}, p(1), r, p_val));
    set(gca, 'FontSize', 12);
end



%% For CS+ movement analysis : 1. Plots heatmaps aligned to first lick from CS+ move onset  2. mPFC activity aligned to first lick

% Parameters
ROI_to_extract = left_mPFC_ROI_mask;
curr_components = n_components;

% Response window
time_window_start = 0;
time_window_end = 0.35;
 
% Extract ROI trace for ALL trials at once
roi_trace_all_trials = ap.wf_roi(wf_U(:,:,1:curr_components), widefield_V_start_to_move, [], [], ROI_to_extract);
roi_trace_all_trials = squeeze(roi_trace_all_trials);  % time × trials

fprintf('ROI traces extracted: %d timepoints × %d trials\n', size(roi_trace_all_trials));

% Calculate response windows
response_window_idx = added_time >= time_window_start & added_time <= time_window_end;
baseline_window_idx = added_time >= -0.1 & added_time < 0;

% Calculate baseline and response for all trials (vectorized)
baseline = mean(roi_trace_all_trials(baseline_window_idx, :), 1);  % 1 × trials
response = max(roi_trace_all_trials(response_window_idx, :), [], 1);  % 1 × trials

% Baseline-corrected mPFC response
mPFC_response_all_trials = response - baseline;  % 1 × trials

% % Create a mask to get just the valid stage 1 trials and filter mPFC
% stage_1_mask_full = (workflows == 1) & has_behavioral_data;
% valid_mPFC_responses_all_trials= mPFC_response_all_trials(stage_1_mask_full)';
% 

% Stage 1 only
stage_1_valid = (all_workflows == 1) & has_behavioral_data;


% Cut the latencies >0.9s (before reward delivery)
latency_cutoff = 0.9;
valid_latency_mask = latency_from_stim_move < latency_cutoff;

% Pre-learning mPFC+ with valid latency only
stage_1_valid_plus_pre = stage_1_valid & is_learner & is_pre_learning & valid_latency_mask;
mPFC_plus_s1_pre = roi_trace_all_trials(:, stage_1_valid_plus_pre);
latencies_plus_pre = latency_from_stim_move(stage_1_valid_plus_pre);

% Sort by valid latencies only
[sorted_latency_plus_pre, sort_idx_plus_pre] = sort(latencies_plus_pre);
sorted_mPFC_plus_s1_pre = mPFC_plus_s1_pre(:, sort_idx_plus_pre);

% Post-learning mPFC+ with valid latency only
stage_1_valid_plus_post = stage_1_valid & is_learner & is_post_learning & valid_latency_mask;
mPFC_plus_s1_post = roi_trace_all_trials(:, stage_1_valid_plus_post);
latencies_plus_post = latency_from_stim_move(stage_1_valid_plus_post);

[sorted_latency_plus_post, sort_idx_plus_post] = sort(latencies_plus_post);
sorted_mPFC_plus_s1_post = mPFC_plus_s1_post(:, sort_idx_plus_post);

% Pre-learning mPFC- with valid latency only
stage_1_valid_minus_pre = stage_1_valid & ~is_learner & is_pre_learning & valid_latency_mask & animal_stage_1_idx ~= 3 & animal_stage_1_idx ~= 10;
mPFC_minus_s1_pre = roi_trace_all_trials(:, stage_1_valid_minus_pre);
latencies_minus_pre = latency_from_stim_move(stage_1_valid_minus_pre);

[sorted_latency_minus_pre, sort_idx_minus_pre] = sort(latencies_minus_pre);
sorted_mPFC_minus_s1_pre = mPFC_minus_s1_pre(:, sort_idx_minus_pre);

% Post-learning mPFC- with valid latency only
stage_1_valid_minus_post = stage_1_valid & ~is_learner & is_post_learning & valid_latency_mask;
mPFC_minus_s1_post = roi_trace_all_trials(:, stage_1_valid_minus_post);
latencies_minus_post = latency_from_stim_move(stage_1_valid_minus_post);

[sorted_latency_minus_post, sort_idx_minus_post] = sort(latencies_minus_post);
sorted_mPFC_minus_s1_post = mPFC_minus_s1_post(:, sort_idx_minus_post);


% Plot pre-post for mPFC+ group

% Calculate shared color limits across all groups
all_data = [sorted_mPFC_plus_s1_pre(:); sorted_mPFC_plus_s1_post(:); ...
            sorted_mPFC_minus_s1_pre(:); sorted_mPFC_minus_s1_post(:)];
clim_max = max(all_data(~isnan(all_data)));
shared_clim = [0, clim_max];

% Apply 5-trial moving average smoothing
smooth_window = 5;
smoothed_mPFC_plus_pre = movmean(sorted_mPFC_plus_s1_pre, smooth_window, 2, 'omitnan');
smoothed_mPFC_plus_post = movmean(sorted_mPFC_plus_s1_post, smooth_window, 2, 'omitnan');
smoothed_mPFC_minus_pre = movmean(sorted_mPFC_minus_s1_pre, smooth_window, 2, 'omitnan');
smoothed_mPFC_minus_post = movmean(sorted_mPFC_minus_s1_post, smooth_window, 2, 'omitnan');

% Plot mPFC+ group
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
imagesc(added_time, 1:size(smoothed_mPFC_plus_pre, 2), smoothed_mPFC_plus_pre');
colormap(gca, ap.colormap('WG'));
clim(shared_clim);
cb = colorbar;
cb.Label.String = 'ΔF/F';
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trials (sorted by lick latency)', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC+ Pre-Learning', 'FontSize', 14, 'FontWeight', 'bold');
xlim([added_time(1), added_time(end)]);
xline(0,'--')

subplot(1,2,2);
imagesc(added_time, 1:size(smoothed_mPFC_plus_post, 2), smoothed_mPFC_plus_post');
colormap(gca, ap.colormap('WG'));
clim(shared_clim);
cb = colorbar;
cb.Label.String = 'ΔF/F';
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trials (sorted by lick latency)', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC+ Post-Learning', 'FontSize', 14, 'FontWeight', 'bold');
xlim([added_time(1), added_time(end)]);
xline(0,'--')

sgtitle('mPFC+ Group: Movement-Aligned Activity', 'FontSize', 16, 'FontWeight', 'bold');

% Plot mPFC- group
figure('Position', [100, 600, 1200, 500]);

subplot(1,2,1);
imagesc(added_time, 1:size(smoothed_mPFC_minus_pre, 2), smoothed_mPFC_minus_pre');
colormap(gca, ap.colormap('WP'));
clim(shared_clim);
cb = colorbar;
cb.Label.String = 'ΔF/F';
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trials (sorted by lick latency)', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC- Pre-Learning', 'FontSize', 14, 'FontWeight', 'bold');
xlim([added_time(1), added_time(end)]);
xline(0,'--')

subplot(1,2,2);
imagesc(added_time, 1:size(smoothed_mPFC_minus_post, 2), smoothed_mPFC_minus_post');
colormap(gca, ap.colormap('WP'));
clim(shared_clim);
cb = colorbar;
cb.Label.String = 'ΔF/F';
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Trials (sorted by lick latency)', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC- Post-Learning', 'FontSize', 14, 'FontWeight', 'bold');
xlim([added_time(1), added_time(end)]);
xline(0,'--')

sgtitle('mPFC- Group: Movement-Aligned Activity', 'FontSize', 16, 'FontWeight', 'bold');

%% Plot the CS+ movement aligned mPFC activity and lick aligned activity (using interapolation)
 
% Extract animal IDs for each condition (before sorting)
animal_ids_plus_pre = animal_stage_1_idx(stage_1_valid_plus_pre);
animal_ids_plus_post = animal_stage_1_idx(stage_1_valid_plus_post);
animal_ids_minus_pre = animal_stage_1_idx(stage_1_valid_minus_pre);
animal_ids_minus_post = animal_stage_1_idx(stage_1_valid_minus_post);

% Sort animal IDs to match sorted data
animal_ids_plus_pre_sorted = animal_ids_plus_pre(sort_idx_plus_pre);
animal_ids_plus_post_sorted = animal_ids_plus_post(sort_idx_plus_post);
animal_ids_minus_pre_sorted = animal_ids_minus_pre(sort_idx_minus_pre);
animal_ids_minus_post_sorted = animal_ids_minus_post(sort_idx_minus_post);

% Get unique animals for each condition
unique_animals_plus_pre = unique(animal_ids_plus_pre_sorted);
unique_animals_plus_post = unique(animal_ids_plus_post_sorted);
unique_animals_minus_pre = unique(animal_ids_minus_pre_sorted);
unique_animals_minus_post = unique(animal_ids_minus_post_sorted);

n_animals_plus_pre = length(unique_animals_plus_pre);
n_animals_plus_post = length(unique_animals_plus_post);
n_animals_minus_pre = length(unique_animals_minus_pre);
n_animals_minus_post = length(unique_animals_minus_post);

fprintf('=== Animal Counts ===\n');
fprintf('mPFC+ Pre: %d animals, %d trials\n', n_animals_plus_pre, length(latencies_plus_pre));
fprintf('mPFC+ Post: %d animals, %d trials\n', n_animals_plus_post, length(latencies_plus_post));
fprintf('mPFC- Pre: %d animals, %d trials\n', n_animals_minus_pre, length(latencies_minus_pre));
fprintf('mPFC- Post: %d animals, %d trials\n', n_animals_minus_post, length(latencies_minus_post));



% Define a restricted time window that ensures all data is within bounds
max_latency = max([sorted_latency_plus_post(:); sorted_latency_minus_post(:)]);
valid_time_range = (added_time >= -2) & (added_time <= (max(added_time) - max_latency));
restricted_time = added_time(valid_time_range);



% Compute per-animal averages for mPFC+ POST (movement and lick-aligned)

animal_means_move_plus_post = nan(length(added_time), n_animals_plus_post);
animal_means_lick_plus_post = nan(length(restricted_time), n_animals_plus_post);

for a = 1:n_animals_plus_post
    animal_id = unique_animals_plus_post(a);
    animal_mask = animal_ids_plus_post_sorted == animal_id;
    
    % Movement-aligned: average across this animal's trials
    animal_means_move_plus_post(:, a) = mean(sorted_mPFC_plus_s1_post(:, animal_mask), 2, 'omitnan');
    
    % Lick-aligned: realign this animal's trials then average
    animal_trials = sorted_mPFC_plus_s1_post(:, animal_mask);
    animal_latencies = sorted_latency_plus_post(animal_mask);
    
    n_trials_animal = size(animal_trials, 2);
    lick_aligned_animal = nan(length(restricted_time), n_trials_animal);
    
    for t = 1:n_trials_animal
        lick_lat = animal_latencies(t);
        original_time_relick = added_time - lick_lat;
        lick_aligned_animal(:, t) = interp1(original_time_relick, animal_trials(:, t), ...
                                            restricted_time, 'linear', NaN);
    end
    
    animal_means_lick_plus_post(:, a) = mean(lick_aligned_animal, 2, 'omitnan');
end

% Compute per-animal averages for mPFC- POST (movement and lick-aligned)

animal_means_move_minus_post = nan(length(added_time), n_animals_minus_post);
animal_means_lick_minus_post = nan(length(restricted_time), n_animals_minus_post);

for a = 1:n_animals_minus_post
    animal_id = unique_animals_minus_post(a);
    animal_mask = animal_ids_minus_post_sorted == animal_id;
    
    % Movement-aligned: average across this animal's trials
    animal_means_move_minus_post(:, a) = mean(sorted_mPFC_minus_s1_post(:, animal_mask), 2, 'omitnan');
    
    % Lick-aligned: realign this animal's trials then average
    animal_trials = sorted_mPFC_minus_s1_post(:, animal_mask);
    animal_latencies = sorted_latency_minus_post(animal_mask);
    
    n_trials_animal = size(animal_trials, 2);
    lick_aligned_animal = nan(length(restricted_time), n_trials_animal);
    
    for t = 1:n_trials_animal
        lick_lat = animal_latencies(t);
        original_time_relick = added_time - lick_lat;
        lick_aligned_animal(:, t) = interp1(original_time_relick, animal_trials(:, t), ...
                                            restricted_time, 'linear', NaN);
    end
    
    animal_means_lick_minus_post(:, a) = mean(lick_aligned_animal, 2, 'omitnan');
end

% Compute grand means and SEMs across animals (POST-LEARNING)

% Movement-aligned
grand_mean_move_plus = mean(animal_means_move_plus_post, 2, 'omitnan') * 100;
grand_sem_move_plus = std(animal_means_move_plus_post, 0, 2, 'omitnan') / sqrt(n_animals_plus_post) * 100;

grand_mean_move_minus = mean(animal_means_move_minus_post, 2, 'omitnan') * 100;
grand_sem_move_minus = std(animal_means_move_minus_post, 0, 2, 'omitnan') / sqrt(n_animals_minus_post) * 100;

% Lick-aligned
grand_mean_lick_plus = mean(animal_means_lick_plus_post, 2, 'omitnan') * 100;
grand_sem_lick_plus = std(animal_means_lick_plus_post, 0, 2, 'omitnan') / sqrt(n_animals_plus_post) * 100;

grand_mean_lick_minus = mean(animal_means_lick_minus_post, 2, 'omitnan') * 100;
grand_sem_lick_minus = std(animal_means_lick_minus_post, 0, 2, 'omitnan') / sqrt(n_animals_minus_post) * 100;

% Plot with SEM shading (POST-LEARNING)
figure('Position', [100, 100, 1400, 600]);

% Movement-aligned
subplot(1,2,1); hold on;

% mPFC+ with SEM
fill([added_time'; flipud(added_time')], ...
     [grand_mean_move_plus + grand_sem_move_plus; flipud(grand_mean_move_plus - grand_sem_move_plus)], ...
     cLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time, grand_mean_move_plus, 'Color', cLearner, 'LineWidth', 3, ...
     'DisplayName', sprintf('mPFC+ (n=%d)', n_animals_plus_post));

% mPFC- with SEM
fill([added_time'; flipud(added_time')], ...
     [grand_mean_move_minus + grand_sem_move_minus; flipud(grand_mean_move_minus - grand_sem_move_minus)], ...
     cNonLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time, grand_mean_move_minus, 'Color', cNonLearner, 'LineWidth', 3, ...
     'DisplayName', sprintf('mPFC- (n=%d)', n_animals_minus_post));

xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement','HandleVisibility', 'off');
xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('Movement-Aligned', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
xlim([added_time(1), added_time(end)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Lick-aligned
subplot(1,2,2); hold on;

% mPFC+ with SEM
fill([restricted_time'; flipud(restricted_time')], ...
     [grand_mean_lick_plus + grand_sem_lick_plus; flipud(grand_mean_lick_plus - grand_sem_lick_plus)], ...
     cLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(restricted_time, grand_mean_lick_plus, 'Color', cLearner, 'LineWidth', 3, ...
     'DisplayName', sprintf('mPFC+ (n=%d)', n_animals_plus_post));

% mPFC- with SEM
fill([restricted_time'; flipud(restricted_time')], ...
     [grand_mean_lick_minus + grand_sem_lick_minus; flipud(grand_mean_lick_minus - grand_sem_lick_minus)], ...
     cNonLearner, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(restricted_time, grand_mean_lick_minus, 'Color', cNonLearner, 'LineWidth', 3, ...
     'DisplayName', sprintf('mPFC- (n=%d)', n_animals_minus_post));

xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'First Lick','HandleVisibility', 'off');
xlabel('Time from First Lick (s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('Lick-Aligned', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
xlim([restricted_time(1), restricted_time(end)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% sgtitle(sprintf('Post-Learning: Movement vs Lick Alignment (latency < %.1fs, mean ± SEM across animals)', latency_cutoff), ...
%         'FontSize', 16, 'FontWeight', 'bold');


%% If you want quartile seperation


% % PRE-LEARNING: Side-by-side comparison
% % Combine all pre-learning latencies to calculate shared quartiles
% all_latencies_pre = [latencies_plus_pre(:); latencies_minus_pre(:)];
% latency_quantiles_pre = quantile(all_latencies_pre, [0, 0.25, 0.5, 0.75, 1]);
% 
% ylim_pre= [-0.1,0.3];
% 
% figure('Position', [100, 100, 1400, 600]);
% colors = lines(4);
% 
% % mPFC+ Pre-Learning
% subplot(1,2,1); hold on;
% 
% for q = 1:4
%     bin_mask_plus = latencies_plus_pre >= latency_quantiles_pre(q) & ...
%                     latencies_plus_pre < latency_quantiles_pre(q+1);
%     bin_mean_plus = mean(sorted_mPFC_plus_s1_pre(:, bin_mask_plus), 2, 'omitnan') * 100;
% 
%     plot(added_time, bin_mean_plus, '-', 'Color', colors(q,:), 'LineWidth', 3, ...
%          'DisplayName', sprintf('Q%d: %.2f-%.2f s', q, latency_quantiles_pre(q), latency_quantiles_pre(q+1)));
% end
% 
% 
% xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement', 'HandleVisibility', 'off');
% xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
% ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
% title('mPFC+ Pre-Learning', 'FontSize', 14, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 10);
% xlim([added_time(1), added_time(end)]);
% ylim(ylim_pre);
% set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% 
% % mPFC- Pre-Learning
% subplot(1,2,2); hold on;
% 
% for q = 1:4
%     bin_mask_minus = latencies_minus_pre >= latency_quantiles_pre(q) & ...
%                      latencies_minus_pre < latency_quantiles_pre(q+1);
%     bin_mean_minus = mean(sorted_mPFC_minus_s1_pre(:, bin_mask_minus), 2, 'omitnan') * 100;
% 
%     plot(added_time, bin_mean_minus, '-', 'Color', colors(q,:), 'LineWidth', 3, ...
%          'DisplayName', sprintf('Q%d: %.2f-%.2f s', q, latency_quantiles_pre(q), latency_quantiles_pre(q+1)));
% end
% 
% 
% xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement', 'HandleVisibility', 'off');
% xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
% ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
% title('mPFC- Pre-Learning', 'FontSize', 14, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 10);
% xlim([added_time(1), added_time(end)]);
% ylim(ylim_pre);
% set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% 
% sgtitle('Pre-Learning: Shared Quartiles', 'FontSize', 16, 'FontWeight', 'bold');
% 
% % POST-LEARNING: Side-by-side comparison
% % Combine all post-learning latencies to calculate shared quartiles
% 
% all_latencies_post = [latencies_plus_post(:); latencies_minus_post(:)];
% latency_quantiles_post = quantile(all_latencies_post, [0, 0.25, 0.5, 0.75, 1]);
% 
% ylim_post= [-0.1,0.5];
% 
% figure('Position', [100, 100, 1400, 600]);
% 
% % mPFC+ Post-Learning
% subplot(1,2,1); hold on;
% 
% for q = 1:4
%     bin_mask_plus = latencies_plus_post >= latency_quantiles_post(q) & ...
%                     latencies_plus_post < latency_quantiles_post(q+1);
%     bin_mean_plus = mean(sorted_mPFC_plus_s1_post(:, bin_mask_plus), 2, 'omitnan') * 100;
% 
%     plot(added_time, bin_mean_plus, '-', 'Color', colors(q,:), 'LineWidth', 3, ...
%          'DisplayName', sprintf('Q%d: %.2f-%.2f s', q, latency_quantiles_post(q), latency_quantiles_post(q+1)));
% end
% 
% 
% xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement', 'HandleVisibility', 'off');
% xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
% ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
% title('mPFC+ Post-Learning', 'FontSize', 14, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 10);
% xlim([added_time(1), added_time(end)]);
% ylim(ylim_post);
% set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% 
% % mPFC- Post-Learning
% subplot(1,2,2); hold on;
% 
% for q = 1:4
%     bin_mask_minus = latencies_minus_post >= latency_quantiles_post(q) & ...
%                      latencies_minus_post < latency_quantiles_post(q+1);
%     bin_mean_minus = mean(sorted_mPFC_minus_s1_post(:, bin_mask_minus), 2, 'omitnan') * 100;
% 
%     plot(added_time, bin_mean_minus, '-', 'Color', colors(q,:), 'LineWidth', 3, ...
%          'DisplayName', sprintf('Q%d: %.2f-%.2f s', q, latency_quantiles_post(q), latency_quantiles_post(q+1)));
% end
% 
% 
% xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement', 'HandleVisibility', 'off');
% xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
% ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
% title('mPFC- Post-Learning', 'FontSize', 14, 'FontWeight', 'bold');
% legend('Location', 'best', 'FontSize', 10);
% xlim([added_time(1), added_time(end)]);
% ylim(ylim_post);
% set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% 
% sgtitle('Post-Learning: Shared Quartiles', 'FontSize', 16, 'FontWeight', 'bold');
% 

%% Plot ROI traces  -> averaged within animal -> averaged across animals 

% Create masks for mPFC+ group based on anticipatory licks
stage_1_valid_plus_pre_no_licks_mask = stage_1_valid & is_learner & is_pre_learning & anticipatory_licks==0;
stage_1_valid_plus_pre_with_licks_mask = stage_1_valid & is_learner & is_pre_learning & anticipatory_licks~=0;

stage_1_valid_plus_post_no_licks_mask = stage_1_valid & is_learner & is_post_learning & anticipatory_licks==0;
stage_1_valid_plus_post_with_licks_mask = stage_1_valid & is_learner & is_post_learning & anticipatory_licks~=0;

% Extract traces and animal IDs for each condition
% Pre-learning
stage_1_valid_plus_pre_no_licks_traces = roi_trace_all_trials(:, stage_1_valid_plus_pre_no_licks_mask);
animal_ids_plus_pre_no_licks = animal_stage_1_idx(stage_1_valid_plus_pre_no_licks_mask);

stage_1_valid_plus_pre_with_licks_traces = roi_trace_all_trials(:, stage_1_valid_plus_pre_with_licks_mask);
animal_ids_plus_pre_with_licks = animal_stage_1_idx(stage_1_valid_plus_pre_with_licks_mask);

% Post-learning
stage_1_valid_plus_post_no_licks_traces = roi_trace_all_trials(:, stage_1_valid_plus_post_no_licks_mask);
animal_ids_plus_post_no_licks = animal_stage_1_idx(stage_1_valid_plus_post_no_licks_mask);

stage_1_valid_plus_post_with_licks_traces = roi_trace_all_trials(:, stage_1_valid_plus_post_with_licks_mask);
animal_ids_plus_post_with_licks = animal_stage_1_idx(stage_1_valid_plus_post_with_licks_mask);

% Get unique animals for each condition (use same animals across conditions)
unique_animals_plus_pre = unique(animal_ids_plus_pre_sorted);
unique_animals_plus_post = unique(animal_ids_plus_post_sorted);

n_animals_plus_pre = length(unique_animals_plus_pre);
n_animals_plus_post = length(unique_animals_plus_post);

% Compute per-animal averages: PRE-LEARNING

% No licks
animal_means_move_plus_pre_no_licks = nan(length(added_time), n_animals_plus_pre);
for a = 1:n_animals_plus_pre
    animal_id = unique_animals_plus_pre(a);
    animal_mask = animal_ids_plus_pre_no_licks == animal_id;
    
    if sum(animal_mask) > 0  % Only average if animal has trials
        animal_means_move_plus_pre_no_licks(:, a) = mean(stage_1_valid_plus_pre_no_licks_traces(:, animal_mask), 2, 'omitnan');
    end
end

% With licks
animal_means_move_plus_pre_with_licks = nan(length(added_time), n_animals_plus_pre);
for a = 1:n_animals_plus_pre
    animal_id = unique_animals_plus_pre(a);
    animal_mask = animal_ids_plus_pre_with_licks == animal_id;
    
    if sum(animal_mask) > 0  % Only average if animal has trials
        animal_means_move_plus_pre_with_licks(:, a) = mean(stage_1_valid_plus_pre_with_licks_traces(:, animal_mask), 2, 'omitnan');
    end
end

%Compute per-animal averages: POST-LEARNING

% No licks
animal_means_move_plus_post_no_licks = nan(length(added_time), n_animals_plus_post);
for a = 1:n_animals_plus_post
    animal_id = unique_animals_plus_post(a);
    animal_mask = animal_ids_plus_post_no_licks == animal_id;

    if sum(animal_mask) > 0  % Only average if animal has trials
        animal_means_move_plus_post_no_licks(:, a) = mean(stage_1_valid_plus_post_no_licks_traces(:, animal_mask), 2, 'omitnan');
    end
end

% With licks
animal_means_move_plus_post_with_licks = nan(length(added_time), n_animals_plus_post);
for a = 1:n_animals_plus_post
    animal_id = unique_animals_plus_post(a);
    animal_mask = animal_ids_plus_post_with_licks == animal_id;

    if sum(animal_mask) > 0  % Only average if animal has trials
        animal_means_move_plus_post_with_licks(:, a) = mean(stage_1_valid_plus_post_with_licks_traces(:, animal_mask), 2, 'omitnan');
    end
end

% Compute grand means and SEMs across animals

% Pre-learning - CORRECTED
grand_mean_pre_no_licks = mean(animal_means_move_plus_pre_no_licks, 2, 'omitnan');  % NO division
grand_sem_pre_no_licks = std(animal_means_move_plus_pre_no_licks, 0, 2, 'omitnan') / sqrt(n_animals_plus_pre);

grand_mean_pre_with_licks = mean(animal_means_move_plus_pre_with_licks, 2, 'omitnan');  % NO division
grand_sem_pre_with_licks = std(animal_means_move_plus_pre_with_licks, 0, 2, 'omitnan') / sqrt(n_animals_plus_pre);

% Post-learning - already correct
grand_mean_post_no_licks = mean(animal_means_move_plus_post_no_licks, 2, 'omitnan');
grand_sem_post_no_licks = std(animal_means_move_plus_post_no_licks, 0, 2, 'omitnan') / sqrt(n_animals_plus_post);

grand_mean_post_with_licks = mean(animal_means_move_plus_post_with_licks, 2, 'omitnan');
grand_sem_post_with_licks = std(animal_means_move_plus_post_with_licks, 0, 2, 'omitnan') / sqrt(n_animals_plus_post);

joint_data= [grand_mean_pre_no_licks;grand_mean_pre_with_licks;grand_mean_post_no_licks;grand_mean_post_with_licks];
mutual_max_y_lim= max(joint_data,[],'all') *1.2;
mutual_min_y_lim= min(joint_data,[],'all');

% Plot: Pre-learning comparison
figure('Position', [100, 100, 900, 600]); hold on;

% No licks
fill([added_time'; flipud(added_time')], ...
     [grand_mean_pre_no_licks + grand_sem_pre_no_licks; flipud(grand_mean_pre_no_licks - grand_sem_pre_no_licks)], ...
     cLearner*0.3, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time', grand_mean_pre_no_licks, 'Color', cLearner*0.3, 'LineWidth', 3, ...
     'DisplayName', sprintf('No Anticipatory Licks (n=%d)', sum(~isnan(animal_means_move_plus_pre_no_licks(1,:)))));

% With licks
fill([added_time'; flipud(added_time')], ...
     [grand_mean_pre_with_licks + grand_sem_pre_with_licks; flipud(grand_mean_pre_with_licks - grand_sem_pre_with_licks)], ...
      cLearner*0.8, 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time', grand_mean_pre_with_licks, 'Color', cLearner*0.8, 'LineWidth', 3, ...
     'DisplayName', sprintf('With Anticipatory Licks (n=%d)', sum(~isnan(animal_means_move_plus_pre_with_licks(1,:)))));

xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement');
xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('mPFC+ Pre-Learning: By Anticipatory Licks', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
ylim([mutual_min_y_lim,mutual_max_y_lim]);
xlim([added_time(1), added_time(end)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Force MATLAB to recompute polygon by explicit vector construction - this
% is done give no lick SEM is an order of magnitude lower than the mean
% which creates a fragmented SEM. This is due to how fill function can have
% rendering problems

figure('Position', [100, 100, 900, 600]); hold on;

% No licks - construct SEM band explicitly
upper_no_licks = grand_mean_post_no_licks + grand_sem_post_no_licks;
lower_no_licks = grand_mean_post_no_licks - grand_sem_post_no_licks;

fill([added_time'; flipud(added_time')], ...
     [upper_no_licks; flipud(lower_no_licks)], ...
     [0.3 0.3 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time, grand_mean_post_no_licks, 'Color', [0.3 0.3 0.8], 'LineWidth', 3, ...
     'DisplayName', sprintf('No Anticipatory Licks (n=%d)', n_filtered));

% With licks - construct SEM band explicitly
upper_with_licks = grand_mean_post_with_licks + grand_sem_post_with_licks;
lower_with_licks = grand_mean_post_with_licks - grand_sem_post_with_licks;

fill([added_time'; flipud(added_time')], ...
     [upper_with_licks; flipud(lower_with_licks)], ...
     [0.8 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(added_time, grand_mean_post_with_licks, 'Color', [0.8 0.3 0.3], 'LineWidth', 3, ...
     'DisplayName', sprintf('With Anticipatory Licks (n=%d)', n_filtered));

xline(0, 'k--', 'LineWidth', 2.5, 'Label', 'Movement');
xlabel('Time from Movement (s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (%ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('mPFC+ Post-Learning: By Anticipatory Licks', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
ylim([mutual_min_y_lim,mutual_max_y_lim]);
xlim([added_time(1), added_time(end)]);
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


% Plot individual animals to see the source of the difference
figure('Position', [100, 100, 1400, 800]);

for a = 1:n_filtered
    subplot(2, 3, a); hold on;
    
    animal_id = filtered_animals(a);
    
    % Plot this animal's traces
    plot(added_time, animal_means_move_plus_pre_no_licks(:, a) * 100, 'b', 'LineWidth', 2, ...
         'DisplayName', 'No Licks');
    plot(added_time, animal_means_move_plus_pre_with_licks(:, a) * 100, 'r', 'LineWidth', 2, ...
         'DisplayName', 'With Licks');
    
    % Report trial counts
    n_no_lick = sum(animal_ids_plus_post_no_licks == animal_id);
    n_with_lick = sum(animal_ids_plus_post_with_licks == animal_id);
    
    xline(0, 'k--');
    xlabel('Time (s)');
    ylabel('Activity (%ΔF/F)');
    title(sprintf('Animal %d\nNo-lick: %d, With-lick: %d', animal_id, n_no_lick, n_with_lick));
    legend('Location', 'best', 'FontSize', 9);
    xlim([added_time(1), added_time(end)]);
end

sgtitle('Individual Animal Traces: mPFC+ Post-Learning', 'FontSize', 14, 'FontWeight', 'bold');


% Plot pre vs post for the SAME animals to see the difference clearly
figure('Position', [100, 100, 1400, 600]);

for a = 1:n_filtered  % First 3 animals
    subplot(1, n_filtered, a); hold on;
    
    % Pre-learning (both conditions averaged)
    pre_avg = mean([animal_means_move_plus_pre_no_licks(:, a), ...
                    animal_means_move_plus_pre_with_licks(:, a)], 2, 'omitnan') * 100;
    
    % Post-learning (both conditions averaged)
    post_avg = mean([animal_means_move_plus_post_no_licks(:, a), ...
                     animal_means_move_plus_post_with_licks(:, a)], 2, 'omitnan') * 100;
    
    plot(added_time, pre_avg, 'g', 'LineWidth', 3, 'DisplayName', 'Pre-learning');
    plot(added_time, post_avg, 'b', 'LineWidth', 3, 'DisplayName', 'Post-learning');
    
    xline(0, 'k--');
    xlabel('Time (s)');
    ylabel('Activity (%ΔF/F)');
    title(sprintf('Animal %d', filtered_animals(a)));
    legend('Location', 'best');
    xlim([added_time(1), added_time(end)]);
end

sgtitle('Pre vs Post Learning Comparison', 'FontSize', 14, 'FontWeight', 'bold');

%% Fixing a potential bug 


% Compare the two averaging methods directly for one animal
animal_id = filtered_animals(3);  % Test with first animal

fprintf('=== Comparing Averaging Methods for Animal %d ===\n', animal_id);

peak_window = added_time >= 0 & added_time <= 0.35;

% Method 1: Current (all trials pooled, then averaged)
pre_trials_mask = is_pre_learning & animal_stage_1_idx == animal_id & stage_1_valid;
post_trials_mask = is_post_learning & animal_stage_1_idx == animal_id & stage_1_valid;

pre_activity_pooled = mean(widefield_V_all_trials(:,:, pre_trials_mask), 3);
post_activity_pooled = mean(widefield_V_all_trials(:,:, post_trials_mask), 3);

pre_peak_pooled = max(pre_activity_pooled(peak_window));
post_peak_pooled = max(post_activity_pooled(peak_window));

fprintf('Method 1 (pooled trials): Pre=%.4f, Post=%.4f\n', pre_peak_pooled, post_peak_pooled);

% Method 2: Day-averaged (average within day, then across days)
animal_sessions = find(widefield_animal_idx == animal_id & workflow_cat == 1);

pre_sessions = animal_sessions(learning_index_animal(animal_sessions) == 0);
post_sessions = animal_sessions(learning_index_animal(animal_sessions) == 1);

% Average within each pre-learning day
pre_day_means = nan(size(widefield_V_all_trials,1),length(added_time), length(pre_sessions));

for i = 1:length(pre_sessions)
    sess = pre_sessions(i);
    sess_trials = session_indices == sess;
    if sum(sess_trials) > 0
        pre_day_means(:,:, i) = mean(widefield_V_all_trials(:, :,sess_trials), 3);
    end
end

% Average within each post-learning day
post_day_means = nan(size(widefield_V_all_trials,1),length(added_time), length(post_sessions));
for i = 1:length(post_sessions)
    sess = post_sessions(i);
    sess_trials = session_indices == sess;
    if sum(sess_trials) > 0
        post_day_means(:, :,i) = mean(widefield_V_all_trials(:, :,sess_trials), 3);
    end
end

% Now average across days
pre_activity_day_avg = mean(pre_day_means, 3, 'omitnan');
post_activity_day_avg = mean(post_day_means, 3, 'omitnan');

pre_peak_day_avg = max(pre_activity_day_avg(peak_window));
post_peak_day_avg = max(post_activity_day_avg(peak_window));

fprintf('Method 2 (day-averaged): Pre=%.4f, Post=%.4f\n', pre_peak_day_avg, post_peak_day_avg);

% Compare
fprintf('\nDifference between methods:\n');
fprintf('  Pre-learning: %.4f (%.1f%% difference)\n', ...
        pre_peak_pooled - pre_peak_day_avg, ...
        100 * (pre_peak_pooled - pre_peak_day_avg) / pre_peak_day_avg);
fprintf('  Post-learning: %.4f (%.1f%% difference)\n', ...
        post_peak_pooled - post_peak_day_avg, ...
        100 * (post_peak_pooled - post_peak_day_avg) / post_peak_day_avg);

ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:n_components),pre_activity_day_avg),plab.wf.svd2px(wf_U(:,:,1:n_components),post_activity_day_avg)],added_time);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;


%% Fixing bug 2

animal_id = 7;

% Method 1: Your visualization (direct from stacked data)
viz_pre_mask = workflow_cat == 1 & learning_index_animal == 0 & widefield_animal_idx == animal_id;
viz_post_mask = workflow_cat == 1 & learning_index_animal == 1 & widefield_animal_idx == animal_id;

pre_V_viz = nanmean(rewarded_stim_v_stacked_data(:, :, viz_pre_mask), 3);
post_V_viz = nanmean(rewarded_stim_v_stacked_data(:, :, viz_post_mask), 3);


% Method 2: Your day-averaging (from widefield_V_all_trials)
animal_sessions = find(widefield_animal_idx == animal_id & workflow_cat == 1);
pre_sessions = animal_sessions(learning_index_animal(animal_sessions) == 0);
post_sessions = animal_sessions(learning_index_animal(animal_sessions) == 1);

% WITHOUT stage_1_valid filter
pre_day_means = nan(size(widefield_V_all_trials,1), size(widefield_V_all_trials,2), length(pre_sessions));
for i = 1:length(pre_sessions)
    sess = pre_sessions(i);
    sess_trials = session_indices == sess;  % No stage_1_valid
    if sum(sess_trials) > 0
        pre_day_means(:, :, i) = mean(widefield_V_all_trials(:, :, sess_trials), 3);
    end
end

post_day_means = nan(size(widefield_V_all_trials,1), size(widefield_V_all_trials,2), length(post_sessions));
for i = 1:length(post_sessions)
    sess = post_sessions(i);
    sess_trials = session_indices == sess;  % No stage_1_valid
    if sum(sess_trials) > 0
        post_day_means(:, :, i) = mean(widefield_V_all_trials(:, :, sess_trials), 3);
    end
end

pre_V_dayavg = mean(pre_day_means, 3, 'omitnan');
post_V_dayavg = mean(post_day_means, 3, 'omitnan');


ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:n_components),post_V_dayavg),plab.wf.svd2px(wf_U(:,:,1:n_components),post_V_viz)],added_time);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;



% Manually reconstruct day averages and compare
animal_id = 7;

% Method 1: Direct from rewarded_stim_v_stacked_data (your visualization)
pre_mask_viz = workflow_cat == 1 & learning_index_animal == 0 & widefield_animal_idx == animal_id;
pre_V_viz = mean(rewarded_stim_v_stacked_data(:, :, pre_mask_viz), 3);

fprintf('=== Reconstruction Test ===\n');
fprintf('Visualization uses %d pre-learning sessions\n', sum(pre_mask_viz));
fprintf('Session indices: %s\n', mat2str(find(pre_mask_viz)'));

% Method 2: Manually average trials for each of those sessions
pre_session_list = find(pre_mask_viz);
pre_V_manual = nan(size(rewarded_stim_v_stacked_data, 1), size(rewarded_stim_v_stacked_data, 2), length(pre_session_list));

for i = 1:length(pre_session_list)
    sess = pre_session_list(i);
    
    % Get trials for this session from widefield_V_all_trials
    sess_trials_mask = session_indices == sess;
    
    fprintf('\nSession %d:\n', sess);
    fprintf('  Trials in widefield_V_all_trials: %d\n', sum(sess_trials_mask));
    
    if sum(sess_trials_mask) > 0
        % Average those trials
        pre_V_manual(:, :, i) = mean(widefield_V_all_trials(:, :, sess_trials_mask), 3);
        
        % Compare to pre-averaged version
        diff = max(abs(pre_V_manual(:, :, i) - rewarded_stim_v_stacked_data(:, :, sess)), [], 'all');
        fprintf('  Max difference from rewarded_stim_v_stacked_data: %.6f\n', diff);
        
        if diff > 1e-6
            fprintf('  *** MISMATCH! These should be identical! ***\n');
        end
    else
        fprintf('  *** NO TRIALS FOUND in widefield_V_all_trials! ***\n');
    end
end

% Average across sessions
pre_V_manual_avg = mean(pre_V_manual, 3, 'omitnan');

% Compare final averages
fprintf('\n=== Final Comparison ===\n');
fprintf('Max difference between methods: %.6f\n', max(abs(pre_V_viz(:) - pre_V_manual_avg(:))));


% For each session, check if ALL trials are in widefield_V_all_trials
fprintf('=== Checking Trial Inclusion ===\n');

total_trials_original = 0;
total_trials_in_V_all = 0;

for idx = 1:length(widefield_cat)
    if workflow_cat(idx) ~= 1
        continue;  % Skip non-stage-1
    end
    
    % Trials in original data
    if ~isempty(widefield_cat(idx).rewarded_stim_on_aligned_V)
        n_orig = size(widefield_cat(idx).rewarded_stim_on_aligned_V, 3);
    else
        n_orig = 0;
    end
    
    % Trials assigned to this session in widefield_V_all_trials
    n_in_V = sum(session_indices == idx);
    
    total_trials_original = total_trials_original + n_orig;
    total_trials_in_V_all = total_trials_in_V_all + n_in_V;
    
    if n_orig ~= n_in_V
        fprintf('Session %d: Original=%d, In V_all=%d, Missing=%d\n', ...
                idx, n_orig, n_in_V, n_orig - n_in_V);
    end
end

fprintf('\n=== TOTALS ===\n');
fprintf('Total trials in widefield_cat: %d\n', total_trials_original);
fprintf('Total trials in widefield_V_all_trials: %d\n', total_trials_in_V_all);
fprintf('Difference: %d trials (%.1f%%)\n', ...
        total_trials_original - total_trials_in_V_all, ...
        100 * (total_trials_original - total_trials_in_V_all) / total_trials_original);


%% Plot ROI Traces Time aligned first averaging within animal and then across with SEMs


% Create masks for mPFC+ group
stage_1_valid_plus_pre_no_licks_mask = stage_1_valid & is_learner & is_pre_learning & anticipatory_licks ==0;
stage_1_valid_plus_pre_with_licks_mask = stage_1_valid & is_learner & is_pre_learning & anticipatory_licks ~=0;
stage_1_valid_plus_post_no_licks_mask = stage_1_valid & is_learner & is_post_learning & anticipatory_licks ==0;
stage_1_valid_plus_post_with_licks_mask = stage_1_valid & is_learner & is_post_learning & anticipatory_licks ~=0;

% Corresponding traces
stage_1_valid_plus_pre_no_licks_traces= roi_trace_all_trials(:,stage_1_valid_plus_pre_no_licks_mask);
stage_1_valid_plus_pre_with_licks_traces= roi_trace_all_trials(:,stage_1_valid_plus_pre_with_licks_mask);
stage_1_valid_plus_post_no_licks_traces= roi_trace_all_trials(:,stage_1_valid_plus_post_no_licks_mask);
stage_1_valid_plus_post_with_licks_traces= roi_trace_all_trials(:,stage_1_valid_plus_post_with_licks_mask);


animal_means_move_plus_pre_no_licks = nan(length(added_time), n_animals_plus_pre);
animal_means_move_plus_pre_with_licks = nan(length(added_time), n_animals_plus_pre);

for a = 1:n_animals_plus_pre
    animal_id = unique_animals_plus_post(a);
    animal_mask = animal_ids_plus_post_sorted == animal_id;
    
    % Movement-aligned: average across this animal's trials
    animal_means_move_plus_pre_no_licks(:, a) = mean(sorted_mPFC_plus_s1_post(:, animal_mask), 2, 'omitnan');
    
    % Lick-aligned: realign this animal's trials then average
    animal_trials = sorted_mPFC_plus_s1_post(:, animal_mask);
    animal_latencies = sorted_latency_plus_post(animal_mask);
    
    n_trials_animal = size(animal_trials, 2);
    lick_aligned_animal = nan(length(restricted_time), n_trials_animal);
    
    for t = 1:n_trials_animal
        lick_lat = animal_latencies(t);
        original_time_relick = added_time - lick_lat;
        lick_aligned_animal(:, t) = interp1(original_time_relick, animal_trials(:, t), ...
                                            restricted_time, 'linear', NaN);
    end
    
    animal_means_lick_plus_post(:, a) = mean(lick_aligned_animal, 2, 'omitnan');
end



% Helper for computing mean and SEM
computeStats = @(data) deal( ...
    mean(data, 2), ... % mean over trials
    std(data, 0, 2) ./ sqrt(sum(~isnan(data), 2)) ... % SEM
);

[m_pre_with, sem_pre_with] = computeStats(mPFC_plus_s1_pre_with_lick);
[m_pre_no, sem_pre_no] = computeStats(mPFC_plus_s1_pre_no_lick);

[m_post_with, sem_post_with] = computeStats(mPFC_plus_s1_post_with_lick);
[m_post_no, sem_post_no] = computeStats(mPFC_plus_s1_post_no_lick);


% First pass: determine global y-axis limits across both groups
all_means_s1 = [m_pre_with, m_pre_no,m_post_with,m_post_no];
all_sems_s1 = [sem_pre_with, sem_pre_no,sem_post_with,sem_post_no];

% y_min_global = min(all_means_s1 + all_sems_s1,[],'all');
% y_max_global = max(all_means_s1 + all_sems_s1,[],'all'); 



% Plot
figure;

% Color definitions
c_with = cLearner * 0.4;  % dark
c_no   = cLearner * 0.8;  % light

% Pre-learning panel
subplot(1,2,1)
hold on


fill([added_time, fliplr(added_time)], ...
     [m_pre_with'+sem_pre_with', fliplr(m_pre_with'-sem_pre_with')], ...
     c_with, 'FaceAlpha',0.25, 'EdgeColor','none')
fill([added_time, fliplr(added_time)], ...
     [m_pre_no'+sem_pre_no', fliplr(m_pre_no'-sem_pre_no')], ...
     c_no, 'FaceAlpha',0.25, 'EdgeColor','none')

plot(added_time, m_pre_with', 'Color', c_with, 'LineWidth', 2)
plot(added_time, m_pre_no', 'Color', c_no, 'LineWidth', 2)

title('Pre-learning')
xlabel('Time (s)')
ylabel('ΔF/F')
ylim([y_min_global y_max_global]);
legend({'With lick','No lick'}, 'Location','best')
xline(0,'--k',HandleVisibility='off')

% Post-learning panel
subplot(1,2,2)
hold on



fill([added_time, fliplr(added_time)], ...
     [m_post_with'+sem_post_with', fliplr(m_post_with'-sem_post_with')], ...
     c_with, 'FaceAlpha',0.25, 'EdgeColor','none')
fill([added_time, fliplr(added_time)], ...
     [m_post_no'+sem_post_no', fliplr(m_post_no'-sem_post_no')], ...
     c_no, 'FaceAlpha',0.25, 'EdgeColor','none')

plot(added_time, m_post_with', 'Color', c_with, 'LineWidth', 2)
plot(added_time, m_post_no', 'Color', c_no, 'LineWidth', 2)

title('Post-learning')
xlabel('Time (s)')
ylim([y_min_global y_max_global]);
legend({'With lick','No lick'}, 'Location','best')
xline(0,'--k',HandleVisibility='off')

sgtitle('mPFC+ Stage 1: Left mPFC activity Seperated Trials with and without Anticipatory Licks')



%% ===== Extract mPFC activity for all trials (vectorized) =====

% widefield_V_all_trials = cat(3, widefield_cat(~cellfun(@isempty, {widefield_cat.rewarded_stim_on_aligned_V})).rewarded_stim_on_aligned_V); % get all the data VxTx n all trials 

widefield_V_all_trials = cat(3, widefield_cat(~cellfun(@isempty, {widefield_cat.rewarded_stim_start_to_move_aligned_V})).rewarded_stim_start_to_move_aligned_V); % get all the data VxTx n all trials 


% Parameters
ROI_to_extract = left_mPFC_ROI_mask;
curr_components = n_components;

% Response window
time_window_start = 0;
time_window_end = 0.35;
 
% Extract ROI trace for ALL trials at once
roi_trace_all_trials = ap.wf_roi(wf_U(:,:,1:curr_components), widefield_V_all_trials, [], [], ROI_to_extract);
roi_trace_all_trials = squeeze(roi_trace_all_trials);  % time × trials

fprintf('ROI traces extracted: %d timepoints × %d trials\n', size(roi_trace_all_trials));

% Calculate response windows
response_window_idx = added_time >= time_window_start & added_time <= time_window_end;
baseline_window_idx = added_time >= -0.1 & added_time < 0;

% Calculate baseline and response for all trials (vectorized)
baseline = mean(roi_trace_all_trials(baseline_window_idx, :), 1);  % 1 × trials
response = max(roi_trace_all_trials(response_window_idx, :), [], 1);  % 1 × trials

% Baseline-corrected mPFC response
mPFC_response_all_trials = response - baseline;  % 1 × trials

% Create a mask to get just the valid stage 1 trials and filter mPFC
stage_1_mask_full = (workflows == 1) & has_behavioral_data;
valid_mPFC_responses_all_trials= mPFC_response_all_trials(stage_1_mask_full)';


% Stage 1 only
stage_1_valid = (workflows == 1) & has_behavioral_data & valid_mPFC_responses_all_trials;

% % Stage 2 only
% stage_2_valid = (workflows == 2) & has_behavioral_data & ~isnan(mPFC_response_all_trials');

fprintf('\n===== Stage-specific Analysis =====\n');
fprintf('Stage 1 valid trials: %d\n', sum(stage_1_valid));
fprintf('Stage 2 valid trials: %d\n', sum(stage_2_valid));

% ===== STAGE 1 ANALYSIS =====

% mPFC+ (Learners) - Stage 1
mPFC_plus_s1_pre_with_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_pre_learning); % 
mPFC_plus_s1_pre_no_lick = mPFC_response_all_trials(stage_1_valid & is_learner & is_pre_learning & slow_trials_median);
mPFC_plus_s1_post_with_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_post_learning);
mPFC_plus_s1_post_no_lick = mPFC_response_all_trials(stage_1_valid & is_learner & is_post_learning & slow_trials_median);


% mPFC- (Non-learners) - Stage 1 also excludes HA006 and HA013 given they
% haven't learned the task 

mPFC_minus_s1_pre_with_lick = roi_trace_all_trials(:,stage_1_valid & is_nonlearner & is_pre_learning);
mPFC_minus_s1_pre_no_lick = roi_trace_all_trials(stage_1_valid & is_nonlearner & is_pre_learning );

mPFC_minus_s1_post_with_lick = roi_trace_all_trials(:,stage_1_valid & is_nonlearner & is_post_learning & anticipatory_licks ~=0);
mPFC_minus_s1_post_no_lick = roi_trace_all_trials(:,stage_1_valid & is_nonlearner & is_post_learning & anticipatory_licks ==0);

% Sort based on latency
[~, mPFC_plus_s1_pre_sorted_idx]= sort(latency_from_stim_move(stage_1_valid & is_learner & is_pre_learning));
sorted_mPFC_plus_s1_pre = mPFC_plus_s1_pre_with_lick(:,mPFC_plus_s1_pre_sorted_idx);

[sorted_latency_mPFC_plus_s1_post, mPFC_plus_s1_post_sorted_idx]= sort(latency_from_stim_move(stage_1_valid & is_learner & is_post_learning));
sorted_mPFC_plus_s1_post = mPFC_plus_s1_post_with_lick(:,mPFC_plus_s1_post_sorted_idx);


% Sort based on latency
[sorted_latency, mPFC_minus_s1_pre_sorted_idx]= sort(latency_from_stim_move(stage_1_valid & is_nonlearner & is_pre_learning));
sorted_mPFC_minus_s1_pre = mPFC_minus_s1_pre_with_lick(:,mPFC_minus_s1_pre_sorted_idx);

% Split based on lick or no lick
[sorted_latency_mPFC_minus_s1_post_no_lick, mPFC_minus_s1_post_sorted_idx_no_lick,]= sort(latency_from_stim_move(stage_1_valid & is_nonlearner & is_post_learning & anticipatory_licks ==0));
sorted_mPFC_minus_s1_post_no_lick = mPFC_minus_s1_post_no_lick(:,mPFC_minus_s1_post_sorted_idx_no_lick);

[sorted_latency_mPFC_minus_s1_post_with_lick, mPFC_minus_s1_post_sorted_idx_with_lick,]= sort(latency_from_stim_move(stage_1_valid & is_nonlearner & is_post_learning & anticipatory_licks ~=0));
sorted_mPFC_minus_s1_post_with_lick = mPFC_minus_s1_post_with_lick(:,mPFC_minus_s1_post_sorted_idx_with_lick);


% Bin by latency quartiles
latency_quantiles = quantile(sorted_latency_mPFC_minus_s1_post_with_lick, [0, 0.25, 0.5, 0.75, 1]);

figure; hold on;
colors = lines(4);

for q = 1:4
    bin_mask = sorted_latency_mPFC_minus_s1_post_with_lick >= latency_quantiles(q) & ...
               sorted_latency_mPFC_minus_s1_post_with_lick < latency_quantiles(q+1);
    
    bin_mean = mean(sorted_mPFC_minus_s1_post_with_lick(:,bin_mask), 2, 'omitnan');
    bin_sem = std(sorted_mPFC_minus_s1_post_with_lick(:,bin_mask), 0, 2, 'omitnan') / sqrt(sum(bin_mask));
    
    plot(added_time, bin_mean*100, 'Color', colors(q,:), 'LineWidth', 2.5, ...
         'DisplayName', sprintf('Q%d: %.2f-%.2f s', q, latency_quantiles(q), latency_quantiles(q+1)));
end

xline(0, 'k--', 'LineWidth', 2, 'Label', 'Movement');
xlabel('Time from Movement (s)');
ylabel('mPFC Activity (ΔF/F)');
title('mPFC Activity by Lick Latency Quartile');
legend('Location', 'best');



%%% Re align the widefield to first lick 

% For each group, realign to first lick
lick_aligned_learner = nan(length(added_time),size(sorted_mPFC_plus_s1_post,2));
lick_aligned_nonlearner = nan(length(added_time),size(sorted_mPFC_minus_s1_post,2));

% new_time = -2:0.01:2;  % Time relative to lick

% Realign learners post learning
for i = 1:length(sorted_latency_mPFC_plus_s1_post)
    lick_lat = sorted_latency_mPFC_plus_s1_post(i);
    
    if isnan(lick_lat)
        continue;
    end

    original_time_relick = added_time - lick_lat;
    lick_aligned_learner(:, i) = interp1(original_time_relick, sorted_mPFC_plus_s1_post(:,i), ...
                                         added_time, 'linear', NaN);
end

% Realign non-learners

for i = 1:length(sorted_latency_mPFC_minus_s1_post)
  
    lick_lat = sorted_latency_mPFC_minus_s1_post(i);

    if isnan(lick_lat)
        continue;
    end

    original_time_relick = added_time - lick_lat;
    lick_aligned_nonlearner(:,i) = interp1(original_time_relick, sorted_mPFC_minus_s1_post(:,i), ...
                                            added_time, 'linear', NaN);
end

% Plot comparison
figure('Position', [100, 100, 1400, 600]);

% Movement-aligned
subplot(1,2,1); hold on;
mean_move_learner = mean(sorted_mPFC_plus_s1_post, 2, 'omitnan');
mean_move_nonlearner = mean(sorted_mPFC_minus_s1_post, 2, 'omitnan');

plot(added_time, mean_move_learner, 'Color', cLearner, 'LineWidth', 3, 'DisplayName', 'mPFC+');
plot(added_time, mean_move_nonlearner, 'Color', cNonLearner, 'LineWidth', 3, 'DisplayName', 'mPFC-');
xline(0, 'k--', 'LineWidth', 2, 'Label', 'Movement');
xlabel('Time from Movement (s)');
ylabel('mPFC Activity (ΔF/F)');
title('Movement-Aligned');
legend('Location', 'best');

% Lick-aligned
subplot(1,2,2); hold on;
mean_lick_learner = mean(lick_aligned_learner, 2, 'omitnan');
mean_lick_nonlearner = mean(lick_aligned_nonlearner, 2, 'omitnan');

plot(added_time, mean_lick_learner, 'Color', cLearner, 'LineWidth', 3, 'DisplayName', 'mPFC+');
plot(added_time, mean_lick_nonlearner, 'Color', cNonLearner, 'LineWidth', 3, 'DisplayName', 'mPFC-');
xline(0, 'k--', 'LineWidth', 2, 'Label', 'First Lick');
xlabel('Time from First Lick (s)');
ylabel('mPFC Activity (ΔF/F)');
title('Lick-Aligned');
legend('Location', 'best');

sgtitle('Movement vs Lick Alignment by Group');


% Plot as heatmap
figure;
imagesc(added_time,1:size(sorted_mPFC_plus_s1_post,2),sorted_mPFC_plus_s1_post');
colormap(ap.colormap('WG'));
colorbar;
clim
xlabel('Time from Movement (s)');
ylabel('Trials (sorted by lick latency)');
title('mPFC Activity Sorted by First Lick Latency');






% Statistical tests - Stage 1
[~, p_plus_s1_pre] = ttest2(mPFC_plus_s1_pre_with_lick, mPFC_plus_s1_pre_no_lick);
[~, p_plus_s1_post] = ttest2(mPFC_plus_s1_post_with_lick, mPFC_plus_s1_post_no_lick);
[~, p_minus_s1_pre] = ttest2(mPFC_minus_s1_pre_with_lick, mPFC_minus_s1_pre_no_lick);
[~, p_minus_s1_post] = ttest2(mPFC_minus_s1_post_with_lick, mPFC_minus_s1_post_no_lick);

fprintf('Stage 1 - mPFC+: Pre p=%.4f, Post p=%.4f\n', p_plus_s1_pre, p_plus_s1_post);
fprintf('Stage 1 - mPFC-: Pre p=%.4f, Post p=%.4f\n', p_minus_s1_pre, p_minus_s1_post);

% ===== STAGE 2 ANALYSIS =====

fprintf('\n===== Stage 2: Extracting data by group and learning stage =====\n');

% mPFC+ (Learners) - Stage 2
mPFC_plus_s2_pre_with_lick = mPFC_response_all_trials(stage_2_valid & is_learner & is_pre_learning & anticipatory_licks == 1);
mPFC_plus_s2_pre_no_lick = mPFC_response_all_trials(stage_2_valid & is_learner & is_pre_learning & anticipatory_licks == 0);
mPFC_plus_s2_post_with_lick = mPFC_response_all_trials(stage_2_valid & is_learner & is_post_learning & anticipatory_licks == 1);
mPFC_plus_s2_post_no_lick = mPFC_response_all_trials(stage_2_valid & is_learner & is_post_learning & anticipatory_licks == 0);

% mPFC- (Non-learners) - Stage 2
mPFC_minus_s2_pre_with_lick = mPFC_response_all_trials(stage_2_valid & is_nonlearner & is_pre_learning & anticipatory_licks == 1);
mPFC_minus_s2_pre_no_lick = mPFC_response_all_trials(stage_2_valid & is_nonlearner & is_pre_learning & anticipatory_licks == 0);
mPFC_minus_s2_post_with_lick = mPFC_response_all_trials(stage_2_valid & is_nonlearner & is_post_learning & anticipatory_licks == 1);
mPFC_minus_s2_post_no_lick = mPFC_response_all_trials(stage_2_valid & is_nonlearner & is_post_learning & anticipatory_licks == 0);

% Statistical tests - Stage 2
[~, p_plus_s2_pre] = ttest2(mPFC_plus_s2_pre_with_lick, mPFC_plus_s2_pre_no_lick);
[~, p_plus_s2_post] = ttest2(mPFC_plus_s2_post_with_lick, mPFC_plus_s2_post_no_lick);
[~, p_minus_s2_pre] = ttest2(mPFC_minus_s2_pre_with_lick, mPFC_minus_s2_pre_no_lick);
[~, p_minus_s2_post] = ttest2(mPFC_minus_s2_post_with_lick, mPFC_minus_s2_post_no_lick);

fprintf('Stage 2 - mPFC+: Pre p=%.4f, Post p=%.4f\n', p_plus_s2_pre, p_plus_s2_post);
fprintf('Stage 2 - mPFC-: Pre p=%.4f, Post p=%.4f\n', p_minus_s2_pre, p_minus_s2_post);

%% Plot ViS ROI Traces Time aligned - to examine whether the increase at post-learning anticipatory stage 1 is due to higher feedback

% mPFC+ (Learners) - Stage 1
mPFC_plus_s1_pre_with_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_pre_learning & anticipatory_licks~=0); % not zero for anticipatory licks
mPFC_plus_s1_pre_no_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_pre_learning & anticipatory_licks==0);
mPFC_plus_s1_post_with_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_post_learning & anticipatory_licks~=0);
mPFC_plus_s1_post_no_lick = roi_trace_all_trials(:,stage_1_valid & is_learner & is_post_learning & anticipatory_licks==0);

% Helper for computing mean and SEM
computeStats = @(data) deal( ...
    mean(data, 2), ... % mean over trials
    std(data, 0, 2) ./ sqrt(sum(~isnan(data), 2)) ... % SEM
);

[m_pre_with, sem_pre_with] = computeStats(mPFC_plus_s1_pre_with_lick);
[m_pre_no, sem_pre_no] = computeStats(mPFC_plus_s1_pre_no_lick);

[m_post_with, sem_post_with] = computeStats(mPFC_plus_s1_post_with_lick);
[m_post_no, sem_post_no] = computeStats(mPFC_plus_s1_post_no_lick);


% First pass: determine global y-axis limits across both groups
all_means_s1 = [m_pre_with, m_pre_no,m_post_with,m_post_no];
all_sems_s1 = [sem_pre_with, sem_pre_no,sem_post_with,sem_post_no];

% y_min_global = min(all_means_s1 + all_sems_s1,[],'all');
% y_max_global = max(all_means_s1 + all_sems_s1,[],'all'); 



% Plot
figure;

% Color definitions
c_with = cLearner * 0.4;  % dark
c_no   = cLearner * 0.8;  % light

% Pre-learning panel
subplot(1,2,1)
hold on


fill([added_time, fliplr(added_time)], ...
     [m_pre_with'+sem_pre_with', fliplr(m_pre_with'-sem_pre_with')], ...
     c_with, 'FaceAlpha',0.25, 'EdgeColor','none')
fill([added_time, fliplr(added_time)], ...
     [m_pre_no'+sem_pre_no', fliplr(m_pre_no'-sem_pre_no')], ...
     c_no, 'FaceAlpha',0.25, 'EdgeColor','none')

plot(added_time, m_pre_with', 'Color', c_with, 'LineWidth', 2)
plot(added_time, m_pre_no', 'Color', c_no, 'LineWidth', 2)

title('Pre-learning')
xlabel('Time (s)')
ylabel('ΔF/F')
ylim([y_min_global y_max_global]);
legend({'With lick','No lick'}, 'Location','best')
xline(0,'--k',HandleVisibility='off')

% Post-learning panel
subplot(1,2,2)
hold on



fill([added_time, fliplr(added_time)], ...
     [m_post_with'+sem_post_with', fliplr(m_post_with'-sem_post_with')], ...
     c_with, 'FaceAlpha',0.25, 'EdgeColor','none')
fill([added_time, fliplr(added_time)], ...
     [m_post_no'+sem_post_no', fliplr(m_post_no'-sem_post_no')], ...
     c_no, 'FaceAlpha',0.25, 'EdgeColor','none')

plot(added_time, m_post_with', 'Color', c_with, 'LineWidth', 2)
plot(added_time, m_post_no', 'Color', c_no, 'LineWidth', 2)

title('Post-learning')
xlabel('Time (s)')
ylim([y_min_global y_max_global]);
legend({'With lick','No lick'}, 'Location','best')
xline(0,'--k',HandleVisibility='off')

sgtitle('mPFC+ Stage 1: Left mPFC activity Seperated Trials with and without Anticipatory Licks')

%% Plot scatter plots of mPFC vs RT 

% Create a cutoff for really long trials- for plotting purposes
slow_trials_to_cut_off = valid_diff_trials & diff_from_optimal <= 10;

% Create the masks
mPFC_plus_pre_stage_1_mask= (stage_1_valid & is_learner & is_pre_learning & slow_trials_to_cut_off);
mPFC_plus_post_stage_1_mask= (stage_1_valid & is_learner & is_post_learning & slow_trials_to_cut_off);

% Get the ROIs for all trials for learners for stage 1 split by pre-post
mPFC_plus_pre_stage_1= mPFC_response_all_trials(:,mPFC_plus_pre_stage_1_mask); 
mPFC_plus_post_stage_1= mPFC_response_all_trials(:,mPFC_plus_post_stage_1_mask); 

% Now get the RTs for the same masks
diff_from_optimal_mPFC_plus_stage_1_pre= diff_from_optimal(mPFC_plus_pre_stage_1_mask);
diff_from_optimal_mPFC_plus_stage_1_post=diff_from_optimal(mPFC_plus_post_stage_1_mask);

% Plot

figure('Position', [100 100 1400 600], 'Color', 'w');

% Subplot 1: mPFC+ (Learners) - Pre-Learning
subplot(2, 2, 1);
hold on;

scatter(diff_from_optimal_mPFC_plus_stage_1_pre', mPFC_plus_pre_stage_1, 80, cLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cLearner*0.8, 'LineWidth', 0.5);

% Add correlation
[r, p_corr] = corr(diff_from_optimal_mPFC_plus_stage_1_pre, mPFC_plus_pre_stage_1', 'rows', 'complete');

% Labels
% xticks([unique(lick_data)]);
% xticklabels('1');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Pre-Learning (n=%d trials)', sum(mPFC_plus_pre_stage_1_mask)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

% Add statistics text
text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


% Subplot 2: mPFC+ (Learners) - Post-Learning
subplot(2, 2, 2);
hold on;

scatter(diff_from_optimal_mPFC_plus_stage_1_post', mPFC_plus_post_stage_1, 80, cLearner, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cLearner*0.8, 'LineWidth', 0.5);

% Statistics
[r, p_corr] = corr(diff_from_optimal_mPFC_plus_stage_1_post, mPFC_plus_post_stage_1', 'rows', 'complete');

% xticks([0 1]);
% xticklabels({'No Lick', 'With Lick'});
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Post-Learning (n=%d trials)', sum(mPFC_plus_post_stage_1_mask)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% xlim([-0.3 1.3]);

%% Plots Trial by Trial mPFC vs RTs (can be split between mPFC+ and mPFC- groups)

% Create a cutoff for really long trials- for plotting purposes
slow_trials_to_cut_off = valid_diff_trials & diff_from_optimal <= 10;

% Create the masks
mPFC_plus_pre_stage_1_mask = (stage_1_valid & is_learner & is_pre_learning & slow_trials_to_cut_off);
mPFC_plus_post_stage_1_mask = (stage_1_valid & is_learner & is_post_learning & slow_trials_to_cut_off);

% Get the ROIs for all trials for learners for stage 1 split by pre-post
mPFC_plus_pre_stage_1 = mPFC_response_all_trials(:, mPFC_plus_pre_stage_1_mask);
mPFC_plus_post_stage_1 = mPFC_response_all_trials(:, mPFC_plus_post_stage_1_mask);

% Now get the RTs for the same masks
diff_from_optimal_mPFC_plus_stage_1_pre = diff_from_optimal(mPFC_plus_pre_stage_1_mask);
diff_from_optimal_mPFC_plus_stage_1_post = diff_from_optimal(mPFC_plus_post_stage_1_mask);

% Define colors with hue variation
color_pre = [0.5 0.5 0.5];   % grey
color_post = cLearner;         % Full color for post-learning

% Create single figure
figure('Position', [100 100 900 700], 'Color', 'w');
hold on;

% Plot pre-learning data
scatter(diff_from_optimal_mPFC_plus_stage_1_pre', mPFC_plus_pre_stage_1, 60, color_pre, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', color_pre*0.7, 'LineWidth', 0.5, ...
    'DisplayName', sprintf('Pre-Learning (n=%d)', sum(mPFC_plus_pre_stage_1_mask)));

% Plot post-learning data
scatter(diff_from_optimal_mPFC_plus_stage_1_post', mPFC_plus_post_stage_1, 60, color_post, 'filled', ...
    'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', color_post*0.7, 'LineWidth', 0.5, ...
    'DisplayName', sprintf('Post-Learning (n=%d)', sum(mPFC_plus_post_stage_1_mask)));

% Add regression lines
% Pre-learning regression
valid_pre = ~isnan(diff_from_optimal_mPFC_plus_stage_1_pre) & ~isnan(mPFC_plus_pre_stage_1');
if sum(valid_pre) > 2
    x_pre = diff_from_optimal_mPFC_plus_stage_1_pre(valid_pre);
    y_pre = mPFC_plus_pre_stage_1(valid_pre)';
    
    [r_pre, p_pre] = corr(x_pre, y_pre, 'rows', 'complete');
    
    % Fit linear model
    p_fit_pre = polyfit(x_pre, y_pre, 1);
    x_fit_pre = linspace(min(x_pre), max(x_pre), 100);
    y_fit_pre = polyval(p_fit_pre, x_fit_pre);
    
    plot(x_fit_pre, y_fit_pre, '-', 'Color', color_pre, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Pre Fit: r=%.3f, p=%.4f', r_pre, p_pre));
end

% Post-learning regression
valid_post = ~isnan(diff_from_optimal_mPFC_plus_stage_1_post) & ~isnan(mPFC_plus_post_stage_1');
if sum(valid_post) > 2
    x_post = diff_from_optimal_mPFC_plus_stage_1_post(valid_post);
    y_post = mPFC_plus_post_stage_1(valid_post)';
    
    [r_post, p_post] = corr(x_post, y_post, 'rows', 'complete');
    
    % Fit linear model
    p_fit_post = polyfit(x_post, y_post, 1);
    x_fit_post = linspace(min(x_post), max(x_post), 100);
    y_fit_post = polyval(p_fit_post, x_fit_post);
    
    plot(x_fit_post, y_fit_post, '-', 'Color', color_post, 'LineWidth', 2.5, ...
        'DisplayName', sprintf('Post Fit: r=%.3f, p=%.4f', r_post, p_post));
end

% Formatting
xlabel('Reaction Time (diff from optimal, s)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('mPFC+ Activity vs Reaction Time: Pre vs Post Learning', ...
    'FontSize', 14, 'FontWeight', 'bold');

% Legend
legend('Location', 'best', 'FontSize', 11);

grid on;
box on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
hold off;

%% Plot a 2D histogram for mPFC vs Anticipatory licks
figure('Position', [100 100 1400 600], 'Color', 'w');

% Subplot 1: mPFC+ (Learners) - Pre-Learning
subplot(2, 2, 1);
hold on;

% Get trial-level data for mPFC+ Pre-learning Stage 2
mask_plus_pre = stage_2_valid & is_learner & is_pre_learning;
lick_data = anticipatory_licks(mask_plus_pre);
mPFC_data = mPFC_response_all_trials(mask_plus_pre);

% Create 2D histogram with edges centered on integer values
x_edges = (min(lick_data)-0.5):1:(max(lick_data)+0.5);  % Centers bins on integers
y_edges = linspace(min(mPFC_data), max(mPFC_data), 30);

h = histogram2(lick_data, mPFC_data', x_edges, y_edges, ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'off', ...
    'EdgeColor', 'none', 'FaceColor', 'flat');

% More intense colormap
colormap(gca, [linspace(1, cLearner(1), 256)', ...
                linspace(1, cLearner(2), 256)', ...
                linspace(1 ,cLearner(3), 256)']);
cb = colorbar;
cb.Label.String = 'Trial Count';
cb.Label.FontSize = 11;

% Add correlation line
[r, p_corr] = corr(lick_data, mPFC_data', 'rows', 'complete');
coeffs = polyfit(lick_data, mPFC_data, 1);
x_line = linspace(min(lick_data), max(lick_data), 100);
y_line = polyval(coeffs, x_line);
plot(x_line, y_line, 'k--', 'LineWidth', 2.5);

% Labels
xlabel('Anticipatory Licks', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Pre-Learning (n=%d trials)', sum(mask_plus_pre)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

% Add statistics text
text(mean([min(lick_data), max(lick_data)]), max(mPFC_data)*0.95, ...
    sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);

% Set x-ticks to integer values for clarity
xticks(min(lick_data):max(lick_data));
set(gca, 'FontSize', 11, 'LineWidth', 1.5);

% Subplot 2: mPFC+ (Learners) - Post-Learning
subplot(2, 2, 2);
hold on;

% Get trial-level data for mPFC+ Post-learning Stage 2
mask_plus_post = stage_2_valid & is_learner & is_post_learning;
lick_data = anticipatory_licks(mask_plus_post);
mPFC_data = mPFC_response_all_trials(mask_plus_post);

% Create 2D histogram with edges centered on integer values
x_edges = (min(lick_data)-0.5):1:(max(lick_data)+0.5);
y_edges = linspace(min(mPFC_data), max(mPFC_data), 30);

h = histogram2(lick_data, mPFC_data', x_edges, y_edges, ...
    'DisplayStyle', 'tile', 'ShowEmptyBins', 'off', ...
    'EdgeColor', 'none', 'FaceColor', 'flat');

% More intense colormap (darker for post-learning)
colormap(gca, [linspace(1, cLearner(1)*0.8, 256)', ...
                linspace(1, cLearner(2)*0.8, 256)', ...
                linspace(1, cLearner(3)*0.8, 256)']);
cb = colorbar;
cb.Label.String = 'Trial Count';
cb.Label.FontSize = 11;

% Add correlation line
[r, p_corr] = corr(lick_data, mPFC_data', 'rows', 'complete');
coeffs = polyfit(lick_data, mPFC_data, 1);
x_line = linspace(min(lick_data), max(lick_data), 100);
y_line = polyval(coeffs, x_line);
plot(x_line, y_line, 'k--', 'LineWidth', 2.5);

% Labels
xlabel('Anticipatory Licks', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Post-Learning (n=%d trials)', sum(mask_plus_post)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

text(mean([min(lick_data), max(lick_data)]), max(mPFC_data)*0.95, ...
    sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold', ...
    'BackgroundColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.5);

xticks(min(lick_data):max(lick_data));
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


%% Stage 1 Bar plots
figure('Position', [100 100 1200 500], 'Color', 'w');

% First pass: determine global y-axis limits across both groups
all_means_s1 = [data_means_plus_s1, data_means_minus_s1];
all_sems_s1 = [data_sems_plus_s1, data_sems_minus_s1];
y_min_global = 0;
y_max_global = max(all_means_s1 + all_sems_s1) * 2; % Extra space for significance markers

% Subplot 1: mPFC+ Stage 2
subplot(1, 2, 1);
hold on;

% Prepare data for mPFC+
data_means_plus_s1 = [mean(mPFC_plus_s1_pre_with_lick), mean(mPFC_plus_s1_pre_no_lick), ...
                      mean(mPFC_plus_s1_post_with_lick), mean(mPFC_plus_s1_post_no_lick)];

data_sems_plus_s1 = [std(mPFC_plus_s1_pre_with_lick)/sqrt(length(mPFC_plus_s1_pre_with_lick)), ...
                     std(mPFC_plus_s1_pre_no_lick)/sqrt(length(mPFC_plus_s1_pre_no_lick)), ...
                     std(mPFC_plus_s1_post_with_lick)/sqrt(length(mPFC_plus_s1_post_with_lick)), ...
                     std(mPFC_plus_s1_post_no_lick)/sqrt(length(mPFC_plus_s1_post_no_lick))];

% Bar positions
x_plus = [1 2 4.5 5.5];

% Pre-learning - With Lick (dark)
b1 = bar(x_plus(1), data_means_plus_s1(1), 'FaceColor', cLearner*0.8, 'EdgeColor', 'k', 'LineWidth', 2);

% Pre-learning - No Lick (light)
b2 = bar(x_plus(2), data_means_plus_s1(2), 'FaceColor', cLearner*0.4, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - With Lick (dark)
bar(x_plus(3), data_means_plus_s1(3), 'FaceColor', cLearner, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - No Lick (light)
bar(x_plus(4), data_means_plus_s1(4), 'FaceColor', cLearner*0.5, 'EdgeColor', 'k', 'LineWidth', 2);

% Error bars
errorbar(x_plus, data_means_plus_s1, data_sems_plus_s1, 'k', 'LineStyle', 'none', ...
    'LineWidth', 2.5, 'CapSize', 12);

% Significance markers - Pre
y_sig_line = max(data_means_plus_s1(1:2) + data_sems_plus_s1(1:2)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_plus_s1_pre < 0.05
    sig_text = ''; 
    if p_plus_s1_pre < 0.001, sig_text = '***';
    elseif p_plus_s1_pre < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_plus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(1:2)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_plus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(1:2)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% Significance markers - Post
y_sig_line = max(data_means_plus_s1(3:4) + data_sems_plus_s1(3:4)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_plus_s1_post < 0.05
    sig_text = ''; 
    if p_plus_s1_post < 0.001, sig_text = '***';
    elseif p_plus_s1_post < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_plus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(3:4)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_plus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(3:4)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% X-axis labels
xticks([mean(x_plus(1:2)), mean(x_plus(3:4))]);
xticklabels({'Pre-Learning', 'Post-Learning'});
set(gca, 'XTickLabelRotation', 0);

ylabel('Left mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('mPFC+ (Learners) - Stage 1', 'FontSize', 15, 'FontWeight', 'bold', 'Color', cLearner);

% Legend
legend([b1, b2], {'With Anticipatory Lick', 'No Anticipatory Lick'}, ...
    'Location', 'northwest', 'FontSize', 12, 'Box', 'on');

grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 7]);
ylim([y_min_global y_max_global]); % Set consistent y-axis

% Subplot 2: mPFC- Stage 2
subplot(1, 2, 2);
hold on;

% Prepare data for mPFC-
data_means_minus_s1 = [mean(mPFC_minus_s1_pre_with_lick), mean(mPFC_minus_s1_pre_no_lick), ...
                       mean(mPFC_minus_s1_post_with_lick), mean(mPFC_minus_s1_post_no_lick)];

data_sems_minus_s1 = [std(mPFC_minus_s1_pre_with_lick)/sqrt(length(mPFC_minus_s1_pre_with_lick)), ...
                      std(mPFC_minus_s1_pre_no_lick)/sqrt(length(mPFC_minus_s1_pre_no_lick)), ...
                      std(mPFC_minus_s1_post_with_lick)/sqrt(length(mPFC_minus_s1_post_with_lick)), ...
                      std(mPFC_minus_s1_post_no_lick)/sqrt(length(mPFC_minus_s1_post_no_lick))];

x_minus = [1 2 4.5 5.5];

% Pre-learning - With Lick (dark)
bar(x_minus(1), data_means_minus_s1(1), 'FaceColor', cNonLearner*0.8, 'EdgeColor', 'k', 'LineWidth', 2);

% Pre-learning - No Lick (light)
bar(x_minus(2), data_means_minus_s1(2), 'FaceColor', cNonLearner*0.4, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - With Lick (dark)
bar(x_minus(3), data_means_minus_s1(3), 'FaceColor', cNonLearner, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - No Lick (light)
bar(x_minus(4), data_means_minus_s1(4), 'FaceColor', cNonLearner*0.5, 'EdgeColor', 'k', 'LineWidth', 2);

% Error bars
errorbar(x_minus, data_means_minus_s1, data_sems_minus_s1, 'k', 'LineStyle', 'none', ...
    'LineWidth', 2.5, 'CapSize', 12);

% Significance markers - Pre
y_sig_line = max(data_means_minus_s1(1:2) + data_sems_minus_s1(1:2)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_minus_s1_pre < 0.05
    sig_text = ''; 
    if p_minus_s1_pre < 0.001, sig_text = '***';
    elseif p_minus_s1_pre < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_minus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(1:2)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_minus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(1:2)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% Significance markers - Post
y_sig_line = max(data_means_minus_s1(3:4) + data_sems_minus_s1(3:4)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_minus_s1_post < 0.05
    sig_text = ''; 
    if p_minus_s1_post < 0.001, sig_text = '***';
    elseif p_minus_s1_post < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_minus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(3:4)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_minus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(3:4)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% X-axis labels
xticks([mean(x_minus(1:2)), mean(x_minus(3:4))]);
xticklabels({'Pre-Learning', 'Post-Learning'});
set(gca, 'XTickLabelRotation', 0);

ylabel('Left mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('mPFC- (Non-learners) - Stage 1', 'FontSize', 15, 'FontWeight', 'bold', 'Color', cNonLearner);

% Legend
legend({'With Anticipatory Lick', 'No Anticipatory Lick'}, ...
    'Location', 'northwest', 'FontSize', 12, 'Box', 'on');

grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 7]);
ylim([y_min_global y_max_global]); % Set consistent y-axis

sgtitle('Stage 1: Widefield Activity by Anticipatory Lick × Learning Stage', 'FontSize', 17, 'FontWeight', 'bold');


%% ===== FIGURE 2: STAGE 2 - Bar plots (CLEAN VERSION) =====

figure('Position', [100 100 1200 500], 'Color', 'w');

% First pass: determine global y-axis limits across both groups
all_means_s2 = [data_means_plus_s2, data_means_minus_s2];
all_sems_s2 = [data_sems_plus_s2, data_sems_minus_s2];
y_min_global = 0;
y_max_global = max(all_means_s2 + all_sems_s2) * 2; % Extra space for significance markers

% Subplot 1: mPFC+ Stage 2
subplot(1, 2, 1);
hold on;

% Prepare data for mPFC+
data_means_plus_s2 = [mean(mPFC_plus_s2_pre_with_lick), mean(mPFC_plus_s2_pre_no_lick), ...
                      mean(mPFC_plus_s2_post_with_lick), mean(mPFC_plus_s2_post_no_lick)];

data_sems_plus_s2 = [std(mPFC_plus_s2_pre_with_lick)/sqrt(length(mPFC_plus_s2_pre_with_lick)), ...
                     std(mPFC_plus_s2_pre_no_lick)/sqrt(length(mPFC_plus_s2_pre_no_lick)), ...
                     std(mPFC_plus_s2_post_with_lick)/sqrt(length(mPFC_plus_s2_post_with_lick)), ...
                     std(mPFC_plus_s2_post_no_lick)/sqrt(length(mPFC_plus_s2_post_no_lick))];

% Bar positions
x_plus = [1 2 4.5 5.5];

% Pre-learning - With Lick (dark)
b1 = bar(x_plus(1), data_means_plus_s2(1), 'FaceColor', cLearner*0.8, 'EdgeColor', 'k', 'LineWidth', 2);

% Pre-learning - No Lick (light)
b2 = bar(x_plus(2), data_means_plus_s2(2), 'FaceColor', cLearner*0.4, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - With Lick (dark)
bar(x_plus(3), data_means_plus_s2(3), 'FaceColor', cLearner, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - No Lick (light)
bar(x_plus(4), data_means_plus_s2(4), 'FaceColor', cLearner*0.5, 'EdgeColor', 'k', 'LineWidth', 2);

% Error bars
errorbar(x_plus, data_means_plus_s2, data_sems_plus_s2, 'k', 'LineStyle', 'none', ...
    'LineWidth', 2.5, 'CapSize', 12);

% Significance markers - Pre
y_sig_line = max(data_means_plus_s2(1:2) + data_sems_plus_s2(1:2)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_plus_s2_pre < 0.05
    sig_text = ''; 
    if p_plus_s2_pre < 0.001, sig_text = '***';
    elseif p_plus_s2_pre < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_plus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(1:2)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_plus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(1:2)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% Significance markers - Post
y_sig_line = max(data_means_plus_s2(3:4) + data_sems_plus_s2(3:4)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_plus_s2_post < 0.05
    sig_text = ''; 
    if p_plus_s2_post < 0.001, sig_text = '***';
    elseif p_plus_s2_post < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_plus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(3:4)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_plus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_plus(3:4)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% X-axis labels
xticks([mean(x_plus(1:2)), mean(x_plus(3:4))]);
xticklabels({'Pre-Learning', 'Post-Learning'});
set(gca, 'XTickLabelRotation', 0);

ylabel('mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('mPFC+ (Learners) - Stage 2', 'FontSize', 15, 'FontWeight', 'bold', 'Color', cLearner);

% Legend
legend([b1, b2], {'With Anticipatory Lick', 'No Anticipatory Lick'}, ...
    'Location', 'northwest', 'FontSize', 12, 'Box', 'on');

grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 7]);
ylim([y_min_global y_max_global]); % Set consistent y-axis

% Subplot 2: mPFC- Stage 2
subplot(1, 2, 2);
hold on;

% Prepare data for mPFC-
data_means_minus_s2 = [mean(mPFC_minus_s2_pre_with_lick), mean(mPFC_minus_s2_pre_no_lick), ...
                       mean(mPFC_minus_s2_post_with_lick), mean(mPFC_minus_s2_post_no_lick)];

data_sems_minus_s2 = [std(mPFC_minus_s2_pre_with_lick)/sqrt(length(mPFC_minus_s2_pre_with_lick)), ...
                      std(mPFC_minus_s2_pre_no_lick)/sqrt(length(mPFC_minus_s2_pre_no_lick)), ...
                      std(mPFC_minus_s2_post_with_lick)/sqrt(length(mPFC_minus_s2_post_with_lick)), ...
                      std(mPFC_minus_s2_post_no_lick)/sqrt(length(mPFC_minus_s2_post_no_lick))];

x_minus = [1 2 4.5 5.5];

% Pre-learning - With Lick (dark)
bar(x_minus(1), data_means_minus_s2(1), 'FaceColor', cNonLearner*0.8, 'EdgeColor', 'k', 'LineWidth', 2);

% Pre-learning - No Lick (light)
bar(x_minus(2), data_means_minus_s2(2), 'FaceColor', cNonLearner*0.4, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - With Lick (dark)
bar(x_minus(3), data_means_minus_s2(3), 'FaceColor', cNonLearner, 'EdgeColor', 'k', 'LineWidth', 2);

% Post-learning - No Lick (light)
bar(x_minus(4), data_means_minus_s2(4), 'FaceColor', cNonLearner*0.5, 'EdgeColor', 'k', 'LineWidth', 2);

% Error bars
errorbar(x_minus, data_means_minus_s2, data_sems_minus_s2, 'k', 'LineStyle', 'none', ...
    'LineWidth', 2.5, 'CapSize', 12);

% Significance markers - Pre
y_sig_line = max(data_means_minus_s2(1:2) + data_sems_minus_s2(1:2)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_minus_s2_pre < 0.05
    sig_text = ''; 
    if p_minus_s2_pre < 0.001, sig_text = '***';
    elseif p_minus_s2_pre < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_minus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(1:2)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_minus(1:2), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(1:2)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% Significance markers - Post
y_sig_line = max(data_means_minus_s2(3:4) + data_sems_minus_s2(3:4)) * 1.08;
y_sig_text = y_sig_line * 1.05;

if p_minus_s2_post < 0.05
    sig_text = ''; 
    if p_minus_s2_post < 0.001, sig_text = '***';
    elseif p_minus_s2_post < 0.01, sig_text = '**'; 
    else, sig_text = '*'; 
    end
    plot(x_minus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(3:4)), y_sig_text, sig_text, 'HorizontalAlignment', 'center', ...
        'FontSize', 20, 'FontWeight', 'bold');
else
    plot(x_minus(3:4), [y_sig_line y_sig_line], 'k-', 'LineWidth', 2.5);
    text(mean(x_minus(3:4)), y_sig_text, 'n.s.', 'HorizontalAlignment', 'center', ...
        'FontSize', 11);
end

% X-axis labels
xticks([mean(x_minus(1:2)), mean(x_minus(3:4))]);
xticklabels({'Pre-Learning', 'Post-Learning'});
set(gca, 'XTickLabelRotation', 0);

ylabel('mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('mPFC- (Non-learners) - Stage 2', 'FontSize', 15, 'FontWeight', 'bold', 'Color', cNonLearner);

% Legend
legend({'With Anticipatory Lick', 'No Anticipatory Lick'}, ...
    'Location', 'northwest', 'FontSize', 12, 'Box', 'on');

grid on;
set(gca, 'FontSize', 12, 'LineWidth', 1.5);
xlim([0 7]);
ylim([y_min_global y_max_global]); % Set consistent y-axis

sgtitle('Stage 2: mPFC Activity by Anticipatory Lick × Learning Stage', 'FontSize', 17, 'FontWeight', 'bold');

%% ===== FIGURE: Scatter Plot - Anticipatory Licks vs mPFC Activity (Trial-level) =====

figure('Position', [100 100 1400 600], 'Color', 'w');

% Subplot 1: mPFC+ (Learners) - Pre-Learning
subplot(2, 2, 1);
hold on;

% Get trial-level data for mPFC+ Pre-learning Stage 2
mask_plus_pre = stage_2_valid & is_learner & is_pre_learning;
lick_data = anticipatory_licks(mask_plus_pre);
mPFC_data = mPFC_response_all_trials(mask_plus_pre);

% Scatter plot with jitter on x-axis
jitter_amount = 0.08;
lick_jittered = lick_data + (rand(size(lick_data)) - 0.5) * jitter_amount;

scatter(lick_jittered, mPFC_data, 80, cLearner*0.6, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cLearner*0.8, 'LineWidth', 0.5);

% Statistical test
[~, p_val] = ttest2(mPFC_data(lick_data == 0), mPFC_data(lick_data == 1));

% Add correlation
[r, p_corr] = corr(lick_data, mPFC_data', 'rows', 'complete');

% Labels
% xticks([unique(lick_data)]);
% xticklabels('1');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Pre-Learning (n=%d trials)', sum(mask_plus_pre)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

% Add statistics text
text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);


% Subplot 2: mPFC+ (Learners) - Post-Learning
subplot(2, 2, 2);
hold on;

% Get trial-level data for mPFC+ Post-learning Stage 2
mask_plus_post = stage_2_valid & is_learner & is_post_learning;
lick_data = anticipatory_licks(mask_plus_post);
mPFC_data = mPFC_response_all_trials(mask_plus_post);

% Scatter plot
lick_jittered = lick_data + (rand(size(lick_data)) - 0.5) * jitter_amount;
scatter(lick_jittered, mPFC_data, 80, cLearner, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cLearner*0.8, 'LineWidth', 0.5);

% % Means and SEMs
% mean_no_lick = mean(mPFC_data(lick_data == 0), 'omitnan');
% sem_no_lick = std(mPFC_data(lick_data == 0), 'omitnan') / sqrt(sum(lick_data == 0));
% 
% mean_with_lick = mean(mPFC_data(lick_data == 1), 'omitnan');
% sem_with_lick = std(mPFC_data(lick_data == 1), 'omitnan') / sqrt(sum(lick_data == 1));
% 
% errorbar(0, mean_no_lick, sem_no_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
%     'LineWidth', 3, 'CapSize', 15);
% errorbar(1, mean_with_lick, sem_with_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
%     'LineWidth', 3, 'CapSize', 15);
% 
% plot([0 1], [mean_no_lick mean_with_lick], 'k--', 'LineWidth', 2);

% Statistics
[~, p_val] = ttest2(mPFC_data(lick_data == 0), mPFC_data(lick_data == 1));
[r, p_corr] = corr(lick_data, mPFC_data', 'rows', 'complete');

% xticks([0 1]);
% xticklabels({'No Lick', 'With Lick'});
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC+ Post-Learning (n=%d trials)', sum(mask_plus_post)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);

text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
% xlim([-0.3 1.3]);

%% Subplot 3: mPFC- (Non-learners) - Pre-Learning
subplot(2, 2, 3);
hold on;

% Get trial-level data for mPFC- Pre-learning Stage 2
mask_minus_pre = stage_2_valid & is_nonlearner & is_pre_learning;
lick_data = anticipatory_licks(mask_minus_pre);
mPFC_data = mPFC_response_all_trials(mask_minus_pre);

% Scatter plot
lick_jittered = lick_data + (rand(size(lick_data)) - 0.5) * jitter_amount;
scatter(lick_jittered, mPFC_data, 80, cNonLeaner*0.6, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cNonLeaner*0.8, 'LineWidth', 0.5);

% Means and SEMs
mean_no_lick = mean(mPFC_data(lick_data == 0), 'omitnan');
sem_no_lick = std(mPFC_data(lick_data == 0), 'omitnan') / sqrt(sum(lick_data == 0));

mean_with_lick = mean(mPFC_data(lick_data == 1), 'omitnan');
sem_with_lick = std(mPFC_data(lick_data == 1), 'omitnan') / sqrt(sum(lick_data == 1));

errorbar(0, mean_no_lick, sem_no_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
    'LineWidth', 3, 'CapSize', 15);
errorbar(1, mean_with_lick, sem_with_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
    'LineWidth', 3, 'CapSize', 15);

plot([0 1], [mean_no_lick mean_with_lick], 'k--', 'LineWidth', 2);

% Statistics
[~, p_val] = ttest2(mPFC_data(lick_data == 0), mPFC_data(lick_data == 1));
[r, p_corr] = corr(lick_data, mPFC_data, 'rows', 'complete');

xticks([0 1]);
xticklabels({'No Lick', 'With Lick'});
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC- Pre-Learning (n=%d trials)', sum(mask_minus_pre)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cNonLeaner);

text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
xlim([-0.3 1.3]);

%% Subplot 4: mPFC- (Non-learners) - Post-Learning
subplot(2, 2, 4);
hold on;

% Get trial-level data for mPFC- Post-learning Stage 2
mask_minus_post = stage_2_valid & is_nonlearner & is_post_learning;
lick_data = anticipatory_licks(mask_minus_post);
mPFC_data = mPFC_response_all_trials(mask_minus_post);

% Scatter plot
lick_jittered = lick_data + (rand(size(lick_data)) - 0.5) * jitter_amount;
scatter(lick_jittered, mPFC_data, 80, cNonLeaner, 'filled', 'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', cNonLeaner*0.8, 'LineWidth', 0.5);

% Means and SEMs
mean_no_lick = mean(mPFC_data(lick_data == 0), 'omitnan');
sem_no_lick = std(mPFC_data(lick_data == 0), 'omitnan') / sqrt(sum(lick_data == 0));

mean_with_lick = mean(mPFC_data(lick_data == 1), 'omitnan');
sem_with_lick = std(mPFC_data(lick_data == 1), 'omitnan') / sqrt(sum(lick_data == 1));

errorbar(0, mean_no_lick, sem_no_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
    'LineWidth', 3, 'CapSize', 15);
errorbar(1, mean_with_lick, sem_with_lick, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'k', ...
    'LineWidth', 3, 'CapSize', 15);

plot([0 1], [mean_no_lick mean_with_lick], 'k--', 'LineWidth', 2);

% Statistics
[~, p_val] = ttest2(mPFC_data(lick_data == 0), mPFC_data(lick_data == 1));
[r, p_corr] = corr(lick_data, mPFC_data, 'rows', 'complete');

xticks([0 1]);
xticklabels({'No Lick', 'With Lick'});
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title(sprintf('mPFC- Post-Learning (n=%d trials)', sum(mask_minus_post)), ...
    'FontSize', 13, 'FontWeight', 'bold', 'Color', cNonLeaner);

text(0.5, max(mPFC_data)*0.95, sprintf('r=%.3f, p=%.4f', r, p_corr), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'BackgroundColor', 'w', 'EdgeColor', 'k');

grid on;
set(gca, 'FontSize', 11, 'LineWidth', 1.5);
xlim([-0.3 1.3]);

sgtitle('Stage 2: Trial-Level Relationship - Anticipatory Lick vs mPFC Activity', ...
    'FontSize', 17, 'FontWeight', 'bold');

fprintf('\n✓ Trial-level scatter plot created!\n');





%% Plots Bimodal Distribution Visualization for the mPFC magnitude at end of stage 1 to illustrate the split between mPFC+ and mPFC- groups

% ===== Extract mPFC Magnitude from Last TWO Days of Stage 1 =====

% Parameters
ROI_to_plot = left_mPFC_ROI_mask; % Left mPFC ROI index
n_components; % Number of SVD components

% Time window for CS+ response
time_window_start = 0;    % Start of response window (s)
time_window_end = 0.35;    % End of response window (s)

% Initialize storage
mPFC_magnitude_last_two_days = nan(numel(behaviour_data), 1);
mPFC_metadata_two_days = struct();


for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    
    if isempty(rd)
        fprintf('Animal %s: No recording days, skipping\n', aid);
        continue;
    end
    
    % Find classical conditioning days (Stage 1)
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'static'), rd);
    
    if ~any(isClassical)
        fprintf('Animal %s: No classical days, skipping\n', aid);
        continue;
    end
    
    % Check if learned Stage 1
    classical_days = find(isClassical);
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');
    
    if isempty(learned_day)
        fprintf('Animal %s: Never learned Stage 1, skipping\n', aid);
        continue;
    end
    
    % Get last TWO classical days
    if length(classical_days) < 2
        fprintf('Animal %s: Only 1 classical day, using just that day\n', aid);
        days_to_use = classical_days(end);
    else
        days_to_use = classical_days(end-1:end); % Last two days
    end

    % Extract mPFC for each of the last two days
    mPFC_responses = nan(1, length(days_to_use)); % Pre-allocate for number of days to use
    
    day_counter=0;

    for day_idx = days_to_use

        day_counter= day_counter +1;
        % Find matching widefield index
        matching_idx = find(widefield_animal_idx == a); % get all days indecies
        matching_idx= matching_idx(day_idx); % get the last two days

        if isempty(matching_idx)
            fprintf('  Day %d: No widefield data\n', day_idx);
            continue;
        end
        
        mPFC_trace_avg = ap.wf_roi(wf_U(:,:,1:kernel_n_components),  rewarded_stim_kernel(:,:,matching_idx(1)), [], [], ROI_to_plot);

        % Find response window indices
        response_window_idx = added_time_Kernel > time_window_start & added_time_Kernel <= time_window_end;
        baseline_window_idx = added_time_Kernel >= -0.1 & added_time_Kernel < 0;

        % Calculate baseline-corrected response
        baseline = mean(mPFC_trace_avg(baseline_window_idx), 'omitnan');
        response = max(mPFC_trace_avg(response_window_idx));
        mPFC_response_corrected = response - baseline;
        
        % Store this day's response
        mPFC_responses(day_counter) = mPFC_response_corrected;
        
    end
    
    % Average across the last two days
    if ~isempty(mPFC_responses)
        mPFC_magnitude_last_two_days(a) = mean(mPFC_responses);
        
        % Store metadata
        mPFC_metadata_two_days(a).animal_id = aid;
        mPFC_metadata_two_days(a).days_used = days_to_use;
        mPFC_metadata_two_days(a).responses_per_day = mPFC_responses;
        mPFC_metadata_two_days(a).mean_response = mean(mPFC_responses);
        
        fprintf('Animal %s: Mean mPFC (last 2 days) = %.4f\n', aid, mean(mPFC_responses));
    else
        fprintf('Animal %s: No valid mPFC data extracted\n', aid);
    end
end



valid_animals = ~isnan(mPFC_magnitude_last_two_days);
n_valid = sum(valid_animals);

learners_group_ID = {'HA005','HA008','HA010','HA011','HA012'}; % Define learners

% create a mask for learners
is_learner = ismember(animal_list, learners_group_ID);

% Get valid data
mPFC_data = mPFC_magnitude_last_two_days(valid_animals);
is_learner_valid = is_learner(valid_animals);

% Prepare learner vs non-learner data
learner_mags = mPFC_data(is_learner_valid);
nonlearner_mags = mPFC_data(~is_learner_valid);

fprintf('\n===== Group Statistics =====\n');
fprintf('Learners: %.4f ± %.4f (n=%d)\n', mean(learner_mags), std(learner_mags), length(learner_mags));
fprintf('Non-learners: %.4f ± %.4f (n=%d)\n', mean(nonlearner_mags), std(nonlearner_mags), length(nonlearner_mags));

% Statistical test
[~, p_ttest] = ttest2(learner_mags, nonlearner_mags);
fprintf('T-test p-value: %.4f\n', p_ttest);

%% ===== Single Figure: Before and After Reveal (Cleaned Up) =====

figure('Color', 'w', 'Position', [100 100 1400 600]);

% Get sorted data
[sorted_mag, sort_idx] = sort(mPFC_data);
valid_indices = find(valid_animals);
sorted_animal_idx = valid_indices(sort_idx);

% Identify learners in sorted order
is_learner_sorted = false(size(sorted_mag));
for i = 1:length(sorted_mag)
    animal_idx = sorted_animal_idx(i);
    aid = string(behaviour_data(animal_idx).animal_id);
    is_learner_sorted(i) = ismember(aid, learners_group_ID);
end

% Find gap
diffs = diff(sorted_mag);
[~, gap_idx] = max(diffs);
gap_y = mean(sorted_mag(gap_idx:gap_idx+1));

% Panel 1: Before Reveal (Black dots, disconnected)
subplot(2, 1, 1);
hold on;

% Plot as individual black dots (no connection)
scatter(1:length(sorted_mag), sorted_mag, 80, 'k', 'filled', 'MarkerFaceAlpha', 0.7);

ylabel('mPFC Activity (ΔF/F)', 'FontSize', 16, 'FontWeight', 'bold');
title('Two Natural Groups in Stage 1 Neural Data', 'FontSize', 18, 'FontWeight', 'bold');
yline(0, 'k--', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'LineWidth', 2, 'XTick', [], 'Box', 'off');
grid on;
xlim([0, length(sorted_mag)+1]);

% Panel 2: After Reveal (Colored circles with your predefined colors)
subplot(2, 1, 2);
hold on;

% Plot colored by outcome using circles
for i = 1:length(sorted_mag)
    if is_learner_sorted(i)
        color = cLearner; % Your predefined color
        marker = 'o'; % Circle for learners
    else
        color = cNonLearner; % Your predefined color
        marker = 'o'; % Circle for non-learners
    end
    
    plot(i, sorted_mag(i), marker, 'MarkerSize', 14, ...
        'MarkerFaceColor', color, 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
end

% Connect within groups with lines
learner_indices = find(is_learner_sorted);
nonlearner_indices = find(~is_learner_sorted);

plot(learner_indices, sorted_mag(learner_indices), '-', ...
    'Color', cLearner, 'LineWidth', 3);
plot(nonlearner_indices, sorted_mag(nonlearner_indices), '-', ...
    'Color', cNonLearner, 'LineWidth', 3);



xlabel('Animals (Sorted by Stage 1 mPFC Activity)', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 16, 'FontWeight', 'bold');
title('Groups Predict Stage 2 Transfer Success!', 'FontSize', 18, 'FontWeight', 'bold');

% Legend with circles
learner_patch = plot(NaN, NaN, 'o', 'MarkerSize', 14, 'MarkerFaceColor', cLearner, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
nonlearner_patch = plot(NaN, NaN, 'o', 'MarkerSize', 14, 'MarkerFaceColor', cNonLearner, ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
legend([learner_patch, nonlearner_patch], ...
    {sprintf('Learners (n=%d)', length(learner_mags)), ...
     sprintf('Non-learners (n=%d)', length(nonlearner_mags))}, ...
    'Location', 'northwest', 'FontSize', 14, 'Box', 'off');

yline(0, 'k--', 'LineWidth', 2);
set(gca, 'FontSize', 14, 'LineWidth', 2, 'Box', 'off');
grid on;
xlim([0, length(sorted_mag)+1]);

%% Plots a Scatter of left vs right mPFC for the last two days of stage 1

% Parameters
left_mPFC_ROI = left_mPFC_ROI_mask;   % Left mPFC ROI index
right_mPFC_ROI = right_mPFC_ROI_mask; % Right mPFC ROI index
n_components; % Number of SVD components

% Time window for CS+ response
time_window_start = 0;    % Start of response window (s)
time_window_end = 0.3;   % End of response window (s)

% Initialize storage
left_mPFC_magnitude = nan(numel(behaviour_data), 1);
right_mPFC_magnitude = nan(numel(behaviour_data), 1);
mPFC_metadata_two_days = struct();

for a = 1:numel(behaviour_data)
    aid = string(behaviour_data(a).animal_id);
    rd = behaviour_data(a).recording_day;
    
    if isempty(rd)
        fprintf('Animal %s: No recording days, skipping\n', aid);
        continue;
    end
    
    % Find classical conditioning days (Stage 1)
    isClassical = arrayfun(@(d) isfield(d,'workflow') && contains(d.workflow,'right_move'), rd);
    
    if ~any(isClassical)
        fprintf('Animal %s: No classical days, skipping\n', aid);
        continue;
    end
    
    % Check if learned Stage 1
    classical_days = find(isClassical);
    sigDays = combined_sig_day_all_protocols{a}(classical_days);
    learned_day = find(sigDays, 1, 'first');
    
    if isempty(learned_day)
        % fprintf('Animal %s: Never learned Stage 1, using last two days\n', aid);
        % days_to_use = classical_days(end-1:end);
        fprintf('Animal %s: Never learned Stage 1, skipping\n', aid);
        continue;
    end
    
    % Get last TWO classical days
    if length(classical_days) < 2
        fprintf('Animal %s: Only 1 classical day, using just that day\n', aid);
        days_to_use = classical_days(end);
    else
        days_to_use = classical_days(end-1:end); % Last two days
    end

    % Extract mPFC for each of the last two days
    left_mPFC_responses = nan(1, length(days_to_use));
    right_mPFC_responses = nan(1, length(days_to_use));
    
    day_counter = 0;

    for day_idx = days_to_use
        day_counter = day_counter + 1;
        
        % Find matching widefield index
        matching_idx = find(widefield_animal_idx == a); % get all days indices
        matching_idx = matching_idx(day_idx); % get the specific day

        if isempty(matching_idx)
            fprintf('  Day %d: No widefield data\n', day_idx);
            continue;
        end
        
        % Extract LEFT mPFC ROI activity
        left_mPFC_trace_avg = ap.wf_roi(wf_U(:,:,1:kernel_n_components), ...
            rewarded_stim_kernel(:,:,matching_idx(1)), [], [], left_mPFC_ROI);
        
        % Extract RIGHT mPFC ROI activity
        right_mPFC_trace_avg = ap.wf_roi(wf_U(:,:,1:kernel_n_components), ...
            rewarded_stim_kernel(:,:,matching_idx(1)), [], [], right_mPFC_ROI);

        % Find response window indices
        response_window_idx = added_time_Kernel > time_window_start & added_time_Kernel <= time_window_end;
        baseline_window_idx = added_time_Kernel >= -0.1 & added_time_Kernel < 0;

        % Calculate baseline-corrected response for LEFT mPFC
        baseline_left = mean(left_mPFC_trace_avg(baseline_window_idx), 'omitnan');
        response_left = mean(left_mPFC_trace_avg(response_window_idx));
        left_mPFC_response_corrected = response_left - baseline_left;
        
        % Calculate baseline-corrected response for RIGHT mPFC
        baseline_right = mean(right_mPFC_trace_avg(baseline_window_idx), 'omitnan');
        response_right = mean(right_mPFC_trace_avg(response_window_idx));
        right_mPFC_response_corrected = response_right - baseline_right;
        
        % Store this day's responses
        left_mPFC_responses(day_counter) = left_mPFC_response_corrected;
        right_mPFC_responses(day_counter) = right_mPFC_response_corrected;
    end
    
    % Average across the last two days
    if ~isempty(left_mPFC_responses) && any(~isnan(left_mPFC_responses))
        left_mPFC_magnitude(a) = mean(left_mPFC_responses, 'omitnan');
        right_mPFC_magnitude(a) = mean(right_mPFC_responses, 'omitnan');
        
        % Store metadata
        mPFC_metadata_two_days(a).animal_id = aid;
        mPFC_metadata_two_days(a).days_used = days_to_use;
        mPFC_metadata_two_days(a).left_responses_per_day = left_mPFC_responses;
        mPFC_metadata_two_days(a).right_responses_per_day = right_mPFC_responses;
        mPFC_metadata_two_days(a).mean_left = mean(left_mPFC_responses, 'omitnan');
        mPFC_metadata_two_days(a).mean_right = mean(right_mPFC_responses, 'omitnan');
        
        fprintf('Animal %s: Left mPFC = %.4f, Right mPFC = %.4f\n', ...
            aid, left_mPFC_magnitude(a), right_mPFC_magnitude(a));
    else
        fprintf('Animal %s: No valid mPFC data extracted\n', aid);
    end
end

% ===== Prepare Data for Plotting =====

valid_animals = ~isnan(left_mPFC_magnitude) & ~isnan(right_mPFC_magnitude);
n_valid = sum(valid_animals);

% Create a mask for learners
is_learner = ismember(animal_list, learners_group_ID);

% Get valid data
left_mPFC_data = left_mPFC_magnitude(valid_animals);
right_mPFC_data = right_mPFC_magnitude(valid_animals);
is_learner_valid = is_learner(valid_animals);

% Separate by group
learner_left = left_mPFC_data(is_learner_valid);
learner_right = right_mPFC_data(is_learner_valid);
nonlearner_left = left_mPFC_data(~is_learner_valid);
nonlearner_right = right_mPFC_data(~is_learner_valid);


% ===== Two-Panel Scatter Plot =====

figure('Color', 'w', 'Position', [100 100 1400 600]);

% Find axis limits (same for both panels)
all_data = [left_mPFC_data; right_mPFC_data];
axis_min = min(all_data) * 1.1;
axis_max = max(all_data) * 1.1;

% Panel 1: All gray (no color coding)
subplot(1, 2, 1);
hold on;

% Plot all animals in gray
scatter(left_mPFC_data, right_mPFC_data, 150, [0.5 0.5 0.5], 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);

% Add unity line
plot([axis_min axis_max], [axis_min axis_max], 'k--', 'LineWidth', 2, ...
    'DisplayName', 'Unity line');

xlabel('Left mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Right mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('Left vs Right mPFC 2 Last Days Stage 1' ,'FontSize', 16, 'FontWeight', 'bold');

% Add diagonal reference lines
% xline(0, 'k:', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1.5);

legend('Location', 'northwest', 'FontSize', 11);
axis equal;
xlim([axis_min axis_max]);
ylim([axis_min axis_max]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');

% Panel 2: Color coded by learner/non-learner
subplot(1, 2, 2);
hold on;

% Plot non-learners
scatter(nonlearner_left, nonlearner_right, 150, cNonLearner, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7, ...
    'DisplayName', sprintf('Non-learners (n=%d)', length(nonlearner_left)));

% Plot learners
scatter(learner_left, learner_right, 150, cLearner, 'filled', ...
    'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7, ...
    'DisplayName', sprintf('Learners (n=%d)', length(learner_left)));

% Add unity line
plot([axis_min axis_max], [axis_min axis_max], 'k--', 'LineWidth', 2, ...
    'DisplayName', 'Unity line');

xlabel('Left mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Right mPFC Activity (ΔF/F)', 'FontSize', 14, 'FontWeight', 'bold');
title('Stage 2 Transfer Groups', 'FontSize', 16, 'FontWeight', 'bold');

% Add diagonal reference lines
% xline(0, 'k:', 'LineWidth', 1.5);
% yline(0, 'k:', 'LineWidth', 1.5);

legend('Location', 'northwest', 'FontSize', 11);
axis equal;
xlim([axis_min axis_max]);
ylim([axis_min axis_max]);
set(gca, 'FontSize', 12, 'LineWidth', 1.5, 'Box', 'on');
% 
% sgtitle('Bilateral mPFC Activity During Stage 1 Predicts Stage 2 Success', ...
%     'FontSize', 18, 'FontWeight', 'bold');


%% Plots a summary plot of ROI activity across the different protocols

% data to index:  V×T×Ndays_total
aligned_all = rewarded_stim_kernel;

% --- 1) Choose time window and extract ROI traces (Ndays × Twin) ---
start_time = 0;  
end_time = 0.35;

% Baseline window for subtraction
baseline_start = -0.1;
baseline_end = 0;

% added_time;
event_time_window_mask = find(added_time_Kernel >= start_time & added_time_Kernel <= end_time);
baseline_window_mask = find(added_time_Kernel >= baseline_start & added_time_Kernel <= baseline_end);

% Define ROIs
ROI_to_plot_1=right_mPFC_ROI_mask;
ROI_to_plot_2=left_mPFC_ROI_mask;

% Extract ROI traces for event window
right_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,event_time_window_mask,:), [], [], ROI_to_plot_1), ...
    [3,2,1]);  % Ndays × Twin

left_mPFC_ROI_trace  = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,event_time_window_mask,:), [], [], ROI_to_plot_2), ...
    [3,2,1]);  % Ndays × Twin

% Extract ROI traces for baseline window
right_mPFC_baseline = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,baseline_window_mask,:), [], [], ROI_to_plot_1), ...
    [3,2,1]);  % Ndays × Tbaseline

left_mPFC_baseline  = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,baseline_window_mask,:), [], [], ROI_to_plot_2), ...
    [3,2,1]);  % Ndays × Tbaseline

% Compute baseline mean for each day
right_baseline_mean = mean(right_mPFC_baseline, 2, 'omitnan');  % Ndays × 1
left_baseline_mean = mean(left_mPFC_baseline, 2, 'omitnan');    % Ndays × 1

% Subtract baseline from event window
right_mPFC_ROI_trace = right_mPFC_ROI_trace - right_baseline_mean;
left_mPFC_ROI_trace = left_mPFC_ROI_trace - left_baseline_mean;

protocols = unique(workflow_cat(:))';          % plot segments in this order
groupColor = {cLearner, cNonLearner};          % 0->Learners, 1->Non-learners
labels     = {'Learners','Non-learners'};

% Minimum number of animals required per day
min_animals_per_day = 3;

% ---------- per-PROTOCOL stats (mean ± SEM vs relative day for each group) ----------
statsP = cell(numel(protocols),1);   % statsP{p}.G(g+1).rel/mR/sR/mL/sL

for ip = 1:numel(protocols)
    p = protocols(ip);

    % Store per-animal data with animal IDs for counting
    rel_all = []; 
    grp_all = []; 
    R_all = []; 
    L_all = [];
    animal_id_all = [];  % NEW: track which animal each data point comes from
    
    for ai = unique(widefield_animal_idx(:))'
        sel = (widefield_animal_idx==ai) & (workflow_cat==p);
        days_idx = find(sel);
        if isempty(days_idx), continue; end

        grp = is_group_animal(days_idx(1)); % skip non-learners for now

        % if grp==1
        %     continue
        % end
        % 

        ld  = find(learning_index_animal(days_idx)==1, 1);   % learning day within this protocol
        if isempty(ld), continue; end 

        n = numel(days_idx);
        rel_day = (1:n) - ld;

        meanR = max(right_mPFC_ROI_trace(days_idx,:),[], 2, 'omitnan');
        meanL = max(left_mPFC_ROI_trace (days_idx,:), [],2, 'omitnan');

        rel_all = [rel_all; rel_day(:)];
        grp_all = [grp_all; repmat(grp,n,1)];
        R_all   = [R_all;   meanR(:)];
        L_all   = [L_all;   meanL(:)];
        animal_id_all = [animal_id_all; repmat(ai, n, 1)];  % NEW: store animal ID
    end

    Gpresent = unique(grp_all(:))';
    stats = struct([]);
    
    for g = Gpresent
        mask = (grp_all==g);
        r    = rel_all(mask);
        Rg   = R_all(mask);
        Lg   = L_all(mask);
        animal_ids = animal_id_all(mask);  % NEW: get animal IDs for this group

        uRel = unique(r,'sorted');
        mR = nan(size(uRel)); 
        sR = nan(size(uRel));
        mL = nan(size(uRel)); 
        sL = nan(size(uRel));
        n_animals = nan(size(uRel));  % NEW: track number of animals per day
        
        for k = 1:numel(uRel)
            idx = (r==uRel(k));
            Rk = Rg(idx); 
            Lk = Lg(idx); 
            animals_this_day = animal_ids(idx);  % NEW: get animals for this day
            
            % Remove NaN values
            valid_R = ~isnan(Rk);
            valid_L = ~isnan(Lk);
            Rk = Rk(valid_R);
            Lk = Lk(valid_L);
            
            % Count unique animals for this day
            n_animals(k) = numel(unique(animals_this_day));  % NEW
            
            % Only compute statistics if we have minimum number of animals
            if n_animals(k) >= min_animals_per_day  % NEW: filter by animal count
                mR(k) = mean(Rk); 
                sR(k) = std(Rk)/max(1,sqrt(numel(Rk)));
                mL(k) = mean(Lk); 
                sL(k) = std(Lk)/max(1,sqrt(numel(Lk)));
            else
                % Set to NaN if insufficient animals
                mR(k) = NaN;
                sR(k) = NaN;
                mL(k) = NaN;
                sL(k) = NaN;
            end
        end
        
        % NEW: Filter out days with insufficient animals
        valid_days = n_animals >= min_animals_per_day;
        
        stats(g+1).rel = uRel(valid_days);
        stats(g+1).mR  = mR(valid_days); 
        stats(g+1).sR  = sR(valid_days);
        stats(g+1).mL  = mL(valid_days); 
        stats(g+1).sL  = sL(valid_days);
        stats(g+1).n_animals = n_animals(valid_days);  % NEW: store animal counts
    end
    statsP{ip} = stats;
end

% ---------- global y-limits across all protocols & groups ----------
yvals = [];
for ip = 1:numel(protocols)
    st = statsP{ip}; if isempty(st), continue; end
    for gi = 1:numel(st)
        if isempty(st(gi)), continue; end
        yvals = [yvals; st(gi).mR+st(gi).sR; st(gi).mR-st(gi).sR; ...
                        st(gi).mL+st(gi).sL; st(gi).mL-st(gi).sL];
    end
end
yvals = yvals(~isnan(yvals));
ymin = max(-3, min(yvals));  
ymax = max(yvals);

% ---------- PLOT: concatenate protocol segments, each aligned internally to its learning day ----------
figure('Color','w','Position', [100, 100, 1400, 800], ...
       'Name', sprintf('Summary by protocol (aligned, baseline-subtracted, min %d animals)', min_animals_per_day));
tiledlayout(2,1,'TileSpacing','compact');   % Row1=Right, Row2=Left

spacer = 1;   % gap (days) between protocol segments

for hemi = 1:2
    nexttile; hold on;
    offset = 0;    % cumulative x-offset (in "days")

    for ip = 1:numel(protocols)
        st = statsP{ip};
        if isempty(st), offset = offset + spacer; continue; end

        % Span within this protocol (relative day range across groups)
        minRel = inf; maxRel = -inf;
        for gi = 1:numel(st)
            if isempty(st(gi)) || isempty(st(gi).rel), continue; end
            minRel = min(minRel, min(st(gi).rel));
            maxRel = max(maxRel, max(st(gi).rel));
        end
        if ~isfinite(minRel) || ~isfinite(maxRel)
            offset = offset + spacer; continue;
        end
        span = maxRel - minRel;

        % Draw each group in its color
        for g = [0 1]   % force order: learners then non-learners
            gi = g+1;
            if gi>numel(st) || isempty(st(gi)) || isempty(st(gi).rel), continue; end
            s = st(gi);
            col = groupColor{gi};

            % choose hemisphere series
            if hemi==1
                m = s.mR; 
                se = s.sR; 
                ylab = 'Right mPFC (ΔF/F, baseline-subtracted)';
            else       
                m = s.mL; 
                se = s.sL; 
                ylab = 'Left mPFC (ΔF/F, baseline-subtracted)';
            end

            % x for this segment: shift so the segment starts at 'offset'
            [rel_sorted, ord] = sort(s.rel);
            m  = m(ord);  
            se = se(ord);
            x  = offset + (rel_sorted - minRel);

            % SEM ribbon
            xf = [x; flipud(x)];
            yf = [m+se; flipud(m-se)];
            fill(xf, yf, col, 'FaceAlpha', 0.2, 'EdgeColor','none', 'HandleVisibility','off');

            % mean line
            plot(x, m, '-o', 'Color', col, 'LineWidth', 2.5, 'MarkerSize', 6, ...
                 'MarkerFaceColor', col, 'DisplayName', labels{gi});
        end

        % learning-day marker within this segment (relative day == 0)
        % Only draw if day 0 has sufficient animals in at least one group
        day_zero_exists = false;
        for gi = 1:numel(st)
            if ~isempty(st(gi)) && any(st(gi).rel == 0)
                day_zero_exists = true;
                break;
            end
        end
        
        if day_zero_exists
            x0 = offset + (0 - minRel);
            xline(x0,'k--','LineWidth', 1.5, 'HandleVisibility','off');
            text(x0, ymin + 0.1*(ymax-ymin), 'Learning', ...
                 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                 'FontSize', 9, 'FontWeight', 'bold');
        end

        % protocol label centered over the segment
        xc = offset + span/2;
        yl = ylim; 
        text(xc, yl(2), sprintf('Protocol %d', protocols(ip)), ...
             'HorizontalAlignment','center','VerticalAlignment','top', ...
             'FontSize', 11, 'FontWeight', 'bold', 'Interpreter','none', 'Color',[0 0 0 0.7]);

        % advance offset for next segment
        offset = offset + span + spacer;
    end

    ylim([ymin ymax]);
    xlabel('Concatenated Relative Day (Learning-Day Aligned)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ylab, 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location','best', 'FontSize', 11);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2);
end

% Print summary
fprintf('\n===== Summary =====\n');
fprintf('Baseline window: %.2f to %.2f s\n', baseline_start, baseline_end);
fprintf('Event window: %.2f to %.2f s\n', start_time, end_time);
fprintf('Minimum animals per day: %d\n', min_animals_per_day);

for ip = 1:numel(protocols)
    fprintf('\nProtocol %d:\n', protocols(ip));
    st = statsP{ip};
    if isempty(st), continue; end
    
    for g = [0 1]
        gi = g+1;
        if gi > numel(st) || isempty(st(gi)), continue; end
        
        fprintf('  %s: %d days plotted\n', labels{gi}, numel(st(gi).rel));
        if isfield(st(gi), 'n_animals')
            fprintf('    Animals per day: min=%d, max=%d, mean=%.1f\n', ...
                    min(st(gi).n_animals), max(st(gi).n_animals), mean(st(gi).n_animals));
        end
    end
end

%% Plots summary plot of ROIs stage 1 aligned to learning day and stage 2 to first day 

% data to index:  V×T×Ndays_total
aligned_all = rewarded_stim_kernel;

% --- 1) Choose time window and extract ROI traces (Ndays × Twin) ---
start_time = 0;  
end_time = 0.35;

% Baseline window for subtraction
baseline_start = -0.1;
baseline_end = 0;

% added_time;
event_time_window_mask = find(added_time_Kernel >= start_time & added_time_Kernel <= end_time);
baseline_window_mask = find(added_time_Kernel >= baseline_start & added_time_Kernel <= baseline_end);

% Define ROIs
ROI_to_plot_1=right_mPFC_ROI_mask;
ROI_to_plot_2=left_mPFC_ROI_mask;

% Extract ROI traces for event window
right_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,event_time_window_mask,:), [], [], ROI_to_plot_1), ...
    [3,2,1]);  % Ndays × Twin

left_mPFC_ROI_trace  = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,event_time_window_mask,:), [], [], ROI_to_plot_2), ...
    [3,2,1]);  % Ndays × Twin

% Extract ROI traces for baseline window
right_mPFC_baseline = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,baseline_window_mask,:), [], [], ROI_to_plot_1), ...
    [3,2,1]);  % Ndays × Tbaseline

left_mPFC_baseline  = permute( ...
    ap.wf_roi(wf_U(:,:,1:kernel_n_components), aligned_all(:,baseline_window_mask,:), [], [], ROI_to_plot_2), ...
    [3,2,1]);  % Ndays × Tbaseline

% Compute baseline mean for each day
right_baseline_mean = mean(right_mPFC_baseline, 2, 'omitnan');  % Ndays × 1
left_baseline_mean = mean(left_mPFC_baseline, 2, 'omitnan');    % Ndays × 1

% Subtract baseline from event window
right_mPFC_ROI_trace = right_mPFC_ROI_trace - right_baseline_mean;
left_mPFC_ROI_trace = left_mPFC_ROI_trace - left_baseline_mean;

protocols = unique(workflow_cat(:))';          % plot segments in this order
groupColor = {cLearner, cNonLearner};          % 0->Learners, 1->Non-learners
labels     = {'Learners','Non-learners'};

% Minimum number of animals required per day
min_animals_per_day = 3;

% ---------- per-PROTOCOL stats (mean ± SEM vs relative day for each group) ----------
statsP = cell(numel(protocols),1);   % statsP{p}.G(g+1).rel/mR/sR/mL/sL

for ip = 1:numel(protocols)
    p = protocols(ip);

    % Store per-animal data with animal IDs for counting
    rel_all = []; 
    grp_all = []; 
    R_all = []; 
    L_all = [];
    animal_id_all = [];  % Track which animal each data point comes from
    
    for ai = unique(widefield_animal_idx(:))'
        sel = (widefield_animal_idx==ai) & (workflow_cat==p);
        days_idx = find(sel);
        if isempty(days_idx), continue; end

        grp = is_group_animal(days_idx(1));
        
        % MODIFIED: Different alignment strategy based on protocol
        if ip == 1  % Protocol 1: align to learning day
            ld = find(learning_index_animal(days_idx)==1, 1);
            if isempty(ld), continue; end  % Skip animals without learning day in Protocol 1
        else  % Protocol 2+: align to first day
            ld = 1;  % First day is the reference point
        end

            if p==3 % skip stage 3 for now
                continue;
            end

        n = numel(days_idx);
        rel_day = (1:n) - ld;

        meanR = max(right_mPFC_ROI_trace(days_idx,:),[], 2, 'omitnan');
        meanL = max(left_mPFC_ROI_trace (days_idx,:), [],2, 'omitnan');

        rel_all = [rel_all; rel_day(:)];
        grp_all = [grp_all; repmat(grp,n,1)];
        R_all   = [R_all;   meanR(:)];
        L_all   = [L_all;   meanL(:)];
        animal_id_all = [animal_id_all; repmat(ai, n, 1)];
    end

    Gpresent = unique(grp_all(:))';
    stats = struct([]);
    
    for g = Gpresent
        mask = (grp_all==g);
        r    = rel_all(mask);
        Rg   = R_all(mask);
        Lg   = L_all(mask);
        animal_ids = animal_id_all(mask);

        uRel = unique(r,'sorted');
        mR = nan(size(uRel)); 
        sR = nan(size(uRel));
        mL = nan(size(uRel)); 
        sL = nan(size(uRel));
        n_animals = nan(size(uRel));
        
        for k = 1:numel(uRel)
            idx = (r==uRel(k));
            Rk = Rg(idx); 
            Lk = Lg(idx); 
            animals_this_day = animal_ids(idx);
            
            % Remove NaN values
            valid_R = ~isnan(Rk);
            valid_L = ~isnan(Lk);
            Rk = Rk(valid_R);
            Lk = Lk(valid_L);
            
            % Count unique animals for this day
            n_animals(k) = numel(unique(animals_this_day));
            
            % Only compute statistics if we have minimum number of animals
            if n_animals(k) >= min_animals_per_day
                mR(k) = mean(Rk); 
                sR(k) = std(Rk)/max(1,sqrt(numel(Rk)));
                mL(k) = mean(Lk); 
                sL(k) = std(Lk)/max(1,sqrt(numel(Lk)));
            else
                mR(k) = NaN;
                sR(k) = NaN;
                mL(k) = NaN;
                sL(k) = NaN;
            end
        end
        
        % Filter out days with insufficient animals
        valid_days = n_animals >= min_animals_per_day;
        
        stats(g+1).rel = uRel(valid_days);
        stats(g+1).mR  = mR(valid_days); 
        stats(g+1).sR  = sR(valid_days);
        stats(g+1).mL  = mL(valid_days); 
        stats(g+1).sL  = sL(valid_days);
        stats(g+1).n_animals = n_animals(valid_days);
    end
    statsP{ip} = stats;
end

% ---------- global y-limits across all protocols & groups ----------
yvals = [];
for ip = 1:numel(protocols)
    st = statsP{ip}; if isempty(st), continue; end
    for gi = 1:numel(st)
        if isempty(st(gi)), continue; end
        yvals = [yvals; st(gi).mR+st(gi).sR; st(gi).mR-st(gi).sR; ...
                        st(gi).mL+st(gi).sL; st(gi).mL-st(gi).sL];
    end
end
yvals = yvals(~isnan(yvals));
ymin = max(-3, min(yvals));  
ymax = max(yvals);

% ---------- PLOT: concatenate protocol segments ----------
figure('Color','w','Position', [100, 100, 1400, 800], ...
       'Name', sprintf('Summary by protocol (P1=learning-aligned, P2+=first-day-aligned, min %d animals)', min_animals_per_day));
tiledlayout(2,1,'TileSpacing','compact');

spacer = 1;

for hemi = 1:2
    nexttile; hold on;
    offset = 0;

    for ip = 1:numel(protocols)
        st = statsP{ip};
        if isempty(st), offset = offset + spacer; continue; end

        % Span within this protocol
        minRel = inf; maxRel = -inf;
        for gi = 1:numel(st)
            if isempty(st(gi)) || isempty(st(gi).rel), continue; end
            minRel = min(minRel, min(st(gi).rel));
            maxRel = max(maxRel, max(st(gi).rel));
        end
        if ~isfinite(minRel) || ~isfinite(maxRel)
            offset = offset + spacer; continue;
        end
        span = maxRel - minRel;

        % Draw each group
        for g = [0 1]
            gi = g+1;
            if gi>numel(st) || isempty(st(gi)) || isempty(st(gi).rel), continue; end
            s = st(gi);
            col = groupColor{gi};

            if hemi==1
                m = s.mR; 
                se = s.sR; 
                ylab = 'Right mPFC (ΔF/F, baseline-subtracted)';
            else       
                m = s.mL; 
                se = s.sL; 
                ylab = 'Left mPFC (ΔF/F, baseline-subtracted)';
            end

            [rel_sorted, ord] = sort(s.rel);
            m  = m(ord);  
            se = se(ord);
            x  = offset + (rel_sorted - minRel);

            % SEM ribbon
            xf = [x; flipud(x)];
            yf = [m+se; flipud(m-se)];
            fill(xf, yf, col, 'FaceAlpha', 0.2, 'EdgeColor','none', 'HandleVisibility','off');

            % Mean line
            plot(x, m, '-o', 'Color', col, 'LineWidth', 2.5, 'MarkerSize', 6, ...
                 'MarkerFaceColor', col, 'DisplayName', labels{gi});
        end

        % Reference day marker (day 0)
        day_zero_exists = false;
        for gi = 1:numel(st)
            if ~isempty(st(gi)) && any(st(gi).rel == 0)
                day_zero_exists = true;
                break;
            end
        end
        
        if day_zero_exists
            x0 = offset + (0 - minRel);
            xline(x0,'k--','LineWidth', 1.5, 'HandleVisibility','off');
            
            % MODIFIED: Different label based on protocol
            if ip == 1
                label_text = 'Learning';
            else
                label_text = 'First Day';
            end
            
            text(x0, ymin + 0.1*(ymax-ymin), label_text, ...
                 'HorizontalAlignment','center', 'VerticalAlignment','bottom', ...
                 'FontSize', 9, 'FontWeight', 'bold');
        end

        % Protocol label
        xc = offset + span/2;
        yl = ylim; 
        text(xc, yl(2), sprintf('Stage %d', protocols(ip)), ...
             'HorizontalAlignment','center','VerticalAlignment','top', ...
             'FontSize', 11, 'FontWeight', 'bold', 'Interpreter','none', 'Color',[0 0 0 0.7]);

        offset = offset + span + spacer;
    end

    ylim([ymin ymax]);
    xlabel('Concatenated Relative Day (P1: Learning-Aligned, P2+: First-Day-Aligned)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ylab, 'FontSize', 12, 'FontWeight', 'bold');
    legend('Location','best', 'FontSize', 11);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2);
end


%% Plot a single animal

% ===== QUICK SINGLE ANIMAL PLOT WITH SEM =====

% Choose which animal to plot (by index or name)
target_animal_idx = 8;  % Change this to the animal index you want
% OR use: target_animal_idx = find(strcmp(animal_list, 'HA005'));  % by name

% Extract data for this animal
animal_mask = (widefield_animal_idx == target_animal_idx);
animal_workflows = workflow_cat(animal_mask);
animal_learning = learning_index_animal(animal_mask);
animal_days_idx = find(animal_mask);

% Get traces for this animal - MEAN and SEM across time window
animal_right_mean = max(right_mPFC_ROI_trace(animal_days_idx, :),[],2, 'omitnan');
animal_right_sem = std(right_mPFC_ROI_trace(animal_days_idx, :), 0, 2, 'omitnan') ./ ...
                   sqrt(sum(~isnan(right_mPFC_ROI_trace(animal_days_idx, :)), 2));

animal_left_mean = max(left_mPFC_ROI_trace(animal_days_idx, :),[],2, 'omitnan');
animal_left_sem = std(left_mPFC_ROI_trace(animal_days_idx, :), 0, 2, 'omitnan') ./ ...
                  sqrt(sum(~isnan(left_mPFC_ROI_trace(animal_days_idx, :)), 2));

global_y_lim_min= min([animal_left_mean;animal_right_mean],[],1);
global_y_lim_max= max([animal_left_mean;animal_right_mean],[],1);

% Determine if learner or non-learner
isLearner = is_group_animal(animal_days_idx(1));
if isLearner == 0
    animal_color = cLearner;
    group_label = 'mPFC+ (Learner)';
else
    animal_color = cNonLearner;
    group_label = 'mPFC- (Non-Learner)';
end

% Get unique protocols for this animal
unique_protocols = unique(animal_workflows);

% Create figure
figure('Color', 'w', 'Position', [100, 100, 1400, 800], ...
       'Name', sprintf('Animal %s - Individual Trajectory with SEM', animal_list{target_animal_idx}));
tiledlayout(2, 1, 'TileSpacing', 'compact');

for hemi = 1:2
    nexttile; hold on;
    
    offset = 0;
    spacer = 2;
    
    for ip = 1:numel(unique_protocols)
        p = unique_protocols(ip);
        
        % Get days for this protocol
        protocol_mask = (animal_workflows == p);
        protocol_days_local = find(protocol_mask);
        
        if isempty(protocol_days_local)
            offset = offset + spacer;
            continue;
        end
        
        % Find learning day within this protocol
        ld = find(animal_learning(protocol_mask) == 1, 1);
        
        if isempty(ld)
            % No learning - just use sequential days
            rel_days = 1:numel(protocol_days_local);
            x = offset + rel_days - 1;
        else
            % Has learning - align to learning day
            n = numel(protocol_days_local);
            rel_days = (1:n) - ld;
            x = offset + rel_days - min(rel_days);  % shift to start at offset
        end
        
        % Get activity for this protocol
        if hemi == 1
            activity_mean = animal_right_mean(protocol_days_local);
            activity_sem = animal_right_sem(protocol_days_local);
            ylab = 'Right mPFC (ΔF/F, baseline-subtracted)';
        else
            activity_mean = animal_left_mean(protocol_days_local);
            activity_sem = animal_left_sem(protocol_days_local);
            ylab = 'Left mPFC (ΔF/F, baseline-subtracted)';
        end
        
        % Plot SEM ribbon
        xf = [x(:); flipud(x(:))];
        yf = [activity_mean + activity_sem; flipud(activity_mean - activity_sem)];
        fill(xf, yf, animal_color, 'FaceAlpha', 0.25, 'EdgeColor', 'none', ...
             'HandleVisibility', 'off');
        
        % Plot mean line
        plot(x, activity_mean, '-o', 'Color', animal_color, 'LineWidth', 2.5, ...
             'MarkerSize', 8, 'MarkerFaceColor', animal_color, ...
             'DisplayName', sprintf('Protocol %d', p));
        % 
        % % Mark learning day if exists
        % if ~isempty(ld)
        %     learning_x = offset + (0 - min(rel_days));
        %     plot(learning_x, activity_mean(ld), 'p', 'MarkerSize', 15, ...
        %          'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'w', 'LineWidth', 2, ...
        %          'HandleVisibility', 'off');
        % end
        
        % Protocol label
        xc = offset + (max(x) - offset) / 2;
        yl = [global_y_lim_min;global_y_lim_max];        
        text(xc, yl(2), sprintf('Protocol %d', p), ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
             'FontSize', 11, 'FontWeight', 'bold', 'Color', [0 0 0 0.7]);
        
        % Update offset for next protocol
        span = max(x) - offset;
        offset = offset + span + spacer;
    end
    

    xlabel('Day within Protocol (Learning-Day Aligned)', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(ylab, 'FontSize', 12, 'FontWeight', 'bold');
    ylim(yl)
    legend('Location', 'best', 'FontSize', 10);
    set(gca, 'FontSize', 10, 'LineWidth', 1.2);
    hold off;
end

sgtitle(sprintf('Animal %s (%s) - Individual Activity (Mean ± SEM)', ...
        animal_list{target_animal_idx}, group_label), ...
        'FontSize', 15, 'FontWeight', 'bold');

% Print summary with noise assessment
fprintf('\n===== Animal %s Summary =====\n', animal_list{target_animal_idx});
fprintf('Group: %s\n', group_label);
fprintf('Total days recorded: %d\n', sum(animal_mask));
fprintf('Learning days: %d\n', sum(animal_learning));
fprintf('Protocols: %s\n', mat2str(unique_protocols'));

fprintf('\n===== Signal-to-Noise Assessment =====\n');
fprintf('Right mPFC:\n');
fprintf('  Mean activity range: %.3f to %.3f\n', min(animal_right_mean), max(animal_right_mean));
fprintf('  Mean SEM: %.3f (avg noise level)\n', mean(animal_right_sem));
fprintf('  Signal-to-noise ratio: %.2f\n', range(animal_right_mean) / mean(animal_right_sem));

fprintf('\nLeft mPFC:\n');
fprintf('  Mean activity range: %.3f to %.3f\n', min(animal_left_mean), max(animal_left_mean));
fprintf('  Mean SEM: %.3f (avg noise level)\n', mean(animal_left_sem));
fprintf('  Signal-to-noise ratio: %.2f\n', range(animal_left_mean) / mean(animal_left_sem));

for ip = 1:numel(unique_protocols)
    p = unique_protocols(ip);
    protocol_mask = (animal_workflows == p);
    n_days = sum(protocol_mask);
    has_learning = any(animal_learning(protocol_mask));
    
    protocol_days_local = find(protocol_mask);
    
    fprintf('\n  Protocol %d: %d days, Learning: %s\n', ...
            p, n_days, string(has_learning));
    fprintf('    Right mPFC SEM range: %.3f to %.3f\n', ...
            min(animal_right_sem(protocol_days_local)), ...
            max(animal_right_sem(protocol_days_local)));
    fprintf('    Left mPFC SEM range: %.3f to %.3f\n', ...
            min(animal_left_sem(protocol_days_local)), ...
            max(animal_left_sem(protocol_days_local)));
end

%% Plot CCF map of Kernels of both AP014/15 and Stage 3 of full traninig data just post-traninig for now

% Define the time window in seconds (0 to 300ms)
start_time = 0; 
end_time = 0.35; 

% Choose 'max' or 'mean' to plot
statistic = 'max';

added_time_Kernel = fliplr((-10:52)/30);

% Find the corresponding indices in added_time for the time window
time_window_idx = find(added_time_Kernel >= start_time & added_time_Kernel <= end_time);

% ===== FULL TRAINING DATA (Original) =====
% post of right move workflow split by learning days and learners vs non-learners
mean_post_learning_learners = full_training_last_3_wf;

% ===== AP014_15 DATA (Stage 3 Only) =====
% Extract last 3-4 days as "post-learning" equivalent
mean_post_learning_AP014_15 = AP014_AP015_last_3_wf;

% ===== CONVERT TO PIXEL SPACE =====
% Full training - Learners
timewindow_mean_post_learning_learners = plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), ...
    mean_post_learning_learners(:,time_window_idx));


% AP014_15 - Stage 3 Only
timewindow_mean_post_learning_AP014_15 = plab.wf.svd2px(wf_U(:,:,1:kernel_n_components), ...
    mean_post_learning_AP014_15(:,time_window_idx));

% ===== COMPUTE IMAGES =====
switch statistic
    case 'max'
        % Full training
        img_post_learning_learners = max(timewindow_mean_post_learning_learners, [], 3);
        
        % AP014_15
        img_post_learning_AP014_15 = max(timewindow_mean_post_learning_AP014_15, [], 3);
        
    case 'mean'
        % Full training
        img_post_learning_learners = mean(timewindow_mean_post_learning_learners, 3);
        
        % AP014_15
        img_post_learning_AP014_15 = mean(timewindow_mean_post_learning_AP014_15, 3);
        
    otherwise
        error('statistic must be ''max'' or ''mean''');
end

% ===== COMPUTE COMMON COLOR LIMIT ACROSS ALL GROUPS =====
all_vals = [
            timewindow_mean_post_learning_learners(:);
            timewindow_mean_post_learning_AP014_15(:)];

switch statistic
    case 'max'
        ref_vals = max(reshape(all_vals,[],size(timewindow_mean_post_learning_learners,3)),[],2);
    case 'mean'
        ref_vals = mean(reshape(all_vals,[],size(timewindow_mean_post_learning_learners,3)),2);
end
clim_val = 0.5 * max(abs(ref_vals));

% ===== PLOT IN 3×2 LAYOUT =====
figure('Position', [100, 100, 1600, 1200], 'Color', 'w', ...
       'Name', 'Pre vs Post Learning: All Groups');

% Row 1: Full Training mPFC+ (Learners)

subplot(3, 2, 2);
imagesc(img_post_learning_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Post-Learning %s\nmPFC+ (Full Training 1→2→3)', statistic), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', cLearner);
ap.wf_draw('ccf');


% Row 3: Stage 3 Only (AP014_15)
subplot(3, 2, 6);
imagesc(img_post_learning_AP014_15);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Late Days %s\nStage 3 Only (AP014, AP015)', statistic), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.7, 0.4, 0.9]);
ap.wf_draw('ccf');

% Add shared colorbar
cb = colorbar('Position', [0.92 0.3 0.015 0.4]);
ylabel(cb, 'ΔF/F', 'FontSize', 13, 'FontWeight', 'bold');

% Overall title
sgtitle(sprintf('Widefield Activity Comparison (Time: %.0f-%.0f ms, %s)', ...
                start_time*1000, end_time*1000, statistic), ...
        'FontSize', 16, 'FontWeight', 'bold');

% ===== CREATE DIFFERENCE MAPS =====
figure('Position', [150, 150, 1600, 600], 'Color', 'w', ...
       'Name', 'Difference Maps: Post - Pre');

% Compute differences
diff_learners = img_post_learning_learners -img_post_learning_AP014_15 ;
% diff_AP014_15 = img_post_learning_AP014_15 - img_pre_learning_AP014_15;

% Common color limit for differences
all_diffs = [diff_learners(:)];
diff_clim_val = 0.8 * max(abs(all_diffs));

% Plot differences
subplot(1, 3, 1);
imagesc(diff_learners);
axis image off;
colormap(ap.colormap('PWG'));
clim([-clim_val, clim_val]);
title(sprintf('Difference map between mPFC+ (Full Training) & Only Stage 3 \nPost - Pre Learning'), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', 'black');
ap.wf_draw('ccf');

subplot(1, 3, 2);
imagesc(diff_non_learners);
axis image off;
colormap(ap.colormap('BWR'));
clim([-diff_clim_val, diff_clim_val]);
title(sprintf('mPFC- (Full Training)\nPost - Pre Learning'), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', cNonLearner);
ap.wf_draw('ccf');

subplot(1, 3, 3);
imagesc(diff_AP014_15);
axis image off;
colormap(ap.colormap('BWR'));
clim([-diff_clim_val, diff_clim_val]);
title(sprintf('Stage 3 Only (AP014, AP015)\nLate - Early'), ...
      'FontSize', 13, 'FontWeight', 'bold', 'Color', [0.7, 0.4, 0.9]);
ap.wf_draw('ccf');

% Add shared colorbar
cb = colorbar('Position', [0.92 0.3 0.02 0.4]);
ylabel(cb, 'Δ Activity (Post - Pre)', 'FontSize', 13, 'FontWeight', 'bold');

sgtitle('Learning-Related Changes in Widefield Activity', ...
        'FontSize', 16, 'FontWeight', 'bold');


%% Compare CCF maps of mPFC+ group compared to AP014/AP015 (Diss group) only for Vs for now


load('C:\Users\havgana\Desktop\DPhil\packaged_data\AP_014_015_task_widefield_2.mat'); % load AP014+15 task
task_data_AP014_AP015_widefield

% get all the widefield data
widefield_cat_AP014_15= cat(2,task_data_AP014_AP015_widefield.widefield); % concat

% for V - create a V x T x all days variable - rewarded stim
rewarded_stim_v_AP014_15 = cellfun(@(x) mean(x, 3), {widefield_cat_AP014_15.rewarded_stim_on_aligned_V}, 'UniformOutput', false); % average across trials
rewarded_stim_v_stacked_data_AP014_15 = cat(3, rewarded_stim_v_AP014_15{:});

% for Kernel
rewarded_stim_AP014_015_kernel= cat(3,widefield_cat_AP014_15.rewarded_stim_on_aligned_kernel);

% Creates an animal index (n all days x animal number ordered)
widefield_animal_idx_AP014_15 = grp2idx(cell2mat(cellfun(@(animal,wf) repmat(animal,length(wf),1), ...
    {task_data_AP014_AP015_widefield.animal},{task_data_AP014_AP015_widefield.widefield},'uni',false)'));

% Find the last 3 days for each animal
AP014_AP015_last_3_days= [find(widefield_animal_idx_AP014_15==1,4,'last');find(widefield_animal_idx_AP014_15==2,4,'last')];

% Get the widefield activity for those days
AP014_AP015_last_3_wf= nanmean(rewarded_stim_AP014_015_kernel(:,:,AP014_AP015_last_3_days),3);

% Get widefield activity for full training group (last 3 days)
full_training_last_3_wf = nanmean(rewarded_stim_kernel(:,:,full_training_last_3_idx), 3);




% pre post of a selected workflow
pre = nanmean(rewarded_stim_v_stacked_data(:,:,workflow_cat==3 & learning_index_animal==0 & is_group_animal==0),3);
post = nanmean(rewarded_stim_v_stacked_data(:,:,workflow_cat==3 & learning_index_animal==1 & is_group_animal==0),3);

ap.imscroll([plab.wf.svd2px(wf_U(:,:,1:kernel_n_components),AP014_AP015_last_3_wf)],added_time_Kernel);
clim([-max(abs(clim)), max(abs(clim))]);
colormap(ap.colormap( ...
    'PWG'));
axis image;


%%
% ===== Compare Stage Progression: Full Training (1→2→3) vs Stage 3 Only =====

% Get widefield data for both groups
% Group 1: Animals that went through stages 1→2→3 (from your main data)
full_training_animals = widefield_animal_idx(is_group_animal==0);  % Your original learners
full_training_mask = ~ismember(string({behaviour_data.animal_id}), group_animals);
unique_animals= unique(widefield_animal_idx);

% Extract last 3 days for full training group
full_training_last_3_idx = [];
for ai = 1:numel(behaviour_data)
    if ~full_training_mask(ai), continue; end
    
    animal_mask = (is_group_animal==0 & widefield_animal_idx==ai & workflow_cat==3);
    animal_days = find(animal_mask);
    
    if length(animal_days) >= 3
        full_training_last_3_idx = [full_training_last_3_idx; animal_days(end-2:end)];
    elseif ~isempty(animal_days)
        full_training_last_3_idx = [full_training_last_3_idx; animal_days];
    end
end


% Get widefield activity for full training group (last 3 days)
full_training_last_3_wf = nanmean(rewarded_stim_v_stacked_data(:,:,full_training_last_3_idx), 3);

% Ensure we have the stage 3 only data
stage3_only_last_3_wf = AP014_AP015_last_3_wf;

% ===== Prepare for Video Creation =====

output_dir = 'CCF_videos_comparison';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Video parameters
frame_rate = 10;
gif_delay = 0.1;

% Prepare data groups
group_data = {full_training_last_3_wf, stage3_only_last_3_wf};
group_names = {'Full Training (Stages 1→2→3)', 'Stage 3 Only'};
group_colors = {cLearner, [0.7, 0.4, 0.9]};  % Purple for stage 3 only

% Convert to pixel space
group_data_px = cell(1, length(group_data));
for g = 1:length(group_data)
    group_data_px{g} = plab.wf.svd2px(wf_U(:,:,1:n_components), group_data{g});
end

% Determine common color limits
all_data = cat(3, group_data_px{:});
clim_max = max(abs(all_data(:)));
clim_range = [-clim_max, clim_max];

fprintf('\n===== Creating CCF Videos: Training Comparison =====\n');
fprintf('Full Training (1→2→3): n=%d days (last 3 from %d animals)\n', ...
        length(full_training_last_3_idx), sum(full_training_mask));
fprintf('Stage 3 Only: n=8 days (last 4 from 2 animals: AP014, AP015)\n');
fprintf('Color range: [%.3f, %.3f]\n', clim_range(1), clim_range(2));

% ===== Create Individual Videos for Each Group =====

for g = 1:length(group_data)
    
    fprintf('\nProcessing: %s\n', group_names{g});
    
    % Get data
    data_px = group_data_px{g};
    [img_height, img_width, n_frames] = size(data_px);
    
    % Output filename
    video_filename = fullfile(output_dir, sprintf('%s_last3days.mp4', ...
        strrep(group_names{g}, ' ', '_')));
    
    % Create Video Writer
    video_writer = VideoWriter(video_filename, 'MPEG-4');
    video_writer.FrameRate = frame_rate;
    video_writer.Quality = 95;
    open(video_writer);
    
    % Create figure for rendering
    fig = figure('Position', [100 100 900 700], 'Color', 'w', 'Visible', 'off');
    
    for frame_idx =1:n_frames % Goes backwards for kernel data
        
        % Clear figure
        clf(fig);
        
        % Display frame
        imagesc(data_px(:, :, frame_idx));
        axis image off;
        clim(clim_range);
        colormap(ap.colormap('PWG'));
        ap.wf_draw('ccf');
        
        % Add colorbar
        cb = colorbar;
        cb.Label.String = 'ΔF/F';
        cb.Label.FontSize = 14;
        cb.Label.FontWeight = 'bold';
        
        % Add time indicator
        current_time = added_time(frame_idx);
        title(sprintf('%s (Last 3 Days)\nTime: %.2f s', group_names{g}, current_time), ...
            'FontSize', 16, 'FontWeight', 'bold', 'Color', group_colors{g});
        
        % Highlight CS+ onset (t = 0)
        if abs(current_time) < 0.05
            hold on;
            rectangle('Position', [0.5, 0.5, img_width, img_height], ...
                'EdgeColor', 'r', 'LineWidth', 5);
            text(img_width/2, 10, 'CS+ ONSET', ...
                'FontSize', 18, 'FontWeight', 'bold', 'Color', 'r', ...
                'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.7]);
        end
        
        % % Add time progress bar
        % hold on;
        % progress = (n_frames - frame_idx + 1) / n_frames;
        % bar_width = img_width * 0.8;
        % bar_x = img_width * 0.1;
        % bar_y = img_height - 10;
        % rectangle('Position', [bar_x, bar_y, bar_width, 5], ...
        %     'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none');
        % rectangle('Position', [bar_x, bar_y, bar_width * progress, 5], ...
        %     'FaceColor', group_colors{g}, 'EdgeColor', 'none');
        % 
        drawnow;
        
        % Capture frame for video
        frame_img = getframe(fig);
        writeVideo(video_writer, frame_img);
    end
    
    % Close video writer
    close(video_writer);
    close(fig);
    
    fprintf('✓ Saved video: %s\n', video_filename);
end

% ===== Create Side-by-Side Comparison Video =====

fprintf('\nCreating side-by-side comparison video...\n');

comparison_filename = fullfile(output_dir, 'Training_Comparison_SideBySide.mp4');

video_writer = VideoWriter(comparison_filename, 'MPEG-4');
video_writer.FrameRate = frame_rate;
video_writer.Quality = 95;
open(video_writer);

% Create larger figure for side-by-side
fig = figure('Position', [100 100 1800 700], 'Color', 'w', 'Visible', 'off');

n_frames = size(group_data_px{1}, 3);

for frame_idx=1:n_frames
    
    clf(fig);
    
    % Left panel: Full Training
    subplot(1, 2, 1);
    imagesc(group_data_px{1}(:, :, frame_idx));
    axis image off;
    clim(clim_range);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf');
    
    current_time = added_time(frame_idx);
    title(sprintf('%s\nTime: %.2f s', group_names{1}, current_time), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', group_colors{1});
    
    % Highlight CS+ onset
    if abs(current_time) < 0.05
        hold on;
        [img_h, img_w] = size(group_data_px{1}(:,:,1));
        rectangle('Position', [0.5, 0.5, img_w, img_h], ...
            'EdgeColor', 'r', 'LineWidth', 4);
    end
    
    % Right panel: Stage 3 Only
    subplot(1, 2, 2);
    imagesc(group_data_px{2}(:, :, frame_idx));
    axis image off;
    clim(clim_range);
    colormap(ap.colormap('PWG'));
    ap.wf_draw('ccf');
    
    title(sprintf('%s\nTime: %.2f s', group_names{2}, current_time), ...
        'FontSize', 14, 'FontWeight', 'bold', 'Color', group_colors{2});
    
    % Highlight CS+ onset
    if abs(current_time) < 0.05
        hold on;
        [img_h, img_w] = size(group_data_px{2}(:,:,1));
        rectangle('Position', [0.5, 0.5, img_w, img_h], ...
            'EdgeColor', 'r', 'LineWidth', 4);
    end
    
    % Add shared colorbar
    cb = colorbar('Position', [0.92 0.3 0.02 0.4]);
    cb.Label.String = 'ΔF/F';
    cb.Label.FontSize = 14;
    cb.Label.FontWeight = 'bold';
    
    % Overall title
    sgtitle(sprintf('Training Comparison: Last 3 Days (Frame %d/%d)', ...
                    n_frames - frame_idx + 1, n_frames), ...
            'FontSize', 18, 'FontWeight', 'bold');
    
    drawnow;
    
    % Capture frame
    frame_img = getframe(fig);
    writeVideo(video_writer, frame_img);
end

close(video_writer);
close(fig);

fprintf('✓ Saved comparison video: %s\n', comparison_filename);

% ===== Create Difference Map Video =====

fprintf('\nCreating difference map video...\n');

difference_filename = fullfile(output_dir, 'Training_Difference_Map.mp4');

% Calculate difference: Full Training - Stage 3 Only
difference_data = group_data_px{1} - group_data_px{2};

% Set color limits for difference map
diff_clim_max = max(abs(difference_data(:)));
diff_clim_range = [-diff_clim_max, diff_clim_max];

video_writer = VideoWriter(difference_filename, 'MPEG-4');
video_writer.FrameRate = frame_rate;
video_writer.Quality = 95;
open(video_writer);

fig = figure('Position', [100 100 900 700], 'Color', 'w', 'Visible', 'on');

for frame_idx =1: n_frames
    
    clf(fig);
    
    imagesc(difference_data(:, :, frame_idx));
    axis image off;
    clim(diff_clim_range);
    colormap(ap.colormap('PWG'));  % Blue-white-red for difference
    ap.wf_draw('ccf');
    
    cb = colorbar;
    cb.Label.String = 'ΔActivity (Full - Stage3Only)';
    cb.Label.FontSize = 14;
    cb.Label.FontWeight = 'bold';
    
    current_time = added_time(frame_idx);
    title(sprintf('Difference Map: Full Training - Stage 3 Only\nTime: %.2f s', current_time), ...
        'FontSize', 16, 'FontWeight', 'bold');
    
    % Highlight CS+ onset
    if abs(current_time) < 0.05
        hold on;
        [img_h, img_w] = size(difference_data(:,:,1));
        rectangle('Position', [0.5, 0.5, img_w, img_h], ...
            'EdgeColor', 'k', 'LineWidth', 5);
        text(img_w/2, 10, 'CS+ ONSET', ...
            'FontSize', 18, 'FontWeight', 'bold', 'Color', 'k', ...
            'HorizontalAlignment', 'center', 'BackgroundColor', [1 1 1 0.8]);
    end
    
    drawnow;
    
    frame_img = getframe(fig);
    writeVideo(video_writer, frame_img);
end

close(video_writer);
close(fig);

fprintf('✓ Saved difference video: %s\n', difference_filename);

fprintf('\n===== All Videos Created Successfully =====\n');
fprintf('Output directory: %s\n', output_dir);

%% Plot Pre learning RTs vs Post learning mPFC

% define the variables to plot
data_to_plot = rewarded_stim_v_stacked_data;
ROI_to_plot = left_mPFC_ROI_mask;
time_to_plot = added_time;
curr_components = n_components;
protocol_idx = 1; % workflow to analyze

% Define event and baseline time windows
event_times = [0, 0.35];
baseline_time_windows = [-0.1, 0];

% Create logical masks for both
event_time_window_mask = (time_to_plot > event_times(1) & time_to_plot < event_times(2));
baseline_mask = (time_to_plot >= baseline_time_windows(1) & time_to_plot <= baseline_time_windows(2));

left_mPFC_ROI_trace = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,event_time_window_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % same dims

left_mPFC_baseline = permute( ...
    ap.wf_roi(wf_U(:,:,1:curr_components), data_to_plot(:,baseline_mask,:),[],[], ROI_to_plot), ...
    [3,2,1] ); % Ndays_total × T_baseline

left_mPFC_ROI_trace = left_mPFC_ROI_trace - mean(left_mPFC_baseline, 2); %substract baseline

% Find mean/max activity for each day
left_mPFC_ROI_trace_max = mean(left_mPFC_ROI_trace,2); % n×1


% 1. EXTRACT PRE-LEARNING RT (days -3 to -1) PER ANIMAL
pre_window = -3:-1;
post_window = 1:3;

protocol_idx= 1;
n_animals = length(animal_list);
animal_list_idx= grp2idx(animal_list);

pre_RT_animal = nan(n_animals, 1);
post_mPFC_animal = nan(n_animals, 1);
group_label = nan(n_animals, 1);


% 2. COLLECT RT DATA
for ai = 1:numel(animal_list)
    curr_animal = animal_list_idx(ai);
    curr_id= animal_list{ai};

    if curr_animal ==1
        continue;
    end

    sel = (widefield_animal_idx == curr_animal) & (workflow_cat == protocol_idx);
    days_idx = find(sel);

    if isempty(days_idx)
        continue;
    end

    % first learning day 
    ld = find(learning_index_animal(days_idx) == 1, 1, 'first');

    % Skip if animal hasn't learnt
    if isempty(ld)
        warning('Animal %s: no learning day found, skipping', curr_id);
        continue;
    end

    n_days = length(days_idx);
    rel_days = (1:n_days) - ld;
  

    local_day_idx= 1:length(days_idx);
     
    % Extract median RT for each valid day
    rts_per_day = cellfun(@(d) median(d.all_stim_diff_from_optimal_reward(d.cs_labels), 'omitnan'), ...
        num2cell(behaviour_data(ai).recording_day(local_day_idx)));

    % Average pre-learning RT
    pre_mask = ismember(rel_days, pre_window);
    if sum(pre_mask) > 0
        pre_RT_animal(ai) = mean(rts_per_day(pre_mask), 'omitnan');
    end
    
    % Group label
    group_label(ai) = ~ismember(curr_id, learner_ids);
end

% 3. COLLECT mPFC DATA (using widefield structure)

for ai = 1:numel(animal_list)

 curr_animal = animal_list_idx(ai);    
curr_id= animal_list{ai};
    % Select this animal & protocol
    sel = (widefield_animal_idx == curr_animal) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    if isempty(days_idx), continue; end
    
    % Find learning day
    learn_local = learning_index_animal(days_idx) == 1;
    ld = find(learn_local, 1, 'first');
    if isempty(ld)
        warning('Animal %s: no learning day found, skipping', curr_id);
        continue;
    end
    
    % Relative days
    n_days = length(days_idx);
    rel_days = (1:n_days) - ld;
    
    % mPFC activity
    mPFC_trace = left_mPFC_ROI_trace_max(days_idx, :);  % n_days × T
  
    
    % Average post-learning mPFC
    post_mask = ismember(rel_days, post_window);
    if sum(post_mask) > 0
        post_mPFC_animal(curr_animal) = mean(mPFC_trace(post_mask), 'omitnan');
    end
end

% 4. CORRELATION ANALYSIS
valid = ~isnan(pre_RT_animal) & ~isnan(post_mPFC_animal);
[r, p] = corr(pre_RT_animal(valid), post_mPFC_animal(valid), 'Type', 'Spearman');

% 5. PLOT
figure;
gscatter(pre_RT_animal(valid), post_mPFC_animal(valid), group_label(valid), [cLearner; cNonLearner], 'o', 8, 'filled');
xlabel('Pre-learning RT (median, days -3:-1)');
ylabel('Post-learning mPFC activity (mean, days 1:3)');
title(sprintf('ρ = %.3f, p = %.3f', r, p));
legend({'Learners', 'Non-learners'});




%% Creates both raw and time-warpped aligned to stim onset with heirarchcal SEM and mean calculations

mPFC_plus_grp = {'HA005','HA008','HA010','HA011','HA012'};
mPFC_minus_grp = {'DS017','HA007','HA009','HA014','HA015'};
non_learners_grp = {'HA006','HA013','AP030','AP031','AP032'};

target_workflow = 'visual_operant_lick_two_stim_right_move';
all_ids = {behaviour_data.animal_id};

n_norm_points = 100;
norm_time = linspace(-0.2, 1, n_norm_points);

% Define time window for raw traces (e.g., -0.5 to 3 seconds from stim onset)
raw_time_window = [-0.5, 3];
raw_time_bins = linspace(raw_time_window(1), raw_time_window(2), 200);

% Two data structures: warped and raw
animal_data_warped = cell(numel(all_ids), 1);
animal_data_raw = cell(numel(all_ids), 1);
animal_groups = nan(numel(all_ids), 1);

for ai = 1:numel(all_ids)
    animal_id = all_ids{ai};
    
    % Determine group
    if ismember(animal_id, mPFC_plus_grp)
        group_label = 1;
    elseif ismember(animal_id, mPFC_minus_grp)
        group_label = 2;
    elseif ismember(animal_id, non_learners_grp)
        group_label = 3;
    else
        continue;
    end
    
    animal_groups(ai) = group_label;
    
    % Get valid days
    days_all = behaviour_data(ai).recording_day;
    isValid = arrayfun(@(d) isfield(d,'workflow') && strcmp(d.workflow, target_workflow), days_all);
    validDays = days_all(isValid);
    
    if isempty(validDays)
        continue;
    end
    
    % Get learning status
    sig_days_all = combined_sig_day_all_protocols{ai};
    sigDays = sig_days_all(isValid);
    
    % Store day-level data
    animal_days_data_warped = {};
    animal_days_data_raw = {};
    
    % Process only PRE-learning days (you had sigDays(d), continue - which skips POST)
    for d = 1:length(validDays)
        % if ~sigDays(d), continue; end  % Skip post-learning (change to ~sigDays(d) for pre-learning)
        
        day_data = validDays(d);
        
        % Check required fields
        if ~isfield(day_data, 'all_static_times') || ~isfield(day_data, 'cs_labels')
            continue;
        end
        
        % Get stim onset
        stim_onset_V = task_data(ai).widefield(d).rewarded_stim_on_aligned_V;
        
        % Extract CS+ trials only
        cs_labels = day_data.cs_labels;
        cs_plus_static_times = day_data.all_static_times(cs_labels == 1);
        
        % Extract mPFC ROI
        mPFC_trace = permute( ...
            ap.wf_roi(wf_U(:,:,1:curr_components), stim_onset_V,[],[], left_mPFC_ROI_mask), ...
            [3,2,1]); % [n_trials × n_frames]
        
        
        % Time-warp trials for this day
        day_warped_trials = [];
        day_raw_trials = [];
        
        for trial_idx = 1:length(cs_plus_static_times)
            duration = cs_plus_static_times(trial_idx);
            
            if isnan(duration) || duration <= 0
                continue;
            end
            
            % Find frames within static period
            valid_frames = added_time >= -0.2 & added_time <= duration;
            
            if sum(valid_frames) < 2
                continue;
            end
            
            % Extract activity and time
            activity_segment = mPFC_trace(trial_idx, valid_frames);
            time_segment = added_time(valid_frames);
            
            % === WARPED VERSION ===
            normalized_time_trial = time_segment / duration;
            
            try
                warped_trial = interp1(normalized_time_trial, activity_segment, ...
                                       norm_time, 'linear', 'extrap');
                day_warped_trials = [day_warped_trials; warped_trial];
            catch
                continue;
            end
            
            % === RAW VERSION (aligned to stim onset) ===
            % Interpolate to fixed time bins
            try
                raw_trial = mPFC_trace(trial_idx, :);
                day_raw_trials = [day_raw_trials; raw_trial];
            catch
                continue;
            end
        end
        
        % Store this day's data (all trials)
        if ~isempty(day_warped_trials)
            animal_days_data_warped{end+1} = day_warped_trials;
            animal_days_data_raw{end+1} = day_raw_trials;
        end
    end
    
    animal_data_warped{ai} = animal_days_data_warped;
    animal_data_raw{ai} = animal_days_data_raw;
end


% Hierarchical aggregation: Trials → Days → Animals → Groups

group_names = {'mPFC+', 'mPFC-', 'Non-learners'};
group_colors = {cLearner, cNonLearner, [0.5 0.5 0.5]};

group_means = cell(3, 1);
group_sems = cell(3, 1);

for g = 1:3
    % Find animals in this group
    group_animals = find(animal_groups == g);
    
    if isempty(group_animals)
        continue;
    end
    
    % Storage for animal-level means
    animal_means = [];
    animal_sems = [];
    
    for ai = group_animals'
        if isempty(animal_data_warped{ai})
            continue;
        end
        
        % Step 1: Average across trials within each day
        day_means = [];
        
        for day_idx = 1:length(animal_data_warped{ai})
            day_trials = animal_data_warped{ai}{day_idx};
            
            if ~isempty(day_trials)
                % Mean across trials for this day
                day_mean = mean(day_trials, 1, 'omitnan');
                day_means = [day_means; day_mean];
            end
        end
        
        if isempty(day_means)
            continue;
        end

        % Step 2: Average across days for this animal
        animal_mean = mean(day_means, 1, 'omitnan');
        animal_sem = std(day_means, 0, 1, 'omitnan') / sqrt(size(day_means, 1));
        
        animal_means = [animal_means; animal_mean];
        animal_sems = [animal_sems; animal_sem];
    end
    
    if isempty(animal_means)
        continue;
    end
    
    % Step 3: Average across animals
    group_means{g} = mean(animal_means, 1, 'omitnan');
    
    % Step 4: Propagate SEM across animals
    pooled_var = mean(animal_sems.^2, 1, 'omitnan');
    group_sems{g} = sqrt(pooled_var) / sqrt(size(animal_means, 1));
end



% Plot with error bands

figure('Color', 'w', 'Position', [100, 100, 900, 600]);
hold on;

for g = 1:3  % Only plot mPFC+ and mPFC-
    if ~isempty(group_means{g})
        means = group_means{g};
        sems = group_sems{g};
        
        % Error band (SEM)
        fill([norm_time, fliplr(norm_time)], ...
             [means + sems, fliplr(means - sems)], ...
             group_colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        
        % Mean line
        plot(norm_time, means, '-', 'Color', group_colors{g}, ...
             'LineWidth', 2.5, 'DisplayName', group_names{g});
    end
end

% Mark events
xline(0, 'k--', 'Stim Onset', 'LineWidth', 1.5, 'HandleVisibility', 'off', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 11);
xline(1, 'k--', 'Movement', 'LineWidth', 1.5, 'HandleVisibility', 'off', ...
      'LabelVerticalAlignment', 'bottom', 'FontSize', 11);

xlabel('Normalized Time (0=Onset, 1=Movement)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity (ΔF/F)', 'FontSize', 13, 'FontWeight', 'bold');
title('Sustained mPFC Activity: Time-Warped Analysis (Post-Learning)', ...
      'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', 11);
set(gca, 'FontSize', 11, 'LineWidth', 1.2);



%% Extract mean sustained activity for first 2 days per animal


first_days_data = [];

% Define event and baseline time windows
event_times = [0, 0.5];

% Find the corresponding indices in added_time for the time window
time_window_idx = (added_time >= event_times(1) & added_time <= event_times(2));

for ai = 1:numel(all_ids)
    if isempty(animal_data_raw{ai}) || isnan(animal_groups(ai))
        continue;
    end
    
    group_label = animal_groups(ai);
    animal_id = all_ids{ai};
    
    % Get first 2 days (or fewer if not available)
    n_days = min(3, length(animal_data_raw{ai}));
    
    if n_days == 0
        continue;
    end
    
    
    day_means = [];
    
    for day_idx = 1:n_days
        day_trials = animal_data_raw{ai}{day_idx};
        
        if isempty(day_trials)
            continue;
        end
        
        % Average across trials and sustained window for this day
        day_mean = max(day_trials(:, time_window_idx), [],'all');
        day_means = [day_means; day_mean];
    end
    
    % Store: [day1_mean, day2_mean (or NaN), group, animal_id]
    if length(day_means) == 1
        first_days_data = [first_days_data; day_means(1), NaN, group_label, ai];
    else
        first_days_data = [first_days_data; day_means(1), day_means(2), group_label, ai];
    end
end

% Plot scatter: Day 1 vs Day 2

group_names = {'mPFC+', 'mPFC-', 'Non-learners'};
group_colors = {cLearner, cNonLearner, [0.5 0.5 0.5]};

figure('Color', 'w', 'Position', [100, 100, 800, 800]);

% Overall scatter
subplot(2, 2, 1);
hold on;

for g = 1:3
    group_mask = first_days_data(:, 3) == g;
    day1 = first_days_data(group_mask, 1);
    day2 = first_days_data(group_mask, 2);
    
    scatter(day1, day2, 100, group_colors{g}, 'filled', 'MarkerFaceAlpha', 0.6, ...
            'DisplayName', group_names{g});
end

% Unity line
plot([0, max(first_days_data(:, 1:2), [], 'all')], ...
     [0, max(first_days_data(:, 1:2), [], 'all')], 'k--', 'LineWidth', 1.5, ...
     'HandleVisibility', 'off');

xlabel('Day 1 Mean Activity', 'FontSize', 12);
ylabel('Day 2 Mean Activity', 'FontSize', 12);
title('Within-Animal Consistency', 'FontSize', 13, 'FontWeight', 'bold');
legend('Location', 'best');
grid on;
axis equal;


sgtitle('First 2 Days: Within-Animal Variability', 'FontSize', 14, 'FontWeight', 'bold');


%% Extract activity trajectory across all days per animal

max_days = 0;
for ai = 1:numel(all_ids)
    if ~isempty(animal_data_raw{ai})
        max_days = max(max_days, length(animal_data_raw{ai}));
    end
end

% Storage: [n_animals x max_days]
animal_trajectories = nan(numel(all_ids), max_days);
animal_group_labels = nan(numel(all_ids), 1);

for ai = 1:numel(all_ids)
    if isempty(animal_data_raw{ai}) || isnan(animal_groups(ai))
        continue;
    end
    
    animal_group_labels(ai) = animal_groups(ai);
    
    % Define sustained window
    sustained_start_idx = find(norm_time >= 0.1, 1, 'first');
    sustained_end_idx = find(norm_time <= 0.5, 1, 'last');
    sustained_window = sustained_start_idx:sustained_end_idx;
    
    % Extract mean activity for each day
    for day_idx = 1:length(animal_data_raw{ai})
        day_trials = animal_data_raw{ai}{day_idx};
        
        if ~isempty(day_trials)
            day_mean = mean(day_trials(:, sustained_window), 'all', 'omitnan');
            animal_trajectories(ai, day_idx) = day_mean;
        end
    end
end

% Plot trajectories per group

figure('Color', 'w', 'Position', [100, 100, 1400, 400]);

for g = 1:3
    subplot(1, 3, g);
    hold on;
    
    % Get animals in this group
    group_animals = find(animal_group_labels == g);
    
    % Plot each animal's trajectory
    for ai = group_animals'
        trajectory = animal_trajectories(ai, :);
        valid_days = ~isnan(trajectory);
        
        if sum(valid_days) > 0
            plot(find(valid_days), trajectory(valid_days), '-o', ...
                 'Color', group_colors{g}, 'LineWidth', 1.5, 'MarkerSize', 6, ...
                 'MarkerFaceColor', group_colors{g}, 'MarkerEdgeColor', 'none');
        end
    end
    
    % Group mean trajectory
    group_mean_trajectory = mean(animal_trajectories(group_animals, :), 1, 'omitnan');
    valid_days = ~isnan(group_mean_trajectory);
    
    plot(find(valid_days), group_mean_trajectory(valid_days), 'k-', ...
         'LineWidth', 3, 'DisplayName', 'Group Mean');
    
    xlabel('Day', 'FontSize', 12);
    ylabel('Mean mPFC Activity', 'FontSize', 12);
    title(sprintf('%s (n=%d animals)', group_names{g}, length(group_animals)), ...
          'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    xlim([0.5, max_days + 0.5]);
    
    % Add learning day marker if applicable
    xline(1.5, '--', 'First 2 Days', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5], ...
          'LabelVerticalAlignment', 'bottom');
end

sgtitle('Individual Animal Trajectories Across Days', 'FontSize', 14, 'FontWeight', 'bold');

%% Statistics: Within-group variability

fprintf('\n=== Within-Group Variability (First 2 Days) ===\n');

for g = 1:3
    group_animals = find(animal_group_labels == g);
    
    % Get first 2 days data for this group
    first_two_days = animal_trajectories(group_animals, 1:min(2, max_days));
    
    % Flatten to get all values
    all_values = first_two_days(:);
    all_values = all_values(~isnan(all_values));
    
    fprintf('%s:\n', group_names{g});
    fprintf('  Mean = %.6f\n', mean(all_values));
    fprintf('  SD = %.6f\n', std(all_values));
    fprintf('  CV = %.2f%%\n', 100 * std(all_values) / mean(all_values));
    fprintf('  Range = [%.6f, %.6f]\n', min(all_values), max(all_values));
    
    % Within-animal variability (day-to-day)
    within_animal_vars = nan(length(group_animals), 1);
    for i = 1:length(group_animals)
        animal_days = first_two_days(i, :);
        animal_days = animal_days(~isnan(animal_days));
        if length(animal_days) > 1
            within_animal_vars(i) = std(animal_days);
        end
    end
    
    fprintf('  Mean within-animal SD = %.6f\n', mean(within_animal_vars, 'omitnan'));
    fprintf('  Between-animal SD = %.6f\n', std(mean(first_two_days, 2, 'omitnan')));
    fprintf('\n');
end

%% Mixed effect analysis 

% Prepare data with diagnostics

all_data = [];
animal_ids_str = {};  % Track which animal is which


for ai = 1:numel(all_ids)
    if isempty(animal_data_raw{ai}) || isnan(animal_groups(ai))
        continue;
    end
    
    group_label = animal_groups(ai);
    animal_id = all_ids{ai};
    
    % Loop through days
    for day_idx = 1:length(animal_data_raw{ai})
        day_trials = animal_data_raw{ai}{day_idx};  % [n_trials x n_timepoints]
        
        if isempty(day_trials)
            continue;
        end
       
        
        % For each trial
        for trial_idx = 1:size(day_trials, 1)
            sustained_activity = mean(day_trials(trial_idx, time_window_idx), 'omitnan');
            
            % Store: [activity, animal_id, day, group, animal_string]
            all_data = [all_data; sustained_activity, ai, day_idx, group_label];
            animal_ids_str{end+1} = animal_id;
        end
    end
end

% Convert to table
data_table = array2table(all_data, 'VariableNames', {'Activity', 'AnimalID', 'Day', 'Group'});

% Add animal string IDs for tracking
data_table.AnimalIDStr = animal_ids_str';

% Convert to categorical with meaningful labels
data_table.AnimalID = categorical(data_table.AnimalID);


% Diagnostic: Check group means
fprintf('\n=== Diagnostic: Group Means ===\n');
for g = 1:3
    group_names = {'mPFC+', 'mPFC-', 'Non-learners'};
    group_data = data_table.Activity(data_table.Group == categorical(g, [1,2,3], {'mPFC_plus', 'mPFC_minus', 'Non_learners'}));
    fprintf('%s: Mean = %.6f, SD = %.6f, N = %d\n', ...
            group_names{g}, mean(group_data, 'omitnan'), std(group_data, 'omitnan'), length(group_data));
end

% Check per-animal means
fprintf('\n=== Per-Animal Means ===\n');
unique_animals = unique(data_table.AnimalIDStr);
for i = 1:length(unique_animals)
    animal = unique_animals{i};
    animal_data = data_table(strcmp(data_table.AnimalIDStr, animal), :);
    fprintf('%s (Group %s): Mean = %.6f, N_trials = %d\n', ...
            animal, char(animal_data.Group(1)), mean(animal_data.Activity), height(animal_data));
end

%% Fit mixed-effects model

fprintf('\n=== Fitting Mixed-Effects Model ===\n');

% Model with proper contrasts
lme = fitlme(data_table, 'Activity ~ Group + (1|AnimalID) + (1|AnimalID:Day)');

% Display results
disp(lme)

%% ANOVA for overall group effect
fprintf('\n=== ANOVA: Overall Group Effect ===\n');
anova_results = anova(lme);
disp(anova_results)

%% Pairwise comparisons

fprintf('\n=== Pairwise Comparisons ===\n');

% Get coefficient names to build contrasts correctly
coef_names = lme.CoefficientNames;
fprintf('Coefficient names: ');
disp(coef_names);

% Build contrast matrix
% Assuming order: (Intercept), Group_mPFC_minus, Group_Non_learners
% Reference is mPFC_plus

% Contrast 1: mPFC_minus vs mPFC_plus (test if Group_mPFC_minus coef ≠ 0)
contrast1 = [0, 1, 0];  % Tests Group_mPFC_minus coefficient

% Contrast 2: Non_learners vs mPFC_plus (test if Group_Non_learners coef ≠ 0)
contrast2 = [0, 0, 1];  % Tests Group_Non_learners coefficient

% Contrast 3: mPFC_minus vs Non_learners (difference of differences)
contrast3 = [0, 1, -1];  % Tests Group_mPFC_minus - Group_Non_learners

[psi1, psici1, stats1] = coefTest(lme, contrast1);
[psi2, psici2, stats2] = coefTest(lme, contrast2);
[psi3, psici3, stats3] = coefTest(lme, contrast3);

fprintf('mPFC- vs mPFC+: Diff = %.6f, p = %.4f\n', psi1, stats1.pValue);
fprintf('Non-learners vs mPFC+: Diff = %.6f, p = %.4f\n', psi2, stats2.pValue);
fprintf('mPFC- vs Non-learners: Diff = %.6f, p = %.4f\n', psi3, stats3.pValue);

%% Effect sizes (Cohen's d)
fprintf('\n=== Effect Sizes ===\n');

group_means_vec = [
    mean(data_table.Activity(data_table.Group == 'mPFC_plus'));
    mean(data_table.Activity(data_table.Group == 'mPFC_minus'));
    mean(data_table.Activity(data_table.Group == 'Non_learners'))
];

pooled_std = std(data_table.Activity, 'omitnan');

cohens_d_1vs2 = (group_means_vec(1) - group_means_vec(2)) / pooled_std;
cohens_d_1vs3 = (group_means_vec(1) - group_means_vec(3)) / pooled_std;

fprintf('mPFC+ vs mPFC-: Cohen''s d = %.3f\n', cohens_d_1vs2);
fprintf('mPFC+ vs Non-learners: Cohen''s d = %.3f\n', cohens_d_1vs3);

%%

