%% 11/12/2025 Files that I can load to work on cohort until HA015

% new_packaging_test_task_widefield
% behaviour_structure_all_animals_11_12_25

%%
%% Load variables

load('C:\Users\havgana\Desktop\DPhil\packaged_data\task_widefield_08_01_26.mat'); % wf task
load('C:\Users\havgana\Desktop\DPhil\packaged_data\behaviour_structure_all_animals_11_12_25.mat') % behaviour
load('C:\Users\havgana\Desktop\DPhil\packaged_data\combined_sig_day_all_protocols_01_12_26.mat') % sig days

% load master_U
wf_U= plab.wf.load_master_U;

added_time = -1:0.03:3;
added_time_Kernel=  fliplr((-10:100)/30); % kernel

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
cmPFC_plus    = [0 0.6 0];   % dark green
cmPFC_minus = [0.6 0 0.6]; % purple

%% Plot line variables setup basline substracted

event_times= [-0.5,0.75];
baseline_time_windows = [-0.1, 0];

% if you want to examine CS- you need to padd the non_rewarded variables

% Find indices where protocol is 3
protocol_3_mask = workflow_cat == 3;

% Initialize non_rewarded_stim_kernel with same size as rewarded_stim_kernel
if  length(non_rewarded_stim_kernel) < length(rewarded_stim_kernel)
    % Create full-size array initialized with NaNs
    non_rewarded_stim_kernel_padded = nan(size(rewarded_stim_kernel));
    
    % Fill in existing data for non-protocol-3 sessions
    non_protocol_3_idx = find(~protocol_3_mask);
    non_rewarded_stim_kernel_padded(:,:,non_protocol_3_idx) = non_rewarded_stim_kernel;
    
    % Replace original with padded version
    % non_rewarded_stim_kernel = non_rewarded_stim_kernel_padded;
end

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



%% Plot pre and post Scatter for Right vs Left mPFC responses for all animals

protocol_idx = 1;
animals = unique(widefield_animal_idx);

% Define group assignment by actual animal IDs
animal_groups = containers.Map('KeyType', 'double', 'ValueType', 'double');

% mPFC+ animals (use actual animal IDs from widefield_animal_idx)
mPFC_plus_IDs = [2 5 7 8 9]; % Replace with actual IDs
for id = mPFC_plus_IDs
    animal_groups(id) = 1;
end

% mPFC- animals
mPFC_minus_IDs = [1 3 4 6 10 11 12]; % Replace with actual IDs
for id = mPFC_minus_IDs
    animal_groups(id) = 0;
end



% Storage
animal_mPFC_right_pre = []; animal_mPFC_left_pre = [];
animal_mPFC_right_post = []; animal_mPFC_left_post = [];
animal_ids_included = []; % Track which animals are included

counter = 0;
for ai = animals(:)'
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx);
    days_idx = find(sel);
    if isempty(days_idx), continue; end
    
    learn_local = learning_index_animal(days_idx) == 1;

    % Find first occurrence of two consecutive significant days
    consecutive = learn_local(1:end-1) & learn_local(2:end);
    ld = find(consecutive, 1, 'first');

    % if isempty(ld), continue; end % No consecutive days found
    if isempty(ld)
        pre_days= 1:length(days_idx); % if no learning - then all days are pre

    else

        pre_days = find((1:length(days_idx)) < ld);
        post_days = find((1:length(days_idx)) >= ld);
    end
    % if isempty(post_days) || isempty(pre_days), continue; end
    
    counter = counter + 1;
    animal_ids_included(counter) = ai;
    
    response_window = t_for_plot >= 0 & t_for_plot <= 0.35;
    
    % Pre-learning
    TR_pre = right_mPFC_ROI_trace(days_idx(pre_days), :);
    TL_pre = left_mPFC_ROI_trace(days_idx(pre_days), :);
    TR_pre_mean = mean(TR_pre, 1, 'omitnan');
    TL_pre_mean = mean(TL_pre, 1, 'omitnan');
    animal_mPFC_right_pre(counter) = max(TR_pre_mean(response_window));
    animal_mPFC_left_pre(counter) = max(TL_pre_mean(response_window));
    
    % Post-learning
    TR_post = right_mPFC_ROI_trace(days_idx(post_days), :);
    TL_post = left_mPFC_ROI_trace(days_idx(post_days), :);
    TR_post_mean = mean(TR_post, 1, 'omitnan');
    TL_post_mean = mean(TL_post, 1, 'omitnan');
    animal_mPFC_right_post(counter) = max(TR_post_mean(response_window));
    animal_mPFC_left_post(counter) = max(TL_post_mean(response_window));
end


% Assign colors based on group
colors_mat = nan(length(animal_ids_included), 3);
for i = 1:length(animal_ids_included)
    ai = animal_ids_included(i);
    
    if isKey(animal_groups, ai)
        if animal_groups(ai) == 1
            colors_mat(i, :) = cmPFC_plus;
        else
            colors_mat(i, :) = cmPFC_minus;
        end
    else
        colors_mat(i, :) = [0.5 0.5 0.5]; % Unclassified
    end
end

% Plot
figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
max_val = max([animal_mPFC_right_pre, animal_mPFC_left_pre, ...
               animal_mPFC_right_post, animal_mPFC_left_post]);

for sp = 1:2
    subplot(1, 2, sp);
    if sp == 1
        data_x = animal_mPFC_right_pre; 
        data_y = animal_mPFC_left_pre;
        tit = 'Pre-Learning';
    else
        data_x = animal_mPFC_right_post; 
        data_y = animal_mPFC_left_post;
        tit = 'Post-Learning';
    end
    
    % Plot each point with its assigned color
    for i = 1:length(data_x)
        scatter(data_x(i), data_y(i), 150, colors_mat(i, :), 'filled', ...
                'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.7, 'LineWidth', 1.5);
        hold on;
    end
    
    plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 2);
    xlabel('Right mPFC Activity', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('Left mPFC Activity', 'FontSize', 13, 'FontWeight', 'bold');
    title(tit, 'FontSize', 14, 'FontWeight', 'bold');
    axis equal; xlim([0, max_val*1.05]); ylim([0, max_val*1.05]);
    set(gca, 'FontSize', 12, 'LineWidth', 1.2);
end

% Add legend
legend([scatter(NaN, NaN, 150, cmPFC_plus, 'filled', 'MarkerEdgeColor', 'k'), ...
        scatter(NaN, NaN, 150, cmPFC_minus, 'filled', 'MarkerEdgeColor', 'k')], ...
       {'mPFC+', 'mPFC-'}, 'Location', 'best');

%% ===== SESSION-LEVEL mPFC vs BEHAVIOR ANALYSIS =====

protocol_idx = 1;
response_window = t_for_plot >= 0 & t_for_plot <= 0.35;

% Pre-allocate
n_sessions = sum(workflow_cat == protocol_idx);
session_animal = nan(n_sessions, 1);
session_mPFC_activity = nan(n_sessions, 1);
session_ITI_licks = nan(n_sessions, 1);
session_discrimination = nan(n_sessions, 1);
session_reaction_time = nan(n_sessions, 1);
session_significance= nan(n_sessions,1);

% Create mapping: for each widefield_cat(idx), what day number is it within that animal?
day_within_animal = nan(length(widefield_cat), 1);
for animal = unique(widefield_animal_idx)'
    animal_mask = widefield_animal_idx == animal;
    day_within_animal(animal_mask) = 1:sum(animal_mask);
end

counter = 0;

for idx = 1:length(widefield_cat)
    if workflow_cat(idx) ~= protocol_idx || isempty(widefield_cat(idx).rewarded_stim_start_to_move_aligned_V)
        continue;
    end
    
    counter = counter + 1;
    animal_idx = widefield_animal_idx(idx);
    day_num = day_within_animal(idx);
    
    session_animal(counter) = animal_idx;
    
    % mPFC activity
    trace_R = right_mPFC_ROI_trace(idx, :);
    trace_L = left_mPFC_ROI_trace(idx, :);
    peak_R = max(trace_R(response_window));
    peak_L = max(trace_L(response_window));
    % session_mPFC_activity(counter) = mean([peak_R, peak_L]);
    session_mPFC_activity(counter)= peak_L;

    % Behavioral metrics
    day_data= behaviour_data(animal_idx).recording_day(day_num);

    % calculate mean ITI lick rate
    iti_licks = day_data.ITI_lick_counts;
    iti_duration =day_data.ITI_actual_duration;
    session_ITI_licks(counter) =nanmean(iti_licks./iti_duration);

    % calculate CS+ and CS- discrimintion index
    cs_plus_mask= day_data.cs_labels ==1;
    cs_plus_licks= mean(day_data.anticipatory_licks(cs_plus_mask));
    cs_minus_licks= mean(day_data.anticipatory_licks(~cs_plus_mask));
    session_discrimination(counter) = cs_plus_licks- cs_minus_licks ./ cs_plus_licks + cs_minus_licks;

    % Reaction times
    session_reaction_time(counter) = nanmedian(day_data.all_stim_diff_from_optimal_reward(cs_plus_mask));

    % Mask whether significant or not
    sig_day = combined_sig_day_all_protocols{animal_idx}(day_num);
    session_significance(counter)= sig_day;
end

% Trim
session_animal = session_animal(1:counter);
session_mPFC_activity = session_mPFC_activity(1:counter);
session_ITI_licks = session_ITI_licks(1:counter);
session_discrimination = session_discrimination(1:counter);
session_reaction_time = session_reaction_time(1:counter);
session_significance = session_significance(1:counter);

% Separate pre vs post learning
pre_mask = session_significance == 0;
post_mask = session_significance == 1;


%% ===== PLOT =====%%

behavior_metric = log(session_reaction_time); 
behavior_label = 'RT';

unique_animals = unique(session_animal);
colors = distinguishable_colors(length(unique_animals));

figure('Color', 'w', 'Position', [100, 100, 700, 600]);
hold on;

for ai = 1:length(unique_animals)
    sess_mask = session_animal == unique_animals(ai);

    combined_mask= sess_mask & post_mask; % in case you want to filter pre-post learning

    x = behavior_metric(combined_mask);
    y = session_mPFC_activity(combined_mask);

    % Scatter
    scatter(x, y, 100, colors(ai, :), 'filled', 'MarkerEdgeColor', 'k', ...
        'LineWidth', 1, 'MarkerFaceAlpha', 0.7);
  

    % Fit line if >2 points
    if sum(combined_mask) > 2
        p = polyfit(x, y, 1);
        x_fit = linspace(min(x), max(x), 50);
        plot(x_fit, polyval(p, x_fit), '-', 'Color', colors(ai, :), 'LineWidth', 1.5);
    end

    % % Plot color bar
    % colormap(colors);
    % colorbar('northoutside');
end



xlabel(behavior_label, 'FontSize', 13, 'FontWeight', 'bold');
ylabel('mPFC Activity', 'FontSize', 13, 'FontWeight', 'bold');
title('Within-Animal Trajectories', 'FontSize', 14, 'FontWeight', 'bold');
% Create discrete colormap
colormap(colors);
cb = colorbar('Location', 'eastoutside');
cb.Label.String = 'Mouse';
cb.Label.FontSize = 12;
cb.Label.FontWeight = 'bold';

% Set discrete ticks for each animal
cb.Ticks = linspace(0, 1, n_animals);
cb.TickLabels = arrayfun(@num2str, 1:n_animals, 'UniformOutput', false);

% Make colorbar smaller
cb.Position(3) = 0.02; % Width
cb.Position(4) = 0.3;  % Height

% Center vertically
cb.Position(2) = 0.35; % Adjust Y position

set(gca, 'FontSize', 12, 'LineWidth', 1.2); box off;


%% Run mixed-effects model to predict mPFC activity based on behavioural metrics

session_ITI_licks;
session_discrimination;
session_reaction_time;


% Requires Statistics Toolbox
tbl = table(session_animal, behavior_metric, session_mPFC_activity, ...
            'VariableNames', {'Animal', 'Behavior', 'mPFC'});

% Random intercept model
lme = fitlme(tbl, 'mPFC ~ Behavior + (1|Animal)');
disp(lme);

%% for a full model

learning_stage= post_mask; % binary

% Requires Statistics Toolbox
tbl = table(session_animal, session_ITI_licks,session_discrimination,log(session_reaction_time), session_mPFC_activity, learning_stage,...
            'VariableNames', {'Animal', 'ITI_licks','Discrimintation','RT' 'mPFC','learning_stage'});


tbl.learning_stage_cat = categorical(tbl.learning_stage, [0 1], {'Pre', 'Post'}); % convert to categorical
lme = fitlme(tbl, 'mPFC ~ (ITI_licks + Discrimintation + RT) * learning_stage_cat + (1|Animal)');

% Random intercept model without pre post
% lme = fitlme(tbl, 'mPFC ~ (ITI_licks + Discrimintation + RT) + (1|Animal)');
disp(lme);

%% Interaction plots (Split by Pre-Post)

figure('Position', [100, 100, 1400, 400]);
behaviors = {'ITI_licks', 'Discrimintation', 'RT'};
titles = {'ITI Licks', 'CS+ vs CS- Discrimination', 'Reaction Time (log)'};

% Find global y-limits first
y_min = min(tbl.mPFC);
y_max = max(tbl.mPFC);
y_range = [y_min - 0.05*(y_max-y_min), y_max + 0.05*(y_max-y_min)];

for i = 1:3
    subplot(1, 3, i);
    
    % Pre-learning
    pre_mask = tbl.learning_stage == 0;
    scatter(tbl.(behaviors{i})(pre_mask), tbl.mPFC(pre_mask), 50, [0.7 0.7 0.7], 'filled', 'MarkerFaceAlpha', 0.5);
    hold on;
    
    % Post-learning  
    post_mask = tbl.learning_stage == 1;
    scatter(tbl.(behaviors{i})(post_mask), tbl.mPFC(post_mask), 50, [0.2 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.5);
    
    % Fit lines
    p_pre = polyfit(tbl.(behaviors{i})(pre_mask), tbl.mPFC(pre_mask), 1);
    p_post = polyfit(tbl.(behaviors{i})(post_mask), tbl.mPFC(post_mask), 1);
    
    x_range = xlim;
    x_fit = linspace(x_range(1), x_range(2), 100);
    plot(x_fit, polyval(p_pre, x_fit), '--', 'Color', [0.7 0.7 0.7], 'LineWidth', 2);
    plot(x_fit, polyval(p_post, x_fit), '-', 'Color', [0.2 0.5 0.8], 'LineWidth', 2);
    
    ylim(y_range); % Set consistent y-axis
    
    xlabel(titles{i}, 'FontSize', 11, 'FontWeight', 'bold');
    if i == 1, ylabel('mPFC Activity', 'FontSize', 11, 'FontWeight', 'bold'); end
    legend({'Pre', 'Post'}, 'Location', 'best');
    set(gca, 'FontSize', 10);
end

%% Plot beta Coefficents 

% Extract interaction coefficients
interactions = [lme.Coefficients.Estimate(6); lme.Coefficients.Estimate(7); lme.Coefficients.Estimate(8)];
CI_lower = [lme.Coefficients.Lower(6); lme.Coefficients.Lower(7); lme.Coefficients.Lower(8)];
CI_upper = [lme.Coefficients.Upper(6); lme.Coefficients.Upper(7); lme.Coefficients.Upper(8)];
p_values = [lme.Coefficients.pValue(6); lme.Coefficients.pValue(7); lme.Coefficients.pValue(8)];

figure('Position', [100, 100, 500, 400]);
hold on;

for i = 1:3
    errorbar(i, interactions(i), interactions(i)-CI_lower(i), CI_upper(i)-interactions(i), ...
             'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerFaceColor', [0.2 0.5 0.8]);
    
    % Add significance asterisks
    if p_values(i) < 0.001
        sig_text = '***';
    elseif p_values(i) < 0.01
        sig_text = '**';
    elseif p_values(i) < 0.05
        sig_text = '*';
    else
        sig_text = 'ns';
    end
    
    % Position above the upper CI
    y_pos = CI_upper(i) + 0.1 * abs(CI_upper(i));
    text(i, y_pos, sig_text, 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');
end

plot([0 4], [0 0], 'k--', 'LineWidth', 1);
xlim([0.5 3.5]);
xticks(1:3);
xticklabels({'ITI × Post', 'Discrim × Post', 'RT × Post'});
ylabel('Interaction Coefficient', 'FontSize', 12, 'FontWeight', 'bold');
title('Learning-Dependent Changes', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);


%% Plot color coded bars for individual animal slopes 

% First, extract individual animal ITI slopes (post-learning)
lme_slopes = fitlme(tbl, 'mPFC ~ Discrimintation * learning_stage + (Discrimintation|Animal)');
[~, ~, stats] = randomEffects(lme_slopes);
animal_ITI_slopes = stats.Estimate(2:2:end);
unique_animals_in_model = unique(tbl.Animal);

% Map slopes to mPFC groups
slope_colors = nan(length(animal_ITI_slopes), 3);
for i = 1:length(unique_animals_in_model)
    ai = unique_animals_in_model(i);
    if isKey(animal_groups, ai)
        if animal_groups(ai) == 1
            slope_colors(i, :) = cmPFC_plus;
        else
            slope_colors(i, :) = cmPFC_minus;
        end
    else
        slope_colors(i, :) = [0.5 0.5 0.5];
    end
end

% Plot histogram with individual animals colored
figure('Position', [100, 100, 600, 500]);
hold on;

% Plot each animal's slope as a bar
for i = 1:length(animal_ITI_slopes)
    bar(i, animal_ITI_slopes(i), 'FaceColor', slope_colors(i, :), 'EdgeColor', 'k', 'LineWidth', 1.5,'FaceAlpha',0.5);
end

% Add zero line
plot([0 length(animal_ITI_slopes)+1], [0 0], 'k--', 'LineWidth', 2);

xlabel('Animal', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('ITI-mPFC Slope (Post-Learning)', 'FontSize', 12, 'FontWeight', 'bold');
title('Individual Animal Slopes by Group', 'FontSize', 13, 'FontWeight', 'bold');

% % Add legend
% legend([bar(NaN, NaN, 'FaceColor', cmPFC_plus, 'EdgeColor', 'k','FaceAlpha'), ...
%         bar(NaN, NaN, 'FaceColor', cmPFC_minus, 'EdgeColor', 'k')], ...
%        {'mPFC+', 'mPFC-'}, 'Location', 'best', 'FontSize', 11);
% set(gca, 'FontSize', 11);



%% Do animals show mPFC responses when CS+ starts to move?

%% ===== EXTRACT PRE-LICK mPFC ACTIVITY =====

protocol_idx = 1;
animals = unique(widefield_animal_idx);
ROI_to_extract = right_mPFC_ROI_mask;
t_for_plot = added_time;
curr_components= n_components;

% Time buffer before first lick (exclude activity too close to motor act)
pre_lick_buffer = 0.1; % seconds

% Storage
animal_results = struct();
counter = 0;

for ai = animals(:)'

    % Find all post-learning days for this animal
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx) & (learning_index_animal == 1);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    counter = counter + 1;
    
    % Map widefield_cat indices to behaviour_data day numbers
    day_within_animal = nan(length(widefield_cat), 1);
    animal_mask = widefield_animal_idx == ai;
    day_within_animal(animal_mask) = 1:sum(animal_mask);
    
    % Collect all trials across post-learning days
    all_mPFC_traces = [];
    all_mPFC_pre_lick = [];
    all_lick_latencies = [];
    
    for idx = days_idx(:)'
        V_data = widefield_cat(idx).rewarded_stim_start_to_move_aligned_V;
        if isempty(V_data), continue; end
        
        % Extract mPFC ROI traces
        mPFC_traces = ap.wf_roi(wf_U(:,:,1:curr_components), V_data, [], [], ROI_to_extract);
        mPFC_traces = squeeze(mPFC_traces); % Time x Trials
        n_trials = size(mPFC_traces, 2);
        
        % Get lick latencies
        day_num = day_within_animal(idx);
        lick_lat = behaviour_data(ai).recording_day(day_num).cs_plus_latency_lick_to_move;
        
        % Ensure matching trial counts
        if length(lick_lat) ~= n_trials
            warning('Trial mismatch for animal %d day %d', ai, day_num);
            continue;
        end
        
        % Create pre-lick masked traces
        mPFC_pre_lick = nan(size(mPFC_traces));
        
        for trial = 1:n_trials
            first_lick_time = lick_lat(trial);
            if isnan(first_lick_time), continue; end
            
            if first_lick_time >= 1
                continue;
            end

            % Cutoff time: first lick - buffer
            cutoff_time = first_lick_time - pre_lick_buffer;

            % Create mask: 1 before cutoff, NaN after
            pre_lick_mask = t_for_plot <= cutoff_time;
            
            % Apply mask
            mPFC_pre_lick(pre_lick_mask, trial) = mPFC_traces(pre_lick_mask, trial);
        end
        
        % Accumulate
        all_mPFC_traces = [all_mPFC_traces, mPFC_traces];
        all_mPFC_pre_lick = [all_mPFC_pre_lick, mPFC_pre_lick];
        all_lick_latencies = [all_lick_latencies; lick_lat(:)];
    end
    
    % Store results
    animal_results(counter).animal_id = ai;
    animal_results(counter).mPFC_traces_full = all_mPFC_traces;
    animal_results(counter).mPFC_traces_pre_lick = all_mPFC_pre_lick;
    animal_results(counter).lick_latencies = all_lick_latencies;
    
    % Compute averages
    animal_results(counter).mean_full = mean(all_mPFC_traces, 2, 'omitnan');
    animal_results(counter).mean_pre_lick = mean(all_mPFC_pre_lick, 2, 'omitnan');
    
    % Group assignment
    if isKey(animal_groups, ai)
        animal_results(counter).group = animal_groups(ai);
    else
        animal_results(counter).group = NaN;
    end
end

% Trim
animal_results = animal_results(1:counter);

%% ===== VISUALIZE: COMPARE FULL VS PRE-LICK =====
mPFC_plus_idx = [animal_results.group] == 1;
mPFC_minus_idx = [animal_results.group] == 0;

% Compute group averages
mean_full_plus = mean(cat(2, animal_results(mPFC_plus_idx).mean_full), 2, 'omitnan');
mean_full_minus = mean(cat(2, animal_results(mPFC_minus_idx).mean_full), 2, 'omitnan');

mean_pre_lick_plus = mean(cat(2, animal_results(mPFC_plus_idx).mean_pre_lick), 2, 'omitnan');
mean_pre_lick_minus = mean(cat(2, animal_results(mPFC_minus_idx).mean_pre_lick), 2, 'omitnan');

% Plot group comparison
figure('Position', [100, 100, 1200, 500]);

% mPFC+ group
subplot(1, 2, 1);
hold on;
plot(t_for_plot, mean_full_plus, '-', 'Color', cmPFC_plus, 'LineWidth', 2.5);
plot(t_for_plot, mean_pre_lick_plus, '-', 'Color', [cmPFC_plus * 0.5], 'LineWidth', 2.5);
plot([0 0], ylim, 'k--', 'LineWidth', 1.5);
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('mPFC Activity', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC+ Animals', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Full trace', 'Pre-lick only'}, 'Location', 'best', 'FontSize', 10);
xlim([-0.5 1]);
set(gca, 'FontSize', 11);

% mPFC- group
subplot(1, 2, 2);
hold on;
plot(t_for_plot, mean_full_minus, '-', 'Color', cmPFC_minus, 'LineWidth', 2.5);
plot(t_for_plot, mean_pre_lick_minus, '-', 'Color', [cmPFC_minus * 0.5], 'LineWidth', 2.5);
plot([0 0], ylim, 'k--', 'LineWidth', 1.5);
xlabel('Time from Movement (s)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('mPFC Activity', 'FontSize', 12, 'FontWeight', 'bold');
title('mPFC- Animals', 'FontSize', 13, 'FontWeight', 'bold');
legend({'Full trace', 'Pre-lick only'}, 'Location', 'best', 'FontSize', 10);
xlim([-0.5 1]);
set(gca, 'FontSize', 11);

sgtitle('mPFC Activity: Full vs Pre-Lick', 'FontSize', 14, 'FontWeight', 'bold');

%% ===== VISUALIZE: INDIVIDUAL ANIMALS =====
figure('Position', [100, 100, 1400, 800]);

n_animals_plus = sum(mPFC_plus_idx);
n_animals_minus = sum(mPFC_minus_idx);

% mPFC+ animals
for i = 1:n_animals_plus
    idx = find(mPFC_plus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), i);
    hold on;
    plot(t_for_plot, animal_results(idx).mean_full, '-', 'Color', cmPFC_plus, 'LineWidth', 2);
    plot(t_for_plot, animal_results(idx).mean_pre_lick, '-', 'Color', [cmPFC_plus * 0.5], 'LineWidth', 2);
    plot([0 0], ylim, 'k--', 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('mPFC Activity', 'FontSize', 8);
    title(sprintf('M%d', animal_results(idx).animal_id), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

% mPFC- animals
for i = 1:n_animals_minus
    idx = find(mPFC_minus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), max(n_animals_plus, n_animals_minus) + i);
    hold on;
    plot(t_for_plot, animal_results(idx).mean_full, '-', 'Color', cmPFC_minus, 'LineWidth', 2);
    plot(t_for_plot, animal_results(idx).mean_pre_lick, '-', 'Color', [cmPFC_minus * 0.5], 'LineWidth', 2);
    plot([0 0], ylim, 'k--', 'LineWidth', 1);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('mPFC Activity', 'FontSize', 8);
    title(sprintf('M%d', animal_results(idx).animal_id), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

sgtitle('Individual Animals: Full (bright) vs Pre-Lick (dark)', 'FontSize', 12, 'FontWeight', 'bold');

%% ===== TEST: mPFC PEAK TIMING vs LICK TIMING =====

% Create a datastructure that includes mPFC activity sorted by latency and
% the lick latencies on trial by trial basis

protocol_idx = 1;
animals = unique(widefield_animal_idx);
ROI_to_extract = left_mPFC_ROI_mask;
t_for_plot = added_time;

% Analysis window for finding peak
peak_window = t_for_plot >= 0 & t_for_plot <= 1; % 

% Storage
animal_results = struct();
counter = 0;

for ai = animals(:)'
    % Find all post-learning days for this animal
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx) & (learning_index_animal == 1);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    counter = counter + 1;
    
    % Map widefield_cat indices to behaviour_data day numbers
    day_within_animal = nan(length(widefield_cat), 1);
    animal_mask = widefield_animal_idx == ai;
    day_within_animal(animal_mask) = 1:sum(animal_mask);
    
    % Collect all trials across post-learning days
    all_mPFC_traces = [];
    all_lick_latencies = [];
    
    for idx = days_idx(:)'
        V_data = widefield_cat(idx).rewarded_stim_start_to_move_aligned_V;
        if isempty(V_data), continue; end
        
        % Extract mPFC ROI traces
        mPFC_traces = ap.wf_roi(wf_U(:,:,1:curr_components), V_data, [], [], ROI_to_extract);
        mPFC_traces = squeeze(mPFC_traces); % Time x Trials
        
        % Get lick latencies
        day_num = day_within_animal(idx);
        lick_lat = behaviour_data(ai).recording_day(day_num).cs_plus_latency_lick_to_move;
        
        % Ensure matching trial counts
        n_trials = size(mPFC_traces, 2);
        if length(lick_lat) ~= n_trials
            warning('Trial mismatch for animal %d day %d', ai, day_num);
            continue;
        end
        
        % Accumulate
        all_mPFC_traces = [all_mPFC_traces, mPFC_traces];
        all_lick_latencies = [all_lick_latencies; lick_lat(:)];
    end
    
    % Remove invalid trials
    valid = ~any(isnan(all_mPFC_traces), 1) & ~isnan(all_lick_latencies') & (all_lick_latencies'<=1); % cut off when latency was bigger than 1s
    all_mPFC_traces = all_mPFC_traces(:, valid);
    all_lick_latencies = all_lick_latencies(valid);
    
    if length(all_lick_latencies) < 30, continue; end % Need sufficient trials
    
    % Calculate peak mPFC latency for each trial
    mPFC_peak_latencies = nan(length(all_lick_latencies), 1);
    for trial = 1:length(all_lick_latencies)
        trial_trace = all_mPFC_traces(peak_window, trial);
        if all(isnan(trial_trace)), continue; end
        
        [~, peak_idx] = max(trial_trace);
        peak_times = t_for_plot(peak_window);
        mPFC_peak_latencies(trial) = peak_times(peak_idx);
    end
    
    % Sort trials by lick latency
    [sorted_lick_lat, sort_idx] = sort(all_lick_latencies);
    sorted_mPFC_traces = all_mPFC_traces(:, sort_idx);
    sorted_mPFC_peaks = mPFC_peak_latencies(sort_idx);
    
    % Store results
    animal_results(counter).animal_id = ai;
    animal_results(counter).sorted_lick_lat = sorted_lick_lat;
    animal_results(counter).sorted_mPFC_traces = sorted_mPFC_traces;
    animal_results(counter).sorted_mPFC_peaks = sorted_mPFC_peaks;
    
    % Group assignment
    if isKey(animal_groups, ai)
        animal_results(counter).group = animal_groups(ai);
    else
        animal_results(counter).group = NaN;
    end
end

% Trim
animal_results = animal_results(1:counter);


%% Plot individual animal plots of mPFC peak vs lick latency with unity line

n_animals_plus = sum(mPFC_plus_idx);
n_animals_minus = sum(mPFC_minus_idx);

% Storage for correlations
correlations_plus = nan(n_animals_plus, 1);
correlations_minus = nan(n_animals_minus, 1);

figure('Position', [100, 100, 1400, 800]);

% mPFC+ animals
for i = 1:n_animals_plus
    idx = find(mPFC_plus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), i);
    
    % Get valid trials
    valid = ~isnan(animal_results(idx).sorted_mPFC_peaks) & ...
            ~isnan(animal_results(idx).sorted_lick_lat);
    
    lick_lat = animal_results(idx).sorted_lick_lat(valid);
    peak_lat = animal_results(idx).sorted_mPFC_peaks(valid);
    
    % Scatter
    scatter(lick_lat, peak_lat, 20, cmPFC_plus, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    
    % Unity line (perfect behavior-locking)
    plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
    
    % Fit line
    if length(lick_lat) > 5
        p = polyfit(lick_lat, peak_lat, 1);
        x_fit = [min(lick_lat), max(lick_lat)];
        plot(x_fit, polyval(p, x_fit), 'r-', 'LineWidth', 2);
    end
    
    % Calculate correlation
    r = corr(lick_lat, peak_lat, 'Type', 'Spearman');
    correlations_plus(i) = r;
    
    xlabel('Lick Latency (s)', 'FontSize', 8);
    ylabel('mPFC Peak Latency (s)', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f', animal_results(idx).animal_id, r), 'FontSize', 9);
    axis equal; xlim([0 1]); ylim([0 1]);
    set(gca, 'FontSize', 7);
end

% mPFC- animals
for i = 1:n_animals_minus
    idx = find(mPFC_minus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), max(n_animals_plus, n_animals_minus) + i);
    
    % Get valid trials
    valid = ~isnan(animal_results(idx).sorted_mPFC_peaks) & ...
            ~isnan(animal_results(idx).sorted_lick_lat);
    
    lick_lat = animal_results(idx).sorted_lick_lat(valid);
    peak_lat = animal_results(idx).sorted_mPFC_peaks(valid);
    
    % Scatter
    scatter(lick_lat, peak_lat, 20, cmPFC_minus, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    
    % Unity line (perfect behavior-locking)
    plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
    
    % Fit line
    if length(lick_lat) > 5
        p = polyfit(lick_lat, peak_lat, 1);
        x_fit = [min(lick_lat), max(lick_lat)];
        plot(x_fit, polyval(p, x_fit), 'r-', 'LineWidth', 2);
    end
    
    % Calculate correlation
    r = corr(lick_lat, peak_lat, 'Type', 'Spearman');
    correlations_minus(i) = r;
    
    xlabel('Lick Latency (s)', 'FontSize', 8);
    ylabel('mPFC Peak Latency (s)', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f', animal_results(idx).animal_id, r), 'FontSize', 9);
    axis equal; xlim([0 1]); ylim([0 1]);
    set(gca, 'FontSize', 7);
end

sgtitle('mPFC Peak Timing vs Lick Timing (Dashed=Unity, Red=Fit)', ...
        'FontSize', 12, 'FontWeight', 'bold');

%% ===== GROUP COMPARISON =====
figure('Position', [100, 100, 500, 400]);
hold on;

scatter(ones(n_animals_plus, 1), correlations_plus, 100, cmPFC_plus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
scatter(2*ones(n_animals_minus, 1), correlations_minus, 100, cmPFC_minus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);

plot([0.5 2.5], [0 0], 'k--', 'LineWidth', 1.5);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'mPFC+', 'mPFC-'});
ylabel('Correlation (Peak Timing vs Lick Timing)', 'FontSize', 12, 'FontWeight', 'bold');
title('Behavior-Locking of mPFC Timing', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

% Statistical test
[~, p_group] = ttest2(correlations_plus, correlations_minus);

fprintf('\n===== mPFC TIMING vs LICK TIMING =====\n');
fprintf('mPFC+: mean r = %.3f ± %.3f\n', mean(correlations_plus), std(correlations_plus));
fprintf('mPFC-: mean r = %.3f ± %.3f\n', mean(correlations_minus), std(correlations_minus));
fprintf('Group difference: p = %.4f\n', p_group);

% Add significance to plot
if p_group < 0.05
    y_max = max([correlations_plus; correlations_minus]) * 1.1;
    plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1.5);
    if p_group < 0.001
        text(1.5, y_max*1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_group < 0.01
        text(1.5, y_max*1.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end

%%
%% ===== ANALYSIS 1: PRE-MOVEMENT mPFC vs STATIC TIME =====

protocol_idx = 1;
animals = unique(widefield_animal_idx);
ROI_to_extract = left_mPFC_ROI_mask;
t_for_plot = added_time;

% Window for extracting pre-movement activity
pre_movement_window = t_for_plot >= -0.2 & t_for_plot <= 0;

% Storage
animal_results = struct();
counter = 0;

for ai = animals(:)'
    % Find all post-learning days for this animal
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx) & (learning_index_animal == 1);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    counter = counter + 1;
    
    % Map widefield_cat indices to behaviour_data day numbers
    day_within_animal = nan(length(widefield_cat), 1);
    animal_mask = widefield_animal_idx == ai;
    day_within_animal(animal_mask) = 1:sum(animal_mask);
    
    % Collect all trials across post-learning days
    all_pre_movement_activity = [];
    all_static_times = [];
    
    for idx = days_idx(:)'
        V_data = widefield_cat(idx).rewarded_stim_start_to_move_aligned_V;
        if isempty(V_data), continue; end
        
        % Extract mPFC ROI traces
        mPFC_traces = ap.wf_roi(wf_U(:,:,1:curr_components), V_data, [], [], ROI_to_extract);
        mPFC_traces = squeeze(mPFC_traces); % Time x Trials
        n_trials = size(mPFC_traces, 2);
        
        % Extract pre-movement activity (mean in window)
        pre_movement_activity = max(mPFC_traces(pre_movement_window, :), [],1)';
        
        % Get static times
        day_num = day_within_animal(idx);

        cs_plus_mask= behaviour_data(ai).recording_day(day_num).cs_labels;
        % Get static times and lick latencies
        static_times = behaviour_data(ai).recording_day(day_num).all_static_times(cs_plus_mask==1);
    
        
        % Ensure matching trial counts
        if length(static_times) ~= n_trials
            warning('Trial mismatch for animal %d day %d', ai, day_num);
            continue;
        end
        
        % Accumulate
        all_pre_movement_activity = [all_pre_movement_activity; pre_movement_activity];
        all_static_times = [all_static_times; static_times(:)];
    end
    
    % Remove invalid trials
    valid = ~isnan(all_pre_movement_activity) & ~isnan(all_static_times);
    all_pre_movement_activity = all_pre_movement_activity(valid);
    all_static_times = all_static_times(valid);
    
    if length(all_static_times) < 10, continue; end
    
    % Calculate correlation
    [r, p] = corr(all_static_times, all_pre_movement_activity, 'Type', 'Spearman');
    
    % Store results
    animal_results(counter).animal_id = ai;
    animal_results(counter).static_times = all_static_times;
    animal_results(counter).pre_movement_activity = all_pre_movement_activity;
    animal_results(counter).correlation = r;
    animal_results(counter).p_value = p;
    
    % Group assignment
    if isKey(animal_groups, ai)
        animal_results(counter).group = animal_groups(ai);
    else
        animal_results(counter).group = NaN;
    end
end

% Trim
animal_results = animal_results(1:counter);

%% ===== COMPARE GROUPS =====
mPFC_plus_idx = [animal_results.group] == 1;
mPFC_minus_idx = [animal_results.group] == 0;

correlations_plus = [animal_results(mPFC_plus_idx).correlation];
correlations_minus = [animal_results(mPFC_minus_idx).correlation];

% Statistical tests
[~, p_group] = ttest2(correlations_plus, correlations_minus);
[~, p_plus_zero] = ttest(correlations_plus);
[~, p_minus_zero] = ttest(correlations_minus);

fprintf('\n===== PRE-MOVEMENT mPFC vs STATIC TIME =====\n');
fprintf('mPFC+: mean r = %.3f ± %.3f\n', mean(correlations_plus), std(correlations_plus));
fprintf('mPFC-: mean r = %.3f ± %.3f\n', mean(correlations_minus), std(correlations_minus));
fprintf('Group difference: p = %.4f\n', p_group);
fprintf('\nmPFC+ different from zero: p = %.4f\n', p_plus_zero);
fprintf('mPFC- different from zero: p = %.4f\n', p_minus_zero);

%% ===== VISUALIZE: INDIVIDUAL SCATTER PLOTS =====
n_animals_plus = sum(mPFC_plus_idx);
n_animals_minus = sum(mPFC_minus_idx);

figure('Position', [100, 100, 1400, 800]);

% mPFC+ animals
for i = 1:n_animals_plus
    idx = find(mPFC_plus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), i);
    scatter(animal_results(idx).static_times, animal_results(idx).pre_movement_activity, ...
            30, cmPFC_plus, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    
    % Fit line
    p_coef = polyfit(animal_results(idx).static_times, animal_results(idx).pre_movement_activity, 1);
    x_fit = [min(animal_results(idx).static_times), max(animal_results(idx).static_times)];
    
    if animal_results(idx).p_value < 0.05
        plot(x_fit, polyval(p_coef, x_fit), 'r-', 'LineWidth', 2);
    else
        plot(x_fit, polyval(p_coef, x_fit), 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Static Time (s)', 'FontSize', 8);
    ylabel('Pre-Movement mPFC', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f, p=%.3f', animal_results(idx).animal_id, ...
          animal_results(idx).correlation, animal_results(idx).p_value), 'FontSize', 9);
    set(gca, 'FontSize', 7);
end

% mPFC- animals
for i = 1:n_animals_minus
    idx = find(mPFC_minus_idx, i, 'first');
    idx = idx(end);
    
    subplot(2, max(n_animals_plus, n_animals_minus), max(n_animals_plus, n_animals_minus) + i);
    scatter(animal_results(idx).static_times, animal_results(idx).pre_movement_activity, ...
            30, cmPFC_minus, 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    
    % Fit line
    p_coef = polyfit(animal_results(idx).static_times, animal_results(idx).pre_movement_activity, 1);
    x_fit = [min(animal_results(idx).static_times), max(animal_results(idx).static_times)];
    
    if animal_results(idx).p_value < 0.05
        plot(x_fit, polyval(p_coef, x_fit), 'r-', 'LineWidth', 2);
    else
        plot(x_fit, polyval(p_coef, x_fit), 'k--', 'LineWidth', 1.5);
    end
    
    xlabel('Static Time (s)', 'FontSize', 8);
    ylabel('Pre-Movement mPFC', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f, p=%.3f', animal_results(idx).animal_id, ...
          animal_results(idx).correlation, animal_results(idx).p_value), 'FontSize', 9);
    set(gca, 'FontSize', 7);
end

sgtitle('Pre-Movement mPFC Activity vs Static Time (-200 to 0ms)', ...
        'FontSize', 12, 'FontWeight', 'bold');

%% ===== VISUALIZE: GROUP COMPARISON =====
figure('Position', [100, 100, 500, 400]);
hold on;

scatter(ones(n_animals_plus, 1), correlations_plus, 100, cmPFC_plus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
scatter(2*ones(n_animals_minus, 1), correlations_minus, 100, cmPFC_minus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);

plot([0.5 2.5], [0 0], 'k--', 'LineWidth', 1.5);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'mPFC+', 'mPFC-'});
ylabel('Correlation (Static Time vs Pre-Movement mPFC)', 'FontSize', 12, 'FontWeight', 'bold');
title('Temporal Expectation Encoding', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

% Add significance
if p_group < 0.05
    y_max = max([correlations_plus; correlations_minus]) * 1.1;
    plot([1 2], [y_max y_max], 'k-', 'LineWidth', 1.5);
    if p_group < 0.001
        text(1.5, y_max*1.05, '***', 'HorizontalAlignment', 'center', 'FontSize', 14);
    elseif p_group < 0.01
        text(1.5, y_max*1.05, '**', 'HorizontalAlignment', 'center', 'FontSize', 14);
    else
        text(1.5, y_max*1.05, '*', 'HorizontalAlignment', 'center', 'FontSize', 14);
    end
end
%% archived







%% ===== TEST TEMPORAL SHIFT ACROSS LICK BINS =====
protocol_idx = 1;
animals = unique(widefield_animal_idx);
ROI_to_extract = left_mPFC_ROI_mask;
t_for_plot = added_time;

% Analysis window for finding peak
peak_window = t_for_plot >= 0 & t_for_plot <= 0.9;

% Number of bins for lick latency
n_bins = 3; % Terciles: fast, medium, slow lickers

% Storage
animal_results = struct();
counter = 0;

for ai = animals(:)'
    % Find all post-learning days for this animal
    sel = (widefield_animal_idx == ai) & (workflow_cat == protocol_idx) & (learning_index_animal == 1);
    days_idx = find(sel);
    
    if isempty(days_idx), continue; end
    
    counter = counter + 1;
    
    % Map widefield_cat indices to behaviour_data day numbers
    day_within_animal = nan(length(widefield_cat), 1);
    animal_mask = widefield_animal_idx == ai;
    day_within_animal(animal_mask) = 1:sum(animal_mask);
    
    % Collect all trials across post-learning days
    all_mPFC_traces = [];
    all_lick_latencies = [];
    
    for idx = days_idx(:)'
        V_data = widefield_cat(idx).rewarded_stim_start_to_move_aligned_V;
        if isempty(V_data), continue; end
        
        % Extract mPFC ROI traces
        mPFC_traces = ap.wf_roi(wf_U(:,:,1:curr_components), V_data, [], [], ROI_to_extract);
        mPFC_traces = squeeze(mPFC_traces); % Time x Trials
        
        % Get lick latencies
        day_num = day_within_animal(idx);
        lick_lat = behaviour_data(ai).recording_day(day_num).cs_plus_latency_lick_to_move;
        
        % Ensure matching trial counts
        n_trials = size(mPFC_traces, 2);
        if length(lick_lat) ~= n_trials
            warning('Trial mismatch for animal %d day %d', ai, day_num);
            continue;
        end
        
        % Accumulate
        all_mPFC_traces = [all_mPFC_traces, mPFC_traces];
        all_lick_latencies = [all_lick_latencies; lick_lat(:)];
    end
    
    % Remove invalid trials
    valid = ~any(isnan(all_mPFC_traces), 1) & ~isnan(all_lick_latencies');
    all_mPFC_traces = all_mPFC_traces(:, valid);
    all_lick_latencies = all_lick_latencies(valid);
    
    if length(all_lick_latencies) < 30, continue; end % Need sufficient trials
    
    % Sort trials by lick latency and bin
    [sorted_lick_lat, sort_idx] = sort(all_lick_latencies);
    sorted_mPFC_traces = all_mPFC_traces(:, sort_idx);
    
    % Divide into bins
    trials_per_bin = floor(length(sorted_lick_lat) / n_bins);
    
    bin_peak_latencies = nan(n_bins, 1);
    bin_traces = cell(n_bins, 1);
    bin_lick_means = nan(n_bins, 1);
    
    for bin = 1:n_bins
        if bin < n_bins
            bin_trials = ((bin-1)*trials_per_bin + 1):(bin*trials_per_bin);
        else
            bin_trials = ((bin-1)*trials_per_bin + 1):length(sorted_lick_lat);
        end
        
        % Average trace for this bin
        bin_trace = mean(sorted_mPFC_traces(:, bin_trials), 2, 'omitnan');
        bin_traces{bin} = bin_trace;
        
        % Find peak latency
        [~, peak_idx] = max(bin_trace(peak_window));
        peak_times = t_for_plot(peak_window);
        bin_peak_latencies(bin) = peak_times(peak_idx);
        
        % Mean lick latency for this bin
        bin_lick_means(bin) = mean(sorted_lick_lat(bin_trials));
    end
    
    % Quantify temporal shift: correlation between bin lick latency and bin peak latency
    temporal_shift = corr(bin_lick_means, bin_peak_latencies, 'Type', 'Spearman');
    
    % Also compute range of peak latencies across bins
    peak_latency_range = max(bin_peak_latencies) - min(bin_peak_latencies);
    
    % Store results
    animal_results(counter).animal_id = ai;
    animal_results(counter).bin_lick_means = bin_lick_means;
    animal_results(counter).bin_peak_latencies = bin_peak_latencies;
    animal_results(counter).bin_traces = bin_traces;
    animal_results(counter).temporal_shift = temporal_shift;
    animal_results(counter).peak_latency_range = peak_latency_range;
    animal_results(counter).sorted_traces = sorted_mPFC_traces;
    
    % Group assignment
    if isKey(animal_groups, ai)
        animal_results(counter).group = animal_groups(ai);
    else
        animal_results(counter).group = NaN;
    end
end

% Trim
animal_results = animal_results(1:counter);

%% ===== COMPARE GROUPS =====
mPFC_plus_idx = [animal_results.group] == 1;
mPFC_minus_idx = [animal_results.group] == 0;

temporal_shift_plus = [animal_results(mPFC_plus_idx).temporal_shift];
temporal_shift_minus = [animal_results(mPFC_minus_idx).temporal_shift];

peak_range_plus = [animal_results(mPFC_plus_idx).peak_latency_range];
peak_range_minus = [animal_results(mPFC_minus_idx).peak_latency_range];

% Statistical tests
[~, p_shift] = ttest2(temporal_shift_plus, temporal_shift_minus);
[~, p_range] = ttest2(peak_range_plus, peak_range_minus);

fprintf('\n===== TEMPORAL SHIFT ANALYSIS =====\n');
fprintf('Temporal shift correlation:\n');
fprintf('  mPFC+: mean = %.3f ± %.3f\n', mean(temporal_shift_plus), std(temporal_shift_plus));
fprintf('  mPFC-: mean = %.3f ± %.3f\n', mean(temporal_shift_minus), std(temporal_shift_minus));
fprintf('  Group difference: p = %.4f\n', p_shift);
fprintf('\nPeak latency range (s):\n');
fprintf('  mPFC+: mean = %.3f ± %.3f\n', mean(peak_range_plus), std(peak_range_plus));
fprintf('  mPFC-: mean = %.3f ± %.3f\n', mean(peak_range_minus), std(peak_range_minus));
fprintf('  Group difference: p = %.4f\n', p_range);

%% ===== VISUALIZE: INDIVIDUAL ANIMAL HEATMAPS =====
n_animals_plus = sum(mPFC_plus_idx);
n_animals_minus = sum(mPFC_minus_idx);

figure('Position', [100, 100, 1400, 800]);

% mPFC+ animals
for i = 1:n_animals_plus
    idx = find(mPFC_plus_idx, i, 'first');
    idx = idx(end);
    
    ax = subplot(2, max(n_animals_plus, n_animals_minus), i);
    imagesc(t_for_plot, 1:size(animal_results(idx).sorted_traces, 2), ...
            animal_results(idx).sorted_traces');
    colormap(ax, ap.colormap('KWG')); % Set colormap for this specific axes
    clim([-max(abs(animal_results(idx).sorted_traces(:))), max(abs(animal_results(idx).sorted_traces(:)))]);
    hold on;
    plot([0 0], ylim, 'k--', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Trials (sorted)', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f', animal_results(idx).animal_id, animal_results(idx).temporal_shift), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

% mPFC- animals
for i = 1:n_animals_minus
    idx = find(mPFC_minus_idx, i, 'first');
    idx = idx(end);
    
    ax = subplot(2, max(n_animals_plus, n_animals_minus), max(n_animals_plus, n_animals_minus) + i);
    imagesc(t_for_plot, 1:size(animal_results(idx).sorted_traces, 2), ...
            animal_results(idx).sorted_traces');
    colormap(ax, ap.colormap('KWP')); % Set colormap for this specific axes
    clim([-max(abs(animal_results(idx).sorted_traces(:))), max(abs(animal_results(idx).sorted_traces(:)))]);
    hold on;
    plot([0 0], ylim, 'k--', 'LineWidth', 2);
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Trials (sorted)', 'FontSize', 8);
    title(sprintf('M%d: r=%.2f', animal_results(idx).animal_id, animal_results(idx).temporal_shift), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

sgtitle('mPFC Activity Sorted by Lick Latency', 'FontSize', 12, 'FontWeight', 'bold');

%% ===== VISUALIZE: GROUP COMPARISON =====
figure('Position', [100, 100, 1000, 400]);

subplot(1, 2, 1);
hold on;
scatter(ones(n_animals_plus, 1), temporal_shift_plus, 100, cmPFC_plus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
scatter(2*ones(n_animals_minus, 1), temporal_shift_minus, 100, cmPFC_minus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
plot([0.5 2.5], [0 0], 'k--', 'LineWidth', 1.5);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'mPFC+', 'mPFC-'});
ylabel('Temporal Shift', 'FontSize', 12, 'FontWeight', 'bold');
title('Peak Shifts with Lick Timing?', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);

subplot(1, 2, 2);
hold on;
scatter(ones(n_animals_plus, 1), peak_range_plus, 100, cmPFC_plus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
scatter(2*ones(n_animals_minus, 1), peak_range_minus, 100, cmPFC_minus, 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1.5, 'MarkerFaceAlpha', 0.7);
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'mPFC+', 'mPFC-'});
ylabel('Peak Latency Range (s)', 'FontSize', 12, 'FontWeight', 'bold');
title('Temporal Variability', 'FontSize', 13, 'FontWeight', 'bold');
set(gca, 'FontSize', 11);



%% ===== VISUALIZE: HEATMAPS WITH PEAK OVERLAY =====
mPFC_plus_idx = [animal_results.group] == 1;
mPFC_minus_idx = [animal_results.group] == 0;

n_animals_plus = sum(mPFC_plus_idx);
n_animals_minus = sum(mPFC_minus_idx);

figure('Position', [100, 100, 1400, 800]);


% mPFC+ animals
for i = 1:n_animals_plus
    idx = find(mPFC_plus_idx, i, 'first');
    idx = idx(end);
    
    ax = subplot(2, max(n_animals_plus, n_animals_minus), i);
    imagesc(t_for_plot, 1:size(animal_results(idx).sorted_mPFC_traces, 2), ...
            animal_results(idx).sorted_mPFC_traces');
    colormap(ax, ap.colormap('KWG'));
    % clim([-max(abs(animal_results(idx).sorted_mPFC_traces(:))), ...
    %       max(abs(animal_results(idx).sorted_mPFC_traces(:)))]);
    clim([-0.026,0.026]);
    hold on;
    
    % % Overlay peak latencies
    % plot(animal_results(idx).sorted_mPFC_peaks, 1:length(animal_results(idx).sorted_mPFC_peaks), ...
    %      'r-', 'LineWidth', 2);
    % 
    % % Overlay lick latencies
    % plot(animal_results(idx).sorted_lick_lat, 1:length(animal_results(idx).sorted_lick_lat), ...
    %      'k--', 'LineWidth', 2);
    
    % Movement onset line
    plot([0 0], ylim, 'k-', 'LineWidth', 2);
    
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Trials (sorted by lick)', 'FontSize', 8);
    title(sprintf('M%d', animal_results(idx).animal_id), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

% mPFC- animals
for i = 1:n_animals_minus
    idx = find(mPFC_minus_idx, i, 'first');
    idx = idx(end);
    
    ax = subplot(2, max(n_animals_plus, n_animals_minus), max(n_animals_plus, n_animals_minus) + i);
    imagesc(t_for_plot, 1:size(animal_results(idx).sorted_mPFC_traces, 2), ...
            animal_results(idx).sorted_mPFC_traces');
    colormap(ax, ap.colormap('KWP'));
    % clim([-max(abs(animal_results(idx).sorted_mPFC_traces(:))), ...
    %       max(abs(animal_results(idx).sorted_mPFC_traces(:)))]);
    clim([-0.026,0.026]);
    hold on;
    
    % % Overlay peak latencies
    % plot(animal_results(idx).sorted_mPFC_peaks, 1:length(animal_results(idx).sorted_mPFC_peaks), ...
    %      'r-', 'LineWidth', 2);
    % 
    % % Overlay lick latencies
    % plot(animal_results(idx).sorted_lick_lat, 1:length(animal_results(idx).sorted_lick_lat), ...
    %      'k--', 'LineWidth', 2);
    % 
    % Movement onset line
    plot([0 0], ylim, 'k-', 'LineWidth', 2);
    
    xlabel('Time (s)', 'FontSize', 8);
    ylabel('Trials (sorted by lick)', 'FontSize', 8);
    title(sprintf('M%d', animal_results(idx).animal_id), 'FontSize', 9);
    xlim([-0.5 1]);
    set(gca, 'FontSize', 7);
end

sgtitle('mPFC Activity Sorted by Lick Latency (Red=mPFC peak, Black=Lick)', ...
        'FontSize', 12, 'FontWeight', 'bold');