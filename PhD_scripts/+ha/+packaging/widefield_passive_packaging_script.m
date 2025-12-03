% Saves a widefield structure the for the lcr_passive across all differnet
% protocols



%% Saves a passive data structure that includes V and Kernel activity for the lcr_passive protocol for all the different protocols

% Define protocols
protocols_to_run = {'visual_operant_lick_two_stim_right_move','visual_operant_lick_two_stim_static', 'visual_operant_lick_two_stim_right_move_big_stim',...
    'visual_operant_lick_two_stim_static_big_stim','stim_wheel*'};

animal_list = {'DS017','HA005','HA006','HA007','HA008','HA009','HA010','HA011','HA012','HA013','HA014','HA015'};


passive_data = struct();

% n SVD components
n_components = 200;

% Configuration for the passive task
stim_conditions = {90, -90, 0};
stim_labels = {'right_stim', 'left_stim', 'center_stim'};
added_time = -0.5:0.03:1;


for animal_idx = 1:numel(animal_list)

    animal = animal_list{animal_idx};

    % Assign fields: name & recordings
    passive_data(animal_idx).animal = animal;
    passive_data(animal_idx).widefield = struct();

    global_day= 1; % define global day

    for protocol_idx = 1:length(protocols_to_run)

        % get all passives (including big stim)
        current_protocol = 'lcr_passive*';

        % Load all recordings
        all_recordings = plab.find_recordings(animal, [],current_protocol);

        all_recording_dates= {all_recordings.day};
        all_recording_dates_flat= datetime(vertcat(all_recording_dates{:}), 'InputFormat','yyyy-MM-dd'); % flatten and format

        delta_recording_days = abs(days(diff(all_recording_dates_flat))); % find the day difference

        % this code removes recordings that were added much later (>month)
        threshold = 30;
        gapIdx= find(delta_recording_days > threshold); % the index where the gap exceedes

        if isempty(gapIdx) % if no gap detected; continue as usual
            all_recordings_filtered= all_recordings;

        else % if gap detected, filter based on it

            removeMask = false(size(all_recording_dates_flat));
            removeMask(1:gapIdx) = true; % create a mask

            all_recordings_filtered= all_recordings(removeMask); % filter

        end

        % Get the just the passive recordings that overlap with current
        % task

        task_recordings = plab.find_recordings(animal, [], protocols_to_run{protocol_idx});
        days_task = {task_recordings.day};
        days_passive = {all_recordings_filtered.day};
        overlap = ismember(days_passive, days_task);
        current_recordings_in_use = all_recordings_filtered(overlap);

        if isempty(current_recordings_in_use)
            fprintf('No task recording overlapping were found');
            continue
        end

        for day_idx = 1: numel(current_recordings_in_use)

            rec = current_recordings_in_use(day_idx);
            rec_day = rec.day;
            workflows = rec.workflow;  % cell array of protocol names

            % Decide which indices to process:
            if numel(workflows) > 1 && numel(unique(workflows)) == 1
                % if duplicate lcr_passive take last recording
                idxs_to_run = numel(workflows);
            else
                % if includes big_stim workflow then run both
                idxs_to_run = 1:numel(workflows);
            end

            for rec_idx = idxs_to_run
                rec_time      = rec.recording{rec_idx};
                workflow_name = workflows{rec_idx};

                fprintf( ...
                    'Animal: %s\nProtocol: %s\nDay: %d/%d (run %d of %d: %s)\n', ...
                    animal, ...
                    protocols_to_run{protocol_idx}, ...
                    day_idx, numel(current_recordings_in_use), ...
                    rec_idx, numel(workflows), workflow_name ...
                    );


                % Snapshot workspace
                vars_to_keep = who;


                % Do your loading & processing
                load_parts.widefield_master = true;
                load_parts.widefield_align= true;
                load_parts.widefield = true;
                ap.load_recording;

                if ~load_parts.widefield_align
                    error('No alignment for animal %s in %s recording day',animal,rec_day)
                end

                timelite;
                bonsai_workflow;
                ha.load_bonsai;

                % Get lick events
                lick_event_times = detect_lick_events(timelite, 'lick_spout', 4, 1);

                % Get X positions
                stimulusX_position = cat(1, trial_events.values.TrialStimX);

                for cond_idx = 1:numel(stim_conditions)
                    stim_position = stim_conditions{cond_idx};
                    stim_label = lower(strrep(stim_labels{cond_idx}, '_stim', ''));  % e.g., 'rewarded'

                    % Get stim on times and wheel movement data
                    stim_timepoints = stimOn_times(stimulusX_position == stim_position);
                    wheel_movement_times = timelite.timestamps(wheel_move == 1);

                    % Filter out trials in which there was a lick or wheel
                    % movement detected

                    added_time_for_artifact_check = 0.75;
                    movement_during_stim = false(size(stim_timepoints));
                    lick_during_stim = false(size(stim_timepoints));

                    for stim_idx = 1:length(stim_timepoints)
                        stim_start = stim_timepoints(stim_idx);
                        stim_end = stim_start + added_time_for_artifact_check;
                        movement_during_stim(stim_idx) = any(wheel_movement_times >= stim_start & wheel_movement_times <= stim_end);
                        lick_during_stim(stim_idx) = any(lick_event_times >= stim_start & lick_event_times <= stim_end);
                    end

                    stim_timepoints_filtered = stim_timepoints(~movement_during_stim & ~lick_during_stim);
                    stim_timepoints_added = stim_timepoints_filtered + added_time;
                    stim_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', stim_timepoints_added), [3,2,1]);


                    %%% regression %%%

                    % define bin edges
                    bin_edges = [wf_t', wf_t(end)];

                    % Bin stimulus onsets
                    stim_regressor = histcounts(stim_timepoints, bin_edges);

                    % define time shifts
                    time_shifts = {-30:10};

                    % predict the stimulus presentation based on the wf_V
                    [kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
                        ap.regresskernel(wf_V(1:n_components,:),stim_regressor,time_shifts,100);

                    % Assign data for the current day

                    % Assign the data into structure seperate
                    % whether big stim or not - have empty fields for
                    % field size consistancy across animals

                    if ~contains(workflow_name,'size60')
                        passive_data(animal_idx).widefield(global_day).([stim_label '_stim_aligned_V'])= stim_aligned_V;
                        passive_data(animal_idx).widefield(global_day).([stim_label '_stim_aligned_kernel'])= kernels;

                        if numel(unique(workflows)) == 1 % if only one recording, then assign empty field to the other condition

                            passive_data(animal_idx).widefield(global_day).([stim_label '_big_stim_aligned_V'])= [];
                            passive_data(animal_idx).widefield(global_day).([stim_label '_big_stim_aligned_kernel'])= [];
                        end


                    else
                        passive_data(animal_idx).widefield(global_day).([stim_label '_big_stim_aligned_V'])= stim_aligned_V;
                        passive_data(animal_idx).widefield(global_day).([stim_label '_big_stim_aligned_kernel'])= kernels;

                        if numel(unique(workflows)) == 1 %if only one recording, then assign empty field to the other condition
                            passive_data(animal_idx).widefield(global_day).([stim_label '_stim_aligned_V'])= [];
                            passive_data(animal_idx).widefield(global_day).([stim_label '_stim_aligned_kernel'])= [];
                        end

                    end

                end

            end

            global_day= global_day+1; % increment
            % clear the day specific variables
            clearvars('-except',vars_to_keep{:});
        end

    end
end


% Construct the dynamic filename based on the current protocol
filename = sprintf('passive_data_all_animals.mat');

%Save the data to the dynamically generated filename
save(filename, 'passive_data','-v7.3');

%% Script to quickly test the packaging data

load('passive_data_all_animals.mat');

% load master_U
wf_U= plab.wf.load_master_U;

single_animal_data= cat(2,passive_data(10).widefield);
right_stim_data= cat(3,single_animal_data(1:10).right_stim_aligned_kernel);

mean_right_stim_data= nanmean(right_stim_data,3);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:n_components),mean_right_stim_data));
axis equal;

