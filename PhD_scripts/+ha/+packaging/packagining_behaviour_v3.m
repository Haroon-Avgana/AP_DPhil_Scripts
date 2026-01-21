% Script that:
% 1. Package the behaviour across all animals and protocols and saves it

% 2. Creates a variable that contains all recording days for all protocols
% and whether they were significant [0,1], and saves it as combined_sig_day_all_protocols

test
%% Packages behaviour as a data structure:

% animal ID
% Workflow
% CS labels (contains a logical map where CS+ =1)
% wheel velocity aligned to stim onset for now

% For right move protocol:
% Rewarded stim on/off times and center times
% Non rewarded stim on/off and putative center times
% Lick events
% ITI licks (defined based on lick bouts)
% Reward times (excluding manuals)
% average PSTH aligned to stim on/stim center
% latency to lick from stim onset
% diff from optimal reward

% Latency to lick from stim movement


% For static protocol:
% Rewarded stim on/off times
% Non rewarded stim on/off
% Lick events
% ITI licks
% Reward times (excluding manuals)
% average PSTH aligned to stim on/off
% average PSTH aligned to first since stim off (for CS-)
% latency to lick from stim onset
% diff from optimal reward

% For wheel protocol:
% Stim on times
% time to move from stim onset

% Also saves reward rates



%% Runs and saves a behaviour structure

% Define animals and protocols
animal_list = {'DS017','HA005','HA006','HA007','HA008','HA009','HA010','HA011','HA012','HA013','HA014','HA015',...
    'AP030','AP031','AP032'};


% animal_list= {'HA016','HA017','HA018','HA019','HA020'};



protocol_list={'visual_operant_lick_two_stim_right_move','visual_operant_lick_two_stim_static','visual_operant_lick_two_stim_right_move_big_stim'...
    ,'visual_operant_lick_two_stim_static_big_stim' ,'stim_wheel'};


% Initialize the nested structure
behaviour_data = struct();
for animal_idx = 1:numel(animal_list)
    behaviour_data(animal_idx).animal_id = animal_list{animal_idx};
    behaviour_data(animal_idx).day_counter = 0; % Initialize day counter
end

for animal_idx = 1:numel(animal_list)

    animal = animal_list{animal_idx};
    % Get all recording days to filter out when there are big gaps in
    % recording
    all_recordings = plab.find_recordings(animal);

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

    for protocol_idx=1:numel(protocol_list)

        % assign name to the data
        behaviour_data(animal_idx).animal_id=animal;

 
        % define current protocol
        current_protocol = protocol_list{protocol_idx};

        % Create a protocol max to extract only the days the correspond to
        % current protocol

        % If protocol is wheel go to a cell function with contain to
        % extract all protocols with that string - to view all stim_wheel
        % protocols as one
        if contains(current_protocol,'stim_wheel')

            protocol_mask = false(length(all_recordings_filtered), 1);

            for i = 1:length(all_recordings_filtered)
                workflow = all_recordings_filtered(i).workflow;

                % Flatten cell array and check if any element contains the string
                if iscell(workflow)
                    % Convert all cells to strings and check
                    workflow_strings = workflow(cellfun(@(x) ischar(x) || isstring(x), workflow));
                    if any(contains(workflow_strings, current_protocol))
                        protocol_mask(i) = true;
                    end
                elseif ischar(workflow) || isstring(workflow)
                    if contains(workflow, current_protocol)
                        protocol_mask(i) = true;
                    end
                end
            end

            filtered_recordings = all_recordings_filtered(protocol_mask);

            % if not - then only look for specific matching string (to
            % avoid big stim protocols being joint with regular ones)
        else

            protocol_mask = cellfun(@(w) any(strcmp(w, current_protocol)), {all_recordings_filtered.workflow}, 'UniformOutput', true);
            filtered_recordings = all_recordings_filtered(protocol_mask);
        end

        % continue if no recording is found
        if isempty(filtered_recordings)
            fprintf('No task recording for protocol %s',current_protocol);
            continue
        end


        % get all current task recordings
        all_recordings_curr_protocol = filtered_recordings;


        for day_idx =1:size(all_recordings_curr_protocol, 2)

            rec = all_recordings_curr_protocol(day_idx); % get all recordings including passive
            rec_day = rec.day;

            % Filter out empty cells and then get only matching protocl
            % index
            valid_mask = cellfun(@(x) ~isempty(x) && contains(x, current_protocol), rec.workflow);

            if ~any(valid_mask)
                warning('Protocol %s not found for day %d', current_protocol, rec_day);
                continue;
            end

            % Get last matching index
            curr_protocol_idx = find(valid_mask, 1, 'last');
            rec_time = rec.recording{curr_protocol_idx};

            fprintf('Animal: %s \nCurrent Potocol :%s \nDay: %d/%d\n', animal, current_protocol, day_idx, numel(all_recordings_curr_protocol));


            % Load recording
            load_parts.mousecam = false; % don't load camera
            ap.load_recording;
            timelite;
            bonsai_workflow;
            ha.load_bonsai;

            behaviour_data(animal_idx).day_counter = behaviour_data(animal_idx).day_counter + 1;
            day_counter = behaviour_data(animal_idx).day_counter;

            % If protocol includes two stim extract the relevant session
            % parameters
            if contains(current_protocol,{'two_stim'})

                % Extract the rewarded X position
                rewardedX = trial_events.parameters.RewardedX;

                % Extract possible X positions
                possible_X_positions = trial_events.parameters.TrialXSpace;

                % Define non-rewarded X positions (excluding the rewardedX)
                non_rewarded_X = setdiff(possible_X_positions, rewardedX);
                
                % Get all trial X positions
                all_trials_x_position= [trial_events.values.TrialX]';

                % Create a logical for CS+
                is_cs_plus= all_trials_x_position == rewardedX;
                
                % Make sure the logical mask is equivilent to the actual n
                % valid trials 
                if length(is_cs_plus) ~= n_valid
                    is_cs_plus= is_cs_plus(1:min(length(is_cs_plus),n_valid));
                end

     
                % Calculate CS+ center and CS- putative center times and
                % combine them

                % Calculate putative center times for both trial types
                cs_plus_center = stimOff_times(is_cs_plus);
                cs_minus_putative_center = stimOff_times(~is_cs_plus) + trial_events.parameters.StimMoveTime;

                % Combine maintaining trial order
                all_stim_center_times = nan(size(is_cs_plus));
                all_stim_center_times(is_cs_plus) = cs_plus_center;
                all_stim_center_times(~is_cs_plus) = cs_minus_putative_center;

                % find rewarded and non rewarded times
                cs_plus_on= stimOn_times(is_cs_plus);
                cs_minus_on= stimOn_times(~is_cs_plus);

                % Final the final position times
                cs_plus_final_position= stimOff_times(is_cs_plus);
                cs_minus_final_position= stimOff_times(~is_cs_plus);

               

                % Calculate the putative left stim in center times (stim off + time to move to center)
                non_rewarded_stim_putative_center_times= cs_minus_final_position + trial_events.parameters.StimMoveTime;

                % Determine lick spout channel and find timelite values
                lick_spout_channel_index = find(contains({timelite.daq_info.channel_name}, 'lick_spout'), 1);
                lick_spout_values = timelite.data(:, lick_spout_channel_index);

                % Define thresholds for detecting licking behavior
                high_lick_detection_threshold = 4;
                low_lick_detection_threshold = 1;

                % Find lick events based on the thresholds
                lick_events_indices = find(lick_spout_values(2:end) > high_lick_detection_threshold & ...
                    lick_spout_values(1:end-1) < low_lick_detection_threshold);

                % licks in timelite
                lick_event_times= timelite.timestamps(lick_events_indices);


                lick_detection_time_window = [-5, 5]; % Time window in seconds
                bin_size = 0.15; % Size of time bins in seconds

                % Define time bins for PSTH
                time_bins = lick_detection_time_window(1):bin_size:lick_detection_time_window(2);

                % define window
                anticipatory_window = 1;

                % Initalize
                anticipatory_licks_all_trials= nan(length(stimOn_times), 1);

                if contains(current_protocol,{'right_move'})

                    % find times aligned to CS+ movement and CS- stim off
                    cs_plus_move_times= rewarded_stim_start_move_times; % for bonsai script

                    rewarded_stim_on_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_plus_on, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    rewarded_stim_in_center_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_plus_center, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    non_rewarded_stim_on_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_minus_on, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    non_rewarded_stim_putative_center_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_minus_putative_center, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    % Find the anticipatory licks aligned to stim center
                    % times for both CS+ and CS- (putative)

                    % Shift the end of the anticipatory window back by 100 ms
                    anticipatory_offset = 0.1;  % 100 ms in seconds


                    % for now calculate anticipatory licks for stim off for left and
                    % not putative center time

                    for trial=1:n_valid
                        window_start= all_stim_center_times(trial)- anticipatory_window;
                        window_end= all_stim_center_times(trial)-anticipatory_offset;

                        % Count licks in the anticipatory window
                        licks_in_window = lick_event_times >= window_start & ...
                            lick_event_times < window_end;
                        anticipatory_licks_all_trials(trial) = sum(licks_in_window);

                    end
                    
                    % Compute the first lick after CS- off  - accidently
                    % used putative center but this currently executes
                    % based on CS- stim off

                    % Initalize variable to store times
                    first_lick_post_cs_putative_center_times= nan(1,numel(cs_minus_final_position));  %non_rewarded_stim_putative_center_times
  
                    for non_rewarded_idx=1:length(cs_minus_final_position) 

                        curr_non_rewarded_off= cs_minus_final_position(non_rewarded_idx); % grab current putative center time
                        first_lick_from_putative_center_idx= find((lick_event_times - curr_non_rewarded_off)>0,1,'first'); % get first lick from that time forward

                        if isempty(first_lick_from_putative_center_idx) % if no licks detected keep it nan
                            first_lick_post_cs_putative_center_times(non_rewarded_idx)= nan;

                        else
                            first_lick_post_cs_putative_center_times(non_rewarded_idx)= lick_event_times(first_lick_from_putative_center_idx);
                        end
                    end

                   
                    % Initialize combined variable for CS+ and CS- reward
                    % and putative reward times

                    all_reward_times = nan(size(is_cs_plus));
                    cs_plus_indices = find(is_cs_plus);
                    num_cs_plus_trials = length(cs_plus_indices);
                    num_actual_rewards = length(reward_times_task);


                    if num_actual_rewards < num_cs_plus_trials
                        warning('Session terminated early: %d CS+ trials missing rewards', ...
                            num_cs_plus_trials - num_actual_rewards);

                        % Only fill the CS+ trials that have corresponding rewards
                        cs_plus_with_rewards = cs_plus_indices(1:num_actual_rewards);
                        all_reward_times(cs_plus_with_rewards) = reward_times_task;

                        % The remaining CS+ trials stay as NaN
                    else
                        % Normal case: all CS+ trials have rewards
                        all_reward_times(is_cs_plus) = reward_times_task;
                    end

                    % Fill in CS- trials (first lick after putative center)
                    all_reward_times(~is_cs_plus) = first_lick_post_cs_putative_center_times';

               % Cacluclate PSTH aligned to reward times

               reward_time_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                   lick_event_times, ...
                   all_reward_times(is_cs_plus), ...
                   time_bins, ...
                   'SmoothSigma', 0, ...     % no smoothing
                   'Return','rate');         % return avg as rate (Hz or 1/s)

                    % Compute the difference from optimal reward for
                    % rewarded and non-rewarded stims 

                    all_stim_latency_structure= ha.helper_func.compute_latency_vs_optimal (...
                        lick_event_times,...
                        stimOn_times,...
                        all_reward_times,...
                        all_stim_center_times);

                        % Compute the first lick after stim in center for CS+
                    first_lick_post_cs_move_times= nan(1,numel(cs_plus_move_times))';

                    for cs_plus_move_idx=1:length(cs_plus_move_times)
                        curr_cs_plus_move_time= cs_plus_move_times(cs_plus_move_idx);
                        first_lick_from_center_idx= find((lick_event_times - curr_cs_plus_move_time)>0,1,'first'); % get first lick from that time forward

                        if isempty(first_lick_from_center_idx) % if no licks detected keep it nan
                            first_lick_post_cs_move_times(cs_plus_move_idx)= nan;

                        else
                            first_lick_post_cs_move_times(cs_plus_move_idx)= lick_event_times(first_lick_from_center_idx);
                        end
                    end
                    
                    % Calculate the latency
                    first_lick_post_cs_move_latency= first_lick_post_cs_move_times-cs_plus_move_times;
     
                    % Calculate #n licks during ITI 

                    % Store data under the incremented counter:

                    %PSTH data:
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_rewarded_stim_on = rewarded_stim_on_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_rewarded_stim_final_position= rewarded_stim_in_center_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_non_rewarded_stim_on= non_rewarded_stim_on_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_non_rewarded_stim_final_position= non_rewarded_stim_putative_center_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_reward_times= reward_time_PSTH_out.avgCount;

                    % Anticipatory licks:
                    behaviour_data(animal_idx).recording_day(day_counter).anticipatory_licks= anticipatory_licks_all_trials;

                    % RT:
                    behaviour_data(animal_idx).recording_day(day_counter).all_stim_latency_to_first_lick= all_stim_latency_structure.latency_to_first_lick;
                    behaviour_data(animal_idx).recording_day(day_counter).all_stim_diff_from_optimal_reward= all_stim_latency_structure.diff_from_optimal_reward;
                    behaviour_data(animal_idx).recording_day(day_counter).cs_plus_latency_lick_to_move= first_lick_post_cs_move_latency;

                    % behaviour_data(animal_idx).recording_day(day_counter).rewarded_stim_latency_lick_stim_onset= cs_plus_diff_latency_structure.latency_filtered;
                    % behaviour_data(animal_idx).recording_day(day_counter).rewarded_stim_diff_from_optimal_reward= cs_plus_diff_latency_structure.diff_filtered;
                    % behaviour_data(animal_idx).recording_day(day_counter).non_rewarded_stim_latency_lick_stim_onset= cs_minus_diff_latency_structure.latency_filtered;
                    % behaviour_data(animal_idx).recording_day(day_counter).non_rewarded_stim_diff_from_optimal_reward= cs_minus_diff_latency_structure.diff_filtered;


                elseif contains(current_protocol,{'static'})
                    % For static compute the psth aligned to reward times

                    % 
                    rewarded_stim_on_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_plus_on, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    rewarded_stim_off_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        reward_times_task, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                
                    non_rewarded_stim_on_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        cs_minus_on, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

          

                    %%%% Compute the first lick after CS- stim off %%%%

                    % Initalize variable to store
                    first_lick_post_cs_minus_off_times= nan(1,numel(cs_minus_final_position));

                    % From each non-rewarded stim off - calculate the duration until next stim on
                    % Calculate the licks during that period and find the first lick
                    for non_rewarded_idx = 1:length(cs_minus_final_position)
                        curr_non_rewarded_off = cs_minus_final_position(non_rewarded_idx);
                        next_stim_on_idx = find((stimOn_times - curr_non_rewarded_off) > 0, 1, "first");

                        if isempty(next_stim_on_idx) % if no stim on afterwards, skip
                            continue;
                        end

                        next_stim_on_times = stimOn_times(next_stim_on_idx);
                        first_lick_in_duration_idx = find((curr_non_rewarded_off <= lick_event_times) & ...
                            (lick_event_times < next_stim_on_times), 1, 'first');

                        if isempty(first_lick_in_duration_idx) % if no licks detected keep it nan
                            first_lick_post_cs_minus_off_times(non_rewarded_idx) = nan;
                        else
                            first_lick_post_cs_minus_off_times(non_rewarded_idx) = lick_event_times(first_lick_in_duration_idx);
                        end
                    end


                  
                    non_rewarded_first_lick_stim_off_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        first_lick_post_cs_minus_off_times, ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)


                    % combine reward times with first lick post CS- off
                    all_reward_times= nan(size(stimOff_times));

                    % A janky solution for a particular recording of HA08 
                    if sum(is_cs_plus) > length(reward_times_task)
                        reward_times_task(end+1)= NaN;
                    end
                   
                    all_reward_times(is_cs_plus)= reward_times_task;
                    all_reward_times(~is_cs_plus)= first_lick_post_cs_minus_off_times';
           


                    % Calculate anticipatory licks based on stim Off

                    % Shift the end of the anticipatory window back by 100 ms
                    anticipatory_offset = 0.1;  % 100 ms in seconds

                    % for now calculate anticipatory licks for stim off for left and
                    % not putative center time

                    for trial=1:length(stimOff_times)
                        window_start= stimOff_times(trial)- anticipatory_window;
                        window_end= stimOff_times(trial)-anticipatory_offset;

                        % Count licks in the anticipatory window
                        licks_in_window = lick_event_times >= window_start & ...
                            lick_event_times < window_end;
                        anticipatory_licks_all_trials(trial) = sum(licks_in_window);

                    end


                    % Calculate the static for all trials

                    % Get the static times for all trials
                    static_times= [trial_events.values.TrialStimStaticTime]';


                    % make it the size of rewarded stim trials as static
                    % can be drawn before stim is on
                    if length(static_times) ~= length(stimOn_times)
                         actual_n_trials= min(numel(static_times), numel(stimOn_times));
                        static_times= static_times(1:actual_n_trials);
                        stimOn_times= stimOn_times(1:actual_n_trials);
                    end
                    
                     % Add the static to stim onset to find the reward
                    % avalible times
                    reward_avalible_times_all= nan(1,numel(stimOn_times))';

                    for trial_idx=1:length(reward_avalible_times_all)
                        reward_avalible_times_all(trial_idx)= stimOn_times(trial_idx)+ static_times(trial_idx);
                    end


                    % Combine the reward times CS+ and first lick post CS-
                    all_first_lick_post_stim = nan(n_valid, 1);

                    % Fill the CS- values
                    all_first_lick_post_stim(~is_cs_plus) = first_lick_post_cs_minus_off_times;

                    % Fill the CS+ values based on reward times
                    % reward_times_task = CS+ rewards only (excluding manual)
                    % reward_times = all rewards (including manual)
                    n_cs_plus = sum(is_cs_plus);

                    if length(reward_times_task) == n_cs_plus
                        % Preferred: use task rewards (excludes manual rewards)
                        all_first_lick_post_stim(is_cs_plus) = reward_times_task;
                    elseif length(reward_times) == n_cs_plus
                        % Fallback: use all rewards (includes manual rewards)
                        all_first_lick_post_stim(is_cs_plus) = reward_times;
                    elseif length(reward_times_task) == n_cs_plus - 1
                        % Missing one reward in reward_times_task - pad with NaN and align
                        warning('reward_times_task has %d rewards but %d CS+ trials. Padding with NaN.', ...
                            length(reward_times_task), n_cs_plus);
                        padded_rewards = [reward_times_task; nan];
                        all_first_lick_post_stim(is_cs_plus) = padded_rewards;

                    end
                    


                    all_stim_latency_structure= ha.helper_func.compute_latency_vs_optimal (...
                        lick_event_times,...
                        stimOn_times,...
                        all_first_lick_post_stim,...
                        reward_avalible_times_all);

                    % Cacluclate PSTH aligned to reward times

                    reward_time_PSTH_out = ha.helper_func.compute_lick_PSTH( ...
                        lick_event_times, ...
                        all_first_lick_post_stim(is_cs_plus), ...
                        time_bins, ...
                        'SmoothSigma', 0, ...     % no smoothing
                        'Return','rate');         % return avg as rate (Hz or 1/s)

                    % Add to the structure

                    % PSTHs:
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_rewarded_stim_on = rewarded_stim_on_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_rewarded_stim_final_position= rewarded_stim_off_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_non_rewarded_stim_on= non_rewarded_stim_on_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_non_rewarded_stim_final_position= non_rewarded_first_lick_stim_off_PSTH_out.avgCount;
                    behaviour_data(animal_idx).recording_day(day_counter).avg_psth_reward_times= reward_time_PSTH_out.avgCount;
                    
                    % Anticipatory licks:
                    behaviour_data(animal_idx).recording_day(day_counter).anticipatory_licks= anticipatory_licks_all_trials;

                    % RT:
                    behaviour_data(animal_idx).recording_day(day_counter).all_stim_latency_to_first_lick= all_stim_latency_structure.latency_to_first_lick;
                    behaviour_data(animal_idx).recording_day(day_counter).all_stim_diff_from_optimal_reward= all_stim_latency_structure.diff_from_optimal_reward;
                end

                %%% Start of ITI Analysis %%

                % Create a variable that includes both CS- time off and CS+ reward
                % times

                reward_times_task= reward_times_task(~isnan(reward_times_task));

                % Preallocate combined times variable
                combined_times = nan(size(stimOff_times));

                % actual trials
                num_csplus_actual = length(reward_times_task);


                % Assign reward times to the first min_csplus_trials CS+ trials
                csplus_indices = find(is_cs_plus);
                %combined_times(csplus_indices(1:min_csplus_trials)) = reward_times_task(1:min_csplus_trials);

                % For CS− trials (not CS+), use stimOff_times directly
                combined_times(~is_cs_plus) = stimOff_times(~is_cs_plus);


                % Define reward offsets (assume reward duration known or constant)

                % Parameters
                bout_thresh = 0.5;  % 500 ms maximum inter‐lick interval within a bout

                % Preallocate
                n_rewards = length(reward_times_task);
                reward_duration = nan(n_rewards,1);

                for r = 1:n_rewards
                    rt = reward_times_task(r);

                    % 1) Find all licks *after* this reward
                    idx = find(lick_event_times > rt);
                    if isempty(idx)
                        continue;
                    end
                    t_after = lick_event_times(idx);

                    % 2) Compute inter‐lick intervals for that sequence
                    dts = diff(t_after);

                    % 3) Find the first gap > bout_thresh
                    gap_idx = find(dts > bout_thresh, 1);

                    % 4) Determine end of bout
                    if isempty(gap_idx)
                        % no large gap: take the last lick
                        bout_end = t_after(end);
                    else
                        % bout ends at the lick *before* that gap
                        bout_end = t_after(gap_idx);
                    end

                    % 5) Reward‐aligned bout duration
                    reward_duration(r) = bout_end - rt;
                end

                % Add the reward durations
                combined_times(csplus_indices(1:num_csplus_actual)) = reward_times_task + reward_duration;

                % grab the ITI start and end points
                ITI_start = nan(length(combined_times),1);
                ITI_end   = nan(length(combined_times),1);

               ITI_actual_duration=nan(length(combined_times),1);

                for i = 1:length(combined_times)-1
                    ITI_start(i) = combined_times(i);
                    ITI_end(i)   = stimOn_times(i+1);
                    ITI_actual_duration(i)= ITI_end(i)-ITI_start(i);
                end

           

                % initalize lick count
                ITI_lick_counts = nan(length(combined_times),1);

                for i = 1:length(combined_times)
                    if ~isnan(ITI_start(i)) && ~isnan(ITI_end(i))
                        ITI_lick_counts(i) = sum( ...
                            lick_event_times >= ITI_start(i) & ...
                            lick_event_times <  ITI_end(i) );
                    end
                end

                % Calculate the #n licks that occured during the static
                % time

                % Get the static times for all trials
                static_times= [trial_events.values.TrialStimStaticTime]';

                % sanity check that static times and stim on times are
                % equal
                actual_n_trials= min(numel(static_times), numel(stimOn_times));

                if length(static_times) ~= length(stimOn_times)

                    static_times= static_times(1:actual_n_trials);
                    stimOn_times= stimOn_times(1:actual_n_trials);
                end

                static_time_lick_counts= nan(length(static_times),1);

                for trial_idx=1:length(static_time_lick_counts)
                    if ~isnan(ITI_start(trial_idx)) && ~isnan(ITI_end(trial_idx))
                        static_time_lick_counts(trial_idx) = sum( ...
                            lick_event_times >= stimOn_times(trial_idx) & ...
                            lick_event_times <  stimOn_times(trial_idx)+ static_times(trial_idx) );
                    end
                end


                % Find reward time indecies
                [~, reward_indices] = ismember(reward_times_task, timelite.timestamps);
                reward_times_timestamp= timelite.timestamps(reward_indices);

                % From each reward osnset add the reward duration
                reward_duration_in_times= reward_times_task + reward_duration;

                % Find the indecies
                [~, reward_duration_indices] = ismember(reward_duration_in_times, timelite.timestamps);

                % Remove licks during reward (keep the rewarded lick itself)
                lick_spout_values_cleaned = lick_spout_values;

                for i = 1:length(reward_indices)
                    % Zero out from reward_index+1 to reward_duration_index (excludes the rewarded lick)
                    if reward_duration_indices(i) > reward_indices(i)
                        lick_spout_values_cleaned((reward_indices(i)+1):reward_duration_indices(i)) = 0;
                        % Or use NaN: licklickSpout1_cleaned((reward_indices(i)+1):reward_duration_indices(i)) = NaN;
                    end
                end

                % Calculate the lick times after removal of post-reward lick

                % Extract lick times (threshold crossings)
                lick_event_times_cleaned = timelite.timestamps(lick_spout_values_cleaned(2:end) > high_lick_detection_threshold & ...
                    lick_spout_values_cleaned(1:end-1) < low_lick_detection_threshold);

                % Parameters
                ILI_time_bins = 0:0.5:20;

                inter_lick_interval = diff(lick_event_times_cleaned);

                bout_threshold = 0.3;
                bout_starts_idx = [1; find(inter_lick_interval > bout_threshold) + 1];
                bout_start_times = lick_event_times_cleaned(bout_starts_idx);

                % Inter-bout intervals
                inter_bout_intervals = diff(bout_start_times);

                % Compute normalized histogram
                inter_bout_hist= histcounts(inter_bout_intervals, ILI_time_bins, ...
                    'Normalization', 'probability');



                % Calculate the inter-bout interval during ITI

                % Remove NaNs and find ITI start/end indcies
                ITI_start= ITI_start(~isnan(ITI_start));
                ITI_end= ITI_end(~isnan(ITI_end));

                % If ITI_start/end are guaranteed to be in timelite.timestamps
                ITI_lick_values = zeros(size(lick_spout_values_cleaned));

                for ITI_idx = 1:length(ITI_start)
                    % Logical mask for this ITI period
                    in_ITI = timelite.timestamps >= ITI_start(ITI_idx) & ...
                        timelite.timestamps <= ITI_end(ITI_idx);

                    ITI_lick_values(in_ITI) = lick_spout_values_cleaned(in_ITI);
                end

                % Extract lick times (threshold crossings)
                lick_event_times_ITI = timelite.timestamps(ITI_lick_values(2:end) > high_lick_detection_threshold & ...
                    ITI_lick_values(1:end-1) < low_lick_detection_threshold);

           
                inter_lick_interval_ITI = diff(lick_event_times_ITI);

                bout_starts_idx_ITI = [1; find(inter_lick_interval_ITI > bout_threshold) + 1];
                bout_start_times_ITI = lick_event_times_ITI(bout_starts_idx_ITI);

                % Inter-bout intervals
                inter_bout_intervals_ITI = diff(bout_start_times_ITI);

                ILI_time_bins = 0:0.1:20;
                % Compute normalized histogram
                inter_bout_hist= histcounts(inter_bout_intervals, ILI_time_bins, ...
                    'Normalization', 'probability');

                % Autocorrelation of Inter-licks during ITI

                [acf, lags] = xcorr(inter_lick_interval - mean(inter_lick_interval), 15, 'coeff');  % 15 lags
                lags_time = lags * mean(inter_lick_interval);  % Convert lags to time


                % get the drawn ITI time
                ITI_drawn_duration= [trial_events.values.TrialITI]';
                ITI_drawn_duration= ITI_drawn_duration(1:actual_n_trials); % make sure it matches

                % add general fields to the structure
                behaviour_data(animal_idx).recording_day(day_counter).right_stim_on_times= cs_plus_on;
                behaviour_data(animal_idx).recording_day(day_counter).right_stim_off_times= cs_plus_final_position;
                behaviour_data(animal_idx).recording_day(day_counter).reward_times= reward_times_task;

                behaviour_data(animal_idx).recording_day(day_counter).cs_labels= is_cs_plus;
                behaviour_data(animal_idx).recording_day(day_counter).ITI_lick_counts= ITI_lick_counts;
                behaviour_data(animal_idx).recording_day(day_counter).ITI_drawn_duration= ITI_drawn_duration;
                behaviour_data(animal_idx).recording_day(day_counter).ITI_actual_duration= ITI_actual_duration;
                behaviour_data(animal_idx).recording_day(day_counter).all_static_times= static_times;
                behaviour_data(animal_idx).recording_day(day_counter).static_time_lick_counts= static_time_lick_counts;

                behaviour_data(animal_idx).recording_day(day_counter).all_reward_times_including_first_lick_post_cs_minus= all_reward_times;
                behaviour_data(animal_idx).recording_day(day_counter).left_stim_on_times= cs_minus_on;
                behaviour_data(animal_idx).recording_day(day_counter).left_stim_center_times= non_rewarded_stim_putative_center_times;
                behaviour_data(animal_idx).recording_day(day_counter).left_stim_off_times= cs_minus_final_position;
                behaviour_data(animal_idx).recording_day(day_counter).lick_event_times= lick_event_times;

                behaviour_data(animal_idx).recording_day(day_counter).inter_lick_bout_histogram= inter_bout_hist;
                behaviour_data(animal_idx).recording_day(day_counter).inter_lick_interval_auto_correlation= acf;
                behaviour_data(animal_idx).recording_day(day_counter).inter_lick_interval_auto_correlation_lag_times= lags_time;

            elseif contains(current_protocol,{'stim_wheel'})

                cs_plus_on= stimOn_times;

                behaviour_data(animal_idx).recording_day(day_counter).right_stim_on_times= cs_plus_on;
                behaviour_data(animal_idx).recording_day(day_counter).right_stim_to_move= stim_to_move;
                behaviour_data(animal_idx).recording_day(day_counter).reward_times= reward_times_task;
            end


            % Calculate the reward rate 
            total_task_time_minutes= length(timelite.timestamps)./ timelite.daq_info(1).rate ./60; % get in minutes
            total_rewards= length([trial_events.values.Outcome]); % n rewarded trials

            reward_rate = total_rewards / total_task_time_minutes; % rewards per minute


            % Save the wheel velocity aligned to stim onset -5 +5s
            surround_time = [-5,5];
            surround_sample_rate = 100;
            surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

            % Align wheel movement to stim onset & reward times 
            stim_on_aligned_wheel_vel = interp1(timelite.timestamps, ...
                wheel_velocity,stimOn_times + surround_time_points);

            reward_aligned_wheel_vel = interp1(timelite.timestamps, ...
                wheel_velocity,reward_times_task + surround_time_points);

            % save general
            behaviour_data(animal_idx).recording_day(day_counter).day= rec_day;
            behaviour_data(animal_idx).recording_day(day_counter).workflow= current_protocol;
            behaviour_data(animal_idx).recording_day(day_counter).reward_rate= reward_rate;

            behaviour_data(animal_idx).recording_day(day_counter).stim_on_aligned_wheel_vel= stim_on_aligned_wheel_vel;
            behaviour_data(animal_idx).recording_day(day_counter).reward_aligned_wheel_vel= reward_aligned_wheel_vel;
        end
    end
end

% Construct the dynamic filename based on the current protocol
filename = sprintf('behaviour_structure_all_animals.mat');

% Save the data to the dynamically generated filename
save(filename, 'behaviour_data','-v7.3');

%% Filter out big_stim protocols from the behaviour, and the regular sized two stim protocols that occured after big stim

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
% opts.following_names = { 'stim_wheel*' }; % specific for HA014 for now


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

    fprintf('Animal %s: kept %d/%d sessions (%s)\n', ...
        animal_list{animal_idx}, numel(keptLocalIdx), numel(rd));
end

%% There is a specific problem with the organization of protocols with HA014 - look into it next time you package the data

% Define the animal to reorganize
target_animal_id = 'HA014';  % Change this to your target animal

% Find the animal index
all_ids = {behaviour_data.animal_id};
animal_idx = find(strcmp(all_ids, target_animal_id));

if isempty(animal_idx)
    error('Animal %s not found in behaviour_data', target_animal_id);
end

% Define the desired workflow order
workflow_order = {
    'visual_operant_lick_two_stim_right_move';
    'visual_operant_lick_two_stim_static';
    'visual_operant_lick_two_stim_right_move_big_stim',
    'visual_operant_lick_two_stim_static_big_stim',
    'stim_wheel*';
};


% Get all recording days for this animal
recording_days = behaviour_data(animal_idx).recording_day;

% Create cell array to hold reorganized days
reorganized_days_cell = {};

% Loop through desired workflow order
for w = 1:length(workflow_order)
    target_workflow = workflow_order{w};
    
    % Find all days matching this workflow
    matching_days_idx = [];
    for d = 1:length(recording_days)
        if isfield(recording_days(d), 'workflow') && ...
           strcmp(recording_days(d).workflow, target_workflow)
            matching_days_idx = [matching_days_idx; d];
        end
    end
    
    % Append matching days to cell array
    if ~isempty(matching_days_idx)
        for i = 1:length(matching_days_idx)
            reorganized_days_cell{end+1} = recording_days(matching_days_idx(i));
        end
        fprintf('Workflow "%s": Found %d days\n', target_workflow, length(matching_days_idx));
    else
        warning('Workflow "%s": No days found', target_workflow);
    end
end

% Convert cell array to struct array
if ~isempty(reorganized_days_cell)
    reorganized_days = [reorganized_days_cell{:}];
else
    reorganized_days = [];
end

% Check if any days were left out
n_original = length(recording_days);
n_reorganized = length(reorganized_days);

if n_original ~= n_reorganized
    warning('Original had %d days, reorganized has %d days. Some days may have been excluded.', ...
            n_original, n_reorganized);
    
    % Find which days were excluded
    original_workflows = arrayfun(@(d) d.workflow, recording_days, 'UniformOutput', false);
    reorganized_workflows = arrayfun(@(d) d.workflow, reorganized_days, 'UniformOutput', false);
    
    fprintf('\nOriginal workflows:\n');
    disp(unique(original_workflows));
    fprintf('\nReorganized workflows:\n');
    disp(unique(reorganized_workflows));
end

% Replace the recording_day field with reorganized version
behaviour_data(animal_idx).recording_day = reorganized_days;

fprintf('\n✓ Successfully reorganized animal %s\n', target_animal_id);
fprintf('  Total days: %d\n', length(reorganized_days));

% Display summary
fprintf('\nReorganized day sequence:\n');
for d = 1:length(reorganized_days)
    fprintf('  Day %d: %s\n', d, reorganized_days(d).workflow);
end

%% Permutation test to save the significant days (one-sided)

% Compares CS minus- CS- in anticipatory licks and then plots the
% observed difference - bonferroni correction for multiple comparisons

% Creates a variable packed_sig_days which is cell array (n animal x 1)
% containing the two stim workflow significant days

n_animals = length(behaviour_data);
two_stim_data = struct([]);

% Step 1: Extract recordings with 'two_stim' in workflow name
for animal_idx = 1:n_animals
    animal_data = behaviour_data(animal_idx);
    valid_days = arrayfun(@(rd) contains(rd.workflow, 'two_stim'), animal_data.recording_day);

    if any(valid_days)
        filtered_days = animal_data.recording_day(valid_days);
        two_stim_data(animal_idx).recording_day = filtered_days;
        two_stim_data(animal_idx).animal_id = animal_data.animal_id;
    else
        two_stim_data(animal_idx).recording_day = [];
        two_stim_data(animal_idx).animal_id = animal_data.animal_id;
    end
end

% Step 2: Run permutation tests
packed_sig_days = cell(n_animals, 1);
num_perms = 5000;

for animal_idx = 1:length(two_stim_data)
    for day_idx = 1:length(two_stim_data(animal_idx).recording_day)

        % relevant metrics from data strucutre
        licks = two_stim_data(animal_idx).recording_day(day_idx).anticipatory_licks;
        labels = two_stim_data(animal_idx).recording_day(day_idx).cs_labels;
        % reaction_times= two_stim_data(animal_idx).recording_day(day_idx).all_stim_diff_from_optimal_reward;
       

        licks_CS_plus = licks(labels);
        licks_CS_minus = licks(~labels);

        observed_diff = mean(licks_CS_plus) - mean(licks_CS_minus);
        n_CS_plus = numel(licks_CS_plus);

        % For RT

        null_diffs = zeros(num_perms, 1);
        for i = 1:num_perms
            perm = licks(randperm(numel(licks)));
            null_diffs(i) = mean(perm(1:n_CS_plus)) - mean(perm(n_CS_plus+1:end));
        end

        p_value = mean(null_diffs >= abs(observed_diff));
        effect_size = observed_diff / std(licks); % Cohen's d approximation

        two_stim_data(animal_idx).recording_day(day_idx).perm_test.observed_diff = observed_diff;
        two_stim_data(animal_idx).recording_day(day_idx).perm_test.p_value = p_value;
        two_stim_data(animal_idx).recording_day(day_idx).perm_test.null_diffs = null_diffs;
        two_stim_data(animal_idx).recording_day(day_idx).perm_test.effect_size= effect_size;
    end
end


% Step 3: Plotting for each workflow separately

% First collect the full set of workflows across all animals
all_workflows = cellfun(@(x) {x.workflow},{two_stim_data.recording_day},'uni',false);
unique_workflows = unique(horzcat(all_workflows{:}),"stable");

for wf_i = 1:numel(unique_workflows)
    wf = unique_workflows{wf_i};
    wf_plot = strrep(wf, '_', ' '); % remove underline

    % Figure for this workflow
    figure('Name',wf,'Color','w','Position',[100 100 1200 800]);

    % Count how many animals actually have this workflow
    has_wf = false(n_animals,1);
    for ai = 1:n_animals
        two_stim_days = two_stim_data(ai).recording_day;
        has_wf(ai) = any(strcmp({two_stim_days.workflow}, wf));
    end

    animals_idx = find(has_wf);
    n_plots = numel(animals_idx);
    if n_plots==0, continue; end

    % Layout: try a square-ish grid
    n_cols = ceil(sqrt(n_plots));
    n_rows =  ceil(n_plots/n_cols);

    for subplot_i = 1:n_plots
        ai = animals_idx(subplot_i);
        two_stim_days = two_stim_data(ai).recording_day;
        seg  = strcmp({two_stim_days.workflow}, wf);
        recs = two_stim_days(seg);

        % compute diffs & p‑values & effect sizes
        n_days = numel(recs);
        diffs   = zeros(1,n_days);
        p_vals  = zeros(1,n_days);
        % effect_sizes= zeros(1,n_days);
        for d = 1:n_days
            diffs(d)  = recs(d).perm_test.observed_diff;
            p_vals(d) = recs(d).perm_test.p_value;
            % effect_sizes(d)= recs(d).perm_test.effect_size;  
        end

      
        sig = (p_vals<0.05) & (diffs>0); 

         % Find first occurrence of two consecutive significant days
        consecutive = sig(1:end-1) & sig(2:end);
        ld = find(consecutive, 1, 'first');

        sig_consecutive= false(size(sig)); % create new array for the consecutive days
        sig_consecutive(ld:end)=sig(ld:end);

        packed_sig_days{ai}.(wf) = sig_consecutive;

        % Subplot
        ax = subplot(n_rows, n_cols, subplot_i);
        plot(1:n_days, diffs, '-ok','LineWidth',1.5,'Parent',ax); hold(ax,'on');
        plot(find(sig_consecutive), diffs(sig_consecutive), 'r*','MarkerSize',8,'Parent',ax);
        % xlim([1 n_days]);
        xlabel(ax,'Day'); ylabel(ax,'CS+−CS−');
        title(ax, two_stim_data(ai).animal_id, 'Interpreter','none');
        grid(ax,'on');
        hold(ax,'off');
    end

    sgtitle(sprintf('Workflow: %s — Anticipatory‐lick diffs', wf_plot), 'FontSize',14);
end

%% Run the wheel analysis to extract signficant days

% Set reaction statistic to use
use_stat = 'median';

% Create master tiled layout
figure;
t = tiledlayout(1,length(animal_list),'TileSpacing','tight');

% create field for stroring learning day
bhv(numel(animal_list)) = struct();  

for curr_animal_idx = 1:length(animal_list)

    animal = animal_list{curr_animal_idx};
    use_workflow = 'stim_wheel*';
    recordings = plab.find_recordings(animal,[],use_workflow);

    if strcmp(animal,'HA014')
        continue;
    end

    % In case no wheel recording for the animal
    if isempty(recordings)
        continue;
    end

    all_recording_dates= {recordings.day};
    all_recording_dates_flat= datetime(vertcat(all_recording_dates{:}), 'InputFormat','yyyy-MM-dd'); % flatten and format

    delta_recording_days = abs(days(diff(all_recording_dates_flat))); % find the day difference
    % this code removes recordings that were added much later (>month)

    threshold = 30;
    gapIdx= find(delta_recording_days > threshold); % the index where the gap exceedes


    if isempty(gapIdx) % if no gap detected; continue as usual
        all_recordings_filtered= recordings;

    else % if gap detected, filter based on it

        removeMask = false(size(all_recording_dates_flat));
        removeMask(1:gapIdx) = true; % create a mask

        all_recordings_filtered= recordings(removeMask); % filter

    end

    % get all current task recordings
    recordings = all_recordings_filtered;
    surround_time = [-5,5];
    surround_sample_rate = 100;
    surround_time_points = surround_time(1):1/surround_sample_rate:surround_time(2);

    n_trials_success = nan(length(recordings),2);
    frac_move_day = nan(length(recordings),1);
    frac_move_stimalign = nan(length(recordings),length(surround_time_points));

    rxn_stat_p = nan(length(recordings),1);
    rxn_stat = nan(length(recordings),1);
    rxn_null_stat = nan(length(recordings),1);

    for curr_recording = 1:length(recordings)

        % Grab pre-load vars
        preload_vars = who;

        % Load data
        rec_day = recordings(curr_recording).day;
        rec_time = recordings(curr_recording).recording{end};
        load_parts = struct;
        load_parts.behavior = true;
        ap.load_recording;

        % Get total trials/water
        n_trials_success(curr_recording,:) = ...
            [length([trial_events.values.Outcome]), ...
            sum([trial_events.values.Outcome])];

        % Align wheel movement to stim onset
        align_times = stimOn_times;
        pull_times = align_times + surround_time_points;

        frac_move_day(curr_recording) = nanmean(wheel_move);

        event_aligned_wheel_vel = interp1(timelite.timestamps, ...
            wheel_velocity,pull_times);
        event_aligned_wheel_move = interp1(timelite.timestamps, ...
            +wheel_move,pull_times,'previous');

        frac_move_stimalign(curr_recording,:) = nanmean(event_aligned_wheel_move,1);

        % Get association stat
        % (skip if only a few trials)
        if n_trials < 10
            continue
        end

        [rxn_stat_p(curr_recording), ...
            rxn_stat(curr_recording),rxn_null_stat(curr_recording)] = ...
            AP_stimwheel_association_pvalue( ...
            stimOn_times,trial_events,stim_to_move,use_stat);

        %%%% WORKING HERE:
        % stim_to_move > stim_to_lastmove - change in stat?
        %%%%

        % Clear vars except pre-load for next loop
        clearvars('-except',preload_vars{:});

    end

    % Define learned day from reaction stat p-value and reaction time
    learned_day = rxn_stat_p < 0.05;

    relative_day = days(datetime({recordings.day}) - datetime({recordings(1).day}))+1;
    nonrecorded_day = setdiff(1:length(recordings),relative_day);

    % Draw in tiled layout nested in master
    t_animal = tiledlayout(t,4,1);
    t_animal.Layout.Tile = curr_animal_idx;
    title(t_animal,animal);

    nexttile(t_animal);
    yyaxis left; plot(relative_day,n_trials_success);
    ylabel('# trials');
    yyaxis right; plot(relative_day,frac_move_day);
    ylabel('Fraction time moving');
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end

    nexttile(t_animal);
    yyaxis left
    plot(relative_day,rxn_stat)
    set(gca,'YScale','log');
    ylabel(sprintf('Rxn stat: %s',use_stat));
    xlabel('Day');
    if any(nonrecorded_day)
        xline(nonrecorded_day,'--k');
    end
    if any(learned_day)
        xline(relative_day(learned_day),'g');
    end


    yyaxis right
    plot(relative_day,(rxn_stat-rxn_null_stat)./(rxn_stat+rxn_null_stat));
    yline(0);
    ylabel(sprintf('Rxn stat idx: %s',use_stat));
    xlabel('Day');

    nexttile(t_animal);
    imagesc(surround_time_points,[],frac_move_stimalign); hold on;
    clim([0,1]);
    colormap(gca,ap.colormap('WK'));
    set(gca,'YTick',1:length(recordings),'YTickLabel', ...
        cellfun(@(day,num) sprintf('%d (%s)',num,day(6:end)), ...
        {recordings.day},num2cell(1:length(recordings)),'uni',false));
    xlabel('Time from stim');
    if any(learned_day)
        plot(0,find(learned_day),'.g')
    end

    nexttile(t_animal); hold on
    set(gca,'ColorOrder',copper(length(recordings)));
    plot(surround_time_points,frac_move_stimalign','linewidth',2);
    xline(0,'color','k');
    ylabel('Fraction moving');
    xlabel('Time from stim');
    if any(learned_day)
        ap.errorfill(surround_time_points,frac_move_stimalign(learned_day,:)', ...
            0.02,[0,1,0],0.1,false);
    end

    drawnow;

    % Store behavior across animals
    bhv(curr_animal_idx).rxn_stat = rxn_stat;
    bhv(curr_animal_idx).rxn_stat = rxn_null_stat;
    bhv(curr_animal_idx).learned_day = learned_day;

end


% Package the structure into a cell array
wheel_move_learned_days= {bhv.learned_day}';

%% Combine the significant days for all tasks

% Combine two stim and wheel tasks and saves the variable

n_animals = numel(packed_sig_days);
combined_sig_day_all_protocols = cell(n_animals, 1);

for animal_idx = 1:n_animals

    % Extract struct of logical arrays
    struct_fields = fieldnames(packed_sig_days{animal_idx});

    % Collect all logical arrays and flatten them
    logicals = cellfun(@(fname) packed_sig_days{animal_idx}.(fname)(:), ...
        struct_fields, 'UniformOutput', false);

    % Add wheel movement learned days (also flattened)
    logicals{end+1} = wheel_move_learned_days{animal_idx}(:);

    % Concatenate into a single column vector
    combined_sig_day_all_protocols{animal_idx} = vertcat(logicals{:});
end

% Construct the dynamic filename based on the current protocol
filename = sprintf('combined_sig_day_all_protocols.mat');

% Save the data to the dynamically generated filename
save(filename, 'combined_sig_day_all_protocols','-v7.3');

%% Exploratory for plotting ITI licks for all animals


% Gather animal IDs
animal_ids = {behaviour_data.animal_id};
n_animals = numel(animal_ids);

figure; hold on;
colors = lines(n_animals);

for a = 1:n_animals
    rec_days = behaviour_data(a).recording_day;
    n_days   = numel(rec_days);

    % Preallocate
    mean_ITI = nan(1, n_days);

    % Compute mean ITI licks per day
    for d = 1:n_days
        counts = rec_days(d).ITI_lick_counts;
        if ~isempty(counts)
            mean_ITI(d) = nanmedian(counts);
        else
            mean_ITI(d) = NaN;
        end
    end

    % Plot
    plot(1:n_days, mean_ITI, '-o', ...
        'Color', colors(a,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', animal_ids{a});
end

xlabel('Recording Day');
ylabel('Mean ITI Licks');
title('Mean Number of Licks During ITI Across Days');
legend('Location','bestoutside');
grid on;
hold off;

%% Explorartory ITI analysis

% Example group definitions
group1_ids = {'HA005','HA008','HA010','HA012'};    % e.g. "learners”
group2_ids = {'DS017','HA007','HA009','HA011'};  % "non-learners"

% Preallocate containers

% Find maximum number of days across *all* animals so we can align days
%all_max_days = max(cellfun(@(x) numel(x.recording_day), {behaviour_data}));

% Containers: group × day
g1_means = nan( numel(group1_ids), 10 );
g2_means = nan( numel(group2_ids), 10 );

% Fill in per-animal, per-day means
for i = 1:numel(behaviour_data)
    id = behaviour_data(i).animal_id;
    days = behaviour_data(i).recording_day;
    n_days = numel(days);
    % extract mean ITI licks per day
    m = nan(1, n_days);
    for d = 1:n_days
        c = days(d).ITI_lick_counts;
        m(d) =nanmean(c);
    end

    % find which group
    idx1 = find(strcmp(id, group1_ids),1);
    if ~isempty(idx1)
        g1_means(idx1,1:n_days) = m;
    end

    idx2 = find(strcmp(id, group2_ids),1);
    if ~isempty(idx2)
        g2_means(idx2,1:n_days) = m;
    end
end

% Now compute group averages and SEM
g1_avg = nanmean(g1_means,1);
g1_sem = nanstd(g1_means,0,1) ./ sqrt(sum(~isnan(g1_means),1));

g2_avg = nanmean(g2_means,1);
g2_sem = nanstd(g2_means,0,1) ./ sqrt(sum(~isnan(g2_means),1));

%last_day = find(~isnan(g1_avg), 1, 'last');
last_day =7;

% Trim all arrays to that day
g1_avg = g1_avg(1:last_day);
g1_sem = g1_sem(1:last_day);
g2_avg = g2_avg(1:last_day);
g2_sem = g2_sem(1:last_day);

% Plot
days = 1:last_day;

figure; hold on;
% Group 1 shaded
fill([days, fliplr(days)], ...
    [g1_avg+g1_sem, fliplr(g1_avg-g1_sem)], ...
    cLearner, 'FaceAlpha',0.2,'EdgeColor','none');
plot(days, g1_avg, 'color',cLearner, 'LineWidth',2);

% Group 2 shaded
fill([days, fliplr(days)], ...
    [g2_avg+g2_sem, fliplr(g2_avg-g2_sem)], ...
    cNonLearner, 'FaceAlpha',0.2,'EdgeColor','none');
plot(days, g2_avg, 'color',cNonLearner, 'LineWidth',2);

xlabel('Recording Day');
ylabel('Mean ITI Licks');
title('Group Comparison of ITI Licking Across Days');
legend({'Learners ± SEM','Learners mean','Non-Learners ± SEM','Non-Learners mean'},'Location','best');
xlim([1,last_day]);
grid on;
hold off;


%% Plots the mean ITI licks split between learners vs non-learners relative to learning day



% get the average ITI licks for all animals across all days
ITI_licks_all_trials= cellfun(@(x) {x.ITI_lick_counts},{behaviour_data.recording_day},'uni',false);
ITI_licks_all_trials_flat = [ITI_licks_all_trials{:}];
mean_ITI_all = cellfun(@nanmean, ITI_licks_all_trials_flat);

% Do it just for right move for now
selected_workflow= workflow_cat==1;

widefield_animal_idx_two_stim_prot= widefield_animal_idx(selected_workflow);
workflow_cat_two_stim_prot= workflow_cat(selected_workflow);
%allITI= allITI(~isnan(allITI));
mean_ITI_all=mean_ITI_all(selected_workflow);
learner_vs_non_learner_grp_idx= is_group_animal(selected_workflow);
learning_index_animal_two_stim_prot= learning_index_animal(selected_workflow);


% 1) find each animal’s learning‐day and compute rel‐day for all sessions
animalIDs = unique(widefield_animal_idx_two_stim_prot);
nSessions = numel(mean_ITI_all);
relDay    = nan(nSessions,1);

for ai = animalIDs(:)'
    % which rows belong to this animal
    idx = find(widefield_animal_idx_two_stim_prot==ai);

    % find within‐animal index of the first learning session
    localLearn = learning_index_animal_two_stim_prot(idx);
    ld_local   = find(localLearn==1,1,'first');
    if isempty(ld_local)
        continue;   % never learned
    end

    % assign relative day: session k → (k – ld_local)
    for k = 1:numel(idx)
        relDay(idx(k)) = k - ld_local;
    end
end

% 2) build common rel-day axis
minD = min(relDay(~isnan(relDay)));
maxD = max(relDay(~isnan(relDay)));
days = minD:maxD;

% 3) preallocate group‐mean/SEM
nDays = numel(days);
mean_L = nan(1,nDays); sem_L = nan(1,nDays);   % learners
mean_N = nan(1,nDays); sem_N = nan(1,nDays);   % non-learners



for i = 1:nDays
    d = days(i);

    % learners
    ix = (learner_vs_non_learner_grp_idx==0) & (relDay==d);
    mean_vals = mean_ITI_all(ix);
    mean_L(i) = mean(mean_vals,'omitnan');
    sem_L(i)  = std(mean_vals,[], 'omitnan')/sqrt(sum(~isnan(mean_vals)));

    % non-learners
    ix = (learner_vs_non_learner_grp_idx==1) & (relDay==d);
    mean_vals = mean_ITI_all(ix);
    mean_N(i) = mean(mean_vals,'omitnan');
    sem_N(i)  = std(mean_vals,[], 'omitnan')/sqrt(sum(~isnan(mean_vals)));
end

% 4) plot with shaded SEM
figure('Color','w','Position',[300 200 600 400]); hold on;

% — Learners —
validL = ~isnan(mean_L) & ~isnan(sem_L);
dL      = days(validL);
mL      = mean_L(validL);
sL      = sem_L(validL);

fill( [dL, fliplr(dL)], ...
    [mL + sL, fliplr(mL - sL)], ...
    cLearner, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(dL, mL, '-','Color',cLearner,'LineWidth',2);

% — Non‐learners —
validN = ~isnan(mean_N) & ~isnan(sem_N);
dN      = days(validN);
mN      = mean_N(validN);
sN      = sem_N(validN);

fill( [dN, fliplr(dN)], ...
    [mN + sN, fliplr(mN - sN)], ...
    cNonLearner, 'FaceAlpha',0.2, 'EdgeColor','none');
plot(dN, mN, '--','Color',cNonLearner,'LineWidth',2);

% mark learning day
xline(0,'k--','LineWidth',1);

xlabel('Days relative to learning');
ylabel('Mean ITI lick count');
title('ITI Licking: Learners vs Non-learners (right move)');
legend({'Learner ±SEM','Learner','Non-learner ±SEM','Non-learner'}, ...
    'Location','best');
xlim([min(days) max(days)]);
grid on;
hold off;
