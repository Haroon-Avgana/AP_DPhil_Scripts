% Load Bonsai data

if verbose; fprintf('Loading Bonsai...'); end

%% Load general Bonsai events

% Get Bonsai workflow
bonsai_dir = dir(plab.locations.filename('server',animal,rec_day,rec_time,'bonsai'));
bonsai_workflow = bonsai_dir([bonsai_dir.isdir] & ~contains({bonsai_dir.name},'.')).name;

% Load Bonsai events (should be included in every workflow)
bonsai_events_fn = plab.locations.filename('server', ...
    animal,rec_day,rec_time,'bonsai','bonsai_events.csv');

if exist(bonsai_events_fn,'file') && ~isempty(readtable(bonsai_events_fn))

    % Check if Bonsai file is bad
    % (when this happens, last line is not written in full, so look for a
    % timestamp entry that is the wrong length)
    bonsai_readtimestamps_opts = detectImportOptions(bonsai_events_fn,'delimiter',',');
    bonsai_readtimestamps_opts.SelectedVariableNames = 'Timestamp';
    bonsai_timestamps_raw = readtable(bonsai_events_fn,bonsai_readtimestamps_opts);
    if length(unique(cellfun(@length,bonsai_timestamps_raw.Timestamp))) ~= 1
        bad_bonsai_csv = true;
        warning('Bonsai file ends improperly: %s',bonsai_events_fn);
    else
        bad_bonsai_csv = false;
    end

    % Set Bonsai timestamp format
    bonsai_table_opts = detectImportOptions(bonsai_events_fn,'delimiter',',');
    bonsai_table_opts = setvaropts(bonsai_table_opts, ...
        'Timestamp','Type','datetime', ...
        'InputFormat','yyyy-MM-dd''T''HH:mm:ss.SSSSSSSZ','TimeZone','local', ...
        'DatetimeFormat','yyyy-MM-dd HH:mm:ss.SSS');

    % Load Bonsai CSV file
    bonsai_events_raw = readtable(bonsai_events_fn,bonsai_table_opts);

    % 'Values' column should be doubles: if this was mixed (e.g. numbers
    % and "true/false" for booleans), convert strings to numbers
    if iscell(bonsai_events_raw.Value)
        bonsai_events_raw.Value = ...
            cellfun(@(x) double(str2num(lower(x))),bonsai_events_raw.Value);
        disp(rec_day);
    end

    % Create nested structure for trial events
    trial_events = struct('parameters',cell(1),'values',cell(1),'timestamps',cell(1));

    % Save anything in "Trial 0" as a parameter
    parameter_idx = bonsai_events_raw.Trial == 0;
    unique_parameters = unique(bonsai_events_raw.Event(parameter_idx));
    for curr_parameter = unique_parameters'
        curr_parameter_idx = parameter_idx & strcmp(bonsai_events_raw.Event,curr_parameter);
        trial_events.parameters.(cell2mat(curr_parameter)) = bonsai_events_raw.Value(curr_parameter_idx);
    end

    % Loop through trials (excluding 0), save values and timestamps for all events
    % (exclude entries with empty events - happens sometimes on last entry)
    empty_events = cellfun(@isempty,bonsai_events_raw.Event);
    n_trials = max(bonsai_events_raw.Trial);
    unique_events = unique(bonsai_events_raw.Event(~parameter_idx & ~empty_events));
    for curr_trial = 1:n_trials
        for curr_event = unique_events'
            curr_event_idx = bonsai_events_raw.Trial == curr_trial & strcmp(bonsai_events_raw.Event,curr_event);
            trial_events.values(curr_trial).(cell2mat(curr_event)) = bonsai_events_raw.Value(curr_event_idx);
            trial_events.timestamps(curr_trial).(cell2mat(curr_event)) = bonsai_events_raw.Timestamp(curr_event_idx);
        end
    end
end


%% --- Workflow-specific loading ---

%% Stim wheel (task)
if contains(bonsai_workflow,'stim_wheel')
    % Task: stim and response times

    % Photodiode bug (old, now fixed): screen could flick to black briefly
    % when clicking on another window. This are always brief, and no way to
    % tell when it happened, so compensate by removing all flips that
    % happen with short duration
    photodiode_flicker_timethresh = 0.2; % off for this duration = flicker
    photodiode_on_times = photodiode_times(photodiode_values == 1);
    photodiode_off_times = photodiode_times(photodiode_values == 0);
    photodiode_offtime = photodiode_on_times(2:end) - photodiode_off_times(1:length(photodiode_on_times)-1);
    photodiode_flicker = find(photodiode_offtime < photodiode_flicker_timethresh);

    % Stim times: when photodiode flips to 1/0
    stimOn_times = photodiode_on_times(setdiff(1:length(photodiode_on_times),photodiode_flicker+1));
    stimOff_times = photodiode_off_times(setdiff(1:length(photodiode_off_times),photodiode_flicker));

    % Use only trials with outcome and stim off
    n_trials = min(length(stimOff_times),length([trial_events.timestamps.Outcome]));

    % If fewer stim times than trials, assume something was broken at the
    % end and cut off
    if length(stimOn_times) < n_trials
        % If fewer stims than trials, assume broken at end and cut off
        n_trials = length(stimOn_times);
        warning('%s %s %s: fewer stim onsets than bonsai trials, truncating stims', ...
            animal,rec_day,rec_time);
    end

    % Find the last move stop before stim on
    % (sometimes this isn't after the stimulus: Bonsai's quiescence clock
    % isn't very accurate)
    last_prestim_move_stop = cellfun(@(x) ...
        find(~wheel_move(timelite.timestamps <= x),1,'last') + 1, ...
        num2cell(stimOn_times(1:n_trials)));

    % Find the first move start after the pre-stim move stop
    % (usually this is the first post-stim move, sometimes this is before
    % stim if quiescence clock didn't work)
    stim_move_time = cellfun(@(x) ...
        timelite.timestamps((x-1) + find(wheel_move(x:end),1,'first')), ...
        num2cell(last_prestim_move_stop));

    % Find outcome movement start (last start before stim off)
    stim_lastmove_time = timelite.timestamps(cellfun(@(x) ...
        find(~wheel_move(timelite.timestamps <= x),1,'last') + 1, ...
        num2cell(stimOff_times(1:n_trials))));

    % Store trial values
    stim_to_move = stim_move_time - stimOn_times(1:n_trials);
    stim_to_lastmove = stim_lastmove_time - stimOn_times(1:n_trials);
    trial_outcome = vertcat(trial_events.values(1:n_trials).Outcome);

    % Get task/manual reward times
    if isfield(trial_events.timestamps,'ManualReward')
        % After 'ManualReward' event was added
        % (check that total rewards = task + manual)
        bonsai_task_reward = vertcat(trial_events.timestamps( ...
            vertcat(trial_events.values.Outcome) == 1).Outcome);
        bonsai_manual_reward = vertcat(trial_events.timestamps.ManualReward);

        if length(reward_times) ~= ...
                (length(bonsai_task_reward) + length(bonsai_manual_reward))
            warning('Total rewards ~= task + manual');
        else
            [~,bonsai_reward_sortidx] = sort(vertcat(bonsai_task_reward,bonsai_manual_reward));
            bonsai_reward_grp = vertcat(1*ones(size(bonsai_task_reward)),2*ones(size(bonsai_manual_reward)));
            reward_times_task = reward_times(bonsai_reward_grp(bonsai_reward_sortidx) == 1);
            reward_times_manual = reward_times(bonsai_reward_grp(bonsai_reward_sortidx) == 2);
        end
    else
        % Without Bonsai event: closest reward to stim off on rewarded trials
        if n_trials > 1
            reward_times_task = interp1(reward_times,reward_times, ...
                stimOff_times(trial_outcome(1:min(n_trials,length(stimOff_times)))==1),'nearest','extrap');
        end
    end

    %     % Get all outcome times (interleaved reward/punish)
    %     % (NOT DONE YET)
    %     bonsai_reward_datetime = vertcat(trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome);
    %     bonsai_reward_reltime = seconds(bonsai_reward_datetime - bonsai_reward_datetime(1));
    %
    %     bonsai_timelite_reward_idx = interp1(bonsai_reward_reltime, ...
    %         1:length(bonsai_reward_reltime),reward_times-reward_times(1),'nearest','extrap');
    %
    %     timelite2bonsai = reward_times;
    %     bonsai2timelite_datetime = vertcat(trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome);
    %
    %     bonsai_relative_t = bonsai2timelite_datetime(1);
    %     bonsai2timelite = seconds(bonsai2timelite_datetime - bonsai_relative_t);
    %
    %     event_datetime = vertcat(trial_events.timestamps.Outcome);
    %
    %     event_t_bonsai = seconds(event_datetime - bonsai_relative_t);
    %
    %     event_t_tl = interp1(bonsai2timelite,timelite2bonsai,event_t_bonsai,'linear','extrap');

%% LCR passive, visual conditioning
elseif contains(bonsai_workflow,{'lcr_passive'})
    % Passive protocol stim on times

    % Photodiode bug (old, now fixed): screen could flick to black briefly
    % when clicking on another window. This are always brief, and no way to
    % tell when it happened, so compensate by removing all flips that
    % happen with short duration
    photodiode_flicker_timethresh = 0.2; % off for this duration = flicker
    photodiode_on_times = photodiode_times(photodiode_values == 1);
    photodiode_off_times = photodiode_times(photodiode_values == 0);
    photodiode_offtime = photodiode_on_times(2:end) - photodiode_off_times(1:length(photodiode_on_times)-1);
    photodiode_flicker = find(photodiode_offtime < photodiode_flicker_timethresh);

    % Stim times: when photodiode flips to 1/0
    stimOn_times = photodiode_on_times(setdiff(1:length(photodiode_on_times),photodiode_flicker+1));
    stimOff_times = photodiode_off_times(setdiff(1:length(photodiode_off_times),photodiode_flicker));
    
    % If bad Bonsai CSV: truncate stim times to what was recorded
    if bad_bonsai_csv
        stimOn_times = stimOn_times(1:sum(vertcat(trial_events.values.StimOn) == 1));
        stimOff_times = stimOff_times(1:sum(vertcat(trial_events.values.StimOn) == 0));
    end

    %% visual static connditioning
elseif contains(bonsai_workflow,{'static'})
    % Passive protocol stim on times

    % Photodiode bug (old, now fixed): screen could flick to black briefly
    % when clicking on another window. This are always brief, and no way to
    % tell when it happened, so compensate by removing all flips that
    % happen with short duration
    photodiode_flicker_timethresh = 0.2; % off for this duration = flicker
    photodiode_on_times = photodiode_times(photodiode_values == 1);
    photodiode_off_times = photodiode_times(photodiode_values == 0);
    photodiode_offtime = photodiode_on_times(2:end) - photodiode_off_times(1:length(photodiode_on_times)-1);
    photodiode_flicker = find(photodiode_offtime < photodiode_flicker_timethresh);

    % Stim times: when photodiode flips to 1/0
    stimOn_times = photodiode_on_times(setdiff(1:length(photodiode_on_times),photodiode_flicker+1));
    stimOff_times = photodiode_off_times(setdiff(1:length(photodiode_off_times),photodiode_flicker));

    % If bad Bonsai CSV: truncate stim times to what was recorded
    if bad_bonsai_csv
        stimOn_times = stimOn_times(1:sum(vertcat(trial_events.values.StimOn) == 1));
        stimOff_times = stimOff_times(1:sum(vertcat(trial_events.values.StimOn) == 0));
    end

    trial_outcome = vertcat(trial_events.values(1:n_trials).Outcome);

    trials_x_pos = [trial_events.values.TrialX]'; % X positions

    % Determine the minimum number of valid trials from all involved arrays
    n_valid = min([n_trials, length(trials_x_pos), length(stimOff_times)]);

    % Truncate arrays to a common length
    stimOff_times = stimOff_times(1:n_valid);
    stimOn_times= stimOn_times(1:n_valid);
    trials_x_pos  = trials_x_pos(1:n_valid);

    % Get task/manual reward times
    if isfield(trial_events.timestamps,'ManualReward')
        % After 'ManualReward' event was added
        % (check that total rewards = task + manual)
        bonsai_task_reward = vertcat(trial_events.timestamps.Outcome);
        bonsai_manual_reward = vertcat(trial_events.timestamps.ManualReward);

        if length(reward_times) ~= ...
                (length(bonsai_task_reward) + length(bonsai_manual_reward))
            warning('Total rewards ~= task + manual');
            
           % reward_times(find(diff(reward_times)<0.5)+1)= []; % remove the second reward task signal
            rewarded_stimOff_times = stimOff_times(trials_x_pos == 90);

            valid_outcomes = trial_outcome(1:min(length(rewarded_stimOff_times), length(trial_outcome))) == 1;

            reward_times_task = interp1( ...
                reward_times, reward_times, ...
                rewarded_stimOff_times(valid_outcomes), ...
                'nearest', 'extrap'  ...
                );
        else
            [~,bonsai_reward_sortidx] = sort(vertcat(bonsai_task_reward,bonsai_manual_reward));
            bonsai_reward_grp = vertcat(1*ones(size(bonsai_task_reward)),2*ones(size(bonsai_manual_reward)));
            reward_times_task = reward_times(bonsai_reward_grp(bonsai_reward_sortidx) == 1);
            reward_times_manual = reward_times(bonsai_reward_grp(bonsai_reward_sortidx) == 2);
        end
    else


        % Without Bonsai event: closest reward to stim off on rewarded trials

        if n_trials > 1

            rewarded_stimOff_times= stimOff_times(trials_x_pos==90);
            % ensure rewarded stim off is the same legnth as reward times
            rewarded_stimOff_times= rewarded_stimOff_times(1:length(reward_times));

            reward_times_task = interp1(reward_times,reward_times, ...
                rewarded_stimOff_times(trial_outcome(1:min(n_trials,length(rewarded_stimOff_times)))==1),'nearest','extrap');
        end
    end


    %% Right move
elseif contains(bonsai_workflow,'right_move')
    % right move protocol stim on/off times

    % Photodiode on and off times
    photodiode_on_times = photodiode_times(photodiode_values == 1);
    photodiode_off_times = photodiode_times(photodiode_values == 0);

    % Stimulus movement time (includes +0.3s for noise)
    time_for_stim_to_move = trial_events.parameters.StimMoveTime + 0.3;

    % Find indices where difference between consecutive off times < move time
    rewarded_stim_first_off_indecies = find(diff(photodiode_off_times) < time_for_stim_to_move);

    % First off indices correspond to first on indices
    rewarded_stim_first_on_indecies = rewarded_stim_first_off_indecies;

    % Second on = first + 1
    rewarded_stim_second_on_indecies = rewarded_stim_first_on_indecies + 1;
    rewarded_stim_second_off_indecies = rewarded_stim_second_on_indecies;

    % Full index list
    stim_on_array_indecies = 1:length(photodiode_on_times);

    % Left stim on = all indices not in right stim (rewarded) on indices
    left_stim_on_indecies = setdiff(stim_on_array_indecies, ...
        [rewarded_stim_first_on_indecies, rewarded_stim_second_on_indecies])';

    % Make sure left stim off indices do not exceed photodiode_off_times
    left_stim_off_indecies = left_stim_on_indecies(left_stim_on_indecies <= length(photodiode_off_times));

    % Combine rewarded and left stim indices
    all_stim_on_indecies = [rewarded_stim_first_on_indecies', left_stim_on_indecies']';
    all_stim_off_indecies = [rewarded_stim_second_off_indecies', left_stim_off_indecies']';
    all_stim_start_move_indecies = rewarded_stim_first_off_indecies;

    % Return corresponding times
    stimOn_times = sort(photodiode_on_times(all_stim_on_indecies));
    stimOff_times = sort(photodiode_off_times(all_stim_off_indecies));
    rewarded_stim_start_move_times = sort(photodiode_off_times(all_stim_start_move_indecies));


    trial_outcome = vertcat(trial_events.values(1:n_trials).Outcome);

    trials_x_pos = [trial_events.values.TrialX]'; % X positions

    % Determine the minimum number of valid trials from all involved arrays
    n_valid = min([n_trials, length(trials_x_pos), length(stimOff_times)]);

    % Truncate arrays to a common length
    stimOff_times = stimOff_times(1:n_valid);
    stimOn_times= stimOn_times(1:n_valid);
    trials_x_pos  = trials_x_pos(1:n_valid);

    % Handle reward time separation
    if isfield(trial_events.timestamps, 'ManualReward')

        % Separate task and manual rewards using original order
        task_reward_times = vertcat(trial_events.timestamps.Outcome);
        manual_reward_times = vertcat(trial_events.timestamps.ManualReward);

        % Sanity check - in case of mismatch just interpolate based on the
        % outcome 

        if length(reward_times) ~= length(task_reward_times) + length(manual_reward_times)
            warning('Total rewards ~= task + manual');

            % reward_times(find(diff(reward_times)<0.5)+1)= []; % remove the second reward task signal
            rewarded_stimOff_times = stimOff_times(trials_x_pos == 90);

            valid_outcomes = trial_outcome(1:min(length(rewarded_stimOff_times), length(trial_outcome))) == 1;

            reward_times_task = interp1( ...
                reward_times, reward_times, ...
                rewarded_stimOff_times(valid_outcomes), ...
                'nearest', 'extrap' ...
                );

        else

            % Build grouping vector and sort
            reward_grp = [ones(size(task_reward_times)); 2 * ones(size(manual_reward_times))];
            [~, sort_idx] = sort([task_reward_times; manual_reward_times]);
            reward_grp_sorted = reward_grp(sort_idx);

            % Separate final reward time arrays
            reward_times_task = reward_times(reward_grp_sorted == 1);
            reward_times_manual = reward_times(reward_grp_sorted == 2);
        end
    else
        % Infer task reward times via proximity to stimOff
        if n_valid > 1
            rewarded_stimOff_times = stimOff_times(trials_x_pos == 90);
            
            valid_outcomes = trial_outcome(1:min(length(rewarded_stimOff_times), length(trial_outcome))) == 1;

            reward_times_task = interp1( ...
                reward_times, reward_times, ...
                rewarded_stimOff_times(valid_outcomes), ...
                'nearest', 'extrap' ...
                );
        end
    end

    %% Sparse noise
elseif strcmp(bonsai_workflow,'sparse_noise')

    % Sparse noise: get noise locations and times

    bonsai_noise_fn = plab.locations.filename('server', ...
        animal,rec_day,rec_time,'bonsai','NoiseLocations.bin');
    fid = fopen(bonsai_noise_fn);

    n_x_squares = trial_events.parameters.ScreenExtentX./trial_events.parameters.StimSize;
    n_y_squares = trial_events.parameters.ScreenExtentY./trial_events.parameters.StimSize;

    noise_locations = reshape(fread(fid),n_y_squares,n_x_squares,[]);
    fclose(fid);

    % Get stim times from photodiode (extrapolate: sparse noise photodiode
    % flips every N stim to give a more robust signal)
    if ~isfield(trial_events.parameters,'NthPhotodiodeFlip')
        % (if it wasn't defined, default to flipping on every stim)
        trial_events.parameters.NthPhotodiodeFlip = 1;
    end
    photodiode_stim_idx = 1:trial_events.parameters.NthPhotodiodeFlip:size(noise_locations,3);
    % (check that the number of photodiode flips is expected)
    if length(photodiode_stim_idx) ~= length(photodiode_times)
        if length(photodiode_stim_idx) < length(photodiode_times)
            % (rarely: Bonsai square to black temporarily, don't know why)
            % (fix? find when time differences on either side are less than a
            % threshold and remove those flips)
            photodiode_diff = diff(photodiode_times);
            photodiode_diff_thresh = mean(photodiode_diff)/2;
            bad_photodiode_idx = ...
                find(photodiode_diff(1:end-1) < photodiode_diff_thresh & ...
                photodiode_diff(2:end) < photodiode_diff_thresh) + 1;

            if (length(photodiode_times) - length(photodiode_stim_idx)) == ...
                    length(bad_photodiode_idx)
                % (if detected bad flips even out the numbers, apply fix)
                photodiode_times(bad_photodiode_idx) = [];
            else
                % (otherwise, error)
                error('Sparse noise: photodiode > stim, unfixable')
            end
        else
            error('Sparse noise: photodiode < stim')
        end
    end

    stim_times = interp1(photodiode_stim_idx,photodiode_times, ...
        1:size(noise_locations,3),'linear','extrap')';

    %% Operant
elseif contains(bonsai_workflow,'operant')

    %% End

end

if verbose; fprintf('%s\n',bonsai_workflow); end


%% For testing: convert bonsai times to timelite times

% %%%% testing: bonsai times into timelite times
% % % (reward)
% % timelite2bonsai = reward_times;
% % bonsai2timelite_datetime = vertcat(trial_events.timestamps([trial_events.values.Outcome] == 1).Outcome);
% % (stim)
% timelite2bonsai = stimOn_times;
% bonsai2timelite_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
%
% bonsai_relative_t = bonsai2timelite_datetime(1);
% bonsai2timelite = seconds(bonsai2timelite_datetime - bonsai_relative_t);
%
% % event_datetime = cellfun(@(x) x(1),{trial_events.timestamps.StimOn}');
% event_datetime = vertcat(trial_events.timestamps.QuiescenceStart);
% % event_datetime = vertcat(trial_events.timestamps.QuiescenceReset);
%
% event_t_bonsai = seconds(event_datetime - bonsai_relative_t);
%
% event_t_tl = interp1(bonsai2timelite,timelite2bonsai,event_t_bonsai,'linear','extrap');




