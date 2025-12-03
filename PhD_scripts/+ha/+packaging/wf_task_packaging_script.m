
% Save task widefield aligned to the following event:

% right move protocol:

% static protocol:


% wheel protocol:



%% Saves a task data structure that includes V and Kernel activity for the lcr_passive protocol for all the different protocols

% Define protocols
protocols_to_run = {'visual_operant_lick_two_stim_right_move','visual_operant_lick_two_stim_static', 'visual_operant_lick_two_stim_right_move_big_stim',...
    'visual_operant_lick_two_stim_static_big_stim','stim_wheel*'};

% animal_list = {'DS017','HA005','HA006','HA007','HA008',...
%     'HA009','HA010','HA011','HA012','HA013','HA014','HA015'};

animal_list= {'AP014','AP015'};

% Initialize structure
task_data = struct();

% n SVD components
n_components = 200;

% Configuration for the two stim tasks
stim_conditions = {90, -90};
added_time = -2:0.03:2;


for animal_idx = 1:numel(animal_list)

    animal = animal_list{animal_idx};

    % Assign fields: name & recordings
    task_data(animal_idx).animal = animal;
    task_data(animal_idx).widefield = struct();

    global_day= 0;

    for protocol_idx = 1:length(protocols_to_run)

        % define current protocol
        current_protocol = protocols_to_run{protocol_idx};

        % get all recordings and re-arrange
        all_recordings = plab.find_recordings(animal, [],current_protocol);

         % continue if no recording is found
        if isempty(all_recordings)
            fprintf('No task recording for protocol %s',current_protocol);
            continue
        end

        all_recording_dates= {all_recordings.day};
        all_recording_dates_flat= datetime(vertcat(all_recording_dates{:}), 'InputFormat','yyyy-MM-dd'); % flatten and format

        delta_recording_days = abs(days(diff(all_recording_dates_flat))); % find the day difference
        % this code removes recordings that were added much later (>month)

        threshold = 10;
        gapIdx= find(delta_recording_days > threshold); % the index where the gap exceedes

        if isempty(gapIdx) % if no gap detected; continue as usual
            all_recordings_filtered= all_recordings;

        else % if gap detected, filter based on it

            removeMask = false(size(all_recording_dates_flat));
            removeMask(1:gapIdx) = true; % create a mask

            all_recordings_filtered= all_recordings(removeMask); % filter

        end
        
        % get all current task recordings
        all_recordings_curr_protocol = all_recordings_filtered;
       
        
        for day_idx = 1: numel(all_recordings_curr_protocol)

            rec = all_recordings_curr_protocol(day_idx);
            rec_day = rec.day;
            rec_time = rec.recording{end};

            fprintf('Animal: %s \nCurrent Potocol :%s \nDay: %d/%d\n', animal, current_protocol, day_idx, numel(all_recordings_curr_protocol));

            % Update global day
            global_day= global_day+1;

            % keep the variables up to this point
            vars_to_keep= who;

            % Load components
            load_parts.widefield_master = true;
            load_parts.widefield_align= true;
            load_parts.widefield = true;
            ap.load_recording;


            if ~load_parts.widefield_align % flag if there is no alignment 
                error('No alignment for animal %s in %s recording day',animal,rec_day)
            end

            timelite;
            bonsai_workflow;
            ha.load_bonsai;

            % If protocol includes two stim extract the relevant session
            % parameters - the key difference is that in right_move the
            % final position is the centre ; whereas for static it is
            % equivilent to reward time

            if contains(current_protocol,{'two_stim'})

                % Extract the rewarded X position
                rewardedX = trial_events.parameters.RewardedX;

                % Extract possible X positions
                possible_X_positions = trial_events.parameters.TrialXSpace;

                % Define non-rewarded X positions (excluding the rewardedX)
                non_rewarded_X = setdiff(possible_X_positions, rewardedX);

                % find rewarded and non rewarded times
                rewarded_stim_on_times= stimOn_times(trials_x_pos==rewardedX);

                rewarded_stim_in_final_position_times= stimOff_times(trials_x_pos==rewardedX);
                non_rewarded_stim_on_times= stimOn_times(trials_x_pos==non_rewarded_X);
                non_rewarded_stim_off_times= stimOff_times(trials_x_pos==non_rewarded_X);

                % Calculate the putative left stim in center times (stim off + time to move to center)
                non_rewarded_stim_putative_center_times= non_rewarded_stim_off_times + trial_events.parameters.StimMoveTime;

                % If it's static protocol then there is no stim move times;
                % therefore keep those empty             

                if contains(current_protocol,'static')

                    event_aligned_labels= {'rewarded_stim_on','rewarded_stim_final_position','non_rewarded_stim_on','reward_times'};
                    % Interpolate V across relevant event times

                    rewarded_stim_onset_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_on_times + added_time), [3,2,1]);
                    rewarded_stim_in_final_position_aligned_V= permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_in_final_position_times + added_time), [3,2,1]);
                    non_rewarded_stim_onset_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', non_rewarded_stim_on_times + added_time), [3,2,1]);
                    reward_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', reward_times_task + added_time), [3,2,1]);

                    % Append V data into the structure
                    task_data(animal_idx).widefield(global_day).rewarded_stim_on_aligned_V= rewarded_stim_onset_aligned_V;
                    task_data(animal_idx).widefield(global_day).rewarded_stim_in_final_position_aligned_V= rewarded_stim_in_final_position_aligned_V;
                    task_data(animal_idx).widefield(global_day).non_rewarded_stim_onset_aligned_V= non_rewarded_stim_onset_aligned_V;
                    task_data(animal_idx).widefield(global_day).rewarded_stim_start_to_move_aligned_V= []; % keep it empty for field consistancy 
                    task_data(animal_idx).widefield(global_day).reward_times_aligned_V= reward_aligned_V;


                    %%% regression analysis %%%%

                    % define bin edges
                    bin_edges = [wf_t', wf_t(end)];

                    % Bin the different regressors

                    rewarded_stim_on_regressor = histcounts(rewarded_stim_on_times, bin_edges);
                    rewarded_stim_final_position_regressor= histcounts(rewarded_stim_in_final_position_times, bin_edges);
                    non_rewarded_stim_regressor = histcounts(non_rewarded_stim_on_times, bin_edges);
                    reward_times_regressor= histcounts(reward_times_task, bin_edges);

                    % join regressors and stack them to dim x time
                    regressors= vertcat(rewarded_stim_on_regressor,rewarded_stim_final_position_regressor,non_rewarded_stim_regressor, ...
                        reward_times_regressor);


                    % define time shifts converted into frames
                    % time_shifts = {-1.5 *wf_framerate:1*wf_framerate};
                    % wf_framerate; % approx 34.5

                    time_shifts= {-52:10};


                    % predict the stimului presentation and reward times based
                    % on wf_V (componented reduced to 100 to speed up)
                    [kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
                        ap.regresskernel(wf_V(1:n_components/2,:),regressors,time_shifts,100);

                    % append Kernels to the data structure
                    for event_idx=1:size(kernels,3)
                        curr_event_label= event_aligned_labels{event_idx};
                        task_data(animal_idx).widefield(global_day).([curr_event_label '_aligned_kernel'])= kernels(:,:,event_idx);
                    end
                      task_data(animal_idx).widefield(global_day).rewarded_stim_start_to_move_aligned_kernel= [];

                else
                    event_aligned_labels= {'rewarded_stim_on','rewarded_stim_final_position','non_rewarded_stim_on','rewarded_stim_start_to_move','reward_times'};

                rewarded_stim_start_move_times; % from bonsai script

                % Interpolate V across relevant event times
                rewarded_stim_onset_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_on_times + added_time), [3,2,1]);
                rewarded_stim_in_final_position_aligned_V= permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_in_final_position_times + added_time), [3,2,1]);
                non_rewarded_stim_onset_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', non_rewarded_stim_on_times + added_time), [3,2,1]);
                rewarded_stim_start_to_move_aligned_V= permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_start_move_times + added_time), [3,2,1]);
                reward_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', reward_times_task + added_time), [3,2,1]);

                % Append V data into the structure
                task_data(animal_idx).widefield(global_day).rewarded_stim_on_aligned_V= rewarded_stim_onset_aligned_V;
                task_data(animal_idx).widefield(global_day).rewarded_stim_in_final_position_aligned_V= rewarded_stim_in_final_position_aligned_V;
                task_data(animal_idx).widefield(global_day).non_rewarded_stim_onset_aligned_V= non_rewarded_stim_onset_aligned_V;
                task_data(animal_idx).widefield(global_day).rewarded_stim_start_to_move_aligned_V= rewarded_stim_start_to_move_aligned_V;
                task_data(animal_idx).widefield(global_day).reward_times_aligned_V= reward_aligned_V;

                %%% regression analysis %%%%

                % define bin edges
                bin_edges = [wf_t', wf_t(end)];

                % Bin the different regressors
               
                rewarded_stim_on_regressor = histcounts(rewarded_stim_on_times, bin_edges);
                rewarded_stim_final_position_regressor= histcounts(rewarded_stim_in_final_position_times, bin_edges);
                non_rewarded_stim_regressor = histcounts(non_rewarded_stim_on_times, bin_edges);
                rewarded_stim_start_to_move_regressor= histcounts(rewarded_stim_start_move_times, bin_edges);
                reward_times_regressor= histcounts(reward_times_task, bin_edges);

                % join regressors and stack them to dim x time
                regressors= vertcat(rewarded_stim_on_regressor,rewarded_stim_final_position_regressor,non_rewarded_stim_regressor,rewarded_stim_start_to_move_regressor, ...
                    reward_times_regressor);


                % define time shifts converted into frames
                % time_shifts = {-1.5 *wf_framerate:1*wf_framerate};
                % wf_framerate; % approx 34.5
                
                time_shifts= {-52:10};
              

                % predict the stimului presentation and reward times based
                % on wf_V (componented reduced to 100 to speed up)
                [kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
                    ap.regresskernel(wf_V(1:n_components/2,:),regressors,time_shifts,100);

                % append Kernels to the data structure
                for event_idx=1:size(kernels,3)
                    curr_event_label= event_aligned_labels{event_idx};
                    task_data(animal_idx).widefield(global_day).([curr_event_label '_aligned_kernel'])= kernels(:,:,event_idx);
                end

                end

                % For wheel currently only extrapolates V aligned to stim
                % onset

            elseif contains(current_protocol,{'stim_wheel'})
                rewarded_stim_on_times= stimOn_times;
                rewarded_stim_onset_aligned_V = permute(interp1(wf_t, wf_V(1:n_components,:)', rewarded_stim_on_times + added_time), [3,2,1]);

                % Append V data into the structure
                task_data(animal_idx).widefield(global_day).rewarded_stim_on_aligned_V= rewarded_stim_onset_aligned_V;

                %%% regression analysis %%%%

                % define bin edges
                bin_edges = [wf_t', wf_t(end)];

                % Bin the different regressors
                rewarded_stim_on_regressor = histcounts(rewarded_stim_on_times, bin_edges);

                regressors= vertcat(rewarded_stim_on_regressor);

                % define time shifts
                time_shifts= {-52:10};

                % predict the stimului presentation and reward times based on wf_V
                [kernels,predicted_signals,explained_var,predicted_signals_reduced] = ...
                    ap.regresskernel(wf_V(1:n_components/2,:),regressors,time_shifts,100);

                % Append kernel to structure
                task_data(animal_idx).widefield(global_day).rewarded_stim_on_aligned_kernel= kernels;


            end


            % clear the day specific variables
            clearvars('-except',vars_to_keep{:});
        end
    end
end


% % Construct the dynamic filename based on the current protocol
filename = sprintf('AP_014_015_task_widefield.mat');

% Save the data to the dynamically generated filename
save(filename, 'task_data','-v7.3');


%% Just to examine the data

% load master_U
wf_U= plab.wf.load_master_U;

kernel_added_time= fliplr((-10:52));
kernel_added_time= kernel_added_time/34;

single_animal_data= cat(2,task_data(1).widefield);

right_stim_data= cat(3,single_animal_data.rewarded_stim_on_aligned_kernel);

mean_right_stim_data= mean(right_stim_data,3);

ap.imscroll(plab.wf.svd2px(wf_U(:,:,1:100),mean_right_stim_data),kernel_added_time);
axis equal;




