animal = 'AP022';
verbose = true; % this turns on/off progress display in command line
protocol= 'visual_conditioning*';
recordings = plab.find_recordings(animal,[],protocol);

for day_indx= 1:length({recordings.day})
    rec_day= recordings(day_indx).day;
    rec_time= char(recordings(day_indx).recording);

    if size(rec_time) >1 % check if there is more than 1 recording, if so then take the second trial
        rec_time= rec_time(2,:)
    end


    load_parts.mousecam= true; % to not load the widefield data as well
    ap.load_recording;
    timelite
    bonsai_workflow

    % Bonsai can change how it saves data, and loading scripts can specify
    % how that data is loaded and parsed. In this demo, ap.load_recording
    % loads data from Bonsai as 'trial_events':
    trial_events



    % find the timepoints in which the stimulus was presented based on
    % changes in photodiode
    photodiode_channel= find(contains({timelite.daq_info.channel_name},'photodiode'),1);
    photodiode_values= timelite.data(:,photodiode_channel);

    all_stim_on_time_indecies= find(photodiode_values>3);
    % gives me all the indecies in which the stimulus was on

    all_stim_on_time_timestamps= timelite.timestamps(all_stim_on_time_indecies);
    % get the timestamps of stimulus onset
    stim_time_thresh = 1; % minimum time between stimuli
    diff_stim_on_indecies= diff(all_stim_on_time_timestamps);
    diff_stim_on_indecies= find(diff_stim_on_indecies > stim_time_thresh) +1 ;
    % get the indecies in which the stimulus was presented by using the
    % timepoint difference

    stim_on_timepoints= all_stim_on_time_timestamps(diff_stim_on_indecies);
    stim_on_timepoints= vertcat(all_stim_on_time_timestamps(1),stim_on_timepoints);
    % get all the timepoints in which stimulus was presented and then add
    % the initial one

    stim_on_timepoints_all= nan(size(timelite.timestamps));
    stim_on_indices = find(ismember(timelite.timestamps, stim_on_timepoints));
    stim_on_timepoints_all(stim_on_indices)=4;
    % intialize a variable with all the timepoints and then change only
    % stimulus onset indecies

    % figure;plot(timelite.timestamps,photodiode_values); hold on; plot
    % (timelite.timestamps,stim_on_timepoints_all,'o','MarkerSize',5);
    %

    % load mousecam
    mousecam_fn;
    mousecam_times;
    mousecam_vr = VideoReader(mousecam_fn);

    stim_onset_timepoints= timelite.timestamps(find(~isnan (stim_on_timepoints_all)));
    % get the timepoints of stimulus onset

    % using interp, find the nearst frame that corresponds to stimulus
    % onset time
    extrap_min_indicies_passive_stim= interp1(mousecam_times,1:length(mousecam_times),stim_onset_timepoints,"nearest");

    %frame band to look into then
    frame_band=[-20,100];
    added_frames_passive_conditioning= extrap_min_indicies_passive_stim +frame_band;



    % initalize the dimensions based on first passive trial
    first_passive_trial= squeeze(read(mousecam_vr,added_frames_passive_conditioning(1,:)));
    mouscam_passive_conditioning_dimensions= [size(first_passive_trial),size(extrap_min_indicies_passive_stim)];
    mouscam_im_passive_conditioning= nan(mouscam_passive_conditioning_dimensions,'single');

    %loop through each trial
    for i=1:length(stim_onset_timepoints)
        if added_frames_passive_conditioning(i,2)> mousecam_vr.NumFrames
            continue
        end

        mouscam_im_passive_conditioning(:,:,:,i)= read(mousecam_vr,added_frames_passive_conditioning(i,:));
        disp(i)
    end

    %frame numbers for labeling
    frame_numbers = frame_band(1):frame_band(2);
    %average across 4th dimension
    avg_mouscam_im_passive_condiotnining= mean(mouscam_im_passive_conditioning,4,"omitnan");


    % plot the differences between frames
    avg_mouscam_im_passive_diff= abs(diff(avg_mouscam_im_passive_condiotnining,1,3));

    % find the indecies of the center point
    center_point_indecies= [round(size(avg_mouscam_im_passive_diff,1)*2/3),round(size(avg_mouscam_im_passive_diff,2)*0.56)];

    % index through the center to calculate the right quadrant
    avg_mouscam_im_passive_diff_right_quadrant= avg_mouscam_im_passive_diff(1:center_point_indecies(1),center_point_indecies(2):size(avg_mouscam_im_passive_diff,2),:);


    % plot the mean differences
    mean_diff_per_frame_right_quadrant{day_indx} = smooth(mean(avg_mouscam_im_passive_diff_right_quadrant,[1,2]));

end

% check if the cell dimensions are eqaul, if not padd with NaN
if checkSimilarDimensions(mean_diff_per_frame_right_quadrant)==0
    mean_diff_per_frame_right_quadrant=padDimensions(mean_diff_per_frame_right_quadrant);
end


% create a vector in which the mean differences for each day is stored in
% each row
mean_diff_per_frame_matrix_right_quadrant = cell2mat(mean_diff_per_frame_right_quadrant);
mean_diff_per_frame_matrix_right_quadrant = mean_diff_per_frame_matrix_right_quadrant'; %transpose for plotting



% find the timepoints in which stim and reward were delivered
frame_rate= 30;
reward_onset_after_stimulus= 1; %reward appears 1 second after stim onset
time_seconds= (frame_numbers * (frame_rate))/1000;
stimulus_onset_frame= find(time_seconds==0); %find stimulus onset frame
[~,reward_onset_frame]= min(abs(time_seconds-reward_onset_after_stimulus)); %find the timepoint of reward

% find the mean pixel difference for the min quiscent period
minimum_quiscent_time= 0.5; % the min quiscent time
[~,quiscent_onset_time]= min(abs(time_seconds+minimum_quiscent_time)); %plus as we know the time is negative

% find the corresponding frame index
quiscent_onset_frame_index = find(frame_numbers== frame_numbers(quiscent_onset_time));

% calculate the mean from quiscent onset to stimulus onset -1 for each day
mean_quiscent_baseline= mean(mean_diff_per_frame_matrix_right_quadrant(:,quiscent_onset_frame_index:(stimulus_onset_frame-1)),2);

% substract baseline from the overall mean difference matrix
baseline_subtracted_mean_diff_frame_matrix= mean_diff_per_frame_matrix_right_quadrant-mean_quiscent_baseline;

%normalized mean difference for each trial
normalized_baseline_subtracted_mean_diff_frame_matrix = nan(size(baseline_subtracted_mean_diff_frame_matrix));
soft_normalization_factor = 0.1;

for i=1:size(baseline_subtracted_mean_diff_frame_matrix,1)
    normalized_baseline_subtracted_mean_diff_frame_matrix(i,:)= baseline_subtracted_mean_diff_frame_matrix(i,:)/mean_quiscent_baseline(i);

end

% find the onset time for stimulus and reward
stim_onset_time= time_seconds(stimulus_onset_frame);
reward_onset_time= time_seconds(reward_onset_frame);

%%
% plot the normalized means across days for all mice- this was done
% manually

% for each loop save the normalized means for each mice

%calculate the mean in the timewindow of stimulus presentation
num_days =length({recordings.day});
stim_presentation_mean= mean(normalized_baseline_subtracted_mean_diff_frame_matrix(:,stimulus_onset_frame:reward_onset_frame),2);

% stim_presentation_means_across_all_animals = cell(3,1);
stim_presentation_means_across_all_animals

% Prepare data
num_animals = length(stim_presentation_means_across_all_animals);
num_days = length(stim_presentation_means_across_all_animals{1})-1; % Assuming all animals have the same number of days

animal_names= {'AP016','AP017','AP018'};
% Plotting
figure;
hold on;
for animal_idx = 1:num_animals
    animal_means = stim_presentation_means_across_all_animals{animal_idx};
    animal_name = animal_names{animal_idx};
    plot(1:num_days, animal_means(1:end-1), 'o-', 'DisplayName', animal_name, 'LineWidth', 2, 'MarkerFaceColor', 'none')
   
end
box on;
xlabel('Day');
ylabel('Mean Pixel Change');
title('Licking Change Across Days for All Animals in Conditioning Task 1 Second Prior to Reward');
legend('Location', 'best', 'Interpreter', 'none', 'FontSize', 18);

%%
% heatmap visualization of days
imagesc(time_seconds,1:length({recordings.day}),normalized_baseline_subtracted_mean_diff_frame_matrix);
axis xy;
colormap gray
colorbar;
% add the recording day labels
yticks(1:length({recordings.day}));
yticklabels({recordings.day});
line([stim_onset_time,stim_onset_time], ylim, 'Color', 'b', 'LineStyle', '--')
line([reward_onset_time,reward_onset_time], ylim, 'Color', 'r', 'LineStyle', '--')
xlabel('Time (seconds)');
ylabel('# Day');
title(['Mean Pixel Difference Substraced from Baseline for ' animal ' Across ' num2str(length({recordings.day})) ' Days']);

% Add grid lines for better visualization
grid on;

%%

% plot all figures individually
h=figure;
for i=1:length({recordings.day})
    subplot(3, 4, i);
    plot(time_seconds(1:end-1),normalized_baseline_subtracted_mean_diff_frame_matrix(i,:));
    h1= xline(stim_onset_time,Color='b');
    xline(reward_onset_time,Color='r');
    title(sprintf('Plot for Day %d', i));  % Use sprintf to format the title string
end

linkaxes(h.Children); % makes the y and x axes alligned across all the plots


%%

% Plot training days 2, 5, and 9 together
plot_days = [1,4,6,9];
h = figure;

% define the axes handles
axes_handle= axes(h);
axes_handle.ColorOrder = copper(length(plot_days));

hold on; % Hold the plot for multiple lines
for i = 1:length(plot_days)
    plot(time_seconds(1:end-1), normalized_baseline_subtracted_mean_diff_frame_matrix(plot_days(i), :), ...
        'DisplayName', sprintf('Day %d', plot_days(i)), 'LineWidth', 2); % Increase line width
end

grey_color = [0.5, 0.5, 0.5]; % RGB values for grey
xline(stim_onset_time, 'Color', grey_color, 'LineWidth', 2, 'HandleVisibility', 'off');
xline(reward_onset_time, 'Color', 'r', 'LineWidth', 2, 'HandleVisibility', 'off'); % Use consistent line colors
xlim(axes_handle,[-1,2.5]);
xlabel('Time (seconds)', 'FontSize', 14); % Increase font size
ylabel('Mean Pixel Difference (a.u)', 'FontSize', 14); % Increase font size
title_str = sprintf('Mean Pixel Change Across Learning Days %s', strjoin(string(plot_days), ', '));
title(title_str, 'FontSize', 16);


% Specify the position of the legend
legend('Location', 'best', 'Interpreter', 'none', 'FontSize', 18);
box off;

h=gca;
[h.XAxis.FontSize, h.YAxis.FontSize] = deal(20);
h.Title.FontSize=15;

%%
% plot all figures together
handle_1=figure;
axes_handle= axes(handle_1);
axes_handle.ColorOrder = copper(length({recordings.day}));

hold on;
for i=1:length({recordings.day})
    plot(time_seconds(1:end-1),normalized_baseline_subtracted_mean_diff_frame_matrix(i,:));
end
h1= xline(stim_onset_time,Color='b');
xline(reward_onset_time,Color='r');
title(sprintf('Plot for %d Days for Animal %s', length({recordings.day}),animal));











