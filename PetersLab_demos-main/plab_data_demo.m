%% Data demo
%
% A demo on how to use the server and load data
% Whenever you see a line of code, run it and continue reading
% [EXERCISE] marks a section with instructions to complete

%% Github repositories to download

% Clone these repos and add to the Matlab path (Home > Set Path > Add with
% Subfolders...)
% https://github.com/PetersNeuroLab/PetersLab_analysis
% https://github.com/petersaj/AP_scripts_peterslab

%% File structure and finding recordings

% To connect to the server: 
% - Open up the file explorer, select This PC in the left hand menu 
% - Right click on this PC and select map network drive 
% - Enter the following address: \\qnap-ap001.dpag.ox.ac.uk\APlab\ 
% - It should prompt you for credentials, please enter your MSD details in the following format: 
% - Username – MSD\(MSD username) 
% - Password – normal MSD password 

% The server contains these top folders:
% \Data: for all raw and preprocessed data (all users can write)
% \Users: all users have a folder for their own use (given users can write)
% \Lab: files for the general lab and rigs (only Andy can write)

% Data is organized on the server in these folders: 
% \Data
% |- animal name (initials & number, e.g. AP012)
% | |- recording date (yyyy-mm-dd, e.g. 2023-03-16)
% |  |- recording time (Recording_HHMM, e.g. Recording_1021 started at 10:21am)
% |  | |- recording data by modality (e.g. mousecam, widefield)
% |  |- day data spanning multiple recordings (e.g. ephys, widefield)
% |- non-recording data (e.g. histology)


%% Finding data from recordings

% The function 'plab.find_recordings' is used to find data from recordings
% The syntax is: plab.find_recordings(animal,recording_day,workflow)
% (type 'help plab.find_recordings' for documentation).
% 'animal' must be filled in, other specifics can be left out to return all
% relevant recordings.

% --- Working with plab.find_recordings

% This line finds all recordings done with animal AP010:
% (day and workflow are not entered)
recordings = plab.find_recordings('AP010');

% 'recordings' is a structure with one entry for each day containing:
% - day: day of recording
% - index: index of recording within the day for that animal (e.g. 3rd recording = [3])
% - recording: time of recording (HHMM)
% - workflow: Bonsai workflow (i.e. stimuli or task protocol)
% - mousecam: if there was mouse camera recording (0/1 = no/yes)
% - widefield: if there was widefield recording (0/1 = no/yes)
% - ephys: if there was electrophysiology recording (0/1 = no/yes)

% [EXERCISE]
% 1) How many days were recordings performed on AP010?

get_unique_days= unique({recordings.day}); %% get all the unique days
num_days= length(get_unique_days); %% return the number of unique days- 16

% 2) How many recordings were performed on AP010 across all days?
all_recordings= {recordings.recording}; %% get a cell that contains the recordings 
num_recordings = sum(cellfun('prodofsize', all_recordings)); %%cellfun(prodofsize) returns the number of elements in each cell- 37

asdf = {[1,2,3],[1,2]};
asdf2 = {1,2};
cellfun(@(x,y) sum(x*y),asdf,asdf2) %pairwaise multiplication

@(x) x+1

all_recordings= length([recordings.recording]);

A= [0,1,2];
~all(A)
% 3) How many unique Bonsai workflows were used, and what were they? %% 
all_workflows= {recordings.workflow}; %% get a cell that contains all workflows for each day
unique_workflows= unique(cat(1,all_workflows{:})); %% concatenates the elements of the array in the first dimension, and then returns the unique values
% num_unique_workflows= length(unique_workflows); %% length of the unique elements-5 


num_unique_workflows= length(unique(cat(1,recordings.workflow)));

% This line finds all recordings of AP010 on 2023-08-24: 
% (workflow is not entered)
recordings2 = plab.find_recordings('AP010','2023-08-24');

% [EXERCISE]
% Which recording modalities were recorded this day? 
counter_mousecam=0;
counter_widefield=0;

for i=1:length(recordings2.recording) %% a loop that goes through the number of recording
    if recordings2.mousecam(i) == 1 %% if mousecam is true then add to the counter
        counter_mousecam = counter_mousecam +1;
    end
    if recordings2.widefield(i) ==1 %% if widefield is true then add to the counter
        counter_widefield= counter_widefield+1;
    end
end

all(recordings2.mousecam)
onset_wheel_move_position= [0,4];
all(~onset_wheel_move_position)

if recordings2.ephys==1 %% check if electrophysiology has been used
   fprintf('Out of %d recordings, mousecam has been used %x times and widefield %d times with electrophysiology \n',length(recordings2.recording),counter_mousecam,counter_widefield)
else 
    fprintf('Out of %d recordings, mousecam has been used %x times and widefield %d times without electrophysiology \n',length(recordings2.recording),counter_mousecam,counter_widefield)
end

% This line finds all recordings of AP010 with workflow 'lcr_passive':
% (date is not entered)
recordings3 = plab.find_recordings('AP010',[],'lcr_passive');

% [EXERCISE]
% 1) How many recordings of this workflow included widefield? 
widefield_recordings= {recordings3.widefield};
number_widefield_recordings= length(find([widefield_recordings{:}]==1));

widefield_recordings= length(cat(2,recordings3.widefield))

% 2) How many recordings of this workflow included electrophysiology? 
ephys_recordings= {recordings3.ephys};
number_ephys_recordings= length(find([ephys_recordings{:}]==1));
% 3) Make a new variable which is a subset of 'recordings' with widefield
% make a new variable that contains only the widefiled recordings from the

store_recordings= {recordings3.recording};
only_widefield_recordings = cell(size(widefield_recordings)); %%initalize the variable
for i=1:length(widefield_recordings) %% loop through the length of the widefield recordings
    for j=1:length(widefield_recordings{i}) %% another loop because we have cells with >1 elements
    if widefield_recordings{i}(j) == 1 %% check if widefield recording took place
        only_widefield_recordings{i}{j}= store_recordings{i}(j); %% take the element from recordings and place it into the new variable
    end
    end
end


widefield_recordings= {recordings3.widefield};

new_var =cellfun(@any,widefield_recordings);


new_var = nan(size(widefield_recordings));
for i=1:length(widefield_recordings)
     new_var(i)=any(widefield_recordings{i});
end



widefield_recordings(10)
widefield_recordings{10}(2)

recordings3(10).widefield(2)

a = rand(0,4,3)

% This line finds recordings of AP010 with workflow 'lcr_passive' on
% 2023-08-16 (all fields are entered)
recordings4 = plab.find_recordings('AP010','2023-08-16','lcr_passive');

% [EXERCISE]
% This day has 2 recordings of the same workflow. This usually is
% because there was an issue with the first one and it needed to be re-run.
% In the number order of that day's recordings, which ones were
% 'lcr_passive'? (hint: 'index')


% The 'workflow' can include * as a wildcard. For example, in the task
% workflow 'stim_wheel_right', there is a 'stage1' and 'stage2' version.
% This returns specifically recordings with stage1:
recordings5 = plab.find_recordings('AP010',[],'stim_wheel_right_stage1');
% And this returns recordings with either stage: 
recordings6 = plab.find_recordings('AP010',[],'stim_wheel_right*');

% [EXERCISE]
% Write code to return the day AP010 switched from stage1 to stage2.
% the second stage is encoded as a string in the workflow variable
store_workflows= {recordings6.workflow}; %store the workflows

for i=1:length({recordings6.workflow}) %loop through all the workflows
    workflow_string=string(store_workflows{i}); %convert the indexed workflow from cell type to string
    if contains(workflow_string,'stage2') == 1 % if the current workflow contains stage2 then true
        day_switch= recordings6(i).day; % store the day in which the animal switched
        break % break from the loop so only the initial day in which the switch occured will be saved
    end
end

onset_wheel_move_position=contains(string([recordings6.workflow]),'stage2');
switch_day_index= find(onset_wheel_move_position,1,"first");
recordings6(switch_day_index).day


recordings6(find(contains(string([recordings6.workflow]),'stage2'), ...
    1,'first')).day

% jfoidpsjafpioasdfjpioajfpioajfiopasdpfiosdapfiojpsiodafjpoafjpoa
% asdfiojjfidosj afiasdofjsdafsdafjsdafiosjdafiosdafj jiofdsafj jiofdsajf
% sidof jiofdsa jiojfjiojiojfiowjefioajaefa jio jioij

onset_wheel_move_position = ...
    3;

% do this to that thingHere's a thing
onset_wheel_move_position = 5; % threshold in mm
onset_wheel_move_position*2; % double this to give extra leeway


%% Constructing server paths

% The class 'plab.locations' is used to find or construct paths on the server.

% --- Standardized generic locations

% This is the general path to the server:
plab.locations.server_path;
% This is the data path on the server:
plab.locations.server_data_path;
% This is where local data should be kept on each computer: 
plab.locations.local_data_path;

% [EXERCISE]
% Use the function 'fullfile' and the above information to generate the
% path: server/Users/(your name)/test_path
fullfile(plab.locations.server_path,"Users/Haron_Avgana/test_path");

% --- Recording locations

% The method 'plab.locations.filename' is used to construct paths/filenames
% for recording data
% (type 'help plab.locations.filename' for documentation). 
% The syntax is:   constructed_filename = filename('server' or 'local',animal,rec_day,rec_time,folder1,...,folderN)

% For example: 
% This constructs the path for animal AP001 on day 2000-01-01 at 12:00:
plab.locations.filename('server','AP001','2000-01-01','1200');
% Input arguments can be added or removed as necessary, for example:
% This constructs the path for the above recording and the subfolder
% 'mousecam'
plab.locations.filename('server','AP001','2000-01-01','1200','mousecam');
% This constructs the path for the above animal and day in the folder
% 'ephys' (not in a recording time folder):
plab.locations.filename('server','AP001','2000-01-01',[],'ephys');
% And this example constructs a path 'histology' in the animal folder (not
% in a day folder): 
plab.locations.filename('server','AP001',[],[],'histology');

% Note that paths with plab.locations.filename are constructed whether or
% not the folder/file exists (e.g. above example paths do not exist).

% [EXERCISE] 
% Use 'plab.find_recording' to find the recording time for AP010 on
% 2023-08-10 with workflow 'lcr_passive', then use
% 'plab.locations.filename' to construct the path to the 'widefield'
% subfolder within the folder for that recording. Write a line to check
% whether that path exists on the server.

recording_time= cell2mat(plab.find_recordings('AP010','2023-08-10','lcr_passive').recording); 
% get as a char the recording time for 'lcr_passive' workflow

widefield_folder= plab.locations.filename('server','AP010','2023-08-10',recording_time,'widefield');
% a string that contains the path to the widefield folder within that recording day

if exist(widefield_folder,'dir')
    % if file exists, argument ==7
     fprintf('The path %s exists on the server\n',widefield_folder);
else
    fprintf('The path %s does not exists on the server\n',widefield_folder);
end


%% Loading data, and data types

% There isn't standardized lab code for loading data, since this is often
% customized for each person depending on their needs. This block demos my
% code, which can be used as a template, or as-is if it works for you (note
% that it's subject to regular changes).

% My loading code is in my repository (petersaj/AP_scripts_peterslab) and
% is called with 'ap.load_recording' after defining animal/day/time. This
% loads an example dataset with passive stimuli, widefield, and ephys. Run
% this data, then each data type will be explained below:
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
verbose = true; % this turns on/off progress display in command line
ap.load_recording;


% --- Timelite

% "Timelite" is our GUI for recording analog signals with a DAQ. These
% include things like wheel movement, stimulus screen photodiode, camera
% frame captures, etc. It is saved as:
% server\Data\animal\day\time\timelite.mat
% as the structure 'timelite':
timelite

% The structure 'timelite' contains these structures: 
% - daq_info: information about the DAQ recording
% - data: the recorded data as N timepoints x M signals
% - timestamps: the timestamps for each data point as N timepoints x 1

% The names of the recorded signals are under
% 'timelite.daq_info(channel number).channel_name'. For example, the
% photodiode is recorded under channel 5, and the name can be found by:
photodiode_channel = 5;
timelite.daq_info(photodiode_channel).channel_name

% The corresponding data is in the 5th column of 'data'. For example, this
% plots photodiode data by time: 
figure;plot(timelite.data(:,photodiode_channel));
xlabel('Time (s)');
ylabel('Voltage');
title('Photodiode signal');


% (The photodiode is over a square on the screen that turns white/black to
% indicate changes to the stimulus).

% [EXERCISE] 
% 1) Write code to identify the channel number for 'widefield_camera'
% (widefield camera exposures) and plot the data for that channel.
widefield_channel= find(contains({timelite.daq_info.channel_name},'widefield_camera'),1)
% get the index of the channel name that corresponds to widefield_camera
figure; plot(timelite.timestamps,timelite.data(:,widefield_channel));
xlabel('Time (s)');
ylabel('Voltage');
title('Widefield Camera');
hold on;


% 2) Find the timestamps when the widefield camera starts each exposure
% (times when the signal flips from low to high). 

timestamp_widefield= timelite.data(:,widefield_channel);
% extract the time stamps from the widefield channel
widefield_threshold= 2;
% a threshold for signal flip
widefield_signal_flip= find(timestamp_widefield(1:end-1)<widefield_threshold...
    & timestamp_widefield(2:end)>widefield_threshold);

signal_flip_timestamp= timelite.timestamps(widefield_signal_flip+1);
% add +1 to the signal flip indicies to get the timepoint in which the
% signal was flipped on

mean_widefield_signal_on= mean(timestamp_widefield(widefield_signal_flip+1));
% calculate the mean signal when flipped on
plot(signal_flip_timestamp,mean_widefield_signal_on,'r.',MarkerSize=10);
% plot the timestamps in which the signal is on at the averged widefield
% signal when on


% 3) Check that the timestamps in 2 match the variable
% 'widefield_expose_times' (this is created in ap.load_recording)
isequal(signal_flip_timestamp,widefield_expose_times);

% One important signal recorded in Timelite is the wheel which the mouse
% can turn left and right. This is in the channel called 'wheel', which
% represents the position of the wheel relative to the start of the
% recording. Besides position, it is also helpful to have wheel velocity,
% and binary classification of whether the wheel is moving or not. This
% information is calculated by 'ap.parse_wheel', and outputs these
% variables - 
% wheel_velocity: velocity of the wheel for each Timelite timestamp
% wheel_move: binary vector representing movement (1) or quiescence (0)

% [EXERCISE] 
% 1) Plot time vs. raw Timelite wheel position and the wheel velocity on
% separate axes. Using 'wheel_move', plot the wheel velocity only when the
% wheel is moving in a separate color (by setting quiescent times to NaN,
% which are not plotted).
figure;plot(timelite.timestamps,wheel_velocity);
hold on;


movement_wheel_move_indices= find(wheel_move);

quiescent_wheel_velocity= wheel_velocity;
quiescent_wheel_velocity(~wheel_move)= NaN;

plot(timelite.timestamps,quiescent_wheel_velocity);

wheel_move= double(wheel_move);
% convert the variable from logical to double
% where wheel_move is zero replace with NaN
figure;plot(timelite.timestamps,wheel_move)
xlabel('Time (ms)');
ylabel('Wheel Movement');
hold on;


% 2) Using 'wheel_move', find the onset and offset times of all movements.
% Plot these as lines (colored differently for onset and offset) on your
% previous velocity plot.

onset_wheel_move_indices= find(~wheel_move(1:end-1) & wheel_move(2:end))+1;

onset_wheel_move_indices = find(diff(wheel_move) == 1) + 1;
offset_wheel_move_indices = find(diff(wheel_move) == -1) + 1;

% gives the indexes of the wheel movement onset, switch from 0 to 1, add
% one to the index as we want the time of the begining of the movement

time_onset_wheel_move= timelite.timestamps(onset_wheel_move_indices);
xline(time_onset_wheel_move,'r');

% 
% onset_wheel_move_position=nan(size(wheel_move));
% % initalize new var for onset position
% 
% onset_wheel_move_position(onset_wheel_move_indices)=1;
% % find the timestamps for onset and make the value 1
% plot(timelite.timestamps,onset_wheel_move_position,'r.');
% plot the onset times as red lines

offset_wheel_move_indicies= find(wheel_move(1:end-1) & ~wheel_move(2:end));
% find when the movement switches from 1 to 0

time_offset_wheel_move= timelite.timestamps(offset_wheel_move_indicies);
xline(time_offset_wheel_move,'g');

% initalize new var for off position

offset_wheel_move_position(offset_wheel_move_indicies)=1;
%find the timestamp for offset and make the value 1
plot(timelite.timestamps,offset_wheel_move_position,'g.');


% --- Bonsai

% We use the program Bonsai to run our stimuli and tasks
% (https://bonsai-rx.org/). Bonsai is an actively supported, highly
% customizable, visual interfaced framework for designing experiments. 

% Bonsai files are called 'workflows'. ap.load_recording loads the name of
% the currently loaded workflow as 'bonsai_workflow'. This one is
% 'lcr_passive', which is passive stimuli presented on the left, center,
% and right screens.
bonsai_workflow

% Bonsai can change how it saves data, and loading scripts can specify how
% that data is loaded and parsed. In this demo, ap.load_recording loads
% data from Bonsai as 'trial_events':
trial_events
% which is a structure containing:
% - parameters: parameter values corresponding to the whole workflow, e.g.
% the duration that each stimulus was displayed was:
trial_events.parameters.StimDuration
% - values: an N trials x 1 array of values of saved events, e.g. the
% order of X-positions for the 3 stimuli presented on trial 4 was:
trial_events.values(4).TrialStimX
% - timestamps: an N trials x 1 array of timestamps of saved events from
% 'values', e.g. the time each of 3 stimuli was turned on/off on trial 4
% was:
trial_events.timestamps(4).StimOn

% Note that timing information in Bonsai is only approximate (when the
% software gave the command), not actual (when the stimulus was physically
% drawn on the screen). All time information should be used from Timelite,
% and only qualitative information should be used from Bonsai (e.g. which
% stimulus was presented on a trial).

% [EXERCISE] 
% The stimulus onset times are loaded by ap.load_recording as
% 'stimOn_times'. Write code to pull out a subset of these timestamps
% corresponding to a stimulus X position (TrialStimX) of 90 (meaning it was
% on the right-hand screen).

stimulusX_position= cat(1,trial_events.values.TrialStimX);
right_stim_timepoint= stimOn_times(stimulusX_position==90);

position90_row_positions = cellfun(@(cell_row) find(cell_row == 90), y, 'UniformOutput', false);
% get the row position where the position was 90

extract_odd_indices = @(cell_values) cell_values(1:2:end);
% anonymous fucntion that extract the odd indicies (i.e onset times)

stim_onset_cells = cellfun(extract_odd_indices, {trial_events.timestamps.StimOn}, 'UniformOutput', false);
% a 1x50 cell that contains the starting points for each presentation




% --- Mousecam

% We record video of the front of the mouse during all experiments. The
% filename is loaded in ap.load_recording as:
mousecam_fn
% and the timestamps of each frame (in Timelite clock) is:
mousecam_times

% Mousecam images can be read into Matlab with a VideoReader object. For
% example:
mousecam_vr = VideoReader(mousecam_fn); % create VideoReader object

load_frame = 1; % define frame to read
mousecam_im = read(mousecam_vr,load_frame); % read frame

figure; % create a figure
imagesc(mousecam_im); % draw the image (sc = scaled colors)
axis image % make the x/y axes have equal aspect ratios
colormap('gray'); % set the colormap to gray

% You can also load in multiple frames at the same time by defining the
% start/end frame. For example:
load_frames = [1,10]; % define frame interval to read
mousecam_im = read(mousecam_vr,load_frames); % read frame
% (the VideoReader by default loads in multiple images in the 4th
% dimension, allowing for colors in the 3rd dimension. Our images are
% grayscale, so we can use 'squeeze' to remove the singleton 3rd dimension.
% Look at the size of the natively loaded data: 
size(mousecam_im)
% and then the size of the 'squeezed' data:)
mousecam_im = squeeze(mousecam_im); 
size(mousecam_im)
% This is an example of how to view a 3D matrix using my function
% ap.imscroll, which plots each 'page' (dim 1+2) which can be scrolled
% through in dim 3:
ap.imscroll(mousecam_im);

% [EXERCISE] 
% Create an average mousecam movie for -20:+20 frames around stimulus X =
% 90 presentations, and separately for stimulus X = -90. Is there a
% difference in behavior when these two stimuli are presented?

frame_band= [-20:20];

left_stim_timepoint= stimOn_times(stimulusX_position==-90);
% find the timepoints in which stimulus was presented to the left

min_indicies_right_stim= nan(size(right_stim_timepoint));
min_indicies_left_stim= nan(size(left_stim_timepoint));

for i=1:length(right_stim_timepoint)
    % assuming that number of right and left presentations is equal

    [~,min_indicies_right_stim(i)]= min(abs(mousecam_times - right_stim_timepoint(i)));
    [~,min_indicies_left_stim(i)]= min(abs(mousecam_times - left_stim_timepoint(i)));
    % find the indices of the minimum difference values for both
    % right and left stimulus presentations

   

end

added_frames_right_indices = min_indicies_right_stim+ frame_band([1,end]);
added_frames_left_indices = min_indicies_left_stim + frame_band([1,end]);
% add +-20 timepoints to the timepoint of stimulus presentation


mouscam_im_avg_right= cell(size(right_stim_timepoint));
mouscam_im_avg_left= cell(size(right_stim_timepoint));

for i=1:length(right_stim_timepoint)
    mouscam_im_right= read(mousecam_vr,added_frames_right_indices(i,:));
    mouscam_im_left= read(mousecam_vr,added_frames_left_indices(i,:));

    mouscam_im_avg_right {i}= mouscam_im_right;
    mouscam_im_avg_left {i}= mouscam_im_left;
end


mouscam_im_avg_right= cellfun(@squeeze,mouscam_im_avg_right,'UniformOutput',false);
mouscam_im_avg_left= cellfun(@squeeze,mouscam_im_avg_left,'UniformOutput',false);

% above you have a cell array - you want to convert this into a ND matrix,
% and then average across the dimension you care about

% matrix_mouscam_im_avg_right= [];
matrix_mouscam_im_avg_right = cat(4, mouscam_im_avg_right{:});
average_mouscam_im_right = mean(matrix_mouscam_im_avg_right, 4);


matrix_mouscam_im_avg_left= [];
matrix_mouscam_im_avg_left = cat(4, mouscam_im_avg_left{:});
average_mouscam_im_left = mean(matrix_mouscam_im_avg_left, 4);

% once you've done that - try writing your for loop storing your frames
% without a cell array, which will let you skip a step


figure;imagesc(average_mouscam_im_right(:,:,1));
colormap gray;

ap.imscroll(average_mouscam_im_right,frame_band);
axis image;

average_mouscam_im_right_diff= abs(diff(average_mouscam_im_right,1,3));
ap.imscroll(average_mouscam_im_right_diff);
axis image;
colorbar;



%%
% convince Andy that AP016 has a stimulation association but not AP018 look
% into AP017 too
% Super H.W; when AP016 developed this association?

% load the first day of training

% My loading code is in my repository (petersaj/AP_scripts_peterslab) and
% is called with 'ap.load_recording' after defining animal/day/time. This
% loads an example dataset with passive stimuli, widefield, and ephys. Run
% this data, then each data type will be explained below:
animal = 'AP016';
verbose = true; % this turns on/off progress display in command line
rec_days= {'2024-02-12','2024-02-14','2024-02-15','2024-02-16','2024-02-19','2024-02-21','2024-02-22','2024-02-23'}; %define recording days
rec_times= {'1345','1615','1453','1349','1329','1420','1401','1435'};% define recording times

for day_indx= 1:length(rec_days)
    rec_day=rec_days{day_indx};
    rec_time= rec_times{day_indx};

    load_parts.mousecam= true;
    ap.load_recording;
    
    
    
    % --- Timelite
    
    % "Timelite" is our GUI for recording analog signals with a DAQ. These
    % include things like wheel movement, stimulus screen photodiode, camera
    % frame captures, etc. It is saved as:
    % server\Data\animal\day\time\timelite.mat
    % as the structure 'timelite':
    timelite
    
    % --- Bonsai
    
    % We use the program Bonsai to run our stimuli and tasks
    % (https://bonsai-rx.org/). Bonsai is an actively supported, highly
    % customizable, visual interfaced framework for designing experiments. 
    
    % Bonsai files are called 'workflows'. ap.load_recording loads the name of
    % the currently loaded workflow as 'bonsai_workflow'. This one is
    % 'lcr_passive', which is passive stimuli presented on the left, center,
    % and right screens.
    bonsai_workflow
    
    % Bonsai can change how it saves data, and loading scripts can specify how
    % that data is loaded and parsed. In this demo, ap.load_recording loads
    % data from Bonsai as 'trial_events':
    trial_events
    
    % how to calculate whether the animal is doing anticipatory licking or not?
    % I can perhaps highlight the pixels surrounding the tongue and tube area
    % and see how much change in pixels there is before stimulus presentation
    
    
    % find the timepoints in which the stimulus was presented based on changes
    % in photodiode
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
    % get all the timepoints in which stimulus was presented and then add the
    % initial one
    
    stim_on_timepoints_all= nan(size(timelite.timestamps));
    stim_on_indices = find(ismember(timelite.timestamps, stim_on_timepoints));
    stim_on_timepoints_all(stim_on_indices)=4;
    % intialize a variable with all the timepoints and then change only
    % stimulus onset indecies
    
    % figure;plot(timelite.timestamps,photodiode_values);
    % hold on;
    % plot (timelite.timestamps,stim_on_timepoints_all,'o','MarkerSize',5);
    % 
    
    % load mousecam
    mousecam_fn
    mousecam_times
    mousecam_vr = VideoReader(mousecam_fn);
     
    stim_onset_timepoints= timelite.timestamps(find(~isnan (stim_on_timepoints_all)));
     % get the timepoints of stimulus onset
    
    % min_indicies_passive_stim=NaN(size(stim_onset_timepoints));
    % for i=1:length(stim_onset_timepoints)
    %     [~,min_indicies_passive_stim(i)]= min(abs(mousecam_times -stim_onset_timepoints(i)));
    % 
    % end
    
    % using interp, find the nearst frame that corresponds to stimulus onset
    % time
    extrap_min_indicies_passive_stim= interp1(mousecam_times,1:length(mousecam_times),stim_onset_timepoints,"nearest");
    
    %frame band to look into then
    frame_band=[-20,100];
    added_frames_passive_conditioning= extrap_min_indicies_passive_stim +frame_band;
    
     % initalize the dimensions based on first passive trial
    first_passive_trial= squeeze(read(mousecam_vr,added_frames_passive_conditioning(1,:)));
    mouscam_passive_conditioning_dimensions= [size(first_passive_trial),150];
    mouscam_im_passive_conditioning= zeros(mouscam_passive_conditioning_dimensions,'single');
    
    %loop through each trial 
    for i=1:length(stim_onset_timepoints)
        mouscam_im_passive_conditioning(:,:,:,i)= read(mousecam_vr,added_frames_passive_conditioning(i,:));
        disp(i)
    end
    
    %frame numbers for labeling 
    frame_numbers = frame_band(1):frame_band(2);
    %average across 4th dimension  
    avg_mouscam_im_passive_condiotnining= mean(mouscam_im_passive_conditioning,4);
    % ap.imscroll(avg_mouscam_im_passive_condiotnining,frame_numbers);
    % axis equal;
    
    % plot the differences between frames
    avg_mouscam_im_passive_diff= abs(diff(avg_mouscam_im_passive_condiotnining,1,3));
    % ap.imscroll(avg_mouscam_im_passive_diff,frame_band(1):frame_band(2)-1);
    % colorbar
    % axis equal;
    
    % plot the mean differences 
    mean_diff_per_frame{day_indx} = smooth(mean(avg_mouscam_im_passive_diff,[1,2]));
    % figure;plot(frame_band(1):frame_band(2)-1,squeeze(mean_diff_per_frame));
    % 
    % % add the stimulus onset and sucrose release timepoints to the plot
    % hold on
    % line([0, 0], ylim, 'Color', 'r', 'LineStyle', '--')
    % line([46, 46], ylim, 'Color', 'b', 'LineStyle', '--')
    % hold off
    % xlabel('Frame');
    % ylabel('Mean Pixel Difference');
    % title('Mean Pixel Difference Across Frames');
    % legend('Mean Pixel Difference', 'Stimulus Onset', 'Location', 'best')

end

% create a vector in which the mean differences for each day is stored in
% each row
mean_diff_per_frame_matrix = cell2mat(mean_diff_per_frame);
mean_diff_per_frame_matrix = mean_diff_per_frame_matrix'; %transpose for plotting



% find the timepoints in which stim and reward were delivered
frame_rate= 30;
reward_onset_after_stimulus= 1; %reward appears 1 second after stim onset
time_seconds= (frame_numbers * (frame_rate))/1000;
stimulus_onset_frame= find(time_seconds==0); %find stimulus onset frame
[~,reward_onset_frame]= min(abs(time_values-reward_onset_after_stimulus)); %find the timepoint of reward

% find the mean pixel difference for the min quiscent period
minimum_quiscent_time= 0.5; % the min quiscent time
[~,quiscent_onset_time]= min(abs(time_values+minimum_quiscent_time)); %plus as we know the time is negative

% find the corresponding frame index
quiscent_onset_frame_index = find(frame_numbers== frame_numbers(quiscent_onset_time));

% calculate the mean from quiscent onset to stimulus onset -1 for each day 
mean_quiscent_baseline= mean(mean_diff_per_frame_matrix(:,quiscent_onset_frame_index:(stimulus_onset_frame-1)),2);

% substract baseline from the overall mean difference matrix
baseline_subtracted_mean_diff_frame_matrix= mean_diff_per_frame_matrix-mean_quiscent_baseline;


% heatmap visualization of days
imagesc(time_values,1:length(rec_days),baseline_subtracted_mean_diff_frame_matrix);
axis xy;
colormap gray
colorbar;
% add the recording day labels 
yticks(1:length(rec_days));
yticklabels(rec_days);
line([time_values(stimulus_onset_frame),time_values(stimulus_onset_frame)], ylim, 'Color', 'b', 'LineStyle', '--')
line([time_values(reward_onset_frame),time_values(reward_onset_frame)], ylim, 'Color', 'r', 'LineStyle', '--')
xlabel('Time (seconds)');
ylabel('# Day');
title(['Mean Pixel Difference Substraced from Baseline for ' animal ' Across ' num2str(length(rec_days)) ' Days']);

% Add grid lines for better visualization
grid on;


% figure for mean pixel difference across all days
figure; plot(time_values(1:end-1),mean((baseline_subtracted_mean_diff_frame_matrix),1));
line([time_values(stimulus_onset_frame),time_values(stimulus_onset_frame)], ylim, 'Color', 'b', 'LineStyle', '--')
line([time_values(reward_onset_frame),time_values(reward_onset_frame)], ylim, 'Color', 'r', 'LineStyle', '--')
xlabel('Time (seconds)');
title(['Mean Pixel Difference Substraced from Baseline for ' animal ' Across ' num2str(length(rec_days)) ' Days']);


% plot all figures 
h=figure;
for i=1:length(rec_days)
    subplot(2, 4, i);
    plot(time_values(1:end-1),baseline_subtracted_mean_diff_frame_matrix(i,:));
    h1= xline(time_values(stimulus_onset_frame),Color='b');
    xline(time_values(reward_onset_frame),Color='r');
    title(sprintf('Plot for Day %d', i));  % Use sprintf to format the title string
end

linkaxes(h.Children); % makes the y and x axes alligned across all the plots



handle_1=figure;
axes_handle= axes(handle_1);
axes_handle.ColorOrder = copper(length(rec_days));

hold on;
for i=1:length(rec_days)
    plot(time_values(1:end-1),baseline_subtracted_mean_diff_frame_matrix(i,:));

end
h1= xline(time_values(stimulus_onset_frame),Color='b');
xline(time_values(reward_onset_frame),Color='r');
title(sprintf('Plot for %d Days', length(rec_days)));

% --- Widefield and electrophysiology

% If available in a recording, ap.load_recording will load and prepare
% widefield and electrophysiology data. For demos on working with that
% data, see: 
% - plab_widefield_demo
% - plab_ephys_demo

% I have a function to scroll through experiment data for exploratory
% purposes, which displays the mousecam (left), widefield (center), and
% multiunit ephys (right), which can be scrolled through time:
ap.expscroll

% [EXERCISE] 
% There is a drop of sucrose available to the mouse at the beginning of
% this recording (left over from the task done previously). Which part of
% the ephys probe has increased activity when the mouse consumes this?

% --- Loading partial datasets
 
% Recordings can be part-loaded by ap.load_recording if only some data is
% necessary by creating a 'load_parts' structure, which can toggle loading
% with fields 'widefield','ephys','mousecam'. Timelite and behavior are
% always loaded. If 'load_parts' exists, any field not defined will default
% to not loading. For example, in this dataset, you can turn off loading of
% widefield and ephys data by doing:
clear all % (clears data loaded above)
animal = 'AP005';
rec_day = '2023-06-21';
rec_time = '1859';
load_parts.mousecam = true; % (this is the new line)
verbose = true;
ap.load_recording; % (this is much faster now and omits widefield and ephys)




















