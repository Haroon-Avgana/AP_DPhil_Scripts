%% Test scripts
%% Widefield preprocessing function
[U,Vrec,im_avg_color,frame_info] = AP_preprocess_widefield_hc('D:\LocalData');

% svd reconstruction
example_fluorescence_blue = AP_svdFrameReconstruct(cell2mat(U(1)), cell2mat(Vrec(1)));
example_fluorescence_violet = AP_svdFrameReconstruct(cell2mat(U(2)), cell2mat(Vrec(2)));

AP_imscroll(example_fluorescence_blue);
axis image;

AP_imscroll(example_fluorescence_violet);
axis image;

%% Save preprocessed widefield data on server

experiment_path = 'D:\LocalData';

% Set number of components to save
max_components_save = 2000;
n_components_save = min(max_components_save,size(U{1},3));

% Assume 2 colors in order of blue/purple
color_names = {'blue','violet'};

% Save frame information in experiment folder
frame_info_fn = fullfile(experiment_path,'widefield_frame_info');
save(frame_info_fn,'frame_info','-v7.3');

% Save mean images in experiment folder by color
for curr_color = 1:length(color_names)
    curr_mean_im_fn = fullfile(experiment_path, ...
        sprintf('meanImage_%s.npy',color_names{curr_color}));
    writeNPY(im_avg_color(:,:,curr_color),curr_mean_im_fn);
end

% Save spatial components in experiment (animal/day) folder by color
for curr_color = 1:length(color_names)
    curr_U_fn = fullfile(experiment_path, ...
        sprintf('svdSpatialComponents_%s.npy',color_names{curr_color}));
    writeNPY(U{curr_color}(:,:,1:n_components_save),curr_U_fn);
end

% Save temporal components in associated recording folders
for curr_recording = 1:size(Vrec,1)
    for curr_color = 1:length(color_names)
        curr_V_fn = fullfile(experiment_path, ...
            sprintf('svdTemporalComponents_%s.npy',color_names{curr_color}));
        writeNPY(Vrec{curr_recording,curr_color}(1:n_components_save,:)',curr_V_fn);
    end
end

%% - experiment split - works! - changed boundary to 20,000
[U,Vrec,im_avg_color,frame_info] = AP_preprocess_widefield_hc('C:\Users\peterslab\Documents\MATLAB\pause_video');

%% - check different ROIs
% centre
[U,Vrec,im_avg_color,frame_info] = AP_preprocess_widefield_hc('C:\Users\peterslab\Documents\MATLAB\centre_ROI');

%% Face camera preprocessing function
[timestamps, frame_num, flipper_signal] = AP_preprocess_face_camera("C:\Users\peterslab\Documents\MATLAB\face_camera\mmmgui_test\face.avi", 2);

%% Test wave save function
AM_save_h_file(trp_wave, 'sq_wave')
AM_save_h_file(trp_wave, 'trp_wave')

