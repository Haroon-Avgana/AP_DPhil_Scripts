% NOTE ON INSTALLING SOFTWARE: 
% install these: https://dcam-api.com/hamamatsu-software/
% in matlab: imaqregister('C:\Users\peterslab\AppData\Roaming\MathWorks\MATLAB Add-Ons\Toolboxes\Hamamatsu Image Acquisition\hamamatsu.dll')

function widefield

% Clear existing matlab image aquisitions objects
imaqreset

% Connect and set up camera
% NOTE: bit depth, binning, and speed mode, set by format
% To see list of formats:
% imaqhwinfo('hamamatsu').DeviceInfo.SupportedFormats';
cam_DeviceName = imaqhwinfo('hamamatsu').DeviceInfo.DeviceName;
video_object = videoinput('hamamatsu',cam_DeviceName,'MONO16_BIN2x2_1024x1024_Std');
src = getselectedsource(video_object);

% Set input trigger (expose while trigger high)
src.TriggerPolarity = "positive";
src.TriggerSource = "external";
src.TriggerActive = "level";

% Set trigger mode (global reset: all lines start together, end asynch)
src.TriggerGlobalExposure = "globalreset";

% Set outputs
% - Trigger ready: no pixels expose, turn on screens
src.OutputTriggerKindOpt1 = "triggerready";
src.OutputTriggerPolarityOpt1 = "positive";
% - Exposure: all lines exposing, turn on illumination
src.OutputTriggerKindOpt2 = "exposure";
src.OutputTriggerPolarityOpt2 = "positive";

%% Set up GUI

gui_fig = figure('MenuBar','none','Units','Normalized', ...
    'Position',[0.01,0.1,0.5,0.5],'color','w','Colormap',gray);

% Set up preview image (two colors + moving difference)
n_axes = 2;
im_axes = gobjects(n_axes,1);
im_preview = gobjects(n_axes,1);
for curr_axes = 1:n_axes
    im_axes(curr_axes) = axes(gui_fig,'Position', ...
        [(curr_axes-1)*(1/n_axes),0.05,1/n_axes,0.8]);
    im_preview(curr_axes) = imshow(zeros(video_object.VideoResolution,'uint16'),Border="tight");
end

% Set up saturation display axes
n_axes = 2;
im_saturation_axes = gobjects(n_axes,1);
im_saturation_preview = gobjects(n_axes,1);
for curr_axes = 1:n_axes
    im_saturation_axes(curr_axes) = axes(gui_fig,'Position', ...
        [(curr_axes-1)*(1/n_axes),0.05,1/n_axes,0.8]);
    im_saturation = im2uint8(cat(3,true(video_object.VideoResolution), ...
        false([video_object.VideoResolution,2])));
    im_saturation_preview(curr_axes) = imshow(im_saturation,Border="tight");
    % (make it not selectable so ROI can be drawn underneath)
    set(im_saturation_preview(curr_axes),'PickableParts','none');
end

% Link all axes
linkaxes(vertcat(im_axes,im_saturation_axes),'xy')

% Metadata text
metadata_text = uicontrol('Style','text','String','Frame info', ...
    'Units','normalized','Position',[0,0,1,0.05], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontSize',12, ...
    'FontName','Consolas');

% Status text
status_text_h = uicontrol('Parent',gui_fig,'Style','text', ...
    'FontSize',12,'FontName','Courier','HorizontalAlignment','left', ...
    'Units','normalized','BackgroundColor','w','Position',[0,0.9,1,0.1]);

% Control buttons
button_fontsize = 12;
button_position = [0.85,0.9,0.15,0.1];
clear controls_h
controls_h(1) = uicontrol('Parent',gui_fig,'Style','pushbutton','FontSize',button_fontsize, ...
    'Units','normalized','Position',button_position,'BackgroundColor','w', ...
    'String','Set ROI','Callback',{@set_roi,gui_fig},'Enable','on');

% Color limit setting
clim_text_h = uicontrol('Parent',gui_fig,'Style','text','FontSize',button_fontsize, ...
    'Units','normalized','Position',[0.85,0,0.05,0.05],'String','White: ', ...
    'BackgroundColor','w','HorizontalAlignment','right');
clim_h = uicontrol('Parent',gui_fig,'Style','edit','FontSize',button_fontsize, ...
    'Units','normalized','Position',[0.9,0,0.10,0.05],'String',num2str(2^16));

% Start listener for experiment controller
update_status_text(status_text_h,'Connecting to experiment server');
try
    client_expcontroller = tcpclient("163.1.249.17",plab.locations.widefield_port,'ConnectTimeout',2);
    configureCallback(client_expcontroller, "terminator", ...
        @(src,event,x) read_expcontroller_data(src,event,gui_fig));
    update_status_text(status_text_h,'Listening for start');
catch me
    % Error if no connection to experiment controller
    update_status_text(status_text_h,'Error connecting to experiment server');
    warning(me.identifier,'Widefield -- Cannot connect to experiment controller: \n %s',me.message)
end

% Set function when frames acquired
video_object.FramesAcquiredFcn = {@upload_data,gui_fig};
video_object.FramesAcquiredFcnCount = 1;
video_object.FramesPerTrigger = Inf;

% If a previous ROI was used, set that ROI
matlab_settings = settings;
if hasGroup(matlab_settings,'widefield')
    video_object.ROIPosition = matlab_settings.widefield.roi.PersonalValue;
    % Set axes tight to ROI
    for curr_axes = 1:length(im_axes)
        axis(im_axes(curr_axes));
        axis tight;
    end
else
    % (if no previous ROI, create setting)
    addGroup(matlab_settings,'widefield');
    addSetting(matlab_settings.widefield,'roi');
    matlab_settings.widefield.roi.PersonalValue = video_object.ROIPosition;
end

% Start the camera
% (keeping the adapter running benefits in removing abberant expose signals
% on starting and removes need for preview button - the downside is that
% the frame number and duration count up from this point onwards rather
% than being recording-specific)
start(video_object);

% Store gui data
gui_data.video_object = video_object;
gui_data.gui_fig = gui_fig;

gui_data.im_axes = im_axes;
gui_data.im_preview = im_preview;

gui_data.im_saturation_axes = im_saturation_axes;
gui_data.im_saturation_preview = im_saturation_preview;

gui_data.metadata_text = metadata_text;
gui_data.status_text_h = status_text_h;
gui_data.controls_h = controls_h;
gui_data.clim_h = clim_h;
if exist('client_expcontroller','var')
    gui_data.client_expcontroller = client_expcontroller;
end

% Initialize relative frame number / timestamp for display
gui_data.frame_info.start_frame_n = 0;
gui_data.frame_info.start_time = datetime('now');
gui_data.frame_info.last_frame_n = 0;

% Update GUI data
guidata(gui_fig,gui_data);


end

%% Button functions

function set_roi(obj,eventdata,gui_fig)
% Draw and set ROI

% Get GUI data
gui_data = guidata(gui_fig);

% Reset ROI
stop(gui_data.video_object);
gui_data.video_object.roi = [0,0,gui_data.video_object.VideoResolution];
start(gui_data.video_object);

% Draw new ROI and want until finished
update_status_text(gui_data.status_text_h,'Draw ROI');
roi = drawrectangle(gui_data.im_axes(1));
wait(roi);

roi_position = roi.Position;
if any(roi_position(2:3) == 0)
    roi_position = [0,0,gui_data.video_object.VideoResolution];
end
delete(roi);

% Set new ROI position
stop(gui_data.video_object);
warning('off','imaq:hamamatsu:invalidROI')
gui_data.video_object.roi = roi_position;
start(gui_data.video_object);

% Set tight axes
for curr_axes = 1:length(gui_data.im_axes)
    axis(gui_data.im_axes(curr_axes));
    axis tight;
end

% Store ROI for use on next load
matlab_settings = settings;
matlab_settings.widefield.roi.PersonalValue = gui_data.video_object.ROIPosition;

% Reset status
update_status_text(gui_data.status_text_h,'Listening for start');

% Update GUI data
guidata(gui_fig,gui_data);

end

%% Upload data functions

function cam_start(gui_fig,save_path)
% Prepare camera for saving

% Get GUI data
gui_data = guidata(gui_fig);

% Disable all buttons
set(gui_data.controls_h,'Enable','off');

% Set relative frame info for display
gui_data.frame_info.start_frame_n = gui_data.frame_info.last_frame_n;
gui_data.frame_info.start_time = datetime('now');

% Open files for saving (data, metadata)
save_data_filename = sprintf('%s_data.bin',save_path);
save_metadata_filename = sprintf('%s_metadata.bin',save_path);

gui_data.save_data_file = fopen(save_data_filename,'w');
gui_data.save_metadata_file = fopen(save_metadata_filename,'w');

% Update status text
set(gui_fig,'color',[1,0.7,0.7]);
update_status_text(gui_data.status_text_h,sprintf('Recording: %s',save_path));

% Update GUI data
guidata(gui_fig,gui_data);

end

function cam_stop(gui_fig)
% Stop recording

% Get GUI data
gui_data = guidata(gui_fig);

% Update status text
update_status_text(gui_data.status_text_h,'Ending recording');

% Upload any remaining frames
upload_data([],[],gui_fig);

% Close save files
fclose(gui_data.save_data_file);
fclose(gui_data.save_metadata_file);

% Clear save file variables
gui_data = rmfield(gui_data,{'save_data_file','save_metadata_file'});

% Enable buttons
set(gui_data.controls_h,'Enable','on');

% Update status text
set(gui_fig,'color','w');
update_status_text(gui_data.status_text_h,'Listening for start');

% Update GUI data
guidata(gui_fig,gui_data);

end

function upload_data(obj,eventdata,gui_fig)
% Display frames and write frames/timestamps to binary file

% Get GUI data
gui_data = guidata(gui_fig);

% Grab all available data off camera
[frames,frame_timestamps,frame_metadata] = ...
    getdata(gui_data.video_object,gui_data.video_object.FramesAvailable);

% If no data (overlapping upload function so data already grabbed), return
if isempty(frames)
    return
end

% If save file exists, write data/metadata to binary file
if isfield(gui_data,'save_data_file') && ~isempty(gui_data.save_data_file)
    % Write images to disk
    fwrite(gui_data.save_data_file,frames,class(frames));
    % Write metadata to disk
    % Metadata format (each frame): 
    % [image height, image width, frame number, timestamp]
    write_metadata = [repmat(size(frames,[1,2]),length(frame_metadata),1), ...
        vertcat(frame_metadata.FrameNumber), ...
        vertcat(frame_metadata.AbsTime)]';
    fwrite(gui_data.save_metadata_file,write_metadata,'double');
end

% Display metadata text
recording_frame_number = frame_metadata(end).FrameNumber - gui_data.frame_info.start_frame_n;
recording_timestamp = datetime('now') - gui_data.frame_info.start_time;
set(gui_data.metadata_text,'String',sprintf( ...
    'Frame %d, %s', ...
    recording_frame_number,recording_timestamp))

% Update preview image (only last frame, if multiple)
update_preview_idx = mod(1+recording_frame_number,length(gui_data.im_preview))+1;
set(gui_data.im_preview(update_preview_idx),'CData',frames(:,:,:,end));

user_clim = str2double(get(gui_data.clim_h,'String'));
if isnan(user_clim) || user_clim < 1
    user_clim = 2^16;
    set(gui_data.clim_h,'String',num2str(user_clim));
end
if max(max(clim(gui_data.im_axes(update_preview_idx)))) ~= user_clim
    clim(gui_data.im_axes(update_preview_idx),[0,user_clim]);
end

% Turn on saturation pixels
% (make saturation pixel mask same size, if different)
if ~all(size(gui_data.im_saturation_preview(update_preview_idx).CData,[1,2]) == ...
        size(frames(:,:,end)))
    im_saturation = im2uint8(cat(3,true(size(frames,[1,2])), ...
        false([size(frames,[1,2]),2])));
    set(gui_data.im_saturation_preview(update_preview_idx), ...
        'CData',im_saturation);
end
% (set saturated pixels as visible)
set(gui_data.im_saturation_preview(update_preview_idx),'AlphaData',frames(:,:,:,end) >= (2^16-1));

% Store frame number (to display relative to recording start, rather than start trigger)
gui_data.frame_info.last_frame_n = frame_metadata(end).FrameNumber;

% Update GUI data
guidata(gui_fig,gui_data);

end

function update_status_text(status_text_h,status)
% Update status text

curr_text = get(status_text_h,'String');
new_text = [{sprintf('Status: %s',status)};curr_text(2:end)];
set(status_text_h,'String',new_text);

end

%% Cross-computer functions

function read_expcontroller_data(client,event,gui_fig)
% Read message from experiment controller

% Get gui data
gui_data = guidata(gui_fig);

% Get message from experiment controller
expcontroller_message = readline(client);

if strcmp(expcontroller_message, 'stop')
    % If experiment controller sends stop, stop DAQ acquisition
    cam_stop(gui_fig);
else
    % If experiment controller send experiment info
    rec_info = jsondecode(expcontroller_message);

    % Set local save path (top-level folder: SVD run on concatenated data
    % and then parsed into experiment folders in preprocessing)
    save_path = ...
        plab.locations.filename('local', ...
        rec_info.mouse,rec_info.date,[],'widefield',sprintf('widefield_%s',rec_info.time));

    % Make local save directory
    if ~exist(fileparts(save_path),'dir')
        mkdir(fileparts(save_path))
    end
    
    % Start camera recording acquisition
    cam_start(gui_fig,save_path)
end

end


