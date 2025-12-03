function mousecam
% Mouse camera recording GUI
% 
% Currently assumes camera is Flir Chameleon3 CM3-U3-13Y3M
%
% Note about strobes: The camera constantly collects frames when triggering
% is off, and the only trigger mode is by frame rather than framerate. This
% means that subsets of strobes will always correspond to recorded frames.

%% Set up camera (specific to rig - make this breakout later?)

% Clears existing matlab image aquisitions objects
imaqreset

% Set camera properties (from imaqhwinfo)

% Turn on embedding and set imaging properties in PointGrey adaptor
cam_DeviceName = imaqhwinfo('pointgrey').DeviceInfo.DeviceName;
% two options for format: 
% - no binning ('F7_Mono8_1280x1024_Mode0')
% - 2x2 binning ('F7_Mono8_640x512_Mode1')
video_object = videoinput('pointgrey',cam_DeviceName,'F7_Mono8_640x512_Mode1');

% (undocumented feature to turn on all embedding)
% (BUT THIS DOESN'T SEEM TO WORK? Turn on embedding from firmware instead)
imaqmex('feature','-pointgreyEmbedMetadata',true)

src = getselectedsource(video_object);
src.FrameRate = 30;
src.ExposureMode = 'manual';
src.Exposure = -6;
src.GainMode = 'manual';
src.Gain = 10;
src.ShutterMode = 'manual';
src.Shutter = 6;
src.Brightness = 7;
% Turn strobe on for GPIO3
src.Strobe3 = 'On';
src.Strobe3Polarity = 'High';

% Reset adapter (keeps avove settings, but releases CPU resources)
imaqreset

% Set GPIO pin direction with GenTL adaptor
cam_DeviceName = imaqhwinfo('gentl').DeviceInfo.DeviceName;
video_object = videoinput('gentl',cam_DeviceName,'Mono8');
src = getselectedsource(video_object);
% GPIO2: flipper input
src.LineSelector = 'Line2';
src.LineMode = 'Input'; 
src.LineMode = 'Input'; % (need to double this line to work??)

video_object.FramesPerTrigger = Inf;
video_object.LoggingMode = 'disk&memory';

% Turn on trigger mode (only used to pause video at start/end of acq)
triggerconfig(video_object,'hardware');

%% Set up GUI

gui_fig = figure('MenuBar','none','Units','Normalized', ...
    'Position',[0.01,0.2,0.32,0.5],'color','w');

% Preview image and embedded information text
im_axes = axes(gui_fig,'Position',[0,0.05,1,0.8]);
im_preview = image(zeros(video_object.VideoResolution));

embedded_info_text = uicontrol('Style','text','String','Embedded header information', ...
    'Units','normalized','Position',[0,0,1,0.05], ...
    'BackgroundColor','w','HorizontalAlignment','left','FontSize',12, ...
    'FontName','Consolas');

setappdata(im_preview,'UpdatePreviewWindowFcn',@preview_cam);
setappdata(im_preview,'gui_fig',gui_fig);

preview(video_object,im_preview);
axis(im_axes);axis tight equal

% Status text
status_text_h = uicontrol('Parent',gui_fig,'Style','text', ...
    'FontSize',12,'FontName','Courier','HorizontalAlignment','left', ...
    'Units','normalized','BackgroundColor','w','Position',[0,0.9,1,0.1]);

% Control buttons
button_fontsize = 12;
view_button_position = [0,0.85,0.3,0.1];
clear controls_h
controls_h(1) = uicontrol('Parent',gui_fig,'Style','pushbutton','FontSize',button_fontsize, ...
    'Units','normalized','Position',view_button_position,'BackgroundColor','w', ...
    'String','Set ROI','Callback',{@set_roi,gui_fig});
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','togglebutton','FontSize',button_fontsize, ...
    'Units','normalized','Position',view_button_position,'BackgroundColor','w', ...
    'String','Manual record','Callback',{@cam_manual,gui_fig});
controls_h(end+1) = uicontrol('Parent',gui_fig,'Style','togglebutton','FontSize',button_fontsize, ...
    'Units','normalized','Position',view_button_position, ...
    'Callback',{@cam_listen,gui_fig},'value',true);
align(controls_h,'fixed',20,'middle');

% Start listener for experiment controller
update_status_text(status_text_h,'Connecting to experiment server');
try
    client_expcontroller = tcpclient("163.1.249.17",plab.locations.mousecam_port,'ConnectTimeout',2);
    configureCallback(client_expcontroller, "terminator", ...
        @(src,event,x) read_expcontroller_data(src,event,gui_fig));
    update_status_text(status_text_h,'Listening for start');
    set(controls_h(end),'String','Autorecord on','value',true,'BackgroundColor','g');
catch me
    % Error if no connection to experiment controller
    client_expcontroller = [];
    update_status_text(status_text_h,'Not connected to experiment server');
    set(controls_h(end),'String','Autorecord off','value',false,'BackgroundColor','r');
end

% Store gui data
gui_data.video_object = video_object;
gui_data.gui_fig = gui_fig;
gui_data.im_axes = im_axes;
gui_data.im_preview = im_preview;
gui_data.embedded_info_text = embedded_info_text;
gui_data.status_text_h = status_text_h;
gui_data.controls_h = controls_h;
gui_data.client_expcontroller = client_expcontroller;

% Update GUI data
guidata(gui_fig,gui_data);

end

%% Button functions

function set_roi(obj,eventdata,gui_fig)
% Draw and set ROI

% Get GUI data
gui_data = guidata(gui_fig);

% Reset ROI
gui_data.video_object.roi = [0,0,gui_data.video_object.VideoResolution];
axes(gui_data.im_axes); axis tight;

% Draw new ROI and want until finished
roi = drawrectangle(gui_data.im_axes);
wait(roi);

% Set new ROI position
roi_position = roi.Position;
if any(roi_position(2:3) == 0)
    roi_position = [0,0,gui_data.video_object.VideoResolution];
end
delete(roi);
gui_data.video_object.roi = roi_position;

% Set tight axes
axes(gui_data.im_axes); axis tight;

% Update GUI data
guidata(gui_fig,gui_data);

end

function cam_manual(obj,eventdata,gui_fig)

switch obj.Value
    case 1
        % Manual recording is turned on

        % Get gui data
        gui_data = guidata(gui_fig);

        % Change button display and disable other buttons
        obj.String = 'Stop';
        obj.BackgroundColor = [0.8,0,0];
        obj.ForegroundColor = 'w';
        set(gui_data.controls_h(gui_data.controls_h ~= obj),'Enable','off');

        % User choose mouse name
        animal = cell2mat(inputdlg('Mouse name'));
        if isempty(animal)
            % (if no mouse entered, do nothing)
            return
        end

        % Set save filename
        save_dir = plab.locations.local_data_path;
        rec_day = datestr(now,'YYYY-mm-DD');
        rec_time = datestr(now,'HHMM');
        save_path = plab.locations.filename('local', ...
            animal,rec_day,rec_time,'mousecam');

        % Make local data directory
        mkdir(save_path);

        % Update gui data
        guidata(gui_fig,gui_data);

        % Start DAQ acquisition
        cam_start(gui_fig,save_path);

    case 0
        % Manual recording is turned off

        % Stop recording
        cam_stop(gui_fig)

        % Get gui data
        gui_data = guidata(gui_fig);

        % Change button display and disable other buttons
        obj.String = 'Manual';
        obj.BackgroundColor = 'w';
        obj.ForegroundColor = 'k';
        set(gui_data.controls_h,'Enable','on');

end

end

function cam_listen(obj,eventdata,gui_fig)

% Get GUI data
gui_data = guidata(gui_fig);

switch obj.Value
    case 0
        % Listen off: clear TCP client
        obj.BackgroundColor = 'r';
        obj.String = 'Autorecord off';

        gui_data.client_expcontroller = [];
        update_status_text(gui_data.status_text_h,'Doing nothing');

    case 1
        % Listen on: connect TCP client
        obj.BackgroundColor = 'g';
        obj.String = 'Autorecord on';

        update_status_text(gui_data.status_text_h,'Connecting to exp controller...');
        gui_data.client_expcontroller = tcpclient("163.1.249.17",plab.locations.mousecam_port,'ConnectTimeout',2);
        configureCallback(gui_data.client_expcontroller, "terminator", ...
            @(src,event,x) read_expcontroller_data(src,event,gui_fig));
        update_status_text(gui_data.status_text_h,'Listening for start');
end

% Update GUI data
guidata(gui_fig,gui_data);

end

function cam_start(gui_fig,save_path)
% Start recording

% Get GUI data
gui_data = guidata(gui_fig);

% Put camera into trigger mode (pause acquisition)
src = getselectedsource(gui_data.video_object);
src.TriggerMode = 'on';

% Create videowriter
vidWriter = VideoWriter(fullfile(save_path,'mousecam.mj2'), 'Motion JPEG 2000');
vidWriter.CompressionRatio = 10;
gui_data.video_object.DiskLogger = vidWriter;

header_fn = fullfile(save_path,'mousecam_header.bin');
gui_data.header_fileID = fopen(header_fn,'w');

% Set header reader function
gui_data.video_object.FramesAcquiredFcnCount = 1;
gui_data.video_object.FramesAcquiredFcn = {@record_cam_header,gui_data.header_fileID};

% Delay to allow Timelite to start
update_status_text(gui_data.status_text_h,sprintf('Pausing to wait for Timelite...'));
pause(5);

% Start recording, turn on free-running mode
start(gui_data.video_object);
src.TriggerMode = 'off';

% Reset relative display info
gui_data.set_relative_info_flag = true;

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
update_status_text(gui_data.status_text_h,'Stopping recording');

% Put into trigger mode (stop acquisition), stop recording
pause(1); % pause 1s to allow flipper to stop
src = getselectedsource(gui_data.video_object);
src.TriggerMode = 'on';
stop(gui_data.video_object)

% Pause to allow Timelite to stop
% (no point stopping acquisiton: putting trigger mode on gives a few
% exposures that aren't saved, so the number of frames is still different)
update_status_text(gui_data.status_text_h,'Pausing to wait for Timelite...');
pause(4);

% Close header file
fclose(gui_data.header_fileID);

% Move data to server
update_status_text(gui_data.status_text_h,'Moving data to server...');
curr_data_path = get(gui_data.video_object.DiskLogger,'path');
move_data_to_server(curr_data_path,gui_data.status_text_h);

% Put free-running mode back on for preview
src.TriggerMode = 'off';

% Update status text
set(gui_fig,'color','w');
update_status_text(gui_data.status_text_h,'Listening for start');

end

%% Preview/record functions

function preview_embedded_info = preview_cam(obj,eventdata,himage)
% Custom preview function: output header information

% Get GUI fig (grab different - can't input to this function)
gui_fig = getappdata(himage,'gui_fig');

% Get GUI data
gui_data = guidata(gui_fig);

% Decode embedded pixel data
flipper_pin = 2; % GPIO input for flipper
header_pixels = eventdata.Data(1,1:40)';
mousecam_header = plab.mousecam.read_mousecam_header(header_pixels, flipper_pin);

% % Set relative timestamp/frame display
% if ~isfield(gui_data,'set_relative_info_flag')
%     gui_data.set_relative_info_flag = true;
% end
% if gui_data.set_relative_info_flag
%     % Store previous and current value
%     gui_data.relative_info.timestamp = [0,0];
%     gui_data.relative_info.frame_num = [0,0];
%     gui_data.set_relative_info_flag = false;
% else
%     % Camera clock resets every 127s: if timestamp goes backwards, set new
%     % relative time
% 
%     new_timestamp = circshift(gui_data.relative_info.timestamp,-1);
%     new_timestamp(2) = new_timestamp(2) + mousecam_header.timestamps;
%     if diff(new_timestamp) <0
%         keyboard
%     end
%     gui_data.relative_info.timestamp = new_timestamp;
%         
% end

% Print camera info
preview_embedded_info = ...
    sprintf('Framerate: %s, Timestamp: %.1f, Frame: %d, Flipper: %d', ...
    eventdata.FrameRate, ...
    mousecam_header.timestamps, ...
    mousecam_header.frame_num, ...
    mousecam_header.flipper);
set(gui_data.embedded_info_text,'String',preview_embedded_info);

% Update preview image
himage.CData = eventdata.Data;

% Update GUI data
guidata(gui_fig,gui_data);

end


function record_cam_header(video_object,vid_info,header_fileID)
% Pull headers from memory frames and save

if video_object.FramesAvailable > 0

    % Set embedded pixel index
    embedded_pixels_idx = 1:40; % 4 px * 10 items = 40 pixels

    % Grab all frames currently in memory
    curr_im = getdata(video_object, video_object.FramesAvailable);

    % Save raw embedded pixel values in bin file
    embedded_pixels = permute(curr_im(1,embedded_pixels_idx,1,:),[2,4,1,3]);
    fwrite(header_fileID,embedded_pixels);

end

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

    % Set local filename
    save_path = ...
        plab.locations.filename('local', ...
        rec_info.mouse,rec_info.date,rec_info.time,'mousecam');

    % Make local save directory
    mkdir(save_path)
    
    % Start camera recording acquisition
    cam_start(gui_fig,save_path)
end

end

function move_data_to_server(curr_data_path,status_text_h)
% Move data from local to server

% Check if the server is available
if ~exist(plab.locations.server_data_path,'dir')
    warning('Server not accessible at %s',plab.locations.server_data_path)
    return
end

% Move local data directories to server
curr_data_path_server = strrep(curr_data_path, ...
    plab.locations.local_data_path,plab.locations.server_data_path);
update_status_text(status_text_h,'Copying to server')
[status,message] = movefile(curr_data_path,curr_data_path_server);
if ~status
    update_status_text(status_text_h,'Last server copy failed! Listening for start');
    warning('Failed copying to server: %s',message);
else
    update_status_text(status_text_h,'Listening for start');
end

% Delete empty local folders
% (3 hierarchy levels: protocol > day > animal)
try
    curr_hierarchy_path = fileparts(curr_data_path);
    for hierarchy_levels = 1:3
        hierarchy_dir = dir(curr_hierarchy_path);
        if all(contains({hierarchy_dir.name},'.'))
            rmdir(curr_hierarchy_path)
            % Move up one step in hierarchy
            curr_hierarchy_path = fileparts(curr_hierarchy_path);
        end
    end
end

end






