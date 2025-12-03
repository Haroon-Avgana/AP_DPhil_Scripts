function timelite
% Light-weight NI-DAQ acquisition GUI (lite version of Timeline)
%
% FOR RIG-SPECIFIC DAQ SETUP: 
% Make local custom version of <plab.local_rig.timelite_config>
%
% Requires: 
% data acquisition toolbox
% instrument control toolbox

%% Make GUI

% Create figure at left edge of the screen
gui_fig = figure('Name','Timelite','Units','Normalized', ...
    'Position',[0,0.05,0.15,0.9],'color','w','menu','none');

% Initialize text
uicontrol('Parent',gui_fig,'Style','text', ...
    'String','TIMELITE','FontSize',20, ...
    'HorizontalAlignment','center', ...
    'Units','normalized','BackgroundColor',[0.8,0.8,0.8],'Position',[0,0.95,1,0.05]);

status_text_h = uicontrol('Parent',gui_fig,'Style','text', ...
    'FontSize',12,'FontName','Courier','HorizontalAlignment','left', ...
    'Units','normalized','BackgroundColor','w','Position',[0,0.1,1,0.85]);

% Initialize buttons
user_button_h = uicontrol('Parent',gui_fig,'Style','pushbutton', ...
    'String','PREVIEW','Callback',{@user_buttonpress,gui_fig}, ...
    'FontSize',16,'BackgroundColor',[0,0.8,0],'ForegroundColor','w', ...
    'Units','normalized','Position',[0,0,1,0.1],'Enable','off');

% Draw initial controls (so user knows startup is in progress)
drawnow;

% Set up DAQ
update_status_text(status_text_h,'Setting up DAQ');
try
    % Reset DAQ (and turn off on-demand channel warning)
    daqreset;
    warning('off','daq:Session:onDemandOnlyChannelsAdded');
    % Set up DAQ according to local config file
    daq_device = plab.local_rig.timelite_config;
    daq_device.analog.ScansAvailableFcn = @(src,evt,x) daq_upload(src,evt,gui_fig);
catch me
    % Error if DAQ could not be configured
    error(me.identifier,'Timelite -- DAQ device could not be configured: \n %s',me.message)
end

% Start listener for experiment controller
update_status_text(status_text_h,'Connecting to experiment server');
try
client_expcontroller = tcpclient("163.1.249.17",plab.locations.timelite_port,'ConnectTimeout',2);
configureCallback(client_expcontroller, "terminator", ...
    @(src,event,x) read_expcontroller_data(src,event,gui_fig));
catch me
    % Error if no connection to experiment controller
    update_status_text(status_text_h,'Error connecting to experiment server');
    warning(me.identifier,'Timelite -- Cannot connect to experiment controller: \n %s',me.message)
end

% Write text
status_text = [ ...
    {sprintf('Status: \n')}
    {sprintf('Sample rate: %d',daq_device.analog.Rate)}
    {sprintf('Samples/upload: %d \n',daq_device.analog.ScansAvailableFcnCount)}
    {'Analog input channels: '}
    join([{daq_device.analog.Channels.ID}',{daq_device.analog.Channels.Name}',],' | ')
    {sprintf('\nDigital output channels: ')}
    join([{daq_device.digital.Channels.ID}',{daq_device.digital.Channels.Name}',],' | ')];
set(status_text_h,'String',status_text)

% Initialize and upload gui data
gui_data = struct;
gui_data.status_text_h = status_text_h;
gui_data.user_button_h = user_button_h;
gui_data.daq_device = daq_device;
if exist('client_expcontroller','var')
    gui_data.client_expcontroller = client_expcontroller;
end

% Update status
update_status_text(gui_data.status_text_h,'Listening for start');

% Update gui data
guidata(gui_fig,gui_data);

% Enable button
set(gui_data.user_button_h,'Enable','on');

end

function user_buttonpress(h,eventdata,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

if gui_data.daq_device.analog.Running
    % If DAQ is running: stop recording (with confirmation dialog)
    confirm = questdlg('\fontsize{14}Stop timelite?','Confirm timelite stop', ...
        'Yes','No',struct('Default','No','Interpreter','tex'));
    if strcmp(confirm,'Yes')
        daq_stop(gui_fig)
    end
else
    % If DAQ isn't running: start DAQ with no saving (preview mode)
    daq_start(gui_fig,[])
end

end

function daq_start(gui_fig,save_filename)

% Get gui data
gui_data = guidata(gui_fig);

if ~isempty(save_filename)
    % If save filename: 

    % Create mat file (to convert at the end)
    % (daq information)
    daq_info = struct( ...
        'rate',gui_data.daq_device.analog.Rate, ...
        'device',gui_data.daq_device.analog.Channels(1).Device.Model, ...
        'type',{gui_data.daq_device.analog.Channels.Type}, ...
        'channel',{gui_data.daq_device.analog.Channels.ID}, ...
        'measurement_type',{gui_data.daq_device.analog.Channels.MeasurementType}, ...
        'channel_name',{gui_data.daq_device.analog.Channels.Name});
    % (daq data - to be filled during streaming)
    [timestamps,data] = deal([]);
    % (save initial variables and keep open for adding at the end)
    save(save_filename,'daq_info','timestamps','data','-v7.3');
    gui_data.save_file_mat = matfile(save_filename,'Writable',true);

    % Create binary file (for streaming)
    [save_path,save_name] = fileparts(save_filename);
    save_bin_filename = fullfile(save_path,sprintf('%s.bin',save_name));
    gui_data.save_file_bin = fopen(save_bin_filename,'w+');

    % Update status
    update_status_text(gui_data.status_text_h,sprintf('Recording: %s',save_path));
    % Update gui data
    guidata(gui_fig,gui_data);
else
    % If no save filename, empty (preview mode - no recording)
    gui_data.save_file_mat = [];
    % Update status
    update_status_text(gui_data.status_text_h,'Previewing');
    % Update gui data
    guidata(gui_fig,gui_data);
end

% Start DAQ input acquisition, set outputs to HIGH
start(gui_data.daq_device.analog,'continuous');
pause(2); % ensure start before outputs high
write(gui_data.daq_device.digital,true);

% Create live plot (unclosable)
gui_fig_position = gui_fig.Position;
live_fig_position = ...
    [gui_fig_position(1)+gui_fig_position(3),gui_fig_position(2), ...
    1-(gui_fig_position(1)+gui_fig_position(3)),gui_fig_position(4)];

live_plot_fig = figure('CloseRequestFcn',[],'color','w', ...
    'Units','normalized','Position',live_fig_position, ...
    'Name','Timelite live plot','menu','none');
gui_data.live_plot_fig = live_plot_fig;

% Enable stop button
set(gui_data.user_button_h,'String','MANUAL/EMERGENCY STOP','BackgroundColor',[0.8,0,0]);

% Update gui data
guidata(gui_fig,gui_data);

end

function daq_stop(gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Stop DAQ input acquisition, set outputs to LOW
update_status_text(gui_data.status_text_h,'Stopping (final 4s)');
write(gui_data.daq_device.digital,false); 
pause(4); % ensure outputs low and other GUIs finished before stopping
stop(gui_data.daq_device.analog)

% Upload remaining data (and get updated gui data)
if gui_data.daq_device.analog.NumScansAvailable > 0
    daq_upload(gui_data.daq_device.analog,[],gui_fig);
    gui_data = guidata(gui_fig);
end

% Delete live plot
delete(gui_data.live_plot_fig);

% Saved file handling (if recorded, not on preview)
if ~isempty(gui_data.save_file_mat)

    % Move binary streamed data into mat file
    update_status_text(gui_data.status_text_h,'Converting data binary to mat');
    drawnow;

    % --> Rewind and read data (reshaped into t x [timestamps,chan])
    n_channels = length(gui_data.daq_device.analog.Channels);

    frewind(gui_data.save_file_bin);
    recorded_daq_data = reshape(fread(gui_data.save_file_bin,'single'),n_channels+1,[])';

    % --> Write data to mat file
    gui_data.save_file_mat.timestamps = recorded_daq_data(:,1);
    gui_data.save_file_mat.data = recorded_daq_data(:,2:end);

    % --> Close and delete temporary bin file
    bin_filename = fopen(gui_data.save_file_bin);
    fclose(gui_data.save_file_bin);
    delete(bin_filename)

    % Move data to server
    move_data_to_server(gui_data.save_file_mat.Properties.Source,gui_data.status_text_h);

end

% Update gui data
guidata(gui_fig,gui_data);

% Reset status and preview button
update_status_text(gui_data.status_text_h,'Listening for start');
set(gui_data.user_button_h,'String','PREVIEW','BackgroundColor',[0,0.8,0]);

end

function daq_upload(obj,event,gui_fig)

% Get gui data
gui_data = guidata(gui_fig);

% Read all buffered data
[daq_data,daq_timestamps] = ...
    read(obj,'all','OutputFormat','Matrix');

% If data is empty, return
% (happens if overlapping function calls so nothing left in buffer)
if isempty(daq_data)
    return
end

% Counter data: convert from unsigned to signed integer type
% (allow for negative values, rather than overflow)
% (assumes 32-bit counter)
position_channels_idx = strcmp({obj.Channels.MeasurementType},'Position');
daq_data(:,position_channels_idx) = ...
    double(typecast(uint32(daq_data(:,position_channels_idx)),'int32'));

% If save file exists, write to binary file
if isfield(gui_data,'save_file_mat') && ~isempty(gui_data.save_file_mat)
    % Format: [timestamp 1, ch1 t1, ch2 t1...timestamp N, chN tN])
    % Precision: single (timestamps lose ~4us precision, 32-bit counter is
    % preserved, 16-bit AI is preserved but 2x space)
    bin_write_data = reshape(single([daq_timestamps,daq_data])',[],1);
    fwrite(gui_data.save_file_bin,bin_write_data,'single');
end

% Plot data
daq_plot(obj,gui_data,daq_data,daq_timestamps,gui_fig);

end

function daq_plot(obj,gui_data,daq_data,daq_timestamps,gui_fig)

plot_data_t = 5; % seconds of data to plot

if isfield(gui_data,'live_plot_fig') && isvalid(gui_data.live_plot_fig)

    if ~isfield(gui_data,'live_plot_traces') || ~any(isgraphics(gui_data.live_plot_traces))
        % If nothing plotted, create traces and plot data at end
        blank_data = zeros(plot_data_t*obj.Rate,size(daq_data,2));
        gui_data.live_plot_traces = tiledlayout(gui_data.live_plot_fig,size(daq_data,2),1,'TileSpacing','none');

        for curr_channel = 1:size(daq_data,2)
            curr_axes = nexttile(gui_data.live_plot_traces,curr_channel);
            plot(curr_axes,(0:size(blank_data,1)-1)/obj.Rate,blank_data(:,curr_channel));
            ylabel(curr_axes,strrep(gui_data.daq_device.analog.Channels(curr_channel).Name,'_',' '), ...
                'Rotation',0,'VerticalAlignment','middle','HorizontalAlignment','right');
            curr_axes.Box = 'off';
            if curr_channel < size(daq_data,2)
                curr_axes.XAxis.Visible = 'off';
            else
                xlabel(curr_axes,'Time (s)');
            end

        end

        % Update gui data
        guidata(gui_fig,gui_data);
    end

    % For each channel: shift off old data, draw in new data
    live_plot_axes = flipud(get(gui_data.live_plot_traces,'Children'));

    old_t = get(get(live_plot_axes(1),'Children'),'XData');
    new_t = horzcat(old_t(length(daq_timestamps)+1:end), ...
            reshape(daq_timestamps,1,[]));

    for curr_channel = 1:size(daq_data,2)
        old_plot_data = get(get(live_plot_axes(curr_channel),'Children'),'YData');
        new_plot_data = horzcat(old_plot_data(size(daq_data,1)+1:end), ...
            daq_data(:,curr_channel)');
        set(get(live_plot_axes(curr_channel),'Children'), ...
            'XData',new_t,'YData',new_plot_data);
    end

    %  Update gui data
    guidata(gui_fig,gui_data);

end

end

function read_expcontroller_data(client,event,gui_fig)
% Read message from experiment controller

% Get gui data
gui_data = guidata(gui_fig);

% Get message from experiment controller
expcontroller_message = readline(client);

if strcmp(expcontroller_message, 'stop')
    % If experiment controller sends stop, stop DAQ acquisition
    daq_stop(gui_fig);
else
    % Otherwise, assume message is information to start protocol    
    rec_info = jsondecode(expcontroller_message);

    % Set local filename
    save_filename = ...
        plab.locations.filename('local', ...
        rec_info.mouse,rec_info.date,rec_info.time,'timelite.mat');

    % Make local save directory
    if ~exist(fileparts(save_filename),'dir')
        mkdir(fileparts(save_filename))
    end
    
    % Start DAQ acquisition
    daq_start(gui_fig,save_filename)
end

end

function update_status_text(status_text_h,status)
% Update status text

curr_text = get(status_text_h,'String');
new_text = [{sprintf('Status: %s',status)};curr_text(2:end)];
set(status_text_h,'String',new_text);

end

function move_data_to_server(curr_data_filename,status_text_h)
% Move data from local to server

% Check if the server is available
if ~exist(plab.locations.server_data_path,'dir')
    warning('Server not accessible at %s',plab.locations.server_data_path)
    return
end

% Move local data directories to server
curr_data_filename_server = strrep(curr_data_filename, ...
    plab.locations.local_data_path,plab.locations.server_data_path);
update_status_text(status_text_h,'Copying to server')
[status,message] = movefile(curr_data_filename,curr_data_filename_server);
if ~status
    update_status_text(status_text_h,'Last server copy failed! Listening for start');
    warning('Failed copying to server: %s',message);
else
    update_status_text(status_text_h,'Listening for start');
end

% Delete empty local folders
% (3 hierarchy levels: protocol > day > animal)
try
    curr_hierarchy_path = fileparts(curr_data_filename);
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







