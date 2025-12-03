function bonsai_server

% AP NOTE 2024-01-18: made figure invisisble (to ensure no chance of
% overlaying bonsai image).
%
% Streamline this in future by making bonsai_server an object rather than a
% GUI

%% Initialize bonsai server

% Create figure on iPad screens
monitors_pos = get(0, 'MonitorPosition');
% (determine ipad screen by x-width)
ipad_screens_idx = monitors_pos(:,3) == 3840;
ipad_screens = monitors_pos(ipad_screens_idx,:);
bonsai_server_fig = ...
    uifigure('Position',ipad_screens,'color','#828282', ...
    'toolbar','none','menubar','none','CloseRequestFcn',@close_bonsai_server, ...
    'visible','off');

% Set up structure for listeners and receivers
communication_handles = struct;

% Open TCP client for MC
communication_handles.client_mc = tcpclient("163.1.249.17",plab.locations.bonsai_port);
configureCallback(communication_handles.client_mc, "terminator", @(src, event, x) readData (src, event, bonsai_server_fig));

% Open UDP connection for Bonsai
communication_handles.u_bonsai = udpport("IPV4");

% Setup listening to Bonsai workflow
% (load the jar file)
java_path = fullfile(fileparts(which('plab.rig.bonsai_server')),'+bonsai_server_helpers');
javaaddpath(fullfile(java_path,'javaosctomatlab.jar'));

% (import java packages)
import com.illposed.osc.*;
import java.lang.String;

% (set up OSC listener and receiver)
oscport= 20000;
communication_handles.oscreceiver = OSCPortIn(oscport);
communication_handles.osclistener = MatlabOSCListener();
communication_handles.oscreceiver.addListener(String('/stop'),communication_handles.osclistener);
communication_handles.oscreceiver.startListening();

% Upload receiver/listener to figure
guidata(bonsai_server_fig,communication_handles);

% Setup arduino for water rewards
% setup_arduino(bonsai_server_fig);

end

function setup_arduino(bonsai_server_fig)
% connect to arduino and calculate amount of time valve open for reward amount
    
    communication_handles = guidata(bonsai_server_fig);
    
    cal_value = readmatrix("C:\Water_calibration\calibration.csv");
    reward_amount = 6; 
    time_valve_open = num2str(cal_value * reward_amount);

    communication_handles.arduino_device = serialport("COM10", 250000);
    configureTerminator(communication_handles.arduino_device,"CR");
    set(bonsai_server_fig,'KeyPressFcn',@give_water);

    communication_handles.time_valve_open = time_valve_open;

    guidata(bonsai_server_fig,communication_handles);
end

function give_water(bonsai_server_fig,event)
    if strcmp(event.Key, 'w')
        communication_handles = guidata(bonsai_server_fig);
        writeline(communication_handles.arduino_device, communication_handles.time_valve_open);
    end
end

function readData(client, ~, bonsai_server_fig)
    disp('Data received')
    client.UserData = readline(client);
    if strcmp(client.UserData, 'stop')
        % send stop to bonsai
        communication_handles = guidata(bonsai_server_fig);
        plab.rig.bonsai_server_helpers.bonsai_oscsend(communication_handles.u_bonsai,'/stop',"localhost",30000,'s','stop');
    else
        run_bonsai(bonsai_server_fig);
    end
end


function run_bonsai(bonsai_server_fig)

    communication_handles = guidata(bonsai_server_fig);

    % clear arduino device
%     communication_handles = rmfield(communication_handles,'arduino_device');

    % decode json
    data_struct = jsondecode(communication_handles.client_mc.UserData);
    
    % Set local filename for bonsai workflow
    [~, bonsai_folder] = fileparts(data_struct.protocol_path);

    local_worfkflow_path = ...
        plab.locations.filename('local', ...
        data_struct.mouse,data_struct.date,data_struct.time, ...
        'bonsai',bonsai_folder);
    [save_path,~,~] = fileparts(local_worfkflow_path);

    % Make local save directory
    mkdir(fileparts(save_path));

    % get paths for files
    workflowpath = fullfile(plab.locations.local_workflow_path, data_struct.protocol_path);
    local_worfkflow_file = fullfile(local_worfkflow_path, data_struct.protocol_name);

    % copy bonsai workflow in new folder
    copyfile(workflowpath, local_worfkflow_path);

    % start bonsai
    plab.rig.bonsai_server_helpers.runBonsaiWorkflow(local_worfkflow_file, {'SavePath', save_path}, [], 1);

    % Update save path into GUI data
    communication_handles.save_path = save_path;
    guidata(bonsai_server_fig,communication_handles);

    % Start timer function to listen for "stopped" message
    bonsai_timer_fcn = timer('TimerFcn', ...
    {@get_bonsai_message,bonsai_server_fig}, ...
    'Period', 1/10, 'ExecutionMode','fixedDelay', ...
    'TasksToExecute', inf);
    start(bonsai_timer_fcn)


end

function get_bonsai_message(obj, ~, bonsai_server_fig)

    communication_handles = guidata(bonsai_server_fig);

    getlastmessage = communication_handles.osclistener.getMessageArgumentsAsString();    
    if ~isempty(getlastmessage)
        disp('Read something')
        disp(getlastmessage)
        % close Bonsai
        pause(4); % Pause to allow Bonsai to cleanly finish
        system('taskkill /F /IM Bonsai.EXE');
        system('taskkill /F /IM OpenConsole.EXE');
        % send done to exp control
        writeline(communication_handles.client_mc, 'done');
        % delete timer
        stop(obj)
        delete(obj)

        % Move data to server
        move_data_to_server(communication_handles.save_path);

        % Setup arduino for water rewards
%         setup_arduino(bonsai_server_fig);
    end
end

function move_data_to_server(curr_data_path)
% Move data from local to server

% Check if the server is available
if ~exist(plab.locations.server_data_path,'dir')
    warning('Server not accessible at %s',plab.locations.server_data_path)
    return
end

% Move local data directories to server
curr_data_path_server = strrep(curr_data_path, ...
    plab.locations.local_data_path,plab.locations.server_data_path);
[status,message] = movefile(curr_data_path,curr_data_path_server);
if ~status
    warning('Failed copying to server: %s',message);
else
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

function close_bonsai_server(obj, ~)

% Get OSC handler and stop listeners
communication_handles = guidata(obj);
communication_handles.oscreceiver.stopListening();
communication_handles.oscreceiver.close();

% Delete the figure
delete(obj);

% Clear TCP client
delete(communication_handles.client_mc)

% Exit matlab (stops the waitfor and resets ports)
exit

end



