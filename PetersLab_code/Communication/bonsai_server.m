function bonsai_server

%% Initialize bonsai server

% Create figure
bonsai_server_fig = ...
    uifigure('Units','normalized','Position',[0,0,1,1],'color',[0.5,0.5,0.5], ...
    'toolbar','none','menubar','none','CloseRequestFcn',@close_bonsai_server);

% Set up structure for listeners and receivers
communication_handles = struct;

% Open TCP client for MC
communication_handles.client_mc = tcpclient("163.1.249.17",50001);

% Open UDP connection for Bonsai
communication_handles.u_bonsai = udpport("IPV4");

% Setup listening to Bonsai workflow
% (load the jar file)
java_path = 'C:\Users\peterslab\Documents\start_run_stop_bonsai\java_stuff';
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

% Start listener to run Bonsai workflows
root_save = 'C:\Users\peterslab\Documents';
root_workflows = 'C:\Users\peterslab\Documents\GitHub\PetersLab_code\Bonsai stuff';
bonsai_controller(root_save, root_workflows, communication_handles)

end

function bonsai_controller(root_save, root_workflows, communication_handles)

% Loop to listen and run Bonsai workflow
while(true)
    % wait for file info from mc
    waitfor(communication_handles.client_mc, 'NumBytesAvailable')
    exp_data = read(communication_handles.client_mc, communication_handles.client_mc.NumBytesAvailable,"string");

    % decode json
    data_struct = jsondecode(exp_data);

    % make mouse dir
    this_save_path = fullfile(root_save,data_struct.mouse);
    mkdir(this_save_path);

    % get paths for files
    workflowpath = fullfile(root_workflows, [data_struct.protocol '.bonsai']);
    local_worfkflowpath = fullfile(this_save_path, [data_struct.protocol '.bonsai']);
    filename = fullfile(this_save_path, 'test.csv');

    % copy bonsai workflow in new folder
    copyfile(workflowpath, local_worfkflowpath);

    % start bonsai
    runBonsaiWorkflow(local_worfkflowpath, {'FileName', filename}, [], 1)

    % get Stop message
    msg_stop = [];
    while ~strcmp(msg_stop,'stop')
        waitfor(communication_handles.client_mc, 'NumBytesAvailable')
        msg_stop = read(communication_handles.client_mc, communication_handles.client_mc.NumBytesAvailable,"string");
    end

    % send stop to bonsai
    bonsai_oscsend(u_bonsai,'/stop',"localhost",30000,'i',45);

    % get message from Bonsai stopping
    getlastmessage = communcation_handles.osclistener.getMessageArgumentsAsDouble();
    while(isempty(getlastmessage))
        getlastmessage = communcation_handles.osclistener.getMessageArgumentsAsDouble();
    end

    % close Bonsai
    system('taskkill /F /IM Bonsai.EXE');
    system('taskkill /F /IM OpenConsole.EXE');
end

end

function close_bonsai_server(obj, event)

% Confirm close
confirm_close = uiconfirm(obj,'Close Bonsai server?','Confirm close');
if strcmp(confirm_close,'OK')

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

end



