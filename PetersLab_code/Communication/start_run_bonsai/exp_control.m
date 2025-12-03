%% init stuff
save_path = 'C:\Users\peterslab\Documents';
all_workflows = 'C:\Users\peterslab\Documents\GitHub\PetersLab_code\Bonsai stuff';

%% setup listening to Bonsai workflow
cd('C:\Users\peterslab\Documents\start_run_stop_bonsai\java_stuff')

% Load the jar file
javaaddpath('javaosctomatlab.jar');

% go back to code folder
cd('C:\Users\peterslab\Documents\start_run_stop_bonsai')

% To initialize:
% javaaddpath(schedule.javaoscpath);
import com.illposed.osc.*;
import java.lang.String;

% To read:
oscport= 20000;
oscreceiver = OSCPortIn(oscport);
osclistener = MatlabOSCListener();
oscreceiver.addListener(String('/stop'),osclistener);
oscreceiver.startListening();

%% setup connection to mc and wait for file info
client = tcpclient("163.1.249.17",10000);
exp_data = read(client, client.NumBytesAvailable,"string");
while isempty(exp_data)
    exp_data = read(client, client.NumBytesAvailable,"string");
end
data_struct = jsondecode(exp_data);

% make mouse dir 
this_save_path = fullfile(save_path,data_struct.mouse);
mkdir(this_save_path);

% get paths for files
workflowpath = fullfile(all_workflows, [data_struct.protocol '.bonsai']);
filename = fullfile(this_save_path, 'test.csv');

%% run Bonsai workflow when get start message
msg = read(client, client.NumBytesAvailable,"string");
while isempty(msg) || msg ~='start'
    msg = read(client, client.NumBytesAvailable,"string");
end
local_runBonsaiWorkflow(workflowpath, {'FileName', filename}, [], 1)

%% stop Bonsai workflow when get stop message
msg = read(client, client.NumBytesAvailable,"string");
while isempty(msg) || msg ~='stop'
    msg = read(client, client.NumBytesAvailable,"string");
end

%% send stop message
u = udpport("IPV4");
bonsai_oscsend(u,'/stop',"localhost",50004,'i',45);

%% get message from Bonsai stopping
getlastmessage=osclistener.getMessageArgumentsAsDouble();
while(isempty(getlastmessage))
    getlastmessage=osclistener.getMessageArgumentsAsDouble();
end

% To close the port
oscreceiver.stopListening();
oscreceiver.close();
oscreceiver=[];
osclistener=[];

%% close Bonsai
system('taskkill /F /IM Bonsai.EXE');  
system('taskkill /F /IM OpenConsole.EXE');  