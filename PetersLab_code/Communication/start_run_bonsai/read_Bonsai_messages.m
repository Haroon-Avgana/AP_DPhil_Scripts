cd('C:\Users\peterslab\Documents\start_run_stop_bonsai\java_stuff')

% Load the jar file
javaaddpath('javaosctomatlab.jar');

% To initialize:
javaaddpath(schedule.javaoscpath);
import com.illposed.osc.*;
import java.lang.String;

% To read:
oscport= 20000;
oscreceiver = OSCPortIn(oscport);
osclistener = MatlabOSCListener();
oscreceiver.addListener(String('/msg'),osclistener);
oscreceiver.startListening();
getlastmessage=osclistener.getMessageArgumentsAsDouble();
% To close the port
oscreceiver.stopListening();
oscreceiver.close();
oscreceiver=[];
osclistener=[];