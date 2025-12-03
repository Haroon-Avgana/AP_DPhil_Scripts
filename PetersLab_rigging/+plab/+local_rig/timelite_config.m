function daq_device = timelite_config
% Connect and configure DAQ for Timelite recording
% This is a template: modify for each rig and keep in rig local code

%% Set DAQ parameters 

% Set sample rate for DAQ
daq_sample_rate = 1000;

% Set DAQ buffer time (uploads and processes every n seconds)
daq_buffer_time = 1;

% Set daq index to use
use_daq_idx = 1;

%% Connect to DAQ

% Find and connect to DAQ, create daq object structure
daqs_available = daqlist;
daq_device = struct;

% Set up DAQ analog input
daq_device.analog = daq(daqs_available.VendorID(use_daq_idx));
daq_device.analog.Rate = daq_sample_rate;
daq_device.analog.ScansAvailableFcnCount = daq_buffer_time*daq_sample_rate;

% Set up DAQ digital output
% (needs separate object: on-demand instead of clocked)
daq_device.digital = daq(daqs_available.VendorID(use_daq_idx));

%% Configure analog inputs

% Analog inputs
ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ctr0','Position');
ch.EncoderType = 'X4';
ch.Name = 'wheel';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai0','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'flipper';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai1','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'widefield_camera';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai2','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'mouse_camera';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai3','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'photodiode';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai4','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'stim_screen';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai5','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'reward_valve';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai11','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'dummy_valve';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai6','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'violet_light';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai14','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'blue_light';

ch = addinput(daq_device.analog,daqs_available.DeviceID(use_daq_idx),'ai7','Voltage');
ch.TerminalConfig = 'SingleEnded';
ch.Name = 'lick_spout';



%% Configure digital outputs

ch = addoutput(daq_device.digital,daqs_available.DeviceID(use_daq_idx),'port1/line0','Digital');
ch.Name = 'acqLive';















