%% setup servers
server_stim = tcpserver("0.0.0.0",plab.locations.bonsai_port);
server_timelite = tcpserver("0.0.0.0",plab.locations.timelite_port);
server_mousecam = tcpserver("0.0.0.0",plab.locations.mousecam_port);

%% send mouse and protocol info
protocol = 'trial_test';
name = 'MissMouse0801';
date_today = string(datetime('today', 'Format', 'yyyy_MM_dd'));

waitfor(server_stim, 'Connected', 1)
s = struct('mouse', name, 'date', date_today, 'protocol', protocol);
encode_s = jsonencode(s);
write(server_stim,encode_s,'string');

%% send start message
write(server_timelite,  'start', 'string');
write(server_mousecam,  'start', 'string');

%% send stop message
write(server_stim,'stop','string');
write(server_timelite, 'stop', 'string');
write(server_mousecam, 'stop', 'string');