%% Testing scale

% Get COM port
serial_available = serialportlist('available');
if length(serial_available) == 1
    % If one port, select that
    scale_com_port = serial_available{1};
elseif length(serial_available) > 1
    % If multiple ports, prompt selection
    idx = listdlg('ListString',serial_available);
    scale_com_port = serial_available{idx};
elseif isempty(serial_available)
    error('Weighing scale -- no available serial ports detected')
end

baud_rate = 9600;
scale_serial_port = serialport(scale_com_port,baud_rate);
configureTerminator(scale_serial_port,"CR")

write(scale_serial_port,"O","single")
writeline(scale_serial_port,'4FH');









