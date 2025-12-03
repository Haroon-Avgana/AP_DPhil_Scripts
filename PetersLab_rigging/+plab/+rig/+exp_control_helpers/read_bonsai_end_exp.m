function read_bonsai_end_exp(server_stim, ~, server_timelite, server_mousecam, server_widefield, StartButton)
    disp('Data received')
    server_stim.UserData = readline(server_stim);
    disp(server_stim.UserData)

    if strfind(server_stim.UserData, 'done')
        % send stop to timelitem mousecam, widefield
        if server_timelite.Connected
            writeline(server_timelite, 'stop');
        end
        if server_mousecam.Connected
            writeline(server_mousecam, 'stop');
        end
        if server_widefield.Connected
            writeline(server_widefield, 'stop');
        end

        % change button
        StartButton.Text = 'Start';
        StartButton.Value = 0;
    end
end

       