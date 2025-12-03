function AM_save_h_file(in_wave, wave_name)
% Function takes a wave and a wave name
% Can give existing opt of square and trapezoid without defining the wave outside this

    if wave_name=='sq_wave'
        wave = [ones(1, 256)*4095, zeros(1, 256)];
    elseif wave_name=='trp_wave'
        wave = [linspace(0,4095,64), 4095*ones(1,128), linspace(4095,0,64), zeros(1, 256)];
    else
        wave = in_wave;
    end

    hex_wave = dec2hex(wave,3);
    
    div_hex_wave = reshape(hex_wave',3, 8, size(hex_wave,1)/8);
    div_hex_wave = permute(div_hex_wave,[2, 1, 3]);
    
    % add name of waveform in the begining of the file
    waveform_title = ['_Waveforms_' wave_name '_h_'];
    beg_file = ['#ifndef ' waveform_title '\n#define ' waveform_title '\n\n#define maxSamplesNum 512\n\nstatic int waveformsTable_' wave_name '[maxSamplesNum] = {\n'];
    end_file = '};\n#endif';
    
    fn = [wave_name '.h'];
    fid = fopen(fn,'w');
    fprintf(fid, beg_file);
    for rows=1:size(div_hex_wave, 3)
        for cols=1:size(div_hex_wave, 1)
            if cols==size(div_hex_wave, 1) && rows==size(div_hex_wave, 3)
                fprintf(fid,'0x%s ', div_hex_wave(cols,:,rows));
            else
                fprintf(fid,'0x%s , ', div_hex_wave(cols,:,rows));
            end
        end
        fprintf(fid,'\n');
    end
    fprintf(fid, end_file);
    fclose(fid);

end