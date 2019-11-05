function csc2kilosort2(Input_folder,Temp_folder)
Num_EperBundle = 4; % Number of electrode/channels per bunddle (here tetrode)
Num_channels = 16; % Number of channels including inactive ones
Num_tetrodes = Num_channels/Num_EperBundle;

CSCFiles = dir(fullfile(Input_folder, '*CSC*.mat'));
Active_channels = nan(length(CSCFiles),1);
% Find out the active channels
for ff=1:length(CSCFiles)
    Istart = strfind(CSCFiles(ff).name, 'CSC') + 3;
    Iend = strfind(CSCFiles(ff).name, '.mat') -1;
    Active_channels(ff) = str2double(CSCFiles(ff).name(Istart:Iend));
end


for tt=1:Num_tetrodes
    % select the active channels for that bundle
    Active_channelsfiles = find((Active_channels <(tt*Num_EperBundle)) .* (Active_channels >=((tt-1)*Num_EperBundle)));
    DATOUT = cell(length(Active_channelsfiles),1);
    for cc=1:length(Active_channelsfiles)
        load(fullfile(CSCFiles(Active_channelsfiles(cc)).folder, CSCFiles(Active_channelsfiles(cc)).name),'AD_count_int16', 'AD_count_to_uV_factor','Bat_id','Date')
        DATOUT{cc} = single(AD_count_int16).*single(AD_count_to_uV_factor);
    end
    clear AD_cound_int16
    DATOUT = cell2mat(DATOUT);
    DATOUT  =  int16(DATOUT);
    whos DATOUT
    fid = fopen(fullfile(Temp_folder, sprintf('%s_%s_TempTetrode%d.bin',Bat_id,Date,tt)),'w');
    fwrite(fid, DATOUT, 'int16');
    fclose(fid);
end
clear DATOUT
