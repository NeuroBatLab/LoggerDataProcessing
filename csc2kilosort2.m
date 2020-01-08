function csc2kilosort2(Input_folder,Output_folder,PerTetrode,Kilosort2Platform)
% Input_folder: the folder containing the CSC files
% Output_folder: a folder where the output file will be saved (please
% prefer local!! not on the server)
% PerTetrode: logical that indicates whether the data should be saved for
% each tetrode (1) or all tetrodes (0)
% Kilosort2Platform: platform for which the file should be prepared
% (platform on which kilosort2 will be run): 'win' or 'unix'

Num_EperBundle = 4; % Number of electrode/channels per bunddle (here tetrode)
Num_channels = 16; % Number of channels including inactive ones
Num_tetrodes = Num_channels/Num_EperBundle;
if nargin<3
    PerTetrode=0;
end
if nargin<4
    Kilosort2Platform = 'win'; % can also be set to 'unix'
end

if ~exist(Output_folder, 'dir')
    mkdir(Output_folder)
end
% Set the path to a working directory on the computer so logger data are
% transfered there and directly accessible for calculations
HerePath = pwd;
WorkDir = fullfile(HerePath(1:(strfind(HerePath,'Documents')+9)), 'WorkingDirectory');

% Bring data back on the computer
fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', Input_folder, WorkDir);
mkdir(WorkDir)
[s,m,e]=copyfile(fullfile(Input_folder, '*CSC*.mat'), WorkDir, 'f');
if ~s
    m %#ok<NOPRT>
    e %#ok<NOPRT>
    error('File transfer did not occur correctly for %s\n', Input_folder);
end

CSCFiles = dir(fullfile(WorkDir, '*CSC*.mat'));
Active_channels = nan(length(CSCFiles),1);
% Find out the active channels
for ff=1:length(CSCFiles)
    Istart = strfind(CSCFiles(ff).name, 'CSC') + 3;
    Iend = strfind(CSCFiles(ff).name, '.mat') -1;
    Active_channels(ff) = str2double(CSCFiles(ff).name(Istart:Iend));
end
Active_channels_sorted = sort(Active_channels);

if PerTetrode
    for tt=1:Num_tetrodes
        fprintf(1,'Creating tetrode file %d/%d\n', tt, Num_tetrodes);
        % select the active channels for that bundle
        Active_channelsfiles = find((Active_channels <(tt*Num_EperBundle)) .* (Active_channels >=((tt-1)*Num_EperBundle)));
        DATOUT = cell(length(Active_channelsfiles),1);
        for cc=1:length(Active_channelsfiles)
            load(fullfile(CSCFiles(Active_channelsfiles(cc)).folder, CSCFiles(Active_channelsfiles(cc)).name),'AD_count_int16', 'AD_count_to_uV_factor','Bat_id','Date')
            DATOUT{cc} = single(AD_count_int16).*single(AD_count_to_uV_factor);
        end
        clear AD_cound_int16
        DATOUT = int16(cell2mat(DATOUT));
        whos DATOUT
        if strcmp('win',Kilosort2Platform)
            fid = fopen(fullfile(Output_folder, sprintf('%s_%s_TempTetrode%d.bin',Bat_id,Date,tt)),'w','l');
            fwrite(fid, DATOUT, 'int16','l');
        elseif strcmp('unix',Kilosort2Platform)
            fid = fopen(fullfile(Output_folder, sprintf('%s_%s_TempTetrode%d.bin',Bat_id,Date,tt)),'w','b');
            fwrite(fid, DATOUT, 'int16','b');
        end
        fclose(fid);
    end
    clear DATOUT
else
    % keep track of the active channels
    A=num2str(Active_channels_sorted);
    A = [repmat('_',length(Active_channels),1) A];
    ActChannels = reshape(A',1,numel(A));
    ActChannels(strfind(ActChannels, ' ')) = [];
    % find the number of samples of each original recording DAT file
    load(fullfile(CSCFiles(1).folder, CSCFiles(1).name),'Indices_of_first_and_last_samples','Bat_id','Date','AD_count_int16')
    DataLength = length(AD_count_int16);
    Ndatfile = size(Indices_of_first_and_last_samples,1);
    IndicesNdat = [1:round(Ndatfile/10):Ndatfile Ndatfile];
    OutputFileName = fullfile(Output_folder, sprintf('%s_%s_TempTetrode%s.bin',Bat_id,Date,ActChannels));
    % Make sure the file does not already exist and erase it in case it does
    if exist(OutputFileName,'file')
        delete(OutputFileName)
    end
    if strcmp('win',Kilosort2Platform)
        fid = fopen(OutputFileName,'a','l');
    elseif strcmp('unix',Kilosort2Platform)
        fid = fopen(OutputFileName,'a','b');
    else
        error('The machine platform for which the file should be written is not recognize, indicate either unix or win')
    end
    % loop through time sections to fill in a matrix with all active
    % channels OUTDAT
    Count = nan(length(IndicesNdat)-1,1);
    for ndat = 1:(length(IndicesNdat)-1)
        fprintf(1, 'Time section %d/%d for single experiment file\n', ndat,length(IndicesNdat)-1)
        OnIndex = Indices_of_first_and_last_samples(IndicesNdat(ndat),1);
        if ndat == (length(IndicesNdat)-1) % Last section
             OffIndex = Indices_of_first_and_last_samples(IndicesNdat(ndat+1),2);
        else
            OffIndex = Indices_of_first_and_last_samples(IndicesNdat(ndat+1)-1,2);
        end
%         NumSamp = OffIndex-OnIndex+1;
%         OUTDAT = nan(length(Active_channels), NumSamp);
        OUTDAT = cell(length(Active_channels),1);
        parfor cc=1:length(Active_channels)
            fprintf(1, 'Channel%d %d/%d \n',Active_channels_sorted(cc), cc,length(Active_channels))
            ch_local = find(Active_channels == Active_channels_sorted(cc));
            Data=load(fullfile(CSCFiles(ch_local).folder, CSCFiles(ch_local).name),'AD_count_int16', 'AD_count_to_uV_factor','DataDeletionOnsetOffset_sample');
            % Replace zones of identified artefacts by mean AD_count_int16
            % value so that artefact are not detected as spikes
            MeanVolt = mean(Data.AD_count_int16);
            dd=1;
            while (dd<=size(Data.DataDeletionOnsetOffset_sample,1)) && (Data.DataDeletionOnsetOffset_sample(dd,2)<=OffIndex)
                Data.AD_count_int16(Data.DataDeletionOnsetOffset_sample(dd,1) : Data.DataDeletionOnsetOffset_sample(dd,2)) = MeanVolt.*ones(diff(Data.DataDeletionOnsetOffset_sample(dd,:))+1,1);
                dd=dd+1;
            end
            OUTDAT{cc}  = single(-Data.AD_count_int16(OnIndex : OffIndex)).*single(Data.AD_count_to_uV_factor); % kilosort2 is based on negative threshold on the data at least for determining active channels, the voltage has been inverted in extract_logger_data, so inverting it back again here
        end
        OUTDAT = int16(cell2mat(OUTDAT));
        whos OUTDAT
        if strcmp('win',Kilosort2Platform)
            Count(ndat) = fwrite(fid,OUTDAT,'int16','l');
        elseif strcmp('unix',Kilosort2Platform)
            Count(ndat) = fwrite(fid,OUTDAT,'int16','b');
        end
        if Count(ndat)~=numel(OUTDAT)
            error('Samples not correctly written to the output file\n')
        end
    end
    if sum(Count)~=(DataLength*length(Active_channels))
        error('Samples not correctly written to the output file, error 2\n')
    end
    clear AD_cound_int16
    clear OUTDAT
    fclose(fid);
end


if s  %erase local data
    [sdel,mdel,edel]=rmdir(WorkDir, 's');
    if ~sdel
        TicErase = tic;
        while toc(TicErase)<30*60
            [sdel,mdel,edel]=rmdir(WorkDir, 's');
            if sdel
                return
            end
        end
    end
    if ~sdel
        sdel %#ok<NOPRT>
        mdel %#ok<NOPRT>
        edel %#ok<NOPRT>
        error('File erase did not occur correctly for %s\n Although we tried for 30min\n', WorkDir);
    end
end
end  
