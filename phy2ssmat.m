function phy2ssmat(InputPath, CodePath,OutputPath)
% this function expect that the original CSC files are sitting in
% 'extracted_data' folder at the same level as InputPath

%% Massage input
%Hard coded input:
Num_EperBundle = 4; % We work with tetrodes
Num_E = 4; %4 tetrodes
TotChannels = Num_E*Num_EperBundle;
TimeStep = 20*60*10^6;% Time Step of 30 min for the calculation of the time varying spike sorting quality measures
DebugFig = 1; % Set to 1 to see figures of spike sorting quality.
BandPassFilter = [600 6000];% Bandpass parameters for the input raw voltage
Compare2CSCData = 0;

% Deal with input variables
if nargin<1
    error('please provide the path to the files containing\nthe snippets and spike arrival time generated by Phy\n');
end

% Get the folder containing the raw data as extracted by extract_logger_data.m
[Root,~]=fileparts(InputPath);
Data_folder = fullfile(Root, 'extracted_data');
if nargin<2
    CodePath = input('Please give the path to the folder containing all Github code (Kilosort2_Tetrode, npy-matlab)\n', 's');
end
addpath(genpath(fullfile(CodePath,'Kilosort2_Tetrode')))
addpath(genpath(fullfile(CodePath,'npy-matlab')))
pathToYourConfigFile = fullfile(CodePath,'Kilosort2_Tetrode/configFiles');

if nargin<3
    OutputPath = Data_folder;
    fprintf(1, 'No outputpath indicated, data will be saved under the raw data folder:\n%s\n', OutputPath);
else
    fprintf(1, 'Converted files will be saved under:\n%s\n', OutputPath);
end



%% Load data
% Load the data output after currating the data with Phy
params.excludeNoise = true;
params.loadPCs = false;
SpikeStruct = loadKSdir(InputPath, params);

% Load the raw data out from kilosort2
% REZ = load(fullfile(InputPath, 'rez.mat'));

% Get the info of that recording
BatID = SpikeStruct.dat_path(1:5);
Date = SpikeStruct.dat_path(7:14);
NTemplatePoints = size(SpikeStruct.temps,2);
NChannels = size(SpikeStruct.temps,3);

% Get the channels ID and the channel ID per tetrode
Ind_ = strfind(SpikeStruct.dat_path, '_');
Ind_ = Ind_(3:end);
NCh = length(Ind_);

ChannelsID = nan(NCh,1);% zero indexed
ChannelsID_perT = cell(Num_E,1);% zero indexed
for cc=1:NCh
    if cc==NCh
        Indend = strfind(SpikeStruct.dat_path, '.');
        ChannelsID(cc) = str2double(SpikeStruct.dat_path(Ind_(cc)+1 : Indend-1));
        TetrodeID = find(ChannelsID(cc)+1<=Num_EperBundle:Num_EperBundle:TotChannels,1,'first');
        ChannelsID_perT{TetrodeID} = [ChannelsID_perT{TetrodeID} ChannelsID(cc)];
    else
        ChannelsID(cc) = str2double(SpikeStruct.dat_path(Ind_(cc)+1 : Ind_(cc+1)-1));
        TetrodeID = find(ChannelsID(cc)+1<=Num_EperBundle:Num_EperBundle:TotChannels,1,'first');
        ChannelsID_perT{TetrodeID} = [ChannelsID_perT{TetrodeID} ChannelsID(cc)];
    end
end

if NCh~=NChannels
    % Find the good map for tetrodes 
    if NCh==16
        chanMapFile = 'Tetrodex4Default_kilosortChanMap.mat';
        ChannelMap = fullfile(pathToYourConfigFile, chanMapFile);
    elseif NCh==15 && strcmp(BatID, '11689') % This is Hodor, who missed channel 11 (12th channel)
        chanMapFile = 'Tetrodex4Ho_kilosortChanMap.mat';
        ChannelMap = fullfile(pathToYourConfigFile, chanMapFile);
    else
        ChannelMap = input('Indicate the path and name of the matfile for your channel map :\n','s');
    end
    ChanMap = load(ChannelMap);
    KSChanMap = nan(NChannels,1);
    for cc=1:NChannels
        KSChanMap(cc) = ChannelsID(logical((ChanMap.xcoords == SpikeStruct.xcoords(cc)) .* (ChanMap.ycoords == SpikeStruct.ycoords(cc))))+1;
    end
        
else
    KSChanMap = 1:NCh;
end



% BandPass filter for the raw data
[b,a]=butter(6,BandPassFilter/(SpikeStruct.sample_rate/2),'bandpass');

%% Loop through each SU or MU and save data
UnitClusters = SpikeStruct.cids;
UnitClustersQ = SpikeStruct.cgs;
Nunits = length(UnitClusters);
Qlabel = {'M','S','U'}; % M: multi-unit, S: single unit, U: Unsorted or Uncertain between noise or unit
for uu=1:Nunits
    fprintf(1, 'Processing unit %d/%d\n', uu, Nunits);
    ClustID = UnitClusters(uu);
    ClustQ = Qlabel{SpikeStruct.cgs(uu)};
    SpikeTemplatesID = SpikeStruct.spikeTemplates(SpikeStruct.clu == ClustID);
    SpikePosition_local = SpikeStruct.ss(SpikeStruct.clu == ClustID);
    
    % find on which channel this unit was detected
    Utemplates = unique(SpikeTemplatesID);
    fprintf(1, '%d templates for that unit\n', length(Utemplates));
    ChannelID = nan(length(Utemplates),1);
    Templates = cell(length(Utemplates),1);
    if uu==1
        FIG = figure();
        pause()
        ColorCode = get(groot,'DefaultAxesColorOrder');
    end
    if length(Utemplates)>6 && ~exist('FIGbis','var')
        FIGbis = figure();
        pause()
    end
    if length(Utemplates)>12 && ~exist('FIGter','var')
        FIGter = figure();
        pause()
    end
    if length(Utemplates)>18 && ~exist('FIGqua','var')
        FIGqua = figure();
        pause()
    end
    for tt=1:length(Utemplates)
        Templates{tt} = squeeze(SpikeStruct.temps(Utemplates(tt),:,:));
        AmpTemplate = max(Templates{tt}) - min(Templates{tt});
        [~,MaxAmpTempInd] = max(AmpTemplate);
        ChannelID(tt) = KSChanMap(MaxAmpTempInd);
        if tt<7
            set(0,'CurrentFigure', FIG)
            subplot(3,2,tt)
        elseif tt<13
            set(0,'CurrentFigure', FIGbis)
            subplot(3,2,tt-6)
        elseif tt<19
            set(0,'CurrentFigure', FIGter)
            subplot(3,2,tt-12)
        else
            set(0,'CurrentFigure', FIGqua)
            subplot(3,2,tt-18)
        end
        
        
        for cc=1:NChannels
            if KSChanMap(cc)<5
                plot(Templates{tt}(:,cc), 'LineWidth',2, 'Color',ColorCode(KSChanMap(cc),:))
            elseif KSChanMap(cc)<9
                plot(Templates{tt}(:,cc), 'LineWidth',2,'LineStyle','--','Color',ColorCode(KSChanMap(cc)-4,:))
            elseif KSChanMap(cc)<13
                plot(Templates{tt}(:,cc), 'LineWidth',2,'LineStyle','-.', 'Color',ColorCode(KSChanMap(cc)-8,:))
            else
                plot(Templates{tt}(:,cc), 'LineWidth',2,'LineStyle',':', 'Color',ColorCode(KSChanMap(cc)-12,:))
            end
            hold on
        end
        hold off
%         if length(Utemplates)>1 && tt==length(Utemplates)
%             legend('location', 'SouthOutside', 'NumColumns',2)
%         elseif length(Utemplates)==1
%             legend('location', 'SouthOutside', 'NumColumns',2)
%         end
        ChannT = mod(ChannelID(tt),4);
        if ~ChannT
            ChannT=4;
        end
        title(sprintf('Template %d max on Channel %d, TT%dC%d', tt, ChannelsID(ChannelID(tt))+1,ceil((ChannelsID(ChannelID(tt))+1)/4),ChannT));
    end
    
    UChannelID = unique(ChannelID);
    if length(UChannelID)>1
        warning('This unit is formed by spikes that were detected by %d templates belonging to different channels (see below)\n', length(Utemplates))
        UChannelID
        UChannelID = input('Which Channel ID do you want to choose as a reference?\n');
    end
    
    % determine to which tetrode that channel belonged and the associated
    % channels
    TetrodeID = find(UChannelID<=Num_EperBundle:Num_EperBundle:TotChannels,1,'first');
    Channel_rows = sum(cellfun(@length,ChannelsID_perT(1:TetrodeID-1))) + (1:length(ChannelsID_perT{TetrodeID}));
    
    % Calculate spike arrival time in transceiver time from SpikeStruct, in
    % microseconds (Spike_arrival_times)
    FileDir = dir(fullfile(Data_folder,sprintf('*CSC%d.mat', UChannelID-1)));
    Filename=fullfile(FileDir.folder,FileDir.name);
    % Convert the indices of spike arrival times to real
    % time since we can load the info regarding time of file
    % onsets, and the sample frequency.
    load(Filename, 'Indices_of_first_and_last_samples');
    load(Filename, 'Estimated_channelFS_Transceiver');
    FS = nanmean(Estimated_channelFS_Transceiver);
    load(Filename, 'Timestamps_of_first_samples_usec');
    % restrict the data set to only the first files (quick fix when there is a
    % clock jump in the recording, only keep the first part of the recording before the jump)
%     IndJump = 1111;
    IndJump = length(Timestamps_of_first_samples_usec);
    Spike_arrival_times=round(get_timestamps_for_Nlg_voltage_samples(SpikePosition_local,Indices_of_first_and_last_samples(1:IndJump,1)',Timestamps_of_first_samples_usec(1:IndJump),10^6/FS));% the time stamps of all the detected spikes, rounded to integer microseconds; note that these are the time stamps of the last channel on this electrode bundle, which differ from the time stamps on the other channels of this electrode bundle by a few sampling periods of the Nlg AD converter
    fprintf(1,'Spike arrival times of %d spikes were not calculated with End of session identified at %d and were discarded\n',sum(Spike_arrival_times==0),IndJump)
    Spike_arrival_times(Spike_arrival_times==0) =[];
    Num_spikes=length(Spike_arrival_times);
    if Num_spikes<length(SpikePosition_local) % part of the recording was discarded, corect the TemplateID vector to the spikes that are collected
        SpikeTemplatesID = SpikeTemplatesID(1:Num_spikes);
    end
    
    % For all channels of the bundle (tetrode), collect the spike snippets
    fprintf(1, 'Collecting snippets\n')
    Spike_snippets=nan(NTemplatePoints,Num_EperBundle,Num_spikes);
    FileID = fopen(fullfile(InputPath, SpikeStruct.dat_path));
    for spike_i=1:Num_spikes
        if ~mod(spike_i,round(Num_spikes/10))
            fprintf(1, '%d Percent\n',spike_i/round(Num_spikes/10)*10);
        end
        fseek(FileID,(SpikePosition_local(spike_i) - NTemplatePoints)*NCh*2,-1);
        Data = fread(FileID,[NCh NTemplatePoints*2], '*int16','l');
        Voltage_Trace = double(Data(Channel_rows,:));
        for channel_i=1:length(Channel_rows)
            Filtered_voltage_trace = filtfilt(b,a,Voltage_Trace(channel_i,:));
            Spike_snippet = Filtered_voltage_trace(NTemplatePoints/2+(1:NTemplatePoints))';
            Spike_snippets(:,channel_i,spike_i)= Spike_snippet; % save the waveforms of the current spike (units are uV)
        end
    end
    fclose(FileID);
    Mat_Filename = fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d.mat',BatID, Date,TetrodeID,ClustQ,ClustID));
    if any(Spike_arrival_times==0)
        fprintf('Some spike arrival times were not correctly calculated, manual input required\n')
        keyboard
    end
    save(Mat_Filename, 'Spike_arrival_times', 'SpikeTemplatesID', 'Spike_snippets', 'Templates', 'ChannelID','UChannelID')
    
    % Now check that we get the same output using the raw CSC data
    if IndJump == length(Timestamps_of_first_samples_usec) && Compare2CSCData % we can easily check, other wise, don't do it, not implemented yet...
        [Spike_times, Spike_snippets2] = extract_tetrode_snippets(SpikePosition_local, Data_folder, ChannelsID_perT{TetrodeID},[-NTemplatePoints/2+1 NTemplatePoints/2]); % here we use the a Spike_window of the same size as in kilosort2
        Spike_snippets2 = -Spike_snippets2; % we are inverting the signal sign for kilosort2 calculations
        if any((Spike_times-Spike_arrival_times)>1000) % check if spike srrival times differ by more than 1ms
            warning('Error in calculations of spike arrival times')
            keyboard
        end
        % loop through snippets and check how different they are!
        ChannelRef = mod(UChannelID-1,4)+1;
        for ss = 1:size(Spike_snippets,3)
            SpikeCorr = corr(Spike_snippets(:,ChannelRef,ss),Spike_snippets2(:,ChannelRef,ss));
            if SpikeCorr<0.85
                figure(10)
                clf
                plot(Spike_snippets(:,ChannelRef,ss),'b-');
                hold on;plot(Spike_snippets2(:,ChannelRef,ss),'r-');
                title(sprintf('Correlation %.4f',SpikeCorr))
                hold off
                warning('Error in the spike extraction!?')
                keyboard
            end
        end
    end
    
    savefig(FIG,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_template.fig',BatID, Date,TetrodeID,ClustQ,ClustID)))
    print(FIG,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_template.pdf',BatID, Date,TetrodeID,ClustQ,ClustID)),'-dpdf','-fillpage')
    clf(FIG)
    if length(Utemplates)>6
        savefig(FIGbis,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templatebis.fig',BatID, Date,TetrodeID,ClustQ,ClustID)))
        print(FIGbis,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templatebis.pdf',BatID, Date,TetrodeID,ClustQ,ClustID)),'-dpdf','-fillpage')
        clf(FIGbis)
    end
    if length(Utemplates)>12
        savefig(FIGter,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templateter.fig',BatID, Date,TetrodeID,ClustQ,ClustID)))
        print(FIGter,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templateter.pdf',BatID, Date,TetrodeID,ClustQ,ClustID)),'-dpdf','-fillpage')
        clf(FIGter)
    end
    if length(Utemplates)>18
        savefig(FIGqua,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templatequa.fig',BatID, Date,TetrodeID,ClustQ,ClustID)))
        print(FIGqua,fullfile(OutputPath,sprintf('%s_%s_TT%d_SS%s_%d_templatequa.pdf',BatID, Date,TetrodeID,ClustQ,ClustID)),'-dpdf','-fillpage')
        clf(FIGqua)
    end
end


end
