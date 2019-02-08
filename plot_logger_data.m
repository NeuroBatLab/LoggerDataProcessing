function plot_logger_data(InputFolder,Start, LoggerID, Duration, NeuralPlotMode)

%% This function plots the audio loggers data as a waveform and a spectrogram and the neural data as a raster plot or a power plot of the neural signal bandpassed by BandPassFilter
% Input folder: folder containing the "loggers" or "audiologger" folder that itself contains
% one folder for each logger's data (e.g. Logger1) as used by extract_logger_data.
% Start: start point in logger time in seconds for plotting the next 10sec of data. If
% unspecified, start point is choosen as the common first recording point
% for the reqested loggers
% LoggerID: vector of doubles indicating logger #. if specified only these loggers will be plotted, otherwise all
% loggers in the input folder will be treated
% NeuralPlotMode: indicate whether to plot neural data as a rasterplot of
% all tetrodes non-spike sorted ('rasterMU'), spike sorted ('rasterSU')
% or amplitude power color plot of all channels
% ('power'). Default is set as 'power'.
Manual=1;
CodePath = pwd;
[RootCodePath,~] = fileparts(CodePath);
addpath(genpath(fullfile(RootCodePath,'SoundAnalysisBats')))
%% Treat input
if nargin<5
    NeuralPlotMode = 'power';
end
if nargin<4
    Duration=10;
end

LoggerFolders  = dir(fullfile(InputFolder, '*ogger*'));
LoggerFolders=LoggerFolders([LoggerFolders.isdir]);
InputLoggerFolder = InputFolder;
if isempty(LoggerFolders)
    InputLoggerFolder = fullfile(InputFolder, 'audiologgers');
    LoggerFolders  = dir(fullfile(InputLoggerFolder, '*ogger*'));
end
if isempty(LoggerFolders)
    InputLoggerFolder = input('Please indicate the full path to the folder containing each individual logger folder (e.g. logger1)','s');
    LoggerFolders = dir(fullfile(InputLoggerFolder, '*ogger*'));
end

if nargin<3
    NLog = length(LoggerFolders);
    LoggerID=nan(NLog,1);
    for ll=1:NLog
        Ind = strfind(LoggerFolders(ll).name,'r');
        LoggerID(NLog) = str2double(LoggerFolders(ll).name(Ind+1:end));
    end
else
    NLog = length(LoggerID);
end

if nargin<2
    % Loop through the loggers and determine the logger that was started
    % last, that will be the start time point for plotting
    StartFile = nan(NLog,1);
    for ll=1:NLog
        Ind = strfind(LoggerFolders(ll).name,'r');
        DataFolder = fullfile(InputLoggerFolder, sprintf('%s%d',LoggerFolders(ll).name(1:(Ind)),LoggerID(ll)), 'extracted_data');
        CSCFiles = dir(fullfile(DataFolder, '*CSC*.mat'));
        Data = load(CSCFiles(1).name); % Just load one of the data file (one file per channel)
        StartFile(ll) = Data.Timestamps_first_samples_usec(1); % extract the time stamp (transceiver time) of the first sample of the first recorded file.
    end
    Start = round(max(StartFile)/10^6);% in seconds
end

% parameters
BandPassFilter = [300 6000]; % frequency cut-offs of the band-pass filter on voltage signal
Fhigh_power = 20; % frequency value of the low-pass filter on voltage signal power ( Frequency upper bound for calculating the envelope (time running RMS))
Fs_env = 1000; % sampling frequency of the neural signal RMS set to be 1000Hz (1 sample per ms) and also of the sound enveloppe
BandPassFilter_spec = [1000 5000 9900]; % Frequency bands chosen for digital signal processing

Stop = Start + Duration;
%% Loop through the loggers, get the snippet of data of interest (from Start to Start + 10sec)
SnipData = cell(NLog,1);
TimeRef_usec = nan(NLog,1);
EstimatedFS = cell(NLog,1);
fprintf(1, '*** Gather data from each logger ***\n')
for ll=1:NLog
    fprintf(1, '%d/%d: Logger%d\n', ll,NLog, LoggerID(ll))
    Ind = strfind(LoggerFolders(ll).name,'r');
    DataFolder = fullfile(InputLoggerFolder, sprintf('%s%d',LoggerFolders(ll).name(1:(Ind)),LoggerID(ll)), 'extracted_data');
    CSCFiles = dir(fullfile(DataFolder, '*CSC*.mat'));
    NChannels = length(CSCFiles);
    % Identify the onset time of the first TTL pulse (That will be
    % used as a reference for transceiver time to allign the
    % loggers
    Events_file = dir(fullfile(CSCFiles(1).folder ,'*_EVENTS.mat'));
    Events = load(fullfile(Events_file.folder ,Events_file.name));
    FirstRisingInd = find(contains(Events.event_types_and_details, 'rising edge'),1,'first');
    FirstFallingInd = find(contains(Events.event_types_and_details, 'falling edge'),1,'first');
    if ~round((Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd))*10^-3)==6
        error('The first recorded TTL pulse is not of the expected duration of 6ms but is: %d usec\n',Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd));
    end
    TimeRef_usec(ll) = Events.event_timestamps_usec(FirstRisingInd);
    if NChannels>1 % this is a Neural Logger, two choices here: plot the amplitude of the filered voltage signal of each channel or extract spikes for each tetrode and plot a raster
        
        
        % Check if spike files are present and plot either a power
        % graph or a ratser.
        ST_files = dir(fullfile(DataFolder, '*Tetrode_spikes_time*.mat'));
        SU_files = dir(fullfile(DataFolder, '*TT*SS*.mat'));
        if strcmp(NeuralPlotMode, 'power') || (isempty(ST_files) && isempty(SU_files))
            NeuralPlotMode= 'power';
            fprintf(1, 'This is a neural logger. The power of each electrode is plotted.\n');
            SnipData{ll} = cell(NChannels,1);
            EstimatedFS{ll} = nan(NChannels,1);
            ParamDir = dir(fullfile(CSCFiles(1).folder, 'extract_logger_data_p*.mat'));
            Param = load(fullfile(ParamDir(1).folder,ParamDir(1).name));
            for cc=1:NChannels
                fprintf(1, 'Processing channel %d/%d\n', cc, NChannels);
                Data = load(fullfile(CSCFiles(cc).folder ,CSCFiles(cc).name)); % Load one of the data file
                
                % Identify the onset time closest to the start requested
                % and the sample start
                EstimatedFS{ll}(cc) = nanmean(Data.Estimated_channelFS_Transceiver(:,1));
                Timestamps_1st_samples_usec = Data.Timestamps_of_first_samples_usec - TimeRef_usec(ll);
                RecOnInd = find(Timestamps_1st_samples_usec/10^6<Start, 1, 'last');
                RecOnTime = Timestamps_1st_samples_usec(RecOnInd);
                RecOnSamp = Data.Indices_of_first_and_last_samples(RecOnInd,1);
                StartSamp = round(RecOnSamp + (Start-RecOnTime*10^-6)*EstimatedFS{ll}(cc));
                
                % Identify the onset time closest to the stop requested
                % and the sample stop
                RecOffInd = find(Timestamps_1st_samples_usec/10^6<Stop, 1, 'last');
                RecOffTime = Timestamps_1st_samples_usec(RecOffInd);
                RecOffSamp = Data.Indices_of_first_and_last_samples(RecOffInd,1);
                StopSamp = round(RecOffSamp + (Stop-RecOffTime*10^-6)*EstimatedFS{ll}(cc));
                
                % Extract data
                Voltage_Trace = double(int16(Data.AD_count_int16(StartSamp : StopSamp)))*Param.AD_count_to_uV_factor;
                
                % Bandpass the raw signal
                [b_band,a_band]=butter(6,BandPassFilter/(EstimatedFS{ll}(cc)/2),'bandpass'); % a 12th order Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency, as required by MATLAB)
                Filtered_voltage_trace=(filtfilt(b_band,a_band,Voltage_Trace)); % band-pass filter the voltage traces
                
                % Low pass filter the power trace
                % Figure out length of filter
                Nframes = length(Filtered_voltage_trace);
                if ( Nframes > 3*512 )
                    Nfilt = 512;
                elseif ( Nframes > 3*64 )
                    Nfilt = 64;
                    
                elseif ( Nframes > 3*16 )
                    Nfilt = 16;
                else
                    error('Data section is too short for filtering');
                end
                % Generate filter and filter signal power
                Lowpass_filter = fir1(Nfilt, Fhigh_power*2.0/EstimatedFS{ll}(cc));
                Amp_env_voltage = (filtfilt(Lowpass_filter, 1, Filtered_voltage_trace.^2)).^.5;
                % Resample to desired sampling rate
                if Data.Estimated_channelFS_Transceiver(cc) ~= Fs_env
                    SnipData{ll}{cc} = resample(Amp_env_voltage, Fs_env, round(EstimatedFS{ll}(cc)));
                else
                    SnipData{ll}{cc} = Amp_env_voltage;
                end
                figure(1)
                subplot(1,2,1)
                plot(Filtered_voltage_trace, '-k')
                hold on
                plot(Amp_env_voltage, '-r', 'LineWidth',2)
                legend('Filtered raw data', 'Running RMS')
                hold off
                subplot(1,2,2)
                plot(Filtered_voltage_trace.^2, '-k')
                hold on
                plot(Amp_env_voltage.^2, '-r', 'LineWidth',2)
                legend('Filtered data power','Amplitude envelope')
                hold off
                title(sprintf('Neural Logger SN%d Channel %d/%d',LoggerID(ll), cc, NChannels))
                
            end
        elseif strcmp(NeuralPlotMode, 'rasterMU') || isempty(SU_files)
            fprintf(1,'This is a neural logger. No data spike sorted, the detected spikes from the 4 tetrodes will be plotted\n');
            Num_MU = length(ST_files);
            SnipData{ll} = cell(Num_MU,1);
            for uu=1:Num_MU
                fprintf(1, 'Processing tetrode activity %d/%d\n', uu, Num_MU);
                load(fullfile(ST_files(uu).folder, ST_files(uu).name), 'Spike_arrival_times')
                % Convert spike arrival times
                %re-align spike arrival times that are in absolute transceiver time to
                % the present reference time,
                Spike_arrival_times = Spike_arrival_times - TimeRef_usec(ll);
                % Find the spike arrival times that are between the
                % requested times
                SnipData{ll}{uu} = Spike_arrival_times(logical((Spike_arrival_times/10^6>Start) .* (Spike_arrival_times/10^6<Stop)))/10^3-Start*10^3;
            end
        elseif strcmp(NeuralPlotMode, 'rasterSU')
            Num_SU = length(SU_files);
            fprintf(1,'This is a neural logger. %d Spike sorted units will be plotted\n', Num_SU);
            SnipData{ll} = cell(Num_SU,1);
            for uu=1:Num_SU
                fprintf(1, 'Processing single unit %d/%d\n', uu, Num_SU);
                load(fullfile(SU_files(uu).folder, SU_files(uu).name), 'Spike_arrival_times')
                % Convert spike arrival times
                %re-align spike arrival times that are in absolute transceiver time to
                % the present reference time,
                Spike_arrival_times = Spike_arrival_times - TimeRef_usec(ll);
                % Find the spike arrival times that are between the
                % requested times
                SnipData{ll}{uu} = Spike_arrival_times(logical((Spike_arrival_times/10^6>Start) .* (Spike_arrival_times/10^6<Stop)))/10^3-Start*10^3;
            end
        end
        
    else % This is an audio logger. Plot the wave form and a spectrogram
        Data = load(fullfile(CSCFiles.folder ,CSCFiles.name));
        % Identify the onset time closest to the start requested
        % and the sample start
        EstimatedFS{ll} = nanmean(Data.Estimated_channelFS_Transceiver);
        Timestamps_1st_samples_usec = Data.Timestamps_of_first_samples_usec - TimeRef_usec(ll);
        RecOnInd = find(Timestamps_1st_samples_usec/10^6<Start, 1, 'last');
        RecOnTime = Timestamps_1st_samples_usec(RecOnInd);
        RecOnSamp = Data.Indices_of_first_and_last_samples(RecOnInd,1);
        StartSamp = round(RecOnSamp + (Start-RecOnTime*10^-6)*EstimatedFS{ll});
        
        % Identify the onset time closest to the stop requested
        % and the sample stop
        RecOffInd = find(Timestamps_1st_samples_usec/10^6<Stop, 1, 'last');
        RecOffTime = Timestamps_1st_samples_usec(RecOffInd);
        RecOffSamp = Data.Indices_of_first_and_last_samples(RecOffInd,1);
        StopSamp = round(RecOffSamp + (Stop-RecOffTime*10^-6)*EstimatedFS{ll});
        
        % extract the snippet of sound
        SnipData{ll} = double(Data.AD_count_int16(StartSamp:StopSamp));
    end
end

%% Check that the time reference was the same for all loggers
if length(unique(TimeRef_usec))>1
    warning('The time reference, set up to be the rising edge of the first TTL pulse was not received at the same time by all loggers:\n');
    for ll=1:NLog
        fprintf('1st rising edge TTL pulse Logger %d: %.3f usec\n', LoggerID(ll), TimeRef_usec(ll));
    end
end

%% Plot the snippets of data
fprintf(1, '*** Plot data ***\n')
% Raw waveforms first
figure(2)
clf
MaxWaveform = nan(NLog,1);
for ll=1:NLog
    if iscell(SnipData{ll}) && strcmp(NeuralPlotMode, 'power') % This is a neural logger and we plot the power of each electrode
        % Snippets of data don't always have the exact same length due to
        % sample rate estimation variation between channels. the following lines ensure that the
        % snippets have the exact same length.
        LLDD = nan(length(SnipData{ll}),1);
        for dd=1:length(SnipData{ll})
            LLDD(dd) = length(SnipData{ll}{dd});
        end
        for dd=1:length(SnipData{ll})
            SnipData{ll}{dd} = SnipData{ll}{dd}(1:min(LLDD));
        end
        Local_Data = cell2mat(SnipData{ll});
        %Local_Data = Local_Data ./ repmat(max(Local_Data, [], 2), 1, size(Local_Data,2));
        subplot(NLog,1, ll)
        imagesc(Local_Data);
        colormap(hot);
        title(sprintf('NeuralLogger %d',LoggerID(ll)));
        ylabel('RMS 0.6-6KHz')
        xlabel('Time (ms)')
        colorbar
    elseif iscell(SnipData{ll}) && contains(NeuralPlotMode, 'raster') % This is a neural logger and we plot the rasters
        subplot(NLog,1, ll)
        hold on
        U_num = length(SnipData{ll});
        for uu=1:U_num
            for spike=1:length(SnipData{ll}{uu})
                plot(SnipData{ll}{uu}(spike)*ones(2,1), uu-[0.9 0.1], 'k-', 'LineWidth',1)
            end
        end
        xlabel('Time (ms)')
        xlim([0 Duration*10^3])
        ylim([0 U_num])
        if contains(NeuralPlotMode, 'MU')
            ylabel('Tetrodes')
        elseif contains(NeuralPlotMode, 'SU')
            ylabel('Single units')
        end
        hold off
    else % This is an audio logger
        subplot(NLog,1, ll)
        plot(SnipData{ll},'k-','LineWidth',2)
        title(sprintf('AudioLogger %d',LoggerID(ll)));
        ylabel('Amplitude')
        xlabel('Samples')
        MaxWaveform(ll) = max(abs(SnipData{ll}));
    end
end
% Set the axis of the waveforms with the same limits:
MaxW = max(MaxWaveform);
for ll=1:NLog
    if ~iscell(SnipData{ll})
        ss=subplot(NLog,1,ll);
        ss.YLim=[-MaxW MaxW];
        ss.XLim=[0 length(SnipData{ll})];
    end
end

% Then the spectrograms
dBScale = 60;
figure(3)
clf
MaxCmap = nan(NLog,1);
for ll=1:NLog
    if iscell(SnipData{ll})  && strcmp(NeuralPlotMode, 'power') % This is a neural logger and we plot the power of each electrode
        % Snippets of data don't always have the exact same length due to
        % sample rate estimation variation between channels. the following lines ensure that the
        % snippets have the exact same length.
        LLDD = nan(length(SnipData{ll}),1);
        for dd=1:length(SnipData{ll})
            LLDD(dd) = length(SnipData{ll}{dd});
        end
        for dd=1:length(SnipData{ll})
            SnipData{ll}{dd} = SnipData{ll}{dd}(1:min(LLDD));
        end
        Local_Data = cell2mat(SnipData{ll});
        %Local_Data = Local_Data ./ repmat(max(Local_Data,[],2), 1, size(Local_Data,2));
        ss=subplot(NLog,1, ll);
        subplot(NLog,1, ll)
        imagesc(Local_Data);
        colormap(ss, hot)
        title(sprintf('NeuralLogger %d',LoggerID(ll)));
        ylabel('RMS 0.6-6KHz')
        xlabel('Time (ms)')
        colorbar
    elseif iscell(SnipData{ll}) && contains(NeuralPlotMode, 'raster') % This is a neural logger and we plot the rasters
        subplot(NLog,1, ll)
        hold on
        U_num = length(SnipData{ll});
        for uu=1:U_num
            for spike=1:length(SnipData{ll}{uu})
                plot(SnipData{ll}{uu}(spike)*ones(2,1), uu-[0.9 0.1], 'k-', 'LineWidth',1)
            end
        end
        xlabel('Time (ms)')
        xlim([0 Duration*10^3])
        ylim([0 U_num])
        if contains(NeuralPlotMode, 'MU')
            ylabel('Tetrodes')
        elseif contains(NeuralPlotMode, 'SU')
            ylabel('Single units')
        end
        hold off
    else
        subplot(NLog,1,ll)
        % design the filters
        [z,p,k] = butter(6,BandPassFilter_spec(1:2)/(EstimatedFS{ll}/2),'bandpass');
        sos_low = zp2sos(z,p,k);
        [z,p,k] = butter(6,BandPassFilter_spec(2:3)/(EstimatedFS{ll}/2),'bandpass');
        sos_high = zp2sos(z,p,k);
        
        % filter the loggers' signals
        LowPassLogVoc = (filtfilt(sos_low,1,SnipData{ll})); % low-pass filter the voltage trace
        HighPassLogVoc = (filtfilt(sos_high,1,SnipData{ll})); % high-pass filter the voltage trace
        Amp_env_LowPassLogVoc=running_rms(LowPassLogVoc, EstimatedFS{ll}, Fhigh_power, Fs_env);
        Amp_env_HighPassLogVoc=running_rms(HighPassLogVoc, EstimatedFS{ll}, Fhigh_power, Fs_env);
        
        % Plot the low pass filtered signal of each logger
        [~, ~, logB, ~, ~] = spec_only_bats(LowPassLogVoc, EstimatedFS{ll});
        MaxCmap(ll) = max(max(logB));
        title(sprintf('AudioLogger %d',LoggerID(ll)));
        hold on
        yyaxis right
        plot((1:length(Amp_env_LowPassLogVoc))/Fs_env*1000, Amp_env_LowPassLogVoc, 'b-','LineWidth', 2);
        hold on
        plot((1:length(Amp_env_HighPassLogVoc))/Fs_env*1000, Amp_env_HighPassLogVoc, 'r-','LineWidth',2);
        ylabel('Amplitude')
        hold off
        
    end
end
% Set all the spectrograms at the same scale:
MaxB = max(MaxCmap);
for ll=1:NLog
    if iscell(SnipData{ll})
        ss=subplot(NLog,1, ll);
        subplot(NLog,1, ll)
        colormap(ss, hot)
%     else
%         ss=subplot(NLog,1,ll);
%         subplot(NLog,1,ll)
%         axis xy;
%         caxis('manual');
%         caxis([MaxB-dBScale MaxB]);
%         cmap = spec_cmap();
%         colormap(ss, cmap);
%         colorbar
%         if Manual
%             pause(0.1)
%             Player= audioplayer((SnipData{ll}-mean(SnipData{ll}))/std(SnipData{ll}), EstimatedFS{ll}); %#ok<TNMLP>
%             play(Player)
%             pause(1)
%         end
    end
end



end






