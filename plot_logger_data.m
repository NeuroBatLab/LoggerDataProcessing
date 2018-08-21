function plot_logger_data(InputFolder,Start, LoggerID, Duration, NeuralPlotMode)

%% This function plots the audio loggers data as a waveform and a spectrogram and the neural data as a raster plot or a power plot of the neural signal bandpassed by BandPassFilter
% Input folder: folder containing the "loggers" folder that itself contains
% one folder for each logger's data as generated by extract_logger_data.
% Start: start point in logger time in seconds for plotting the next 10sec of data. If
% unspecified, start point is choosen as the common first recording point
% for the reqested loggers
% LoggerID: vector of doubles indicating logger #. if specified only these loggers will be plotted, otherwise all
% loggers in the input folder will be treated
% NeuralPlotMode: indicate whether to plot neural data as a rasterplot of
% all tetrodes ('raster') or amplitude power color plot of all channels
% ('power'). Default is set as 'power'.

addpath(genpath('/Users/elie/Documents/CODE/tlab/src/'))
%% Treat input
if nargin<5
    NeuralPlotMode = 'power';
end
if nargin<4
    Duration=10;
end
if nargin<3
     LoggerFolders  = dir(fullfile(InputFolder, 'Loggers', 'Logger*'));
     NLog = length(LoggerFolders);
     LoggerID=nan(NLog,1);
     for ll=1:NLog
         Ind = strcmp('r', LoggerFolders(ll).name);
         LoggerID(NLog) = str2double(LoggerFolders(ll).name(Ind+1:end));
     end
else
     NLog = length(LoggerID);
end

if nargin<2
    % Loop through the loggers and determine the logger that was started
    % last, that will be the start time point for plotting
    StartFile1 = nan(NLog,1);
    for ll=1:NLog
        DataFolder = fullfile(InputFolder, 'Loggers', sprintf('Logger%d',LoggerID(ll)), 'extracted_data');
        DataFiles = dir(fullfile(DataFolder, 'CSC*.mat'));
        Data = load(DataFiles(1).name); % Just load one of the data file (one file per channel)
        StartFile1(ll) = Data.Timestamps_first_samples_usec(1); % extract the time stamp (trasnceiver time) of the first sample of the first recorded file.
    end
    Start = round(max(StartFile)/10^6);% in seconds
end

BandPassFilter = [300 6000]; % frequency cut-offs of the band-pass filter on voltage signal
Fhigh_power = 20; % frequency value of the low-pass filter on voltage signal power
Fs_env = 1000; % sampling frequency of the neural signal RMS set to be 1000Hz (1 sample per ms)

Stop = Start + Duration;
%% Loop through the loggers, get the snippet of data of interest (from Start to Start + 10sec)
SnipData = cell(NLog,1);
TimeRef_usec = nan(NLog,1);
EstimatedFS = cell(NLog,1);
for ll=1:NLog
        DataFolder = fullfile(InputFolder, 'Loggers', sprintf('Logger%d',LoggerID(ll)), 'extracted_data');
        DataFiles = dir(fullfile(DataFolder, 'CSC*.mat'));
        NChannels = length(DataFiles);
        if NChannels>1 % this is a Neural Logger, two choices here: pot the amplitude of the filered voltage signal of each channel or extract spikes for each tetrode and plot a raster
            if strcmp(NeuralPlotMode, 'power')
                SnipData{ll} = cell(NChannels,1);
                EstimatedFS{ll} = nan(NChannels,1);
                ParamDir = dir(fullfile(DataFiles(1).folder, 'extract_logger_data_p*.mat'));
                Param = load(fullfile(ParamDir(1).folder,ParamDir(1).name));
                for cc=1:NChannels
                    Data = load(fullfile(DataFiles(cc).folder ,DataFiles(cc).name)); % Load one of the data file
                        % Identify the onset time of the first TTL pulse (That will be
                     % used as a reference for transceiver time to allign the
                     % loggers
                     Events = load(fullfile(DataFiles(cc).folder ,'EVENTS.mat'));
                     FirstRisingInd = find(contains(Events.event_types_and_details, 'rising edge'),1,'first');
                     FirstFallingInd = find(contains(Events.event_types_and_details, 'falling edge'),1,'first');
                     if ~round((Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd))*10^-3)==6
                         error('The first recorded TTL pulse is not of the expected duration of 6ms but is: %d usec\n',Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd));
                     end
                     TimeRef_usec(ll) = Events.event_timestamps_usec(FirstRisingInd);

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
            elseif strcmp(NeuralPlotMode, 'raster')
                % HERE goes the code to plot a raster plot of the data
            end
        else % This is an audio logger. Plot the wave form and a spectrogram
             Data = load(fullfile(DataFiles.folder ,DataFiles.name));
             % Identify the onset time of the first TTL pulse (That will be
             % used as a reference for transceiver time to allign the
             % loggers
             Events = load(fullfile(DataFiles.folder ,'EVENTS.mat'));
             FirstRisingInd = find(contains(Events.event_types_and_details, 'rising edge'),1,'first');
             FirstFallingInd = find(contains(Events.event_types_and_details, 'falling edge'),1,'first');
             if ~round((Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd))*10^-3)==6
                 error('The first recorded TTL pulse is not of the expected duration of 6ms but is: %d usec\n',Events.event_timestamps_usec(FirstFallingInd) - Events.event_timestamps_usec(FirstRisingInd));
             end
             TimeRef_usec(ll) = Events.event_timestamps_usec(FirstRisingInd);
             
             
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
                SnipData{ll} = Data.AD_count_int16(StartSamp:StopSamp) - mean(Data.AD_count_int16);
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

% Raw waveforms first
figure(2)
MaxWaveform = nan(NLog,1);
for ll=1:NLog
    if iscell(SnipData{ll})
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
    else
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
fband = 100;
dBScale = 60;
Fhigh = 8000;
figure(3)
MaxCmap = nan(NLog,1);
for ll=1:NLog
    if iscell(SnipData{ll})
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
    else
        subplot(NLog,1,ll)
        [~, ~, logB, ~, ~] = spec_only_bats(SnipData{ll}, fband, EstimatedFS{ll}, dBScale, Fhigh);
        MaxCmap(ll) = max(max(logB));
        title(sprintf('AudioLogger %d',LoggerID(ll)));
    end
end
% Set all the spectrograms at the same scale:
MaxB = max(MaxCmap);
for ll=1:NLog
    if iscell(SnipData{ll})
        ss=subplot(NLog,1, ll);
        subplot(NLog,1, ll)
        colormap(ss, hot)
    else
        ss=subplot(NLog,1,ll);
        subplot(NLog,1,ll)
        axis xy;
        caxis('manual');
        caxis([MaxB-dBScale MaxB]); 
        cmap = spec_cmap();
        colormap(ss, cmap);
        colorbar
    end
end



end


                
            
     

