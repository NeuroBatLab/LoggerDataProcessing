function [Raw, LFP, SpikeT, SpikeSU, SpikeTDeNoiseInd, EventOnset_time, EventOffset_time] = extract_timeslot_LFP_spikes(LData_folder, OnsetOffset_transc_time_refined, Buffer, MaxEventDur, Flag, DenoiseT,Rthreshold)
%% This function extract the neural data of loggers according to onset and offset times of behavioral events
% LData_folder is the folder containing the extracted and spike sorted
% neural logger data (CSC matlab files, Tetrode_spikes_time matlab files and TT SS matlab files)
% OnsetOffset_transc_time_refined is a 2 column matrix with each row
% corresponding to a behavioral event which onset and offset times is ms in
% transceiver time are given by the first and second column.

BandPassFilter = [300 6000];
if nargin<3
    Buffer = 100; % delay in ms to add before and after the requested onset/offset time points for extracting data
end
if nargin<4
    MaxEventDur = NaN; % max duration of events in ms = MaxEventDur + Buffer. Neural response is cut beyond that.
end
if nargin<5
    Flag = [1 1 1 1];
end
if nargin<6
    DenoiseT = 0;
end
if nargin<7
    Rthreshold = [0.92 0.94 0.96 0.98];
end
%% list available data for that logger
% List CSC files (raw data)
LCSCDir = dir(fullfile(LData_folder, '*CSC*.mat'));
% List Spike time files (unsorted tetrode spikes)
ST_files = dir(fullfile(LData_folder, '*Tetrode_spikes_time*.mat'));
% List Spike shape files (unsorted tetrode spikes)
SS_files = dir(fullfile(LData_folder, '*Tetrode_spikes_snippets*.mat'));
% List SU Spike time files (spike-sorted tetrode spikes)
SU_files = dir(fullfile(LData_folder, '*TT*SS*.mat'));

%% Extract the LFP of each channel for each indicated set of onset/offset time points
LFP = [];
Raw = [];
if Flag(1) || Flag(2)
    Nevent = size(OnsetOffset_transc_time_refined,1);
    Nchannels = length(LCSCDir);
    Raw = cell(Nevent,Nchannels);
    if Flag(1)
        fprintf(1,'-- extracting Raw Voltage trace\n')
    end
    if Flag(2)
        fprintf(1,'-- extracting LFP\n')
        LFP = cell(Nevent,Nchannels);
    end
    FS = nan(Nevent,Nchannels);
    EventOnset_time = nan(Nevent,1);
    EventOffset_time = nan(Nevent, 1);
    for vv=1:Nevent
        if prod(~isnan(OnsetOffset_transc_time_refined(vv,:)))
            EventOnset_time(vv) = OnsetOffset_transc_time_refined(vv,1) - Buffer;
            if isnan(MaxEventDur)
                EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
            elseif diff(OnsetOffset_transc_time_refined(vv,:))>(MaxEventDur + Buffer)
                EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
            elseif diff(OnsetOffset_transc_time_refined(vv,:))<MaxEventDur
                EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
            elseif diff(OnsetOffset_transc_time_refined(vv,:))>MaxEventDur
                EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
            end
        end
    end

    parfor cc=1:Nchannels
        LData = load(fullfile(LCSCDir(cc).folder, LCSCDir(cc).name), 'Timestamps_of_first_samples_usec', 'Estimated_channelFS_Transceiver','AD_count_int16','Indices_of_first_and_last_samples','Sampling_period_usec_Logger');
        % loop through events and extract the snippet of Raw data
        NanFSInd = find(isnan(LData.Estimated_channelFS_Transceiver));
        GoodFSInd = find(~isnan(LData.Estimated_channelFS_Transceiver));
        for vv=1:Nevent
            fprintf('Raw data Channel %d Event %d/%d\n', cc, vv, Nevent);
            if prod(~isnan(OnsetOffset_transc_time_refined(vv,:))) %#ok<PFBNS>
                % find the time stamp on the logger that is closest to before
                % the snippet of event onset
                IndTSOn = find(LData.Timestamps_of_first_samples_usec<(EventOnset_time(vv)*10^3), 1, 'Last'); %#ok<PFBNS>
                if ~isempty(IndTSOn) %This event did happen after recording onset
                    % deduct the corresponding onset  sample
                    if ~isempty(intersect(NanFSInd,IndTSOn)) || IndTSOn>length(LData.Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                        OKFs = [find(GoodFSInd<IndTSOn,1, 'Last') find(GoodFSInd>IndTSOn,1, 'First')];
                        if ~isempty(OKFs)
                            Local_FS = LData.Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));
                        else
                            Local_FS = 1/LData.Sampling_period_usec_Logger*10^6;
                        end
                    else
                        Local_FS =LData.Estimated_channelFS_Transceiver(IndTSOn);
                    end
                    IndSampOn = round(LData.Indices_of_first_and_last_samples(IndTSOn,1) + Local_FS*(10^-6)*(EventOnset_time(vv)*10^3 - LData.Timestamps_of_first_samples_usec(IndTSOn)));

                    % find the time stamp on the logger that is closest to after
                    % the snippet of event offset
                    if (EventOffset_time(vv)*10^3)<LData.Timestamps_of_first_samples_usec(end) %#ok<PFBNS> % The end of the event is before the onset of the last recorded file
                        IndTSOff = find(LData.Timestamps_of_first_samples_usec>(EventOffset_time(vv)*10^3), 1, 'First');
                        % deduct the corresponding onset offset samples
                        if ~isempty(intersect(NanFSInd,IndTSOff)) || IndTSOff>length(LData.Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                            OKFs = [find(GoodFSInd<IndTSOff,1, 'Last') find(GoodFSInd>IndTSOff,1, 'First')];
                            if ~isempty(OKFs)
                                Local_FS = LData.Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));  
                                FS(vv,cc) = nanmean(LData.Estimated_channelFS_Transceiver(IndTSOn:GoodFSInd(OKFs(1)))); %#ok<PFOUS>
                            else
                                Local_FS = 1/LData.Sampling_period_usec_Logger*10^6;
                            end
                            
                        else
                            Local_FS =LData.Estimated_channelFS_Transceiver(IndTSOff);
                            FS(vv,cc) = nanmean(LData.Estimated_channelFS_Transceiver(IndTSOn:IndTSOff));
                        end
                        IndSampOff = round(LData.Indices_of_first_and_last_samples(IndTSOff,1) - Local_FS*(10^-6)*(LData.Timestamps_of_first_samples_usec(IndTSOff) - EventOffset_time(vv)*10^3));
                        
                    else % The end of the event is after the onset of the last recorded file
                        IndTSOff = length(LData.Timestamps_of_first_samples_usec);
                        % deduct the corresponding offset samples
                        if ~isempty(intersect(NanFSInd,IndTSOff)) || IndTSOff>length(LData.Estimated_channelFS_Transceiver) % there is no sample frequency estimate for that recorded file, take the previous or following estimate
                            OKFs = [find(GoodFSInd<IndTSOff,1, 'Last') find(GoodFSInd>IndTSOff,1, 'First')];
                            if ~isempty(OKFs)
                                Local_FS = LData.Estimated_channelFS_Transceiver(GoodFSInd(OKFs(1)));   
                            else
                                Local_FS = 1/LData.Sampling_period_usec_Logger*10^6;
                            end
                        else
                            Local_FS =LData.Estimated_channelFS_Transceiver(IndTSOff);
                        end
                        IndSampOff = round(LData.Indices_of_first_and_last_samples(IndTSOff,1) + Local_FS*(10^-6)*(EventOffset_time(vv)*10^3 - LData.Timestamps_of_first_samples_usec(IndTSOff)));
                        FS(vv,cc) = nanmean(LData.Estimated_channelFS_Transceiver(IndTSOn:end));
                    end
                    % extract the voltage snippet
                    IndSampOff = min(length(LData.AD_count_int16), IndSampOff);
                    Raw{vv, cc} = double(LData.AD_count_int16(IndSampOn:IndSampOff) - mean(LData.AD_count_int16));
                    if Flag(2) %#ok<PFBNS>
                        % Bandpass the raw signal
                        [z,p,k] = butter(12,BandPassFilter(1)/(FS(vv,cc)/2),'low'); %#ok<PFBNS>
                        sos = zp2sos(z,p,k);
                        LFP{vv,cc} = (filtfilt(sos,1,Raw{vv,cc})); % % low-pass filter the voltage trace
                    end
                end
            end
        end
    end
    clear LData
end

%% Extract the MU spike of each tetrode for each indicated set of onset/offset time points
if Flag(3)
    fprintf(1,'-- extracting MU tetrode data\n')
    if isempty(ST_files)
        SpikeT=[];
        fprintf(1,'********** WARNING: Spikes of tetrode cannot be extracted, no input files ********\n')
    else
        Nevent = size(OnsetOffset_transc_time_refined,1);
        Num_MU = length(ST_files);
        SpikeT = cell(Nevent,Num_MU);
        SpikeTDeNoiseInd = cell(Nevent,Num_MU);
        EventOnset_time = nan(Nevent, 1);
        EventOffset_time = nan(Nevent, 1);
        for uu=1:Num_MU
             fprintf(1, '--- Processing tetrode %d/%d\n', uu, Num_MU);
            % Load the spike arrival times for that tetrode
            Spikes = load(fullfile(ST_files(uu).folder, ST_files(uu).name), 'Spike_arrival_times');
            Snip = load(fullfile(SS_files(uu).folder, [SS_files(uu).name(1:(end-5)) ST_files(uu).name((end-4):end)]), 'Snippets');
            
            
            % loop through vocalizations and extract spike arrival times
            for vv=1:Nevent
                if uu==1
                    EventOnset_time(vv) = OnsetOffset_transc_time_refined(vv,1) - Buffer;
                    if isnan(MaxEventDur)
                        EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))>(MaxEventDur + Buffer)
                        EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))<MaxEventDur
                        EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))>MaxEventDur
                        EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                    end
                end
                % Find the spike arrival times that are between the
                % requested times and center them to the onset of the
                % behavioral event
                SpikeT_local = logical((Spikes.Spike_arrival_times>(EventOnset_time(vv)*10^3)) .* (Spikes.Spike_arrival_times<(EventOffset_time(vv)*10^3)));
                SpikeT{vv,uu} = Spikes.Spike_arrival_times(SpikeT_local)/10^3 - OnsetOffset_transc_time_refined(vv,1);
                SpikeTDeNoiseInd{vv,uu} = cell(length(Rthreshold),1);
                if DenoiseT % get the indices of denoised spikes according to threshold(s)
                    for rr=1:length(Rthreshold)
                        SpikeTDeNoiseInd{vv,uu}{rr} = sort_spike_from_noise(Snip.Snippets(:,:,SpikeT_local), Rthreshold(rr));
                    end
                end
            end
        end
    end
else
    SpikeT=[];
    SpikeTDeNoiseInd=[];
end
%% Extract the spike sorted spike of each single unit for each indicated set of onset/offset time points
if Flag(4)
    fprintf(1,'-- extracting Spike sorted tetrode data\n')
    if isempty(SU_files)
        SpikeSU=[];
        fprintf(1,'********** WARNING: Spikes of single units cannot be extracted, no input files ********\n')
    else
        Nevent = size(OnsetOffset_transc_time_refined,1);
        Num_SU = length(SU_files);
        SpikeSU = cell(Nevent,Num_SU);
        EventOnset_time = nan(Nevent, 1);
        EventOffset_time = nan(Nevent, 1);
        for uu=1:Num_SU
            fprintf(1, '--- Processing single unit %d/%d\n', uu, Num_SU);
            % loading the single unit spike arrival times
            Spikes = load(fullfile(SU_files(uu).folder, SU_files(uu).name), 'Spike_arrival_times');
            % loop through vocalizations and extract spike arrival times
            for vv=1:Nevent
                
                if uu==1
                    EventOnset_time(vv) = OnsetOffset_transc_time_refined(vv,1) - Buffer;
                    if isnan(MaxEventDur)
                        EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))>(MaxEventDur + Buffer)
                        EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))<MaxEventDur
                        EventOffset_time(vv) = OnsetOffset_transc_time_refined(vv,2) + Buffer;
                    elseif diff(OnsetOffset_transc_time_refined(vv,:))>MaxEventDur
                        EventOffset_time(vv) = EventOnset_time(vv) + MaxEventDur + 2*Buffer;
                    end
                end
                % Find the spike arrival times that are between the
                % requested times and center them to the onset of the
                % behavioral event
                SpikeSU{vv,uu} = Spikes.Spike_arrival_times(logical((Spikes.Spike_arrival_times>(EventOnset_time(vv)*10^3)) .* (Spikes.Spike_arrival_times<(EventOffset_time(vv)*10^3))))/10^3 - OnsetOffset_transc_time_refined(vv,1);
            end
        end
    end
else
    SpikeSU=[];
end
end