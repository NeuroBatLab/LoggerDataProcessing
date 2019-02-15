function [Pulse_idx, Pulse_TimeStamp_Transc, File_number, Pulse_samp_audio, Slope_and_intercept,Mean_std_x, Mean_std_Pulse_TimeStamp_Transc,Mean_std_Pulse_samp_audio] = align_soundmexAudio_2_logger(Audio_dir, Loggers_dir, ExpStartTime,varargin)
%%
% Function to correct for clock drift between audio recordings and
% Deuteron logger recordings. This function extract the positions of pulses
% in raw audio recordings and in the time at which these pulses are
% detected by the Deuteron transceiver. A linear fit is run for each audio
% file such that the transceiver time of audio events can be reatrieved
% from the sample indices of these audio events in the continuous raw
% recordings.
%
% INPUT:
% Audio_dir: folder containing the audio recordings
%
% Loggers_dir: base directory containing loggers data. This script expects
% this directory to contain the sufolders 'loggerx' wher x is the serial
% number of the logger
%
% ExpStartTime: identify the exact session you want to look at by the time
% it started (indicated in the filename
%
% Session_strings: cell of strings used to demarcate start and stop of
% time period to analyze in this script from EVENTLOG file generated by extract_logger_data.
%
% Logger_list: vector of numbers corresponding to the logger ID to which
% the extraction should be limited to. By default all loggers are
% considered
%
% TTL_pulse_generator: Indicate how the TTL pulse where generated either by
% the MOTU sound card ('MOTU') or by Avisoft ('Avisoft'). Default is
% 'MOTU'.
%
% OUTPUT:
% Pulse_idx: common pulse indices identified between soundcard recordings
% and logger recordings.
%
% Pulse_TimeStamp_Transc: times (ms) in transceiver time when TTL pulses arrived on
% the transceiver.
%
% Pulse_samp_audio: sample indices when TTL pulses arrived 
% on the soundcard.
%
% File_number: identity of the audio file to which each pulse belongs
%
% 
%%%
save_options_parameters_CD_figure = 1;
pnames = {'Method', 'Session_strings','TTL_pulse_generator', 'Logger_list'};
dflts  = {'risefall',{'record only start' 'record only stop'},'MOTU',[]};
[Method, Session_strings, TTL_pulse_generator, Logger_list] = internal.stats.parseArgs(pnames,dflts,varargin{:});


%% Find  the audio files and parameter file for TTL pulses
% Wav_files = dir(fullfile(Audio_dir, sprintf('*%s*mic*.wav', ExpStartTime))); % all microphone recordings .WAV files for the requested experiment/session (soundmexpro format)
TTL_files = dir(fullfile(Audio_dir, sprintf('*%s*ttl*.wav',ExpStartTime)));% all ttl recordings .WAV files for the requested experiment/session (soundmexpro format)
ExpDate = TTL_files(1).name(6:11);
if strcmp(TTL_pulse_generator,'MOTU')
    TTL_paramfile = dir(fullfile(Audio_dir,sprintf('*%s*unique_ttl_params.mat',ExpStartTime))); % this is the parameter file generated by generate_zero_playback_file.m (code generating TTL pulse playback sequences)
    TTL_param = load(fullfile(TTL_paramfile.folder, TTL_paramfile.name));
elseif strcmp(TTL_pulse_generator, 'Avisoft')
    TTLFolder = fullfile(Audio_dir(1:(strfind(Audio_dir, 'audio')-1)),'TTLPulses4h');
    TTL_paramfile = dir(fullfile(TTLFolder,'*unique_ttl_params.mat')); % this is the parameter file generated by generate_zero_playback_file.m (code generating TTL pulse playback sequences)
    TTL_param = load(fullfile(TTL_paramfile.folder, TTL_paramfile.name));
end
AudioEvent_file = dir(fullfile(Audio_dir, sprintf('*%s*events.txt',ExpStartTime)));

%% Extract TTL pulses from the audio recordings
fprintf('Extract TTL pulses index and value from audio recordings\n');
Pulse_idx_audio = cell(length(TTL_files),1); % This cell array will contain the pulse train index of each TTL pulse identified in the recordings
File_number = cell(length(TTL_files),1); % this cell array will contain the TTL file index number where the corresponding pulses are found
Pulse_samp_audio = cell(length(TTL_files),1); % This cell array will contain the sample index of each TTL pulse onset in the recordings
Num_Samp_audiofile = nan(length(TTL_files),2); % First column is file ID, second is number of samples

for w = 1:length(TTL_files) % run through all .WAV files and extract audio data and TTL status at each sample
    %     if ~strcmp(Wav_files(w).name(end-4), TTL_files(w).name(end-4))
    %         error('There is an error in audio files. Each Audio file should have a twin brother TTL file and here:\n%s\n%s\n', Wav_files(w).name, TTL_files(w).name);
    %     end
    fprintf('TTL file %d/%d\n', w, length(TTL_files));
    Ind_ = strfind(TTL_files(w).name, '_');
    CurrentFile = dir(fullfile(Audio_dir, sprintf('%s%d.wav',TTL_files(w).name(1:Ind_(end)),w)));
    try
        [Ttl_status, FS] = audioread(fullfile(CurrentFile.folder, CurrentFile.name));
    catch
        fprintf('Error with the TTL file: %s\nNo allignment for these data\n',CurrentFile.name)
        Ttl_status = [];
    end
    
    if ~isempty(Ttl_status)
        Num_Samp_audiofile(w,2)= length(Ttl_status);
        if strcmp(TTL_pulse_generator, 'MOTU')
            TTLHigh = find(diff(Ttl_status)>0.5)+1; % identify the increases in volatge
            TTLHigh((find(diff(TTLHigh)==1))+1)=[]; % eliminate consecutive points that show a large increase in the TTL pulses, they are just the continuity of a single pulse start (voltage going up)
            TTLLow = find(diff(Ttl_status)<-0.80)+1; % identify the decreases in volatge
            TTLLow((find(diff(TTLLow)==1))+1)=[]; % eliminate consecutive points that show a large decrease in the TTL pulses, they are just the continuity of a single pulse start (voltage going up)
            if length(TTLHigh)>length(TTLLow)% Try a lower threshold we are missing some pulses
                TTLLow = find(diff(Ttl_status)<-0.70)+1; % identify the decreases in volatge
                TTLLow((find(diff(TTLLow)==1))+1)=[]; % eliminate consecutive points that show a large decrease in the TTL pulses, they are just the continuity of a single pulse start (voltage going up)
                if length(TTLHigh)~=length(TTLLow)%
                    error('Error in align_soundmexAudio_2_logger: A ttl pulse was truncated\n');
                end
            end
            Pulse_dur_ms = (TTLLow - TTLHigh)/TTL_param.Base_ttl_length; % duration of each pulse in ms
            InterPulse_dur_ms = (TTLHigh(2:end) - TTLLow(1:end-1))/TTL_param.Base_ttl_length; % duration between each pulse in ms (length(InterPulse_dur_ms)=length(Pulse_dur_ms)-1)
            PulseTrainInd = [1; find(round(InterPulse_dur_ms) ~= TTL_param.IPI)+1]; % these are the indices of the elements in Pulse_dur_ms that correspond to the first digit coded by pulse trains
            if unique(round(InterPulse_dur_ms(PulseTrainInd(2:end)-1)/10^3)) ~= TTL_param.IPTI
                fprintf(1, 'Observed inter-train of pulses intervals in seconds:\n')
                unique(round(InterPulse_dur_ms(PulseTrainInd(2:end)-1)/10^3))
                error('trains of pulses are not spaced by the expected InterPulseTrain interval of %d seconds\n', TTL_param.IPTI);
            end
            % Check that PulseTrainInd is what we expect to be
            NPulses = length(PulseTrainInd);
            if any(TTLHigh(PulseTrainInd)' ~= TTLHigh(1) + FS*TTL_param.IPTI*(0:(NPulses-1)))
                error('The TTL pulses are not correctly detected, they are not where we expect them to be!\n')
            end
            % Now extract the pulses indices coded in the trains of pulses'
            % durations
            Digits = cell(NPulses,1);
            for pp=1:NPulses
                if pp==NPulses % special case for the last pulse
                    Digits{pp} = int2str(round(Pulse_dur_ms(PulseTrainInd(pp) :end)' - TTL_param.Min_ttl_length));
                else
                    Digits{pp} = int2str(round(Pulse_dur_ms(PulseTrainInd(pp):PulseTrainInd(pp+1)-1)' - TTL_param.Min_ttl_length));
                end
            end
            Pulse_idx_audio{w} = cell2mat(cellfun(@(X) str2double(regexprep(X, ' ','')), Digits, 'UniformOutput',0));
        elseif strcmp(TTL_pulse_generator, 'Avisoft')
            TTLHigh = find(abs(diff(Ttl_status))>0.2)+1; % identify the increases in volatge
            TTLHigh((find(diff(TTLHigh)==1))+1)=[]; % eliminate consecutive points that show a large increase in the TTL pulses, they are just the continuity of a single pulse start (voltage going up)
            TTLHigh((find(diff(TTLHigh)==2))+1)=[]; % eliminate almost consecutive points that show a large increase in the TTL pulses, they are just the continuity of a single pulse start (voltage going up)
            DiffTTL = round(diff(TTLHigh)/FS*TTL_param.Base_ttl_length);
            TTLHigh(find(DiffTTL<TTL_param.Min_ttl_length))=[]; % eliminate detection artefacts
            DiffTTL = round(diff(TTLHigh)/FS*TTL_param.Base_ttl_length);
            PulseTrainInd = [1; find(round(DiffTTL/TTL_param.Base_ttl_length)>=(TTL_param.IPTI-0.5))+1];% these are the indices of the elements in TTLHigh that correspond to the first digit coded by pulse trains
            
            
%             Pulse_dur_ms = DiffTTL(DiffTTL<TTL_param.IPI);% duration of each pulse in ms
            InterPulse_dur_ms = DiffTTL(DiffTTL>=TTL_param.IPI);% duration between each pulse in ms (length(InterPulse_dur_ms)=length(Pulse_dur_ms)-1)
            Pulse_dur_TrainInd = [1; intersect(find(round(InterPulse_dur_ms) ~= TTL_param.IPI),find(round(InterPulse_dur_ms/TTL_param.Base_ttl_length) == TTL_param.IPTI))+1]; % these are the indices of the elements in Pulse_dur_ms that correspond to the first digit coded by pulse trains
            if unique(round(InterPulse_dur_ms(Pulse_dur_TrainInd(2:end)-1)/10^3)) ~= TTL_param.IPTI
                fprintf(1, 'Observed inter-train of pulses intervals in seconds:\n')
                unique(round(InterPulse_dur_ms(Pulse_dur_TrainInd(2:end)-1)/10^3))
                error('trains of pulses are not spaced by the expected InterPulseTrain interval of %d seconds\n', TTL_param.IPTI);
            end
            % Now extract the pulses indices coded in the trains of pulses'
            % durations
            NPulses = length(PulseTrainInd);
            Digits = cell(NPulses,1);
            for pp=1:NPulses
                if pp==NPulses % special case for the last pulse
                    Pulse_dur_ms = DiffTTL(PulseTrainInd(pp):end);
                else
                    Pulse_dur_ms = DiffTTL(PulseTrainInd(pp):(PulseTrainInd(pp+1)-2));
                end
                Pulse_dur_ms = Pulse_dur_ms(Pulse_dur_ms<TTL_param.IPI);
                Digits{pp} = int2str(round(Pulse_dur_ms)' - TTL_param.Min_ttl_length);
            end
            Pulse_idx_audio{w} = cell2mat(cellfun(@(X) str2double(regexprep(X, ' ','')), Digits, 'UniformOutput',0));
        end
        
        % identify the TTL file ID
        Ind_b = strfind(CurrentFile.name, '_');
        Ind_e = strfind(CurrentFile.name, '.wav');
        Num_Samp_audiofile(w,1)=str2double(CurrentFile.name(Ind_b(end)+1 : Ind_e-1));
        File_number{w} = ones(length(Digits),1)*Num_Samp_audiofile(w,1);
        
        % save the 1st sample of each pulse train
        Pulse_samp_audio{w} = TTLHigh(PulseTrainInd);
    end
end
Pulse_idx_audio = cell2mat(Pulse_idx_audio);
Pulse_samp_audio = cell2mat(Pulse_samp_audio);
File_number = cell2mat(File_number);

% Check if several pulses have the same value, the first one is most
% probably the good one
[~,IP,~] = unique(Pulse_idx_audio);
Duplicates = setdiff(1:length(Pulse_idx_audio),IP);
if sum(Pulse_idx_audio(Pulse_idx_audio(Duplicates))==Pulse_idx_audio(Duplicates))==length(Duplicates)
    Pulse_idx_audio(Duplicates)=[];
    Pulse_samp_audio(Duplicates)=[];
    File_number(Duplicates)=[];
else
    warning('There are some Duplicates in the sound TTL indices that might cause issues!\n');
end

%% extract TTL status change from the event log generated by Deuteron and extracted by extract_logger_data
% find eventfile from all loggers (TTL pulses are logged in
% transceiver time which is common to all loggers but they don't always log the same pulses)
All_loggers = dir(fullfile(Loggers_dir, '*ogger*'));
NLog = length(All_loggers);
Transceiver_time_drise1 = cell(NLog,1);
Transceiver_time_dfall1 = cell(NLog,1);
fprintf('Extract TTL status changes from ')
for ll=1:NLog
    LoggerID = str2double(All_loggers(ll).name(strfind(All_loggers(ll).name, 'r')+1 :end));
    if ~isempty(Logger_list)
        if isempty(intersect(Logger_list, LoggerID))
            fprintf('%s not requested, discarded;', All_loggers(ll).name);
            continue
        end
    end
    fprintf('%s ', All_loggers(ll).name);
    Eventfile = dir(fullfile(All_loggers(ll).folder,All_loggers(ll).name, 'extracted_data', '*_EVENTS.mat')); % load file with TTL status info
    load(fullfile(Eventfile.folder, Eventfile.name), 'event_types_and_details', 'event_timestamps_usec');
    B=find(cellfun(@(x) contains(x,Session_strings{1}),event_types_and_details));
    C=find(cellfun(@(x) contains(x,Session_strings{2}),event_types_and_details));
    if ~isempty(B) || ~isempty(C) % This is the reference logger, extract the informations about the session
        OnsetTime = 1e-3*event_timestamps_usec(B);
        OffsetTime = 1e-3*event_timestamps_usec(C);
        if isempty(OnsetTime) % there was a problem with the onset logging of the experiment let's keep everything
            OnsetTime = 0;
        end
        if isempty(OffsetTime)
            OffsetTime = Inf; % there was a problem with the onset logging of the experiment let's keep everything
        end
    else
        fprintf(1,'No session info ')
    end
        
    Din = find(cellfun(@(x) contains(x,'Digital in'),event_types_and_details));
    Drise1 = find(cellfun(@(x) contains(x,'Digital in rising edge on pin number 1'),event_types_and_details(Din))); % extract which lines in EVENTS correspond to TTL status changes
    Dfall1 = find(cellfun(@(x) contains(x,'Digital in falling edge on pin number 1'),event_types_and_details(Din))); % extract which lines in EVENTS correspond to TTL status changes
%     Ind1=find(diff(Drise1)>2,1,'first');
%     Ind2=find(diff(Drise2)>2,1,'first');
%     Din(Drise1(Ind1))
    Transceiver_time_drise1{ll} = 1e-3*event_timestamps_usec(Din(Drise1)); %#ok<FNDSB> % find times (ms) when TTL status changes to up
    Transceiver_time_dfall1{ll} = 1e-3*event_timestamps_usec(Din(Dfall1)); %#ok<FNDSB> % find times (ms) when TTL status changes to down
end
fprintf('\n')
Transceiver_time_drise1 = unique(cell2mat(Transceiver_time_drise1));
Transceiver_time_dfall1 = unique(cell2mat(Transceiver_time_dfall1));

% restrict the dataset to the time frame we want to look at
if length(OnsetTime)>1 || length(OffsetTime)>1
    if length(OnsetTime)>1
        fprintf(1,'%d Onset times where detected for all voc reward start, indicate which one you want to select:\n', length(OnsetTime))
        for oo=1:length(OnsetTime)
            fprintf(1,'%d. %d\n', oo, OnsetTime(oo));
        end
        OnInd = input('Your choice:');
        OnsetTime = OnsetTime(OnInd);
    end
    if length(OffsetTime)>1
        fprintf(1,'%d Offset times where detected for all voc reward start, indicate which one you want to select:\n', length(OffsetTime))
        for oo=1:length(OffsetTime)
            fprintf(1,'%d. %d\n', oo, OffsetTime(oo));
        end
        OffInd = input('Your choice:');
        OffsetTime = OffsetTime(OffInd);
    end
end    
Transceiver_time_drise1=Transceiver_time_drise1(logical((Transceiver_time_drise1>OnsetTime) .* (Transceiver_time_drise1<OffsetTime)));
Transceiver_time_dfall1=Transceiver_time_dfall1(logical((Transceiver_time_dfall1>OnsetTime) .* (Transceiver_time_dfall1<OffsetTime)));

if strcmp(Method, 'rise')
    % hypothesised indices in Transceiver_time_drise of the onset timestamp of each reported Pulse train
    IndPulse = [1; find((mod(round(diff(Transceiver_time_drise1)/10^3),TTL_param.IPTI)==0).*round(diff(Transceiver_time_drise1)/10^3))+1];
    
    % Identify missing rising edges trains or missing rising edges in pulse trains
    MissingPulsesIdx = find(round(diff(Transceiver_time_drise1(IndPulse)).*10^-3)./TTL_param.IPTI -1); % This is the indices of received rising edge trains after which there are some missing rising edges
    MissingPulsesNum = round(diff(Transceiver_time_drise1(IndPulse)).*10^-3)./TTL_param.IPTI -1; % This is the number of missing rising edge trains
    MissingPulsesNum = MissingPulsesNum(MissingPulsesIdx);
    
    % Initialize vectors
    Pulse_Idx_Transc = cell(length(IndPulse)+sum(MissingPulsesNum),1);
    Pulse_TimeStamp_Transc = nan(length(Pulse_Idx_Transc),1);
    
    % Loop through non-missing pulses and check their value
    Current_pulse = 0;
    for pp=1:length(IndPulse)
        Current_pulse = Current_pulse +1;
        if pp==length(IndPulse)
            PulseTrain_Dur = diff(Transceiver_time_drise1(IndPulse(pp):end));
        else
            PulseTrain_Dur = diff(Transceiver_time_drise1(IndPulse(pp):IndPulse(pp+1)));
        end
        if length(PulseTrain_Dur)==1 % One digit pulse either it's an index below 10 or an error
            Pulse_Idx_Transc{Current_pulse} = '0';
            Pulse_TimeStamp_Transc(Current_pulse) = Transceiver_time_drise1(IndPulse(pp));
        else
            Pulse_Idx_hyp = [int2str(round(PulseTrain_Dur(1:(end-1))-TTL_param.IPI-TTL_param.Min_ttl_length)') '0'];
            if length(Pulse_Idx_hyp)>8 % this is most likely an error
                Pulse_Idx_Transc{Current_pulse} = 'NaN';
            else
                Pulse_Idx_Transc{Current_pulse} = Pulse_Idx_hyp;
            end
            Pulse_TimeStamp_Transc(Current_pulse) = Transceiver_time_drise1(IndPulse(pp));
        end
        % increment the indexing by the number of missing pulses after that
        % one if there are some missing ones
        if any(MissingPulsesIdx == pp)
            Current_pulse = Current_pulse + MissingPulsesNum(find(MissingPulsesIdx==pp));
        end
    end
    EmptyIdx = cellfun('isempty',Pulse_Idx_Transc);
    Pulse_Idx_Transc(EmptyIdx) = {'NaN'};
    Pulse_Idx_Transc = cell2mat(cellfun(@(X) str2double(regexprep(X, ' ','')), Pulse_Idx_Transc, 'UniformOutput',0));
    % Now check that the first pulse of a series of ten is correctly
    % detected at the right spot or replace all the pulses of the series by
    % Nan since the timing would not be picked up correctly
    UniquePulses = unique(Pulse_Idx_Transc);
    for upp = 1:length(UniquePulses)
        IdxfirstIdx = find(Pulse_Idx_Transc == UniquePulses(upp),1,'First');
        if ~isempty(IdxfirstIdx)
            if (UniquePulses(upp) == 0) && (IdxfirstIdx~= 1) % Treat the sppecial case of single digits
                Pulse_Idx_Transc(Pulse_Idx_Transc == UniquePulses(upp)) = NaN;
            elseif (UniquePulses(upp) ~= 0) && (UniquePulses(upp) ~= IdxfirstIdx)
                Pulse_Idx_Transc(Pulse_Idx_Transc == UniquePulses(upp)) = NaN;
            end
        end
    end
    
%     % We might have not catch if some pulse trains were missing right at the begining
%     % Check TTL pulses Idx we're expecting 0 value for the first 9 then 10,
%     % 20...etc
%     TransIndObs = (find(diff(Pulse_Idx_Transc)==10)+1)'; % indices were we reach a new 10
%     TransIndExp = (10:10:(round(length(Pulse_Idx_Transc)/10)*10));
%     if any( TransIndObs - TransIndExp)
%         % some of the consecutive first pulses are missing add some nan or
%         % empty spots according to the number of missing pulse trains
%         NumMissing = 10 - (find(diff(Pulse_Idx_Transc)==10,1,'first')+1);
%         Pulse_Idx_Transc = [nan(NumMissing,1); Pulse_Idx_Transc];
%         Pulse_TimeStamp_Transc = [nan(NumMissing,1); Pulse_TimeStamp_Transc];
%         TransIndObs = (find(diff(Pulse_Idx_Transc)==10)+1)'; % indices were we reach a new 10
%         TransIndExp = (10:10:(round(length(Pulse_Idx_Transc)/10)*10));
%         if any( TransIndObs - TransIndExp) % if this is still not right, error
%             error('In aligh_soundmexAudio_2_logger: Impossible to correct for the missing pulses\n')
%         end
%     end
        
        
    elseif strcmp(Method, 'risefall')
        if length(Transceiver_time_drise1) ~= length(Transceiver_time_dfall1)
            error('The number of rising edge is not the same as the number of falling edges on Deuteron')
        end
        if strcmp(TTL_pulse_generator, 'MOTU')
            Transceiver_time_doff = Transceiver_time_dfall1;
            Transceiver_time_don = Transceiver_time_drise1;
        elseif strcmp(TTL_pulse_generator, 'Avisoft')
            Transceiver_time_doff = Transceiver_time_drise1;
            Transceiver_time_don = Transceiver_time_dfall1;
        end
        Pulse_dur_ms = round(Transceiver_time_doff-Transceiver_time_don); % This is the duration of each pulse
        IP_transc = round(Transceiver_time_don(2:end)-Transceiver_time_doff(1:(end-1))); % This is the interval between pulses
        PulseInd = [1; find(round(IP_transc/1000)>=(TTL_param.IPTI-0.5))+1]; % These are the indices of the first pulse of each pulse train in Transceiver_time_don and Transciever_time_doff
        Pulse_TimeStamp_Transc = Transceiver_time_don(PulseInd);% These are the time onsets of each pulse train
        % Now extract the pulses indices coded in the trains of pulses'
        % durations
        NPulses = length(PulseInd);
        Digits = cell(NPulses,1);
        for pp=1:NPulses
            if pp==NPulses % special case for the last pulse
                Digits{pp} = int2str(round(Pulse_dur_ms(PulseInd(pp) :end)' - TTL_param.Min_ttl_length));
            else
                Digits{pp} = int2str(round(Pulse_dur_ms(PulseInd(pp):PulseInd(pp+1)-1)' - TTL_param.Min_ttl_length));
            end
        end
        Pulse_Idx_Transc = cell2mat(cellfun(@(X) str2double(regexprep(X, ' ','')), Digits, 'UniformOutput',0));
        % Eliminate obvious errors of indices (consecutive indices that
        % have a wrong value)
        ErroneousInd = find(diff(Pulse_Idx_Transc)~=1);
        ErrorInd = ErroneousInd(diff(ErroneousInd)==1)+1;
        Pulse_Idx_Transc(ErrorInd) = [];
        Pulse_TimeStamp_Transc(ErrorInd) = [];
end

%% find common indices and eliminate outlayers
[Pulse_idx, Iaudio, ITransc] = intersect(Pulse_idx_audio, Pulse_Idx_Transc);
Pulse_TimeStamp_Transc = Pulse_TimeStamp_Transc(ITransc);
Pulse_samp_audio = Pulse_samp_audio(Iaudio);
File_number = File_number(Iaudio);
% Outlayers are most likely indices that were not correctly decoded
% resulting in obvious error in the Clock drift report
AbsDiff_PTST = abs(diff(Pulse_TimeStamp_Transc));
Outsider_diff = find((AbsDiff_PTST > (nanmean(AbsDiff_PTST)+ 4*nanstd(AbsDiff_PTST))) + (AbsDiff_PTST < nanmean(AbsDiff_PTST - 4*nanstd(AbsDiff_PTST)))); % identify indices of the derivative of Pulse_TimeStamp_Transc that are 4 standard deviation away from the mean
% consecutive indices of the derivative that are away from the average
% distribution correspond to outsider points that we can eliminate.
Outsider_local = Outsider_diff(find(diff(Outsider_diff)==1)+1); % identify consecutive indices of the derivative that are 4 standard deviation away from the mean derivative
Outsider_idx = Pulse_idx(Outsider_local); % % Indices of TTL pulses that are discarted
Outsider_idx_PTST=Pulse_TimeStamp_Transc(Outsider_local);
Outsider_idx_PSA = Pulse_samp_audio(Outsider_local);
Outsider_FileNumber = File_number(Outsider_local);
% eliminate the wrong TTL pulses from the data
Pulse_idx(Outsider_local) = [];
Pulse_TimeStamp_Transc(Outsider_local)= [];
Pulse_samp_audio(Outsider_local) = [];
File_number(Outsider_local) = [];




%% Perform a linear fit for each file
Slope_and_intercept = cell(length(TTL_files),1);
Mean_std_x = cell(length(TTL_files),1);
Mean_std_Pulse_TimeStamp_Transc = nan(length(TTL_files),2);
Mean_std_Pulse_samp_audio = nan(length(TTL_files),2);
FileNum_u = unique(File_number);
if length(FileNum_u)~=length(TTL_files)
    fprintf('There is no usable TTL pulse for TTL file %d, this file will not be alligned\n', setdiff(1:length(TTL_files), FileNum_u))
end
for ff=1:length(FileNum_u)
    Pulse_TimeStamp_Transc_local = Pulse_TimeStamp_Transc(File_number == FileNum_u(ff));
    Pulse_samp_audio_local = Pulse_samp_audio(File_number == FileNum_u(ff));
    Mean_std_Pulse_TimeStamp_Transc(ff,1) = nanmean(Pulse_TimeStamp_Transc_local);
    Mean_std_Pulse_TimeStamp_Transc(ff,2) = nanstd(Pulse_TimeStamp_Transc_local);
    Pulse_TimeStamp_Transc_localzscore = (Pulse_TimeStamp_Transc_local - Mean_std_Pulse_TimeStamp_Transc(ff,1))/Mean_std_Pulse_TimeStamp_Transc(ff,2);
    Mean_std_Pulse_samp_audio(ff,1) = nanmean(Pulse_samp_audio_local);
    Mean_std_Pulse_samp_audio(ff,2) = nanstd(Pulse_samp_audio_local);
    Pulse_samp_audio_localzscore = (Pulse_samp_audio_local - Mean_std_Pulse_samp_audio(ff,1))/Mean_std_Pulse_samp_audio(ff,2);
    [Slope_and_intercept{ff},~,Mean_std_x{ff}]=polyfit(Pulse_samp_audio_localzscore,Pulse_TimeStamp_Transc_localzscore, 1);
    if ff==1
        F=figure();
    else
        clf(F)
    end
    plot(Pulse_samp_audio_local, Pulse_TimeStamp_Transc_local, 'bo')
    hold on
    xlabel('Audio sample Indices of TTL pulses')
    ylabel('Transceiver time in ms')
    x1 = (min(Pulse_samp_audio_localzscore)-0.5):0.1:(max(Pulse_samp_audio_localzscore)+0.5);
    y1 = polyval(Slope_and_intercept{ff}, x1, [], Mean_std_x{ff});
    x1 = x1 .* Mean_std_Pulse_samp_audio(ff,2) + Mean_std_Pulse_samp_audio(ff,1);
    y1 = y1 .* Mean_std_Pulse_TimeStamp_Transc(ff,2) + Mean_std_Pulse_TimeStamp_Transc(ff,1);
    plot(x1,y1,'r-')
    legend('Observed TTL positions', 'Fit')
    title(sprintf('TTL file %d',FileNum_u(ff)))
    hold off
    pause(1)
    if save_options_parameters_CD_figure
        saveas(F,fullfile(Audio_dir,sprintf('%s_%s_CD_correction_audio_piezo_TTLpositions_file_%d.fig', ExpDate,ExpStartTime,ff)))
    end
end

%% save data to file
save(fullfile(Audio_dir, sprintf('%s_%s_TTLPulseTimes.mat', ExpDate,ExpStartTime)),'Pulse_idx', 'Pulse_TimeStamp_Transc', 'File_number', 'Pulse_samp_audio', 'Slope_and_intercept','Mean_std_x','Mean_std_Pulse_TimeStamp_Transc','Mean_std_Pulse_samp_audio');

%% Plot the drift between clocks
% Get the delay of audio samples due to file changes
FID = fopen(fullfile(AudioEvent_file.folder, AudioEvent_file.name));
Header = textscan(FID, '%s\t%s\t%s\t%s\t%s\t%s\t%s',1);
Data = textscan(FID, '%s\t%f\t%s\t%s\t%f\t%f\t%f');
fclose(FID);
IndFChange = contains(Data{4}(:), 'ChangeFile');
DelayFChange = cumsum(Data{end}(IndFChange));

[~,Isort] = sort(Num_Samp_audiofile(:,1));
CumNum_Samp_audiofile = cumsum(Num_Samp_audiofile(Isort,2));
AudioSamp = Pulse_samp_audio + [zeros(sum(File_number==1),1); CumNum_Samp_audiofile(File_number(File_number>1)-1)];
AudioTime = ((AudioSamp/FS - AudioSamp(1)/FS) + [zeros(sum(File_number==1),1);DelayFChange(File_number(File_number>1)-1)])*10^3;
clock_differences_at_pulses = (Pulse_TimeStamp_Transc - Pulse_TimeStamp_Transc(1)) - AudioTime; % determine difference between transceiver and audio soundcard timestamps when pulses arrived

fprintf('The following %d pulse indices were discarded:\n', length(Outsider_idx))
Outsider_idx
fprintf('Transciever values for these indices:\n')
Outsider_idx_PTST
fprintf('Audio samples for these indices:\n')
Outsider_idx_PSA
fprintf('Audio File numbers for these indices\n')
Outsider_FileNumber
 

figure
hold on
plot(AudioTime,clock_differences_at_pulses,'.-');
xlabel('Incoming Audio Pulse Times')
ylabel('Difference between Piezo clock and soundcard clock in ms');
legend('real clock difference');
legend off
MinY = min(clock_differences_at_pulses);
MaxY = max(clock_differences_at_pulses);
for uu=1:length(FileNum_u)
    line([AudioTime(find(File_number == FileNum_u(uu),1,'first')) AudioTime(find(File_number == FileNum_u(uu),1,'first'))], [MinY MaxY])
end
hold off

if save_options_parameters_CD_figure
    saveas(gcf,fullfile(Audio_dir,sprintf('%s_%s_CD_correction_audio_piezo.fig', ExpDate,ExpStartTime)))
end
end