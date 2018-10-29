function [Spike_times, Spike_snippets] = extract_tetrode_snippets(Sample_indices_of_peaks, Data_folder, Channels, Spike_window)
%% This function extracts the spike arrival times in microseconds and the snippets from all the electrode of a given tetrode

% Sample_indices_of_peaks
%                              Position of detected spikes (unit = sample) 
% Data_folder           Folder containing the Data of each electrode
% Channels               ID of the channels. Vector of 4 digits
% Spike_window        2 element vector indicating how many samples should
%                               be saved before and after the position of
%                               the spike peak. Default are [-7 24]
%                               
FigCheck =0;
if nargin<4
    Spike_window = [-7 24];
end

if isempty(Sample_indices_of_peaks)
    Spike_snippets = [];
    Spike_times = [];
    warning('Warning: Input Sample_indices_of_peaks is empty: No spike!\n');
else
    for channel_i = 1:length(Channels)
        FileDir = dir(fullfile(Data_folder,sprintf('*CSC%d.mat', Channels(channel_i))));
        Filename=fullfile(FileDir.folder,FileDir.name);
        % If this is the first channel of the bundle
        % First delete the peaks that would lead to spike waveform
        % windows that extend beyond the recording time
        % and initialize output variable
        % Second convert the indices of spike arrival times to real
        % time since we can load the info regarding time of file
        % onsets, and the sample frequency.
        if channel_i==1
            load(Filename, 'Indices_of_first_and_last_samples');
            Sample_indices_of_peaks(((Sample_indices_of_peaks+Spike_window(1))<1) | ((Sample_indices_of_peaks+Spike_window(2))>Indices_of_first_and_last_samples(end,2)))=[]; %#ok<NODEF>
            Num_spikes=length(Sample_indices_of_peaks);
            Spike_snippets=zeros(sum(abs(Spike_window))+1,length(Channels),Num_spikes);
            load(Filename, 'Estimated_channelFS_Transceiver');
            FS = nanmean(Estimated_channelFS_Transceiver);
            load(Filename, 'Timestamps_of_first_samples_usec');
            Spike_times=round(get_timestamps_for_Nlg_voltage_samples(Sample_indices_of_peaks,Indices_of_first_and_last_samples(:,1)',Timestamps_of_first_samples_usec,10^6/FS)); %#ok<IDISVAR> % the time stamps of all the spike peaks, rounded to integer microseconds; note that these are the time stamps of the last channel on this electrode bundle, which differ from the time stamps on the other channels of this electrode bundle by a few sampling periods of the Nlg AD converter
        end

        % For all channels, collect the spike snippets
        load(Filename, 'Filtered_voltage_trace');
        for spike_i=1:Num_spikes
            Spike_snippets(:,channel_i,spike_i)=Filtered_voltage_trace(Sample_indices_of_peaks(spike_i)+Spike_window(1):Sample_indices_of_peaks(spike_i)+Spike_window(2))'; % save the waveforms of the current spike (units are uV)
        end
    end

    if FigCheck
        figure() %#ok<UNRCH>
        for spike_i = 1:Num_spikes
            clf
            for channel_i = 1:length(Channels)
                subplot(2,2,channel_i)
                plot(Spike_snippets(:,channel_i,spike_i),'LineWidth',2)
                xlabel('samples')
                ylabel('Voltage (uV)')
                %MM = max(abs(Filtered_voltage_trace(Sample_indices_of_peaks)));
                ylim([-100 100])
            end
            stop_plotting=input('Enter anything to stop plotting: ','s');
            if any(stop_plotting)
                break
            end
         end
    end
end
end