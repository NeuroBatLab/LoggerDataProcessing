function [Sample_indices_of_peaks, Peak_Values, Voltage_Trace, Voltage_dynamic_range]=detect_spikes(Voltage_Trace, DataDeletionOnsetOffset_sample, FS, IndFirstNLastSamples, varargin)
%% This function return the potential spike position in a raw voltage trace as long as the filtered voltage trace
% Voltage_Trace         Raw voltage trace given as a vector.
% FS                            Sample frequency of the voltage trace
% IndFirstNLastSamples
%                               First column corresponds to deuteron original .dat onset sample, second
%                               column to end sample of that file. Raw
%                               voltage should be set to zero for samples in the row indicated by MisingFiles, since
%                               it is missing data.
% SpikeTreshMethod: 'manual' manually set a threshold for detecting spikes
%                               in the raw voltage trace of each channel;
%                               'auto' automatically set a threshold as a multiple of
%                               the estimated standard deviation of the
%                               noise in the voltage trace.(modified from Quian
%                               Quiroga et al., 2004, Neural Computation).
%                               Assume the fluctuations of the filtered voltage during
%                               non-spiking periods (ie. the noise) is normally distributed.
%                               We approximate the median and the 75th percentile of the
%                               noise voltage by the median and 75th percentile of the entire
%                               voltage trace, because spikes constitute a small proportion
%                               of the voltage trace. Then the difference between the 75th
%                               percentile and the median, divided by the 75th percentile of
%                               the standard normal distribution (whose mean is 0 and
%                               standard deviation is 1) is the standard deviation of the
%                               noise distribution. Because the Neurologger voltage
%                               measurements may sometimes have an offset (eg. all voltages
%                               are shifted slightly up from their true values), here we
%                               don't assume the noise distribution is centered around zero,
%                               unlike Quian Quiroga et al. (2004).
%                               Default = 'auto'
%                               
% ManualSpikeThresh: 'Value of the spike threshold in uV if in manual mode.
%                               Default is 40
% AutoSpikeThreshFactor: the threshold is this number times the estimated
%                               noise standard deviation, automatic threshold is selected
%                               in SpikeThreshMethod; Rey et al. (2015, Brain Res Bull)
%                               recommends 3 to 5. Default set to 3.
% BandPassFilter:       lower and upper cut-off frequencies in Hz;
%                               the raw voltage traces will be band-pass filtered to retain
%                               only the frequency components relevant for spike detection and sorting
% MissingFiles           Indices of missing .dat files. Default is empty.
% FigCheck                whether or not to plot some of the filtered voltage traces with
%                               the spike threshold; note: this will pause the processing and
%                               require user input to continue

%% Sorting input arguments
pnames = {'SpikeTreshMethod', 'ManualSpikeThresh', 'AutoSpikeThreshFactor', 'BandPassFilter','MissingFiles','FigCheck','Voltage_dynamic_range'};
dflts  = {'auto', 40,3, [600 6000],[],0,1e4};
[SpikeThreshMeth, ManualSpikeThresh, AutoSpikeThreshFactor, BandPassFilter,  MissingFiles, FigCheck, Voltage_dynamic_range] = internal.stats.parseArgs(pnames,dflts,varargin{:});

%% Bandpass filter the input raw voltage
[b,a]=butter(6,BandPassFilter/(FS/2),'bandpass'); % a 12th order Butterworth band-pass filter; the second input argument is normalized cut-off frequency (ie. normalized to the Nyquist frequency, which is half the sampling frequency, as required by MATLAB)
Voltage_Trace = double(Voltage_Trace);
% Bandpass filtering is applied to continous chunks of recordings.
Chunks = [1 length(Voltage_Trace)];
if ~isempty(MissingFiles)
    Ind_firstNlast_samples_MissingFiles = IndFirstNLastSamples(MissingFiles,:);
    Chunks = [1;Ind_firstNlast_samples_MissingFiles(1:(end-1),2)+1 Ind_firstNlast_samples_MissingFiles(:,1)];
end
% Identify NaNs chunks within the continuous recording
for cc=1:size(Chunks,1)

    DataDeletionOnsetOffset_local = DataDeletionOnsetOffset_sample(DataDeletionOnsetOffset_sample(:,1) > Chunks(cc,1) & DataDeletionOnsetOffset_sample(:,2) < Chunks(cc,2),:);    
    if isempty(DataDeletionOnsetOffset_local)
        SmallChunks = Chunks(cc,:);
    else
        SmallChunks = [[Chunks(cc,1); DataDeletionOnsetOffset_local(:,2)+1], [DataDeletionOnsetOffset_local(:,1)-1; Chunks(cc,2)]];
    end
    
    for ccin=1:size(SmallChunks,1)
        Voltage_Trace(SmallChunks(ccin,1):SmallChunks(ccin,2))=filtfilt(b,a,Voltage_Trace(SmallChunks(ccin,1):SmallChunks(ccin,2))); % band-pass filter the voltage traces
        if ccin > 1
            Voltage_Trace(DataDeletionOnsetOffset_local(ccin-1,1):DataDeletionOnsetOffset_local(ccin-1,2)) = NaN;
        end
    end
end

%% Detect spikes as threshold-crossing by the filtered voltage traces
if strcmp(SpikeThreshMeth, 'manual')
    Spike_threshold=ManualSpikeThresh;
elseif strcmp(SpikeThreshMeth, 'auto')
    Estimate_of_voltage_noise_std=(quantile(Voltage_Trace,0.75,2)-nanmedian(Voltage_Trace,2))/icdf('Normal',0.75,0,1);
    Spike_threshold=AutoSpikeThreshFactor*Estimate_of_voltage_noise_std;
end
 
% If requested, plot the filtered voltage trace along with the calculated
% or chosen threshold
if FigCheck
    figure()
    num_samples_to_plot=500; % number of samples to plot for each plotting window
    for window_i=1:100 % plot some windows sequentially
        clf
        samples_to_plot=(window_i-1)*num_samples_to_plot+1:window_i*num_samples_to_plot;
        times_to_plot=(1:length(samples_to_plot))*(1000/FS); % in ms
        plot(times_to_plot,Voltage_Trace(samples_to_plot),'k')
        hold on
        plot(times_to_plot,repmat(Spike_threshold,1,length(times_to_plot)),'r--')
        ylabel('Channel voltage (uV)')
        xlim(times_to_plot([1 end]))
        xlabel('Time (ms)')
        legend('Filtered voltage trace','Spike threshold', 'Location','northoutside','Orientation','horizontal')
        hold off
        stop_plotting=input('Enter anything to stop plotting: ','s');
        if any(stop_plotting)
            break
        end
    end
end

n_voltage_points = length(Voltage_Trace);
chunkSize = round(FS*1e3);
chunkIdx = [1:chunkSize:n_voltage_points n_voltage_points];

if exist('islocalmax','file')
    peak_idx = false(1,n_voltage_points);
    for k = 1:length(chunkIdx)-1
        local_chunk_idx = chunkIdx(k):chunkIdx(k+1);
        peak_idx(local_chunk_idx) = islocalmax(Voltage_Trace(local_chunk_idx)) & Voltage_Trace(local_chunk_idx) > Spike_threshold;
    end
    
    Sample_indices_of_peaks = find(peak_idx);
    Peak_Values = Voltage_Trace(peak_idx);
    
elseif exist('findpeaks','file')
    Sample_indices_of_peaks = cell(1,length(chunkIdx)-1);
    Peak_Values = cell(1,length(chunkIdx)-1);
    for k = 1:length(chunkIdx)-1
        local_chunk_idx = chunkIdx(k):chunkIdx(k+1);
        [Peak_Values{k},Sample_indices_of_peaks{k}] = findpeaks(Voltage_Trace(local_chunk_idx),'MinPeakHeight',Spike_threshold);
        Sample_indices_of_peaks{k} = Sample_indices_of_peaks{k} + local_chunk_idx(1)-1;
    end
    Sample_indices_of_peaks = [Sample_indices_of_peaks{:}];
    Peak_Values = [Peak_Values{:}];
else
    error('couldn''t find function to find peaks for spike detection')
end

Voltage_Trace = single(Voltage_Trace);

if FigCheck
    figure()
    num_samples_to_plot=500; % number of samples to plot for each plotting window
    for window_i=1:100 % plot some windows sequentially
        clf
        MinSamp = (window_i-1)*num_samples_to_plot+1;
        MaxSamp = (window_i)*num_samples_to_plot;
        samples_to_plot=MinSamp:MaxSamp;
        plot(samples_to_plot,Voltage_Trace(samples_to_plot), 'k-', 'LineWidth',1)
        hold on
        LocalPeaks_Ind = find((Sample_indices_of_peaks>=MinSamp) .* (Sample_indices_of_peaks<=MaxSamp));
        if ~isempty(LocalPeaks_Ind)
            plot(Sample_indices_of_peaks(LocalPeaks_Ind),Peak_Values(LocalPeaks_Ind), 'go', 'MarkerSize',10)
            hold on
        end
        plot(samples_to_plot,Spike_threshold*ones(length(samples_to_plot),1), 'r-')
        if ~isempty(LocalPeaks_Ind)
            legend('Filtered voltage trace', 'Peaks detected',  'Spike threshold', 'Location','northoutside','Orientation','horizontal')
        else
            legend('Filtered voltage trace',  'Spike threshold', 'Location','northoutside','Orientation','horizontal')
        end
        hold off
        xlabel('Samples')
        ylabel('Voltage (uV)')
        stop_plotting=input('Enter anything to stop plotting: ','s');
        if any(stop_plotting)
            break
        end
    end
    
end
