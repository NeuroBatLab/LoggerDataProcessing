% When using extract_Nlg_data.m to extract the voltage trace data from
% Neurologger .DAT file and saving as .mat files, we do not save a time
% stamp for each voltage sample. Instead, we save all information needed to
% calculate the time stamps in the same .mat file as the voltage data, and
% use this function to calculate the time stamps for the samples we're
% interested in.
% 7/12/2016, Wujie Zhang
% Last updated, 9/5/2016, Wujie Zhang
%
% Inputs:
% -sample_indices: the indices of the AD count or voltage samples in a
% single channel, counting from the beginning of the recording, whose time
% stamps we'd like to calculate; can be arrays of any dimensions, and the
% indices don't need to be in increasing order
% -indices_of_first_samples: for a given channel, the indices of the first
% sample of every Nlg .DAT file, counting from the beginning of the
% recording; this is the same for all recording channels; this has been
% saved in the same .mat file as the voltage data by extract_Nlg_data.m
% -timestamps_of_first_samples_usec: the time stamps of the first sample of
% each file for each channel; this has been saved in the same .mat file as
% the voltage data by extract_Nlg_data.m
% -sampling_period_usec: the sampling period in us for the samples in a
% given recording channel; this has been saved in the same .mat file as the
% voltage data by extract_Nlg_data.m
%
% Output:
% -timestamps_usec: an array with the same dimensions as sample_indices,
% whose elements are the time stamps for corresponding samples in
% sample_indices

function timestamps_usec=get_timestamps_for_Nlg_voltage_samples(sample_indices,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)
timestamps_usec=nan(size(sample_indices)); % initialize an array of NaNs the same size as sample_indices
indices_of_first_samples=[indices_of_first_samples inf];
for file_i=1:length(indices_of_first_samples)-1 % for each of the Nlg .DAT files
    logical_indices_samples_in_file=sample_indices>=indices_of_first_samples(file_i) & sample_indices<indices_of_first_samples(file_i+1); % all the samples requested in the current file
    periods_from_first_sample_in_file=sample_indices(logical_indices_samples_in_file)-indices_of_first_samples(file_i); % the number of sampling periods between each requested sample and the first sample in the file
    timestamps_usec(logical_indices_samples_in_file)=timestamps_of_first_samples_usec(file_i)+periods_from_first_sample_in_file*sampling_period_usec; % time stamp of a requested sample = time stamp of the first sample in that file + number of sampling periods since that first sample * sampling period
end
