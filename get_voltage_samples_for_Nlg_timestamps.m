% Once we have an event in logger time and we want to retrieve the logger
% samples corresponding to that time this function converts logger time in
% usec to total samples from beginning of recording in order to index into
% the continuous recording vector that results from extract_logger_data.mat
% 9/20/2018 Maimon Rose
%
% Inputs:
% -requested_timestamps_usec: timestamps in raw logger time (usec from
% midnight).
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
% -requested_sample_idxs: sample numbers corresponding to the requested
% timestamps to index in to the AD_counts vector.

function requested_sample_idxs=get_voltage_samples_for_Nlg_timestamps(requested_timestamps_usec,indices_of_first_samples,timestamps_of_first_samples_usec,sampling_period_usec)

n_requested_timestamps = length(requested_timestamps_usec);
file_idxs = nan(1,n_requested_timestamps);
requested_sample_idxs = nan(1,n_requested_timestamps);

for k = 1:n_requested_timestamps
    file_idxs(k) = find(requested_timestamps_usec(k)>timestamps_of_first_samples_usec,1,'last');
    usec_from_first_sample_in_file = (requested_timestamps_usec(k) - timestamps_of_first_samples_usec(file_idxs(k)));
    requested_sample_idxs(k) = indices_of_first_samples(file_idxs(k)) + round(usec_from_first_sample_in_file /sampling_period_usec);
end


