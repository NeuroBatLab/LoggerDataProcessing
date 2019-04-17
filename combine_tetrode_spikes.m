function [Sample_indices_of_peaks] = combine_tetrode_spikes(Peaks_positions, Peaks_voltage,MinPeakSpacing)
%% This function combine spike positions (sample indices) that were detected on the same electrode bundle...
... and return a unified list of positions (in sample unit) of potential spikes
% MinPeakSpacing  Minimum distance in number of samples between two
%                               detected peaks in the voltage trace to be
%                               retained as 2 distinct spike arrival times,
%                               otherwise the largest is kept. Default = 7
%                               (to be in accordance with the snippet
%                               extraction of 32 samples, 7 before the
%                               peak, 24 after the spike peak).

if nargin<3
    MinPeakSpacing = 7;
end

FIG=0;

Num_channels = length(Peaks_positions);
Num_Peaks_Max = sum(cellfun('length', Peaks_positions));
Peak_positions_out = nan(1,Num_Peaks_Max);
Peak_voltage_out = zeros(1,Num_Peaks_Max);

%% Compile the peak positions on all channels and keep, for each peak position, the voltage value 
% that is the highest of all the channels for that bundle
Num_Peaks_Current = 0;
for channel_i = 1:Num_channels
    [CommonP, I_all, I_i] = intersect(Peak_positions_out, Peaks_positions{channel_i});
    if ~isempty(CommonP)
        Higher_V = Peak_voltage_out(I_all) < Peaks_voltage{channel_i}(I_i);
        Peak_voltage_out(I_all(Higher_V)) = Peaks_voltage{channel_i}(I_i(Higher_V));
    end
    Num_Peaks_New = length(Peaks_positions{channel_i}) - length(CommonP);
    [New_Peaks, I_New_Peaks] = setdiff(Peaks_positions{channel_i}, CommonP);
    Peak_positions_out(Num_Peaks_Current +(1:Num_Peaks_New)) = New_Peaks;
    Peak_voltage_out(Num_Peaks_Current +(1:Num_Peaks_New)) = Peaks_voltage{channel_i}(I_New_Peaks);
    Num_Peaks_Current = Num_Peaks_Current + Num_Peaks_New;
end
Peak_positions_out = Peak_positions_out(1:Num_Peaks_Current);
Peak_voltage_out = Peak_voltage_out(1:Num_Peaks_Current);


%% Find peaks that are separated by a minimum distance; if two peaks are
% within that distance, only find the higher peak they most likely
% correspond to the same spike
% 
% Voltage_peaks = zeros(max(Peak_positions_out)+1,1);
% Voltage_peaks(Peak_positions_out) = Peak_voltage_out;
% [Voltage_of_peaks,Sample_indices_of_peaks]=findpeaks(Voltage_peaks,'MinPeakDistance',MinPeakSpacing);

[Peak_positions_out, peak_sort_idx] = sort(Peak_positions_out);
Peak_voltage_out = Peak_voltage_out(peak_sort_idx);
neighboring_peak_idx = find(diff(Peak_positions_out)<MinPeakSpacing);
round_k = 1;
max_rounds = 10;
while ~isempty(neighboring_peak_idx)
    
    peak_positions_to_keep = true(1,length(Peak_positions_out));
    
    for neighbor_k = 1:length(neighboring_peak_idx)
        neighboring_idxs = [neighboring_peak_idx(neighbor_k) neighboring_peak_idx(neighbor_k)+1];
        [~,min_peak_idx] = min(Peak_voltage_out(neighboring_idxs));
        peak_positions_to_keep(neighboring_idxs(min_peak_idx)) = false;
    end
    if round_k > max_rounds
        disp('max # of rounds reached, choose to continue or not')
        keyboard
        continueFlag = input('continue?');
        if continueFlag
            round_k = 1;
        else
            return
        end
    end
    round_k = round_k + 1;
    
    Peak_positions_out = Peak_positions_out(peak_positions_to_keep);
    Peak_voltage_out = Peak_voltage_out(peak_positions_to_keep);    
    
    neighboring_peak_idx = find(diff(Peak_positions_out)<MinPeakSpacing);
end

Sample_indices_of_peaks = Peak_positions_out;

if FIG
    figure()
    ColorCode = get(groot,'DefaultAxesColorOrder');
    subplot(2,1,1)
    for channel_i = 1:Num_channels
        plot(Peaks_positions{channel_i},Peaks_voltage{channel_i}, '*','Color',ColorCode(channel_i,:))
        hold on
    end
    plot(Peak_positions_out, Peak_voltage_out, 'ko','LineWidth',1, 'MarkerSize',10)
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Highest Voltage')
    hold off
    xlabel('Sample')
    ylabel('Peak voltage')
    
    subplot(2,1,2)
    for channel_i = 1:Num_channels
        plot(Peaks_positions{channel_i},Peaks_voltage{channel_i}, '*','Color',ColorCode(channel_i,:))
        hold on
    end
    plot(Sample_indices_of_peaks, Voltage_of_peaks, 'ko','LineWidth',1, 'MarkerSize',10)
    legend('Channel 1', 'Channel 2', 'Channel 3', 'Channel 4', 'Final Peaks')
    hold off
    xlabel('Sample')
    ylabel('Peak voltage')
end

clear Peaks_positions Peaks_voltage

end