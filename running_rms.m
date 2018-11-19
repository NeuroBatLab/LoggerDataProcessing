function [Amp_env_voltage]=running_rms(Filtered_voltage_trace, FS_in, Fhigh_power, FS_env, FigFlag)
if nargin<5
    FigFlag=0;
end
if nargin<4
    FS_env = 100; % sampling frequency of the output running RMS set to be 100Hz (1 sample per 10ms)
end
if nargin<3
    Fhigh_power = 20; % frequency value of the low-pass filter on voltage signal power
end

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
Lowpass_filter = fir1(Nfilt, Fhigh_power*2.0/FS_in);
Amp_env_voltage = (filtfilt(Lowpass_filter, 1, Filtered_voltage_trace.^2)).^.5;
% Resample to desired sampling rate
if FS_in ~= FS_env
    Amp_env_voltage = resample(Amp_env_voltage, FS_env, round(FS_in));
end
if FigFlag
    figure()
    subplot(1,2,1)
    plot(Filtered_voltage_trace, '-k')
    hold on
    t2 = (0:(length(Amp_env_voltage)-1))*FS_in/(FS_env);
    plot(t2, Amp_env_voltage, '-r', 'LineWidth',2)
    legend('Filtered raw data', 'Running RMS')
    hold off
    subplot(1,2,2)
    plot(Filtered_voltage_trace.^2, '-k')
    hold on
    plot(t2,Amp_env_voltage.^2, '-r', 'LineWidth',2)
    legend('Filtered data power','Amplitude envelope')
    hold off
end
end