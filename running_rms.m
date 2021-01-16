function [Amp_env_voltage, Power_env]=running_rms(Filtered_voltage_trace, FS_in, Fhigh_power, FS_env, Method, FigFlag)
if nargin<6
    FigFlag=0;
end
if nargin<5
    Method='filter'; % can be set either as filter (uses a low pass filtering appproach) or movrms using the function of matlab
end

if nargin<4
    FS_env = 100; % sampling frequency of the output running RMS set to be 100Hz (1 sample per 10ms)
end
if nargin<3
    Fhigh_power = 100; % frequency value of the low-pass filter on voltage signal power
end
Nfilt = 2^ceil(log2(2*FS_in/Fhigh_power));
Nframes = length(Filtered_voltage_trace);
if Nframes<3*Nfilt
    error('The signal is too short to calculate an enveloppe with such a low frequency for the low-pass filter on signal power Nframes=%d and Nfilt=%d!\n',Nframes,Nfilt)
end

if strcmp(Method, 'filter')
    % Generate filter and filter signal power
    Lowpass_filter = fir1(Nfilt, Fhigh_power*2.0/FS_in);
    Power_env = filtfilt(Lowpass_filter, 1, Filtered_voltage_trace.^2);
    Filtered_env = Power_env;
    Filtered_env(Filtered_env<0) = 0;
    Amp_env_voltage = Filtered_env.^.5;
elseif strcmp(Method, 'movrms')
    Movrms = dsp.MovingRMS(FS_in/FS_env);
    Amp_env_voltage = Movrms(Filtered_voltage_trace);
end


% Resample to desired sampling rate
if FS_in ~= FS_env
    Amp_env_voltage = resample(Amp_env_voltage, FS_env, round(FS_in));
    Power_env = resample(Power_env, FS_env, round(FS_in));
end
if FigFlag
    figure()
    title(sprintf('%s', Method))
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