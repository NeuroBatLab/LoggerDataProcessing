function plot_spike_snippets(Spike_arrival_times, Snippets, FS)
% Spike_arrival_times: spike arrival times in microseconds
% Snippets: 3D matirx, dim 1 = extracted snippets of length
% length(WinSpike), dim2 = channel #, dim3 = spike #.
% FS= sample frequency of the voltage trace
Delay = 2;
WinSpike = [-7 24];

% change the unit of Spike_arrival_times from us to ms
Spike_arrival_times = Spike_arrival_times*10^-3;

% Make sure Spike_arrival_times is sorted
[Spike_arrival_times, OrdI] = sort(Spike_arrival_times);
Snippets = Snippets(:,:,OrdI);

NChannels = size(Snippets,2);

% Plot the snippets of consecutive peaks that are spaced by less than Delay
% ms together
LargeDelay = [find(diff(Spike_arrival_times)>Delay) length(Spike_arrival_times)]; % Find the inter-spike interval that are larger than the Delay (these will be the limits of each plot) 
SeqInd = find(diff(LargeDelay>1)); % first focus on only plotting sequences of spikes that are less than Delay ms appart
% randomize SeqInd so we don't plot in the temporal order
SeqIndrand = randperm(length(SeqInd)-1);
for SeqI = length(SeqIndrand)
    ConsSnipInd = LargeDelay(SeqInd(SeqIndrand(SeqI))):LargeDelay(SeqInd(SeqIndrand(SeqI)+1));
    F50=figure(50);
    for ss = 1:length(ConsSnipInd)
        spike = ConsSnipInd(ss);
        for cc=1:length(Channels)
            subplot(NChannels,1,cc)
            hold on
            plot(Spike_arrival_times(spike)*FS*10^-3 + WinSpike, Snippets(:,cc,spike), 'k-','LineWidth',2)
            hold on
            plot(Spike_arrival_times(spike)*FS*10^-3, Snippets(8,cc,spike), 'r.','MarkerSize',5)
        end
    end
    Resp = input(sprintf('Keep plotting sequence of spikes that are less than %d ms apart? y/n\n', Delay), 's');
    if strcmp(Resp, 'n')
        break
    end
    clf(F50)
end
    
% plot isolated spikes
Iso = setdiff(1:length(LargeDelay),SeqInd);
Iso = Iso(ranperm(length(Iso)));
F51=figure(51);
for ss=1:length(Iso)
    spike = Iso(ss);
    for cc=1:length(Channels)
        subplot(NChannels,1,cc)
        hold on
        plot(Spike_arrival_times(spike)*FS*10^-3 + WinSpike, Snippets(:,cc,spike), 'k-','LineWidth',2)
        hold on
        plot(Spike_arrival_times(spike)*FS*10^-3, Snippets(8,cc,spike), 'r.','MarkerSize',5)
    end
    Resp = input(sprintf('Keep plotting isolated spikes that are more than %d ms apart? y/n\n', Delay), 's');
    if strcmp(Resp, 'n')
        break
    end
    clf(F51)
end