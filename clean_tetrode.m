function clean_tetrode(OutputPath)
if nargin<1
    error('please provide the path to the files containing\nthe snippets and spike arrival time generated by extract_logger_data\n');
end


Snippets_file = dir(fullfile(OutputPath, '*Tetrode_spikes_snippets*'));
% SpikeTimes_file = dir(fullfile(OutputPath, '*Tetrode_spikes_time*'));
ADBitVolts=800/(32767*10^6); % 32767 is the largest number that can be represented as a signed 16-bit integer, which SpikeSort3D requires to be mapped to the upper bound of the voltage data, which we have picked to be 500 uV here, so that the voltage axis during preview has reasonable scales
   
Num_tetrodes = length(Snippets_file);
for Tetrode_i=1:Num_tetrodes % for each of the electrode bundles, eg. tetrodes
%     load(fullfile(SpikeTimes_file(Tetrode_i).folder, SpikeTimes_file(Tetrode_i).name), 'Spike_arrival_times', 'BatID', 'Date');
    load(fullfile(Snippets_file(Tetrode_i).folder, Snippets_file(Tetrode_i).name), 'Snippets');
    
    warning('off','all')
    TTFile = '/Volumes/server_home/users/JulieE/LMC/LMC_CoEd/logger/20190607/Logger16/extracted_data/59834_20190607_Tetrode_spikes_snippets_T1.mat'    ; 
    load(TTFile)
    Snippets=round(Snippets./(ADBitVolts*10^6)); % the AD counts here, multiplied by the "ADBitVolts" factor above, equals voltages in V
    NumSnipp = size(Snippets,3);
    GoodSnip = ones(NumSnipp,1);
    for sn=1:NumSnipp
        LocalSnips = squeeze(Snippets(:,:,sn));
        fprintf(1,'Snippet %d/%d\t', sn, NumSnipp)
        % Send to trash if there are NaNs in the snippet value
        if any(any(isnan(LocalSnips)))
            GoodSnip(sn) = 0;
            fprintf(1,'Snippets with NaNs\t')
            continue
        end
        
        % Identify the changes of slopes sign
        LocalSnipsSign = diff(sign(diff(LocalSnips)));
        NumCc = size(LocalSnips,2);
        LocalSnipsSign(1,:) = ones(1,NumCc); % assign a value different than zero for the first element of each snippet to make sure we similarly count the number of slopes
        SlopeNum = nan(NumCc,1);
        for cc=1:NumCc
            SlopeNum(cc) = sum(diff(find(LocalSnipsSign(:,cc)))>=3); % This is the number of slopes
        end
        if sum(SlopeNum>=6)==NumCc
            GoodSnip(sn) = 0;
            fprintf(1,'Wrong Slope Number\t')
        end
        
        % Identify the peak frequency of the snippets with the max
        % amplitude
        FS = 32/1000;
        [P,F] = pspectrum(LocalSnips,FS);
        LocalSnipsAmp = max(LocalSnips) - min(LocalSnips);
        [LargestSnipAmp,LargestSnip] = max(LocalSnipsAmp);
        FmaxHigh = nan(NumCc,1);
        for cc=1:NumCc
            [~,Locs] = findpeaks(P(:,cc), 'MinPeakHeight',max(P(:,cc))/5);
            FmaxHigh(cc) = any(F(Locs)>=0.003);
        end
        % There's no way we can have a spike that shows a maximum power
        % higher than 3 periods per ms so higher than 0.003Hz
        if sum(FmaxHigh)==NumCc
            GoodSnip(sn) = 0;
            fprintf(1,'Fmax too high on all channels\t')
        end
        
        % Estimate the number of large peaks and their heights
        Peaks = findpeaks(LocalSnips(:,LargestSnip), 'MinPeakWidth',3, 'MinPeakHeight', max(LocalSnips(:,LargestSnip))-LargestSnipAmp/3,'MinPeakProminence', LargestSnipAmp/3);
        if length(Peaks)>1
            GoodSnip(sn) = 0;
            fprintf(1,'two many high peaks\t')
        end
        
        if isempty(Peaks)
            GoodSnip(sn) = 0;
            fprintf(1,'No Peak!\t')
        end
        
        
        figure(1);
        clf
        if GoodSnip(sn)
            Col = 'k';
        else
            Col = 'r';
        end
        Ylim = nan(NumCc,2);
        Legend = cell(NumCc,1);
        for cc=1:NumCc
            subplot(2,2,cc)
            if cc == LargestSnip
                plot(LocalSnips(:,cc),sprintf('%s-',Col), 'LineWidth',3)
            else
                plot(LocalSnips(:,cc),sprintf('%s-',Col), 'LineWidth',2)
            end
            Ylim(cc,:) = get(gca, 'YLim');
            Legend{cc} = sprintf('Channel %d',cc);
            xlabel(Legend{cc})
        end
        YLim = [min(Ylim(:,1)) max(Ylim(:,2))];
        for cc=1:NumCc
            subplot(2,2,cc)
            set(gca, 'YLim',YLim);
            if cc==1
                title(sprintf('Snippet %d/%d',sn,NumSnipp))
            end
            
        end
        
        
        
        figure(2)
        clf
        plot(F,P,'LineWidth',2);
        legend(Legend)
        xlabel('Frequency (Hz)')
        ylabel('Power')
        fprintf('\n')
%         if ~GoodSnip(sn)
%             pause()
%         end
    end
    warning('on','all')
    
    
    
    save(fullfile(Snippets_file(Tetrode_i).folder, Snippets_file(Tetrode_i).name), 'GoodSnip','-append');
end