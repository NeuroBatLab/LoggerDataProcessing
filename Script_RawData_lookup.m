% exploring the raw data at each eventlog line
function Script_RawData_lookup(ExtData_dir,Date)
if nargin<1
    ExtData_dir = '/Users/elie/Documents/ManipBats/JulieBatsDrive/180711/loggers/logger5/extracted_data';
end
if nargin<2
    Date = '180711';
end
NExMax = 20;
EventfilePath = dir(fullfile(ExtData_dir, sprintf('*%s*EVENTS.mat', Date)));
Eventfile = load(fullfile(EventfilePath.folder, EventfilePath.name));
EventCat = cell(size(Eventfile.event_types_and_details));
EventCatInd = strfind(Eventfile.event_types_and_details,'.');
EventCatInd = cellfun(@(x) x(1), EventCatInd);
for ee=1:length(EventCat)
    EventCat{ee} = Eventfile.event_types_and_details{ee}(1:(EventCatInd(ee)-1));
end
UEventCat = unique(EventCat);
UEventCat = UEventCat(~cellfun('isempty', UEventCat));
OnsetOffset = Eventfile.event_timestamps_usec*ones(1,2).*10^(-3) + [-100 100];
Raw_E = cell(length(UEventCat),1);
IndicesUE = cell(length(UEventCat),1);
for ue=1:length(UEventCat)
    IndicesUE{ue} = find(contains(Eventfile.event_types_and_details, UEventCat{ue}));
    if ~isnan(NExMax)
        IndicesUE{ue}=IndicesUE{ue}(randperm(length(IndicesUE{ue}), min(length(IndicesUE{ue}),NExMax)));
    end
    fprintf('Extracting data for %s events\n', UEventCat{ue})
    [Raw_E{ue},~] = extract_timeslot_LFP_spikes(ExtData_dir, OnsetOffset(IndicesUE{ue},:), 0, NaN, [1 0 0 0]);
    if ue==1
        save(fullfile(ExtData_dir,sprintf('Raw_Data_%s.mat',Date)), 'Raw_E', 'UEventCat','EventCat','IndicesUE')
    else
        save(fullfile(ExtData_dir,sprintf('Raw_Data_%s.mat',Date)), 'Raw_E','IndicesUE','-append') 
    end
    
    Nchannels = size(Raw_E{ue},2);
    for iue = 1:length(IndicesUE{ue})
        F1 = figure(1);
        if Nchannels>1
            for cc=1:Nchannels
                subplot(Nchannels/2,2,cc)
                plot(Raw_E{ue}{iue,cc}, '-k', 'LineWidth',2)
                title(sprintf('Voltage trace Ch %d event %s #%d %d/%d', cc-1,UEventCat{ue},IndicesUE{ue}(iue),iue, length(IndicesUE{ue})))
            end
        else
            plot(Raw_E{ue}{iue,Nchannels}, '-k', 'LineWidth',2)
            title(sprintf('Voltage trace Ch %d event %s #%d %d/%d', Nchannels-1,UEventCat{ue},IndicesUE{ue}(iue),iue, length(IndicesUE{ue})))
        end
        In = input('Type anything to stop plotting events of the same type\n');
        if ~isempty(In)
            clf(F1)
            break
        else
            clf(F1)
        end
    end
end

% Investigate a little further Clock drift reports
IndicesCD = find(contains(Eventfile.event_types_and_details, 'CD='));
if ~isnan(NExMax)
    IndicesCD=IndicesCD(randperm(length(IndicesCD), min(length(IndicesCD),NExMax)));
end
fprintf('Extracting data for Clock Drift events\n')
OnsetOffset = Eventfile.event_timestamps_usec*ones(1,2).*10^(-3) + [-200 100];
[Raw_CD,~] = extract_timeslot_LFP_spikes(ExtData_dir, OnsetOffset(IndicesCD,:), 0, NaN, [1 0 0 0]);
save(fullfile(ExtData_dir,sprintf('Raw_Data_%s.mat',Date)), 'Raw_CD','IndicesCD','-append')
Nchannels = size(Raw_CD,2);
for iue = 1:length(IndicesCD)
    F1 = figure(1);
    if Nchannels>1
        for cc=1:Nchannels
            subplot(Nchannels/2,2,cc)
            plot(Raw_CD{iue,cc}, '-k', 'LineWidth',2)
            title(sprintf('Voltage trace Ch %d event Clock Drift #%d %d/%d', cc-1,IndicesCD(iue),iue, length(IndicesCD)))
        end
    else
        plot(Raw_CD{iue,Nchannels}, '-k', 'LineWidth',2)
        title(sprintf('Voltage trace Ch %d event Clock Drift #%d %d/%d', Nchannels-1,IndicesCD(iue),iue, length(IndicesCD)))
    end
    In = input('Type anything to stop plotting events of the same type\n');
    if ~isempty(In)
        clf(F1)
        break
    else
        clf(F1)
    end
end
        
        