%% Inputs
Date = '07262018';
% BatIDs = {'65430' '71300' '71335' '71132' '59899' '65430'};
BatIDs = {'71173' '65430' '71300' '71335' '71137' '65430'};
% LoggerNum = [5 6 7 9 10 16];
LoggerNum = [10 5 6 7 9 16];

%% Set paths and dependencies
addpath(genpath('C:\Users\Julie\Documents\GitHub\LoggerDataProcessing'))
Server_path = 'Z:\users\Julie E\GiMo_65430_71300\Loggers';
Local_path = 'C:\Users\Julie\Documents\TempData';
Input_serverfolder=fullfile(Server_path, Date);
Input_localfolder = fullfile(Local_path, Date);

%% Bring data back on the computer
mkdir(Input_localfolder)
[s,m,e]=copyfile(Input_serverfolder, Input_localfolder, 'f');
if ~s
    m
    e
    error('File transfer did not occur correctly for %s\n', Date);
end

%% Run the extraction and transfer data back on the server
Logger_dirs = dir(fullfile(Input_localfolder, 'logger*'));
NLogger = length(Logger_dirs);
Input_dir = cell(1,NLogger);

for logger_k = 1:NLogger
    Input_dir{logger_k} = [Logger_dirs(logger_k).folder filesep Logger_dirs(logger_k).name filesep];
%     extract_logger_data(Input_dir{logger_k}, 'BatID', BatIDs{logger_k}, 'NlxSave',1);
end

%% Bring data back on the server
Remote_dir = cell(1,NLogger);
for logger_k = 1:NLogger
    Remote_dir{logger_k} = fullfile(Input_serverfolder, Logger_dirs(logger_k).name, 'extracted_data');
    mkdir(Remote_dir{logger_k})
    [s,m,e]=copyfile(fullfile(Input_dir{logger_k}, 'extracted_data'), Remote_dir{logger_k}, 'f');
    if ~s
        TicTransfer = tic;
        while toc(TicTransfer)<30*60 || ~s
            [s,m,e]=copyfile(fullfile(Input_dir{logger_k}, 'extracted_data'), Remote_dir{logger_k}, 'f');
        end
        if ~s
            s
            m
            e
            error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', Date);
        end
    end
    if s  %erase local data
        [sdel,mdel,edel]=rmdir(Input_dir{logger_k}, 's');
        if ~sdel
            TicErase = tic;
            while toc(TicErase)<30*60 || ~sdel
                [sdel,mdel,edel]=rmdir(Input_dir{logger_k}, 's');
            end
        end
        if ~sdel
            sdel
            mdel
            edel
            error('File erase did not occur correctly for %s\n Although we tried for 30min\n', Date);
        end
        
        
        
    end
end
