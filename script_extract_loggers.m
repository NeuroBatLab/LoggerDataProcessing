function script_extract_loggers(Date, BatIDs, LoggerNums, TransferData)
if nargin<4
    TransferData = 1;
end
%% Inputs
%Date = '06272018';
% BatIDs = {'65430' '71300' '71335' '71132' '59899' '65430'};
% BatIDs = {'65430'};
% BatIDs = {'59899' '65430' '71300' '71335' '71137' '65430'};
% LoggerNums = 16;
% LoggerNum = [5 6 7 9 10 16];
% LoggerNum = [10 5 6 7 9 16];

%% Set paths and dependencies
addpath(genpath('C:\Users\Julie\Documents\GitHub\LoggerDataProcessing'))
Server_path = 'Z:\users\Julie E\GiMo_65430_71300\Loggers';
if strcmp('BATMAN-PC', getenv('computername'))
    Local_path = 'C:\Users\Batman\Documents\AnalysisFolder\';
elseif strcmp('TADARIDA', getenv('computername'))
    Local_path = 'C:\Users\Julie\Documents\TempData';
else
    Local_path = input('What is the path to the folder containing the data?\n','s');
end
Input_serverfolder=fullfile(Server_path, Date);
Input_localfolder = fullfile(Local_path, Date);

% Bring data back on the computer
if TransferData
    fprintf(1,'Transferring data from the server %s\n on the local computer %s\n', Input_serverfolder, Input_localfolder);
    mkdir(Input_localfolder)
    for logger_k = 1:length(LoggerNums)
        fprintf('-> Logger%d\n',LoggerNums(logger_k))
        Input_localfolder_logger = fullfile(Input_localfolder, sprintf('logger%d',LoggerNums(logger_k)));
        Input_serverfolder_logger = fullfile(Input_serverfolder, sprintf('logger%d',LoggerNums(logger_k)));
        mkdir(Input_localfolder_logger)
        [s,m,e]=copyfile(Input_serverfolder_logger, Input_localfolder_logger, 'f');
        if ~s
            m
            e
            error('File transfer did not occur correctly for %s\n', Date);
        end
    end
end

%% Run the extraction and transfer data back on the server
Logger_dirs = dir(fullfile(Input_localfolder, 'logger*'));
NLogger = length(Logger_dirs);
Input_dir = cell(1,NLogger);

for logger_k = 1:NLogger
    LoggerNumInd = strfind(Logger_dirs(logger_k).name, 'er');
    LoggerNum = str2double(Logger_dirs(logger_k).name(LoggerNumInd+2 : end));
    if ~ isempty(intersect(LoggerNum, LoggerNums))
        LogInd = find(LoggerNums==LoggerNum);
        Input_dir{logger_k} = [Logger_dirs(logger_k).folder filesep Logger_dirs(logger_k).name filesep];
        extract_logger_data(Input_dir{logger_k}, 'BatID', BatIDs{LogInd}, 'NlxSave',1);
    end
end

%% Bring data back on the server
if TransferData
    fprintf(1,'Transferring data from the local computer %s\n back on the server %s\n', Input_localfolder, Input_serverfolder);
    
    Remote_dir = cell(1,NLogger);
    for logger_k = 1:NLogger
        LoggerNumInd = strfind(Logger_dirs(logger_k).name, 'er');
        LoggerNum = str2double(Logger_dirs(logger_k).name(LoggerNumInd+2 : end));
        if ~ isempty(intersect(LoggerNum, LoggerNums))
            Remote_dir{logger_k} = fullfile(Input_serverfolder, Logger_dirs(logger_k).name, 'extracted_data');
            mkdir(Remote_dir{logger_k})
            [s,m,e]=copyfile(fullfile(Input_dir{logger_k}, 'extracted_data'), Remote_dir{logger_k}, 'f');
            if ~s
                TicTransfer = tic;
                while toc(TicTransfer)<30*60
                    [s,m,e]=copyfile(fullfile(Input_dir{logger_k}, 'extracted_data'), Remote_dir{logger_k}, 'f');
                    if s
                        return
                    end
                end
                if ~s
                    s
                    m
                    e
                    error('File transfer did not occur correctly for %s\n Although we tried for 30min\n', Remote_dir{logger_k});
                else
                    fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir{logger_k});
                end
            else
                fprintf('Extracted data transfered back on server in:\n%s\n',  Remote_dir{logger_k});
            end
            if s  %erase local data
                [sdel,mdel,edel]=rmdir(Input_dir{logger_k}, 's');
                if ~sdel
                    TicErase = tic;
                    while toc(TicErase)<30*60
                        [sdel,mdel,edel]=rmdir(Input_dir{logger_k}, 's');
                        if sdel
                            return
                        end
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
    end
end
