%% Inputs
Date = '06252018';
% BatIDs = {'65430' '71300' '71335' '71132' '59899' '65430'};
BatIDs = {'60014' '65430' '71300' '71335' '71137' '65430'};
% LoggerNum = [5 6 7 9 10 16];
LoggerNum = [3 5 6 7 9 16];

%% Set paths and dependencies
addpath(genpath('C:\Users\Julie\Documents\GitHub\LoggerDataProcessing'))
Server_path = 'Z:\users\Julie E\GiMo_65430_71300\Loggers';
Local_path = 'C:\Users\Julie\Documents\TempData';
Input_serverfolder=fullfile(Server_path, Date);
Input_localfolder = fullfile(Local_path, Date);

%% Bring data back on the computer
mkdir(Input_localfolder)
[s,m,e]=copyfile(Input_serverfolder, Input_localfolder);
if ~s
    m
    e
    error('File transfer did not occur correctly for %s\n', Date);
end

%% Run
Logger_dirs = dir(fullfile(Input_localfolder, 'logger*'));
NLogger = length(Logger_dirs);
Input_dir = cell(1,NLogger);

for logger_k = 1:NLogger
    Input_dir{logger_k} = [Logger_dirs(logger_k).folder filesep Logger_dirs(logger_k).name filesep];
    extract_logger_data(Input_dir{logger_k}, 'BatID', BatIDs{logger_k}, 'NlxSave',1);
end