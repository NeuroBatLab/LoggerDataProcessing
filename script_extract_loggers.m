Input_folder='Z:\users\Julie E\GiMo_65430_71300\Loggers\07112018';
BatIDs = {'65430' '71300' '71335' '71132' '59899' '65430'};
LoggerNum = [5 6 7 9 10 16];
addpath(genpath('C:\Users\Julie\Documents\GitHub\LoggerDataProcessing'))
Logger_dirs = dir(fullfile(Input_folder, 'logger*'));
NLogger = length(Logger_dirs);
Input_dir = cell(1,NLogger);

for logger_k = 1:NLogger
    Input_dir{logger_k} = [Logger_dirs(logger_k).folder filesep Logger_dirs(logger_k).name filesep];
    extract_logger_data(Input_dir{logger_k}, 'BatID', BatIDs{logger_k}, 'NlxSave',1);
end