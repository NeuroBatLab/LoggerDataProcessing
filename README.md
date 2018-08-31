# LoggerDataProcessing
This repository contains the code written to extract  and align any data recorded with one or several loggers from Deuteron technology.

extract_logger_data is the code that you can run on any logger to extract data from it. If it is a neural logger it will also extract the spikes and save spike arrival times and snippets in a matfile format and, if you specifically ask for it, in neuralynx format (.ntt file) by calling the function mat2ntt.m. Note thate the conversion from matlab variable to neuralynx format only works on Windows and needs a piece of matlab code (Mat2NlxSpike) written by Neuralynx people and that can be downloaded there:https://neuralynx.com/software/category/matlab-netcom-utilities. mat2ntt.m can also be called independantly after running extract_logger_data.

extract_logger_data_par is the same as extract_logger_data except that it is running a parfor loop for all the channels which is useful for neural loggers. Note that it might make your computer crash if you're short on memory....

Once data have been spike sorted on Neuralynx, you can convert back your spike sorted files ([…]SS_xx.ntt) to malt alb format by calling Nlx2MatSpike()
