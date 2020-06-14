function [data] = extractData(file,fileType,VarRow)
% extractData: Extracts data from MotionMonitor .exp file and saves it to the
% workspace
% *************************************************************************
% file: file name encased in "*"
% filetype = 'text' for .exp or .txt files
% VarRow = row containing variable names
%% Get and set file characteristics
opts = detectImportOptions(file,'FileType',fileType);

% Variable names are located in VarRow
opts.VariableNamesLine = VarRow;

% Preserve variable names as is
opts.PreserveVariableNames = true;

% Ignore extra column at end of .exp file
opts.ExtraColumnsRule = 'ignore';

% Data lines start 1 after VarRow
opts.DataLines = [VarRow+1 inf];

%%Read in file to variable 'data'
data = readtable(file,opts);

%%Erase temp variables
clear opts
end

