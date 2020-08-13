function [data] = extractData(file,fileType,varRow)
% extractData: Extracts data from MotionMonitor file and saves it to 
% the workspace
% *************************************************************************
% Inputs:
%   file: file path encased in "*"
%   filetype = 'text' for .exp (Classic) or .txt (xGen) files
%   VarRow = row containing variable names
%       VarRow = 9 for TMM Classic
%       VarRow = 10 for TMM xGen
%
% Outputs:
%   data: table of imported MotionMonitor variables
%
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-08-13
% *************************************************************************
%% Get and set file characteristics
opts = detectImportOptions(file,'FileType',fileType);

% Variable names are located in VarRow
opts.VariableNamesLine = varRow;

% Preserve variable names as is
opts.PreserveVariableNames = true;

% Ignore extra column at end of .exp file
opts.ExtraColumnsRule = 'ignore';

% Data lines start 1 after VarRow
opts.DataLines = [varRow+1 inf];

%%Read in file to variable 'data'
data = readtable(file,opts);

%%Erase temp variables
clear opts
end

