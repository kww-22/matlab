function [firstRowMaster,aveFirstRowMaster] = grabFirstRow(fileNames,numFiles,numTrials,VarRow,saveFile,saveAveFile)
% grabFirstRow: Extract first data row from MotionMonitor .exp files
% *************************************************************************
% Grabs first row of data from MotionMonitor .exp exports
% and compile them into a master data table
%
% Inputs (provided by firstRowExtract.m):
%   fileNames: table of file names in directory
%   numFiles: total number of files in directory
%   numTrials: number of trials per participant
%   VarRow: table row with variable names
%   saveFile: file path for saved master table
%   saveAveFile: file path for save average master table
%
% Outputs:
%   firstRowMaster: compiled master table of each trial
%   aveFirstRowMaster: compiled mater table of averages within each
%   participant
%
% Saves '*.csv' and 'ave*.csv' to specified directory
% 
% Designed for Bordelon et al manuscript examining mechanical energy
% flow through the kinetic chain during softball hitting
%
% Requires in directory: 
%   1) exported files 
%   2) extractData.m
% 
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-06-14
% *************************************************************************
%% Get file dimensions from first file in directory to create master table

% Extract first file in directory and it's import options
data = extractData(fileNames.fileNames{1},'text',VarRow);
opts = detectImportOptions(fileNames.fileNames{1},'FileType','text');

% Calculate number of columns
numVars = width(data);

% Create master table and give it variable names from first file in
% directory
firstRowMaster = array2table(NaN(numFiles,numVars));
firstRowMaster.Properties.VariableNames = opts.VariableNames;

%% Populate master table with data from individual trials

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',VarRow);
    firstRowMaster{i,:} = data{1,:};
end

%% Adds sorted file names to first column of master

files = sort(fileNames.fileNames);
firstRowMaster = addvars(firstRowMaster,files,'before',1);

%% Create average values across three trials

% Create vector that goes from 1 to numFiles every x
pStartRow = 1:numTrials:numFiles;

% Number of participants
numPeeps = length(pStartRow);

% Initialize average master table and give it appropriate column names
aveFirstRowMaster = array2table(NaN(numPeeps,numVars));
aveFirstRowMaster.Properties.VariableNames = opts.VariableNames;

for i = 1:numPeeps
    aveFirstRowMaster{i,2:end} = mean(firstRowMaster{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Add every third trial name to the beginning of avemaster
aveFirstRowMaster = addvars(aveFirstRowMaster,files(1:3:length(files)),'before',1);

%% Save master and avemaster

% Save hand and no hand files with different file names
writetable(firstRowMaster,saveFile{:});
writetable(aveFirstRowMaster,saveAveFile{:});
end

