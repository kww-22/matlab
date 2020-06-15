function [maxMaster,aveMaxMaster] = maxFinder(fileNames,numFiles,numTrials,varRow,saveFile,saveAveFile)
% maxFinder: find maximaum values between first and last event
% *************************************************************************
% Excracts maximum data values from MotionMonitor .exp exports
% and compiles them into a master data table
%
% saves 'maxMaster.csv' & 'aveMaxMaster.csv' to selected directory
%
% Inputs: 
%   fileNames: list of text files in current directory
%   numFiles: total number of text files in current directory
%   numTrials: number of trials per participant
%   VarRow: row containing variable names
%
% Outputs:
%
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
data = extractData(fileNames.fileNames{1},'text',varRow);
opts = detectImportOptions(fileNames.fileNames{1},'FileType','text');

% Calculate number of columns
numVars = width(data);

% Create master table and give it variable names from first file in
% directory
maxMaster = array2table(NaN(numFiles,numVars));
maxMaster.Properties.VariableNames = opts.VariableNames;

%% Populate master table with data from individual trials

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',varRow);
    maxData = NaN(1,numVars);
    
    trialRange = find(data.VEM_0 == 1);
    trialData = data(trialRange(1):trialRange(end),:);
    maxData = max(trialData{:,:});

    maxMaster{i,:} = maxData;
end
%% Append sorted file names to beginning of table

files = sort(fileNames.fileNames);
maxMaster = addvars(maxMaster,files,'before',1);

%% Create average values across n trials

% Create vector that goes from 1 to numFiles every n
pStartRow = 1:numTrials:numFiles;

% Number of participants
numPeeps = length(pStartRow);

% Initialize average master table and give it appropriate column names
aveMaxMaster = array2table(NaN(numPeeps,numVars));

for i = 1:numPeeps
    aveMaxMaster{i,2:end} = mean(maxMaster{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Add every nth trial name to the beginning of avemaster
aveMaxMaster = addvars(aveMaxMaster,files(1:numTrials:length(files)),'before',1);
aveMaxMaster.Properties.VariableNames = maxMaster.Properties.VariableNames;

%% Save maxMaster and aveMaxMaster

writetable(maxMaster,saveFile{:});
writetable(aveMaxMaster,saveAveFile{:});

end

