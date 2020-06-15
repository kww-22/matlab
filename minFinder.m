function [minMaster,aveMinMaster] = minFinder(path,fileNames,numFiles,numTrials,varRow)
% minFinder: find maximaum values between first and last event
% *************************************************************************
% Excracts maximum data values from MotionMonitor .exp exports
% and compiles them into a master data table
%
% saves 'minMaster.csv' & 'aveMinMaster.csv' to selected directory
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
minMaster = array2table(NaN(numFiles,numVars));
minMaster.Properties.VariableNames = opts.VariableNames;

%% Populate master table with data from individual trials

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',varRow);
    minData = NaN(1,numVars);
    
    trialRange = find(data.VEM_0 == 1);
    trialData = data(trialRange(1):trialRange(end),:);
    minData = min(trialData{:,:});

    minMaster{i,:} = minData;
end
%% Append sorted file names to beginning of table

files = sort(fileNames.fileNames);
minMaster = addvars(minMaster,files,'before',1);

%% Create average values across n trials

% Create vector that goes from 1 to numFiles every n
pStartRow = 1:numTrials:numFiles;

% Number of participants
numPeeps = length(pStartRow);

% Initialize average master table and give it appropriate column names
aveMinMaster = array2table(NaN(numPeeps,numVars));

for i = 1:numPeeps
    aveMinMaster{i,2:end} = mean(minMaster{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Add every nth trial name to the beginning of avemaster
aveMinMaster = addvars(aveMinMaster,files(1:numTrials:length(files)),'before',1);
aveMinMaster.Properties.VariableNames = minMaster.Properties.VariableNames;

%% Save master and avemaster

% Save hand and no hand files with different file names
maxFile = [path 'minMaster.csv'];
aveMaxFile = [path "aveMinMaster.csv"];

saveFile = join(maxFile,'/');
saveAveFile = join(aveMaxFile,'/');

writetable(minMaster,saveFile{:});
writetable(aveMinMaster,saveAveFile{:});

end
