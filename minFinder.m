function [minMaster,aveMinMaster] = minFinder(fileNames,numFiles,...
    numTrials,numEvents,varRow,saveFile,saveAveFile)
% minFinder: find maximaum values between first and last event
% *************************************************************************
% Extracts minimum data values from MotionMonitor .exp exports
% and compiles them into a master data table
%
% Inputs (all provided by smml_gui): 
%   fileNames: list of text files in current directory
%   numFiles: total number of text files in current directory
%   numTrials: number of trials per participant
%   numEvents: number of events in each trial
%   VarRow: row containing variable names
%   saveFile: name of master csv file
%   saveAveFile: name of ave master csv file
% 
% Outputs:
%   minMaster.csv: table containing event data from each individual trial
%   aveMinMaster.csv: table containing averaged event data across all
%   trials within each participan
% 
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-06-25
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
    if length(trialRange) ~= numEvents
        disp('A trial does not have the correct number of event marks');
        disp('You may see which trials by checking minMaster in the workspace')
    else
        trialData = data(trialRange(1):trialRange(end),:);
        minData = min(trialData{:,:});
        minMaster{i,:} = minData;
    end
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

%% Save minMaster and aveMinMaster with provided file names

writetable(minMaster,saveFile{:});
writetable(aveMinMaster,saveAveFile{:});

end
