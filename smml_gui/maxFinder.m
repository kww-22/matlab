function [maxMaster,aveMaxMaster] = maxFinder(fileNames,numFiles,...
    numTrials,numEvents,varRow,saveFile,saveAveFile,sum_stats)
% maxFinder: find maximaum values between first and last event
% *************************************************************************
% Extracts maximum data values from MotionMonitor .exp exports
% and compiles them into a master data table
%
% Inputs (all provided by smml_gui): 
%   fileNames: list of text files in current directory
%   numFiles: total number of text files in current directory
%   numTrials: number of trials per participant
%   numEvents: number of events in trial
%   VarRow: row containing variable names
%   saveFile: name of master csv file
%   saveAveFile: name of ave master csv file
%
% Outputs:
%   maxMaster.csv: table containing event data from each individual trial
%   aveMaxMaster.csv: table containing averaged event data across all
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
maxMaster = array2table(NaN(numFiles,numVars));
maxMaster.Properties.VariableNames = opts.VariableNames;

%% Populate master table with data from individual trials

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',varRow);
    maxData = NaN(1,numVars);
    
    trialRange = find(data.VEM_0 == 1);
    if length(trialRange) ~= numEvents
        disp('A trial does not have the correct number of event marks');
        disp('You may see which trials by checking maxMaster in the workspace')
    else
        trialData = data(trialRange(1):trialRange(end),:);
        maxData = max(trialData{:,:});
        maxMaster{i,:} = maxData;
    end
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

% Give aveMaxMaster the same variable names as maxMaster
aveMaxMaster.Properties.VariableNames = maxMaster.Properties.VariableNames;

%% Summary statistics?
if sum_stats == 1

    %% Write summary statistics table for maxMaster
    
    % I don't know what this does but I found it online and it get the
    % variable types for each table column
    varTypes = varfun(@class,eventMaster,'OutputFormat','cell');
    
    % Find elements of varClasses that match "double"
    numericVars = find(varTypes == "double");

    % Trim eventMaster to only include numeric variables
    myVars = maxMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    summary_stats = array2table(summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', maxMaster.Properties.VariableNames(numericVars));

    %% Write Summary Statistics Table for aveMaxMaster
    
    % Trim eventMaster to only include numeric variables
    myVars = aveMaxMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    ave_summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    ave_summary_stats = array2table(ave_summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', aveMaxMaster.Properties.VariableNames(numericVars));
end

%% Save maxMaster and aveMaxMaster with provided file names

writetable(maxMaster,saveFile{:});
writetable(aveMaxMaster,saveAveFile{:});
writetable(summary_stats,'max_sum_stat.csv','WriteRowNames',true);
writetable(ave_summary_stats,'aveMax_sum_stat.csv','WriteRowNames',true);

end

