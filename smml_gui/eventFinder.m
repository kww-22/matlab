function [eventMaster,aveEventMaster] = eventFinder(fileNames,numFiles,...
    numTrials,numEvents,varRow,eventParams,eventSort,saveFile,saveAveFile,sum_stats)
% eventFinder: find values at events
% *************************************************************************
% Excracts event data from MotionMonitor .exp exports
% and compiles them into a master data table
%
% saves 'eventMaster.csv' & 'aveEventMaster.csv' to selected directory
%
% Inputs (all provided by smml_gui.m): 
%   fileNames: list of text files in current directory
%   numFiles: total number of text files in current directory
%   numTrials: number of trials per participant
%   numEvents: number of events in trial
%   VarRow: row containing variable names
%   eventParams: event abbreviations
%   eventSort: logical to determine how output files are sorted
%   saveFile: name of master csv file
%   saveAveFile: name of ave master csv file
%
% Outputs:
%   eventMaster.csv: table containing event data from each individual trial
%   aveEventMaster.csv: table containing averaged event data across all
%   trials within each participant
% 
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-07-07
% *************************************************************************
%% Get file dimensions from first file in directory to create master table

% Extract first file in directory and it's import options
data = extractData(fileNames.fileNames{1},'text',varRow);
opts = detectImportOptions(fileNames.fileNames{1},'FileType','text');

% Calculate number of columns
numVars = width(data);

% Create master table and give it variable names from first file in
% directory
eventMaster = array2table(NaN(numFiles*numEvents,numVars));
eventMaster.Properties.VariableNames = opts.VariableNames;
varTypes = opts.VariableTypes;

%% Populate master table with data from individual trials

pStartRow = 1:numEvents:numFiles*numEvents;

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',varRow);
    eventData = NaN(numEvents,numVars);
    
    trialRange = find(data.VEM_0 == 1);
    if length(trialRange) ~= numEvents
        disp('A trial does not have the correct number of event marks');
        disp('You may see which trials by checking maxMaster in the workspace')
    else
        eventData = data{trialRange,:};
        eventMaster{pStartRow(i):pStartRow(i)+numEvents-1,:} = eventData;
    end
end
%% Append sorted file names to beginning of table

files = sort(fileNames.fileNames);
efiles = repelem(files,numEvents);
eventMaster = addvars(eventMaster,efiles,'before',1);

%% Append event abbreviations to table

events = string(eventParams);
repevents = repmat(events,numFiles);
repevents = repevents(:,1);

eventMaster = addvars(eventMaster,repevents,'after',1);

%% Create average values across n trials (3-5 events)

% Create vector that goes from 1 to numFiles every n
pStartRow = 1:numTrials:numFiles;
% Number of participants
numPeeps = length(pStartRow);

% Initialize average master table and give it appropriate column names
aveEventMaster = array2table(nan(height(eventMaster)/numTrials,width(eventMaster)));
if numEvents == 3
% Smaller master sheets with only certain events
e1Master = eventMaster(eventMaster.repevents == events(1),:);
e2Master = eventMaster(eventMaster.repevents == events(2),:);
e3Master = eventMaster(eventMaster.repevents == events(3),:);

% event 1
for i = 1:numPeeps
    aveEventMaster{i,3:end} = mean(e1Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 2
for i = 1:numPeeps
    aveEventMaster{numPeeps+i,3:end} = mean(e2Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 3
for i = 1:numPeeps
    aveEventMaster{2*numPeeps+i,3:end} = mean(e3Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Append file names and event abbreviations to beginning of table

efiles = repmat(files(1:numTrials:length(files)),[numEvents 1]);
repevents = repelem(events,numPeeps);
aveEventMaster = addvars(aveEventMaster,efiles,repevents,'after',2);
aveEventMaster = removevars(aveEventMaster,[1 2]);
aveEventMaster.Properties.VariableNames = eventMaster.Properties.VariableNames;
end

if numEvents == 4
% Smaller master sheets with only certain events
e1Master = eventMaster(eventMaster.repevents == events(1),:);
e2Master = eventMaster(eventMaster.repevents == events(2),:);
e3Master = eventMaster(eventMaster.repevents == events(3),:);
e4Master = eventMaster(eventMaster.repevents == events(4),:);

% event 1
for i = 1:numPeeps
    aveEventMaster{i,3:end} = mean(e1Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 2
for i = 1:numPeeps
    aveEventMaster{numPeeps+i,3:end} = mean(e2Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 3
for i = 1:numPeeps
    aveEventMaster{2*numPeeps+i,3:end} = mean(e3Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 4
for i = 1:numPeeps
    aveEventMaster{3*numPeeps+i,3:end} = mean(e4Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Append file names and event abbreviations to beginning of table

efiles = repmat(files(1:numTrials:length(files)),[numEvents 1]);
repevents = repelem(events,numPeeps);
aveEventMaster = addvars(aveEventMaster,efiles,repevents,'after',2);
aveEventMaster = removevars(aveEventMaster,[1 2]);
aveEventMaster.Properties.VariableNames = eventMaster.Properties.VariableNames;
end

if numEvents == 5
% Smaller master sheets with only certain events
e1Master = eventMaster(eventMaster.repevents == events(1),:);
e2Master = eventMaster(eventMaster.repevents == events(2),:);
e3Master = eventMaster(eventMaster.repevents == events(3),:);
e4Master = eventMaster(eventMaster.repevents == events(4),:);
e5Master = eventMaster(eventMaster.repevents == events(5),:);

% event 1
for i = 1:numPeeps
    aveEventMaster{i,3:end} = mean(e1Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 2
for i = 1:numPeeps
    aveEventMaster{numPeeps+i,3:end} = mean(e2Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 3
for i = 1:numPeeps
    aveEventMaster{2*numPeeps+i,3:end} = mean(e3Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 4
for i = 1:numPeeps
    aveEventMaster{3*numPeeps+i,3:end} = mean(e4Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% event 5
for i = 1:numPeeps
    aveEventMaster{4*numPeeps+i,3:end} = mean(e5Master{pStartRow(i):pStartRow(i)+numTrials-1,3:end},1);
end

% Append file names and event abbreviations to beginning of table

efiles = repmat(files(1:numTrials:length(files)),[numEvents 1]);
repevents = repelem(events,numPeeps);
aveEventMaster = addvars(aveEventMaster,efiles,repevents,'after',2);
aveEventMaster = removevars(aveEventMaster,[1 2]);

% Give aveEventMaster the same variable names as eventMaster
aveEventMaster.Properties.VariableNames = eventMaster.Properties.VariableNames;


eventMaster = renamevars(eventMaster,["efiles" "repevents"],["fileName" "event"]);
aveEventMaster = renamevars(aveEventMaster,["efiles" "repevents"],["fileName" "event"]);

if eventSort == 1
    eventMaster = sortrows(eventMaster,'event');
    aveEventMaster = sortrows(aveEventMaster,'event');
elseif eventSort == 2
    eventMaster = sortrows(eventMaster,'fileName');
    aveEventMaster = sortrows(aveEventMaster,'fileName');
end

end

%% Summary statistics?
if sum_stats == 1

    %% Write summary statistics table for eventMaster

    % Find elements of varClasses that match "double"
    numericVars = find(varTypes == "double");

    % Trim eventMaster to only include numeric variables
    myVars = eventMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    summary_stats = array2table(summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', eventMaster.Properties.VariableNames(numericVars));

    %% Write Summary Statistics Table for aveEventMaster

    % Trim eventMaster to only include numeric variables
    myVars = aveEventMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    ave_summary_stats = array2table(summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', aveEventMaster.Properties.VariableNames(numericVars));
end
%% Save eventMaster and aveEventMaster

writetable(eventMaster,saveFile{:});
writetable(aveEventMaster,saveAveFile{:});
writetable(summary_stats,'events_sum_stat.csv','WriteRowNames',true);
writetable(ave_summary_stats,'aveEvents_sum_stat.csv','WriteRowNames',true);
