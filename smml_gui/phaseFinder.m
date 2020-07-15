function [phaseMaster,avePhaseMaster] = phaseFinder(fileNames,numFiles,...
    numTrials,numEvents,varRow,phaseParams,phaseSort,saveFile,saveAveFile,sum_stats)
% phaseFinder: find average values across phases
% *************************************************************************
% Excracts phase data from MotionMonitor .exp exports
% and compiles them into a master data table
%
% saves 'phaseMaster.csv' & 'avePhaseMaster.csv' to selected directory
%
% Inputs (all provided by smml_gui.m): 
%   fileNames: list of text files in current directory
%   numFiles: total number of text files in current directory
%   numTrials: number of trials per participant
%   numPhases: number of phases in trial
%   VarRow: row containing variable names
%   eventParams: event abbreviations
%   saveFile: name of master csv file
%   saveAveFile: name of ave master csv file
%
% Outputs:
%   phaseMaster.csv: table containing phase data from each individual trial
%   avePhaseMaster.csv: table containing averaged phase data across all
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
phaseMaster = array2table(NaN(numFiles*(numEvents-1),numVars));
phaseMaster.Properties.VariableNames = opts.VariableNames;
varTypes = opts.VariableTypes;

%% Populate master table with data from individual trials

pStartRow = 1:numEvents-1:numFiles*(numEvents-1);

for i = 1:numFiles
    data = extractData(fileNames.fileNames{i},'text',varRow);
    phaseData = NaN(numEvents-1,numVars);
    
    trialRange = find(data.VEM_0 == 1);
    if length(trialRange) ~= numEvents
        disp('A trial does not have the correct number of phase marks');
        disp('You may see which trials by checking maxMaster in the workspace')
    else
        for j = 1:numEvents-1
        phaseData(j,:) = nanmean(data{trialRange(j):trialRange(j+1)-1,:});
        end
        phaseMaster{pStartRow(i):pStartRow(i)+numEvents-2,:} = phaseData;
    end
end
%% Append sorted file names to beginning of table

files = sort(fileNames.fileNames);
efiles = repelem(files,numEvents-1);
phaseMaster = addvars(phaseMaster,efiles,'before',1);

%% Append phase abbreviations to table

phases = string(phaseParams);
repphases = repmat(phases,numFiles);
repphases = repphases(:,1);

phaseMaster = addvars(phaseMaster,repphases,'after',1);

%% Sort and separate phaseMaster for easier indexing of avePhaseMaster

phaseMaster = sortrows(phaseMaster,"repphases");

if numEvents == 3
    phase_1_master = phaseMaster(phaseMaster.repphases == phaseParams{1},:);
    phase_2_master = phaseMaster(phaseMaster.repphases == phaseParams{2},:);
    ave_phase_1_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_2_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
elseif numEvents == 4
    phase_1_master = phaseMaster(phaseMaster.repphases == phaseParams{1},:);
    phase_2_master = phaseMaster(phaseMaster.repphases == phaseParams{2},:);
    phase_3_master = phaseMaster(phaseMaster.repphases == phaseParams{3},:);
    ave_phase_1_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_2_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_3_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
elseif numEvents == 5
    phase_1_master = phaseMaster(phaseMaster.repphases == phaseParams{1},:);
    phase_2_master = phaseMaster(phaseMaster.repphases == phaseParams{2},:);
    phase_3_master = phaseMaster(phaseMaster.repphases == phaseParams{3},:);
    phase_4_master = phaseMaster(phaseMaster.repphases == phaseParams{4},:);
    ave_phase_1_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_2_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_3_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
    ave_phase_4_master = array2table(nan(height(phase_1_master)/numTrials,width(phase_1_master)));
end
%% Create average values across n trials (2-4 phases)

% Create vector that goes from 1 to numFiles every n
pStartRow = 1:numTrials:numFiles;
% Number of participants
numPeeps = length(pStartRow);

if numEvents == 3
    for i = 1:numPeeps
    ave_phase_1_master(i,3:end) = array2table(nanmean(phase_1_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_2_master(i,3:end) = array2table(nanmean(phase_2_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end
    
elseif numEvents == 4
    for i = 1:numPeeps
    ave_phase_1_master(i,3:end) = array2table(nanmean(phase_1_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_2_master(i,3:end) = array2table(nanmean(phase_2_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_3_master(i,3:end) = array2table(nanmean(phase_3_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end
    
elseif numEvents == 5
    for i = 1:numPeeps
    ave_phase_1_master(i,3:end) = array2table(nanmean(phase_1_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_2_master(i,3:end) = array2table(nanmean(phase_2_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_3_master(i,3:end) = array2table(nanmean(phase_3_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end

    for i = 1:numPeeps
    ave_phase_4_master(i,3:end) = array2table(nanmean(phase_4_master{pStartRow(i):pStartRow(i)+numTrials-1,3:end}));
    end
end
%% Combine individual phase masters into grand phase master

phaseMaster = renamevars(phaseMaster,["efiles","repphases","Frame_"],...
["participant","phase","frame"]);

aveFileNames = phaseMaster.participant(pStartRow);
aveEventNames = phaseMaster.phase(1:numTrials:height(phaseMaster));

avePhaseMaster = [ave_phase_1_master ; ave_phase_2_master ; ave_phase_3_master ; ave_phase_4_master];


% Give avePhaseMaster the same variable names and properties as phaseMaster
avePhaseMaster.Properties.VariableNames = phaseMaster.Properties.VariableNames;

avePhaseMaster.phase = aveEventNames;
avePhaseMaster.participant = repmat(aveFileNames,4,1);

%% Sort tables

if phaseSort == 1
    phaseMaster = sortrows(phaseMaster,"phase");
    avePhaseMaster = sortrows(avePhaseMaster,"phase");
else
    phaseMaster = sortrows(phaseMaster,"participant");
    avePhaseMaster = sortrows(avePhaseMaster,"participant");
end

%% Summary statistics?
if sum_stats == 1

    %% Write summary statistics table for phaseMaster

    % I don't know what this does but I found it online and it get the
    % variable types for each table column
    varTypes = varfun(@class,phaseMaster,'OutputFormat','cell');
    
    % Find elements of varClasses that match "double"
    numericVars = find(varTypes == "double");

    % Trim eventMaster to only include numeric variables
    myVars = phaseMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    summary_stats = array2table(summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', phaseMaster.Properties.VariableNames(numericVars));

    %% Write Summary Statistics Table for avePhaseMaster

    % Trim eventMaster to only include numeric variables
    myVars = avePhaseMaster(:,numericVars);

    % Compute common descriptive statistics for numeric variables
    means = mean(myVars{:,:})';
    std_dev = std(myVars{:,:})';
    quarts = prctile(myVars{:,:},[25 50 75])';

    % Combine descriptives into one array
    ave_summary_stats = [means std_dev quarts];

    % Convert descriptive array to table
    ave_summary_stats = array2table(ave_summary_stats,...
        'VariableNames',{'mean','std','quart_25','median','quart_75'},...
        'RowNames', avePhaseMaster.Properties.VariableNames(numericVars));
    
    % Write summary stats tables
    writetable(summary_stats,'phase_sum_stat.csv','WriteRowNames',true);
    writetable(ave_summary_stats,'avePhase_sum_stat.csv','WriteRowNames',true);
end
%% Save phaseMaster and avePhaseMaster

writetable(phaseMaster,saveFile{:});
writetable(avePhaseMaster,saveAveFile{:});

end

