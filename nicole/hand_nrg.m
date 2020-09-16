%% Gather directory information

clear
clc

% set directory
% choose folder path containing exported and necessary matlab files
path = uigetdir();
cd(path);

% isolate files with .exp extension
fileDir = dir('*.exp');

% isolate file names from directory structure
fileNames = {fileDir.name}.';

% convert cell array to table
fileNames = cell2table(fileNames);

% get number of files in directory
numFiles = height(fileNames);

% get table properties
opts = detectImportOptions(fileNames.fileNames{1},'filetype','text');
varNames = opts.VariableNames;

websave('extractData.m',...
    'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');

%%

%% read in data, trim trial, and find events
for i = 1:numFiles
    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});
    clear txt
    % load i-th file from directory
    data = extractData(fileNames.fileNames{i},'text',9);
    
    % get height, weight, sampling frequency
    bodyheight(i,1) = data.BodyHeight(1);
    bodymass(i,1) = data.BodyMass(1);
    measRate(i,1) = fs;
    
    % find event indices
    events = find(data.VEM_0 == 1);
    
    data = data(events(1):events(end),:);
    events = find(data.VEM_0 == 1);
    
    bHand_sp = data.Bhand_JfpP + data.Bhand_StpP;
    fHand_sp = data.Fhand_JfpP + data.Fhand_StpP;
    
    bHand_IF_12(i,1) = (trapz(bHand_sp(events(1):events(2)) + abs(bHand_sp(events(1):events(2)))))/ (2*fs);
    bHand_IF_23(i,1) = (trapz(bHand_sp(events(2):events(3)) + abs(bHand_sp(events(2):events(3)))))/ (2*fs);
    bHand_IF_34(i,1) = (trapz(bHand_sp(events(3):events(4)) + abs(bHand_sp(events(3):events(4)))))/ (2*fs);
    bHand_IF_45(i,1) = (trapz(bHand_sp(events(4):events(5)) + abs(bHand_sp(events(4):events(5)))))/ (2*fs);

    fHand_IF_12(i,1) = (trapz(fHand_sp(events(1):events(2)) + abs(fHand_sp(events(1):events(2)))))/ (2*fs);
    fHand_IF_23(i,1) = (trapz(fHand_sp(events(2):events(3)) + abs(fHand_sp(events(2):events(3)))))/ (2*fs);
    fHand_IF_34(i,1) = (trapz(fHand_sp(events(3):events(4)) + abs(fHand_sp(events(3):events(4)))))/ (2*fs);
    fHand_IF_45(i,1) = (trapz(fHand_sp(events(4):events(5)) + abs(fHand_sp(events(4):events(5)))))/ (2*fs);
    
    pbHand_IF_12(i,1) = max(bHand_sp(events(1):events(2)));
    pbHand_IF_23(i,1) = max(bHand_sp(events(2):events(3)));
    pbHand_IF_34(i,1) = max(bHand_sp(events(3):events(4)));
    pbHand_IF_45(i,1) = max(bHand_sp(events(4):events(5)));
    
    pfHand_IF_12(i,1) = max(fHand_sp(events(1):events(2)));
    pfHand_IF_23(i,1) = max(fHand_sp(events(2):events(3)));
    pfHand_IF_34(i,1) = max(fHand_sp(events(3):events(4)));
    pfHand_IF_45(i,1) = max(fHand_sp(events(4):events(5)));
end