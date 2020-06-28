% smml_gui: interactive script to handle basic data cleaning/organizing
% from MotionMonitor
% *************************************************************************
% Gives user a series of interactive prompts to select which data
% cleaning/organizing processes they would like. Options include:
%       1. Max of variables btwn first & last event
%       2. Min of variables btwn first & last event
%       3. Values of variables @ events
%       4. Average value of variables between events (over phases) (in
%       progress)
%
% User may also provide custom event labels and custom outpus file names
% 
% Inputs:
%   NA
%
% Outputs:
%   depending on which processes user selects saves various .csv master and 
%   average master files to selected directory
%
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-06-25
% *************************************************************************
%% Gather directory information

% Set directory
% Choose folder path containing exported and necessary matlab files
path = uigetdir('Navigate to the folder containing MotionMonitor Exports');
cd(path);

% Isolate files with .exp extension
fileDir = dir('*.exp');

% Isolate file names from directory structure
fileNames = {fileDir.name}.';

% Convert cell array to table
fileNames = cell2table(fileNames);

% Get number of files in directory
numFiles = height(fileNames);

%% Interactive prompt to collect trial parameters
prompt = {'How many trials per participant?',...
    'How many events?'...
    'Which row contains the variable names?'};
dlgtitle = 'Define trial parameters';
dims = [1 45];
defParams = {'3','4','9'};
trialParams = inputdlg(prompt,dlgtitle,dims,defParams);

numTrials = str2double(trialParams{1});
numEvents = str2double(trialParams{2});
varRow = str2double(trialParams{3});

%% Values at events?

% Builds yes/no selection box. Yes: valEvent = 1; No: valEvent = 2
valEvent = listdlg('PromptString', 'Value of variables at events?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Yes', 'No'}, ...
    'Name', 'Values at Events', ...
    'ListSize', [200 100]);

%% Maxes?

% Builds yes/no selection box. Yes: valMax = 1; No: valMax = 2
valMax = listdlg('PromptString', 'Max of variables between first and last event?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Yes', 'No'}, ...
    'Name', 'Values at Events', ...
    'ListSize', [225 100]);

%% Mins?

% Builds yes/no selection box. Yes: valMin = 1; No: valMin = 2
valMin = listdlg('PromptString', 'Min of variables between first and last event?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Yes', 'No'}, ...
    'Name', 'Values at Events', ...
    'ListSize', [225 100]);

%% Provide output file names

% Events, maxes, and mins
if valEvent == 1 && valMax == 1 && valMin == 1
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','minMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Events and Maxes
elseif valEvent == 1 && valMax == 1 && valMin == 2
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Events and Mins    
elseif valEvent == 1 && valMax == 2 && valMin == 1
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','minMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Maxes and Mins only
elseif valEvent == 2 && valMax == 1 && valMin == 1
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','minMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Only Events    
elseif valEvent == 1 && valMax == 2 && valMin == 2
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Only Maxes
elseif valEvent == 2 && valMax == 1 && valMin == 2
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);

% Only mins
elseif valEvent == 2 && valMax == 2 && valMin == 1
    prompt = {'Events file name?',...
        'Maxes file name?'...
        'Mins file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','NA (leave alone)','minMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
    eventOutFile = outputParams(1);
    maxOutFile = outputParams(2);
    minOutFile = outputParams(3);


% Nothing selected
elseif valEvent == 2 && valMax == 2 && valMin == 2
    error('You dont want any data?');

end   
%% Event Abbreviations?

if valEvent == 1

    if numEvents == 3
        prompt = {'Event 1 abbreviation'...
            'Event 2 abbreviation'...
            'Event 3 abbreviation'};
        dlgtitle = 'Event Abbreviations';
        dims = [1 45];
        defEvents = {'e1','e2','e3'};
        eventParams = inputdlg(prompt,dlgtitle,dims,defEvents);
        
    elseif numEvents == 4
        prompt = {'Event 1 abbreviation'...
            'Event 2 abbreviation'...
            'Event 3 abbreviation'...
            'Event 4 abbreviation'};
        dlgtitle = 'Event Abbreviations';
        dims = [1 45];
        defEvents = {'e1','e2','e3','e4'};
        eventParams = inputdlg(prompt,dlgtitle,dims,defEvents);
        
    elseif numEvents == 5
        prompt = {'Event 1 abbreviation'...
            'Event 2 abbreviation'...
            'Event 3 abbreviation'...
            'Event 4 abbreviation'...
            'Event 5 abbreviation'};
        dlgtitle = 'Event Abbreviations';
        dims = [1 45];
        defEvents = {'e1','e2','e3','e4','e5'};
        eventParams = inputdlg(prompt,dlgtitle,dims,defEvents);
    end
        % Script stopper if cancel button is selected when defining 
        % event abbreviations
        eventStop = size(eventParams);
        
        if eventStop == 0
            error('event abbreviations were not defined');
        end
end

%% Run eventFinder if selected

if valEvent == 1
    % download maxFinder from online repository
     websave('extractData.m',...
             'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');
     websave('eventFinder.m',...
         'https://raw.githubusercontent.com/kww-22/matlab/master/eventFinder.m');

% Create character strings for saved output files
    avefile = ['ave' outputParams(1)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(1)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');

% Builds selection box. Sort by event: eventSort = 1; Sort by participant: eventSort = 2
eventSort = listdlg('PromptString', 'Average event master grouped by...?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Event', 'Participant'}, ...
    'Name', 'VAverage Event Sorting', ...
    'ListSize', [225 100]);

    % Run maxFinder.m;
    [eventMaster,aveEventMaster] = eventFinder(fileNames,numFiles,...
        numTrials,numEvents,varRow,eventParams,eventSort,saveFile,saveAveFile);

    % Remove downloaded files from selected directory
    delete eventFinder.m
    delete extractData.m
end

%% Run maxFinder if selected

if valMax == 1
    % download extractData and maxFinder from online repository
    websave('extractData.m',...
            'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');
    websave('maxFinder.m',...
        'https://raw.githubusercontent.com/kww-22/matlab/master/maxFinder.m');
    
    % Create character strings for saved output files
    avefile = ['ave' outputParams(2)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(2)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');
    
    % Run maxFinder.m;
    [maxMaster, aveMaxMaster] = maxFinder(fileNames,numFiles,numTrials,...
        numEvents,varRow,saveFile,saveAveFile);
    
    % Remove downloaded files from selected directory
    delete maxFinder.m
    delete extractData.m
end

%% Run minFinder if selected

if valMin == 1
    % download extractData and minFinder from online repository
    websave('extractData.m',...
            'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');
    websave('minFinder.m',...
        'https://raw.githubusercontent.com/kww-22/matlab/master/minFinder.m');
    
    % Create character strings for saved output files
    avefile = ['ave' outputParams(3)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(3)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');
    
    % Run mixFinder.m;
    [minMaster, aveMinMaster] = minFinder(fileNames,numFiles,numTrials,...
        numEvents,varRow,saveFile,saveAveFile);
    
    % Remove downloaded files from selected directory
    delete minFinder.m
    delete extractData.m
end

%% Display when complete   
disp('the script finished')