% smml_gui: interactive script to handle basic data cleaning/organizing
% from MotionMonitor
% *************************************************************************
% Gives user a series of interactive prompts to select which data
% cleaning/organizing processes they would like. Options include:
%       1. Max of variables btwn first & last event
%       2. Min of variables btwn first & last event
%       3. Values of variables @ events
%       4. Average value of variables between events (over phases)
%       5. Plot time series for variables of interest (in progress; works
%       but is suuuuuuper slow; needs work!)
%       6. Basic summary statistics (mean; standard deviation; 25th, 50th,
%       and 75th percentiles)
%
% User may also provide custom phase labels, event labels and custom output
% file names
% 
% Inputs:
%   NA
%
% Outputs:
%   depending on which processes user selects, saves various .csv master and 
%   average master files to selected directory
%
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-07-21
% *************************************************************************
%% Gather directory information

clear
clc

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

%% Phases?

% Builds yes/no selection box. Yes: valPhase = 1; No: valPhase = 2
valPhase = listdlg('PromptString', 'Average of variables across phases?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Yes', 'No'}, ...
    'Name', 'Values at Events', ...
    'ListSize', [225 100]);

%% Summary statistics?

% Builds yes/no selection box. Yes: sum_stats = 1; No: sum_stats = 2
sum_stats = listdlg('PromptString', 'Summary statistics?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Yes', 'No'}, ...
    'Name', 'Values at Events', ...
    'ListSize', [225 100]);
%% Provide output file names

% 1. Events, maxes, mins, & phases
if valEvent == 1 && valMax == 1 && valMin == 1 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','minMaster.csv','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
 % 2. Events, maxes, and mins
elseif valEvent == 1 && valMax == 1 && valMin == 1 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','minMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 3. Events, maxes, and phases
elseif valEvent == 1 && valMax == 1 && valMin == 2 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','NA (leave alone)','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 4. Events, Mins, and Phases
elseif valEvent == 1 && valMax == 2 && valMin == 1 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','minMaster.csv','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
 % 5. Maxes, mins, and phases
    elseif valEvent == 2 && valMax == 1 && valMin == 1 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','minMaster.csv','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 6. Events and Maxes
elseif valEvent == 1 && valMax == 1 && valMin == 2 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','maxMaster.csv','NA (leave alone)','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 7. Events and Mins    
elseif valEvent == 1 && valMax == 2 && valMin == 1 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','minMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 8. Events and Phases
elseif valEvent == 1 && valMax == 2 && valMin == 2 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','NA (leave alone)','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 9. Maxes and Mins
elseif valEvent == 2 && valMax == 1 && valMin == 1 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phase file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','minMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 10. Maxes and Phases
elseif valEvent == 2 && valMax == 1 && valMin == 2 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','NA (leave alone)','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 11. Mins and Phases
elseif valEvent == 2 && valMax == 2 && valMin == 1 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','NA (leave alone)','minMaster.csv','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 12. Only Events    
elseif valEvent == 1 && valMax == 2 && valMin == 2 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'eventMaster.csv','NA (leave alone)','NA (leave alone)','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 13. Only Maxes
elseif valEvent == 2 && valMax == 1 && valMin == 2 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','maxMaster.csv','NA (leave alone)','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);

% 14. Only mins
elseif valEvent == 2 && valMax == 2 && valMin == 1 && valPhase == 2
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','NA (leave alone)','minMaster.csv','NA (leave alone)'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 15. Only Phases
elseif valEvent == 2 && valMax == 2 && valMin == 2 && valPhase == 1
    prompt = {'Events file name?',...
        'Maxes file name?',...
        'Mins file name?',...
        'Phases file name?'};
    dlgtitle = 'Provide output file names (Inlcude extension; .csv or .xls recommended)';
    dims = [1 45];
    defParams = {'NA (leave alone)','NA (leave alone)','NA (leave alone)','phaseMaster.csv'};
    outputParams = inputdlg(prompt,dlgtitle,dims,defParams);
    
% 16 . Nothing selected
elseif valEvent == 2 && valMax == 2 && valMin == 2 && valPhase == 2
    error('You dont want any data?');
end  
%% Event Abbreviations?

if valEvent == 1
    if numEvents == 2
        prompt = {'Event 1 abbreviation'...
            'Event 2 abbreviation'};
        dlgtitle = 'Event Abbreviations';
        dims = [1 45];
        defEvents = {'e1','e2'};
        eventParams = inputdlg(prompt,dlgtitle,dims,defEvents);
        
    elseif numEvents == 3
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

%% Phase abbreviations

if valPhase == 1
    if numEvents == 2
    prompt = {'Phase 1 abbreviation'};
        dlgtitle = 'Phase Abbreviations';
        dims = [1 45];
        defEvents = {'p1'};
        phaseParams = inputdlg(prompt,dlgtitle,dims,defEvents);
    
    elseif numEvents == 3
    prompt = {'Phase 1 abbreviation',...
            'Phase 2 abbreviation'};
        dlgtitle = 'Phase Abbreviations';
        dims = [1 45];
        defEvents = {'p1','p2'};
        phaseParams = inputdlg(prompt,dlgtitle,dims,defEvents);
    
    elseif numEvents == 4
        prompt = {'Phase 1 abbreviation',...
            'Phase 2 abbreviation',...
            'Phase 3 abbreviation'};
        dlgtitle = 'Phase Abbreviations';
        dims = [1 45];
        defEvents = {'p1','p2','p3'};
        phaseParams = inputdlg(prompt,dlgtitle,dims,defEvents);
    
    elseif numEvents == 5
            prompt = {'Phase 1 abbreviation'...
            'Phase 2 abbreviation'...
            'Phase 3 abbreviation',...
            'Phase 4 abbreviation'};
        dlgtitle = 'Phase Abbreviations';
        dims = [1 45];
        defEvents = {'p1','p2','p3','p4'};
        phaseParams = inputdlg(prompt,dlgtitle,dims,defEvents);
    end
        % Script stopper if cancel button is selected when defining 
        % event abbreviations
        phaseStop = size(phaseParams);
        
        if phaseStop == 0
            error('phase abbreviations were not defined');
        end
end

%% Plots?

websave('extractData.m',...
    'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');

% Extract first file in directory and it's import options
data = extractData(fileNames.fileNames{1},'text',varRow);
opts = detectImportOptions(fileNames.fileNames{1},'FileType','text');

% Get variable names for selection dialog
varNames = opts.VariableNames;
clear opts

% create character list for variable selection dialog box
list = string(varNames);

% Dialog box to select which variables user would like plotted
[indx,tf] = listdlg('ListString',list,...
    'PromptString','Would you like any time-series plots?',...
    'CancelString','No thanks',...
    'ListSize',[200 300]);

if valEvent == 2 && tf == 1
    
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
    % download eventFinder from online repository
     websave('eventFinder.m',...
         'https://raw.githubusercontent.com/kww-22/matlab/master/smml_gui/eventFinder.m');

    % Create character strings for saved output files
    avefile = ['ave' outputParams(1)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(1)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');

    % Builds selection box. Sort by event: eventSort = 1; Sort by participant: eventSort = 2
    eventSort = listdlg('PromptString', 'eventMaster grouped by...?', ...
        'SelectionMode', 'single', ...
        'ListString',{'Event', 'Participant'}, ...
        'Name', 'VAverage Event Sorting', ...
        'ListSize', [225 100]);

    % Error trapper 
    if size(eventSort) == 0
    error('events have to be sorted by some way if we are to avoid anarchy')
    end

    % Run eventFinder.m;
    [eventMaster,aveEventMaster] = eventFinder(fileNames,numFiles,...
        numTrials,numEvents,varRow,eventParams,eventSort,saveFile,...
        saveAveFile,sum_stats);

    % Remove downloaded files from selected directory
    delete eventFinder.m
end

%% Run maxFinder if selected

if valMax == 1
    % download maxFinder from online repository
    websave('maxFinder.m',...
        'https://raw.githubusercontent.com/kww-22/matlab/master/smml_gui/maxFinder.m');
    
    % Create character strings for saved output files
    avefile = ['ave' outputParams(2)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(2)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');
    
    % Run maxFinder.m;
    [maxMaster, aveMaxMaster] = maxFinder(fileNames,numFiles,numTrials,...
        numEvents,varRow,saveFile,saveAveFile,sum_stats);
    
    % Remove downloaded files from selected directory
    delete maxFinder.m
end

%% Run minFinder if selected

if valMin == 1
    % download minFinder from online repository
    websave('minFinder.m',...
        'https://raw.githubusercontent.com/kww-22/matlab/master/smml_gui/minFinder.m');
    
    % Create character strings for saved output files
    avefile = ['ave' outputParams(3)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(3)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');
    
    % Run mixFinder.m;
    [minMaster, aveMinMaster] = minFinder(fileNames,numFiles,numTrials,...
        numEvents,varRow,saveFile,saveAveFile,sum_stats);
    
    % Remove downloaded files from selected directory
    delete minFinder.m
end

%% Run if phaseFinder is selected

if valPhase == 1
    % download minFinder from online repository
    websave('phaseFinder.m',...
        'https://raw.githubusercontent.com/kww-22/matlab/master/smml_gui/phaseFinder.m');
    
    % Create character strings for saved output files
    avefile = ['ave' outputParams(4)];
    avefile = join(avefile,'');

    strMaster = [path outputParams(4)];
    strAveMaster = [path avefile];

    saveFile = join(strMaster,'/');

    saveAveFile = join(strAveMaster,'/');
    
    % Builds selection box. Sort by phase: phaseSort = 1; Sort by participant: phaseSort = 2
    phaseSort = listdlg('PromptString', 'phaseMaster grouped by...?', ...
    'SelectionMode', 'single', ...
    'ListString',{'Phase', 'Participant'}, ...
    'Name', 'VAverage Event Sorting', ...
    'ListSize', [225 100]);

    % Error trapper 
    if size(phaseSort) == 0
    error('phases have to be sorted by some way if we are to avoid anarchy')
    end

    % Run phaseFinder.m
    [phaseMaster, avePhaseMaster] = phaseFinder(fileNames,numFiles,...
    numTrials,numEvents,varRow,phaseParams,phaseSort,saveFile,saveAveFile,...
    sum_stats);
    
    % Remove downloaded files from selected directory
    delete phaseFinder.m
end

%% Build plots (if selected) (in progress)

if tf == 1 % tf == 1 when user selected at least one variable to plot
    
    list = string(eventParams);

    % Dialog box to select starting event for time-normalized plots
    e_start = listdlg('ListString',list,...
    'PromptString','Plot start event?',...
    'CancelString','Nevermind',...
    'ListSize',[200 300],...
    'SelectionMode','single');

    % Dialog box to select stopping event for time-normalized plots
    e_stop = listdlg('ListString',list,...
    'PromptString','Plot start event?',...
    'CancelString','Nevermind',...
    'ListSize',[200 300],...
    'SelectionMode','single');
    
    % get number of requested plots from dialog box
    numPlots = length(indx);
    
    % initialize figure window
    figure('color','w');
    
    for i = 1:numFiles
        data = extractData(fileNames.fileNames{i},'text',9);
        
        % get sample frequency from .txt file
        txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
        txt = txt{1};
        fs  = str2double(txt{10});
        clear txt
        
        % find event indices and trim data
        events = find(data.VEM_0 == 1);
        data = data(events(e_start):events(e_stop),:);
        
        % Create original time axis
        time_old = (0:size(data)-1)'/fs;

        % Normalize original time axis
        time_old = time_old/time_old(end);

        % Create new time axis for interpolation
        time = linspace(time_old(1),time_old(end),100);

    
        % sets up and plots selected variables in a 3 x n plotting space
        if numPlots / 3 <= 1
            for j = 1:numPlots
            hold on
            interp = interp1(time_old,data{:,indx(j)},time,'spline');
            subplot(1,numPlots,j)
            plot(time,interp);
            title(varNames{indx(j)})
            yline(0,'LineWidth',2)
            end
        else
            for j = 1:numPlots
            hold on
            interp = interp1(time_old,data{:,indx(j)},time,'spline');
            subplot(ceil(numPlots/3),3,j)
            plot(time,interp);
            title(varNames{indx(j)})
            yline(0,'LineWidth',2);
            end
        end
    end
end


%% Clean up workspace and display complete message
delete extractData.m
clear ans avefile data defEvents defParams dims dlgtitle eventOutFile...
    eventParams eventSort eventStop indx list maxOutFile minOutFile...
    outputParams path phaseOutFile phaseParams phaseSort phaseStop prompt...
    saveAveFile saveFile strAveMaster strMaster sum_stats trialParams tf...
    valEvent valMax valMin valPhase
clc
disp('the script finished; your quest is over')