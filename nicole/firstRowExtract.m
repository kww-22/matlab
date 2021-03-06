%% firstRowExtract.m
% *************************************************************************
% Interactive script to grab the first row of data from MotionMonitor
% Classic .exp exports and compile them into a master data table
%
% Promts user for desired output file names(s) and saves indidual trial 
% and participant average master tables to specified directory
% 
% Designed for Bordelon et al manuscript examining mechanical energy
% flow through the kinetic chain during softball hitting
%
% Functions called: 
%   1) extractData.m
%   2) grabFirstRow.m
%
% Downloads required functions from GitHub repository
%
% Author: Kyle Wasserberger
% Sports Medicine and Movement Lab
% School of Kinesiology; Auburn University
% Auburn, AL, USA
% Last Updated: 2020-06-14
% *************************************************************************
%% Gather directory information

% Set directory
% Choose folder path containing exported and necessary matlab files
path = uigetdir();
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
          'Which row contains the variable names?',...
          'Output file name? (Include extension; .csv or .xls recommended)'};
dlgtitle = 'Define trial parameters';
dims = [1 45];
defParams = {'3','4','9','fileName.csv'};
trialParams = inputdlg(prompt,dlgtitle,dims,defParams);

numTrials = str2double(trialParams{1});
numEvents = str2double(trialParams{2});
VarRow = str2double(trialParams{3});
file = trialParams(4);

%% Error Trapping

% Numer of total files divided by number of trials per participant
% should equal zero
if rem(numFiles,numTrials) ~= 0
    
    disp('ERROR: The total number of files indicates unequal trials per participant')
    
else

    %% Cancel script instead of entering file name
    paramsResponse = size(trialParams);
    
    if paramsResponse(1) == 0
            clear
            clc
            disp('ERROR: You hit cancel instead of entering a file name');
    else      
    %% Create character strings for saved output files
        avefile = ['ave' file];
        avefile = join(avefile,'');
        
        strMaster = [path file];
        strAveMaster = [path avefile];
        
        saveFile = join(strMaster,'/');
        
        saveAveFile = join(strAveMaster,'/');
        
     %% Load required functions from online respository
        websave('extractData.m',...
            'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');
        websave('grabFirstRow.m',...
            'https://raw.githubusercontent.com/kww-22/matlab/master/nicole/grabFirstRow.m');
        
     %% Compile the first row of data from individual exports
      % Arguments: fileNames, numFiles, numTrials, saveFile, saveAveFile 
        [firstRowMaster, aveFirstRowMaster] = grabFirstRow(fileNames,numFiles,numTrials,VarRow,saveFile,saveAveFile);
        
    end
delete extractData.m
delete grabFirstRow.m
end