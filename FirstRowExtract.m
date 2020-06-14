%% FirstRowExtract.m
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
% Required in directory: 
%   1) exported files (.exp or .txt) 
%   2) extractData.m
%   3) grabFirstRow.m
%
% FirstRowExtract does NOT have to be in directory
% Author: Kyle Wasserberger
% Last Updated: 2020-06-12
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


%% Imput trial parameters

prompt = {'How many trials per participant?',...
          'Which row contains the variable names?',...
          'What should your output file be called? (Include .csv extension)'};
dlgtitle = 'Define trial parameters';
dims = [1 35];
defParams = {'3','9','fileName.ext'};
trialParams = inputdlg(prompt,dlgtitle,dims,defParams);

numTrials = str2double(trialParams{1});
VarRow = str2double(trialParams{2});
file = trialParams(3);

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
        avefile = ['ave' file];
        avefile = join(avefile,'');
        
        strMaster = [path file];
        strAveMaster = [path avefile];
        
        saveFile = join(strMaster,'/');
        
        saveAveFile = join(strAveMaster,'/');
        
        %% Compile the first row of data from individual exports
        % Arguments: fileNames, numFiles, numTrials, saveFile, saveAveFile
        
        grabFirstRow(fileNames,numFiles,numTrials,VarRow,saveFile,saveAveFile);
        

    end
end