%% Gather directory information

% clear
% clc

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
list = string(varNames);

websave('extractData.m',...
    'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');
%% chose variables to interpolate

% [indx,tf] = listdlg('ListString',list,...
%     'PromptString','Would you like any time-series plots?',...
%     'CancelString','No thanks',...
%     'ListSize',[200 300]);

%%
knee_flex_mstr_c = [];
elb_flex_mstr_c = [];

thor_rot_mstr_c = [];
thor_flex_mstr_c = [];
thor_lflex_mstr_c = [];
thor_rot_velo_mstr_c = [];

pelv_rot_mstr_c = [];
pelv_antpst_mstr_c = [];
pelv_ltilt_mstr_c = [];
pelv_rot_velo_mstr_c = [];

com_velo_mstr_c = [];
prctRL_mstr_c = [];

clock3_time_c = [];
clock12_time_c = [];
fc_time_c = [];
br_time_c = [];

for i = 1:numFiles
    
    %% extract sampling frequency (fs)
    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});

    %% read in and clean data
    data = extractData(fileNames.fileNames{i},'text',9);
    events = find(data.VEM_0 == 1);
    
    % round to .1 seconds before first and after fourth event
    data = data(events(1):events(4)+round(fs*.1),:);
    % redefine event indicies
    events = find(data.VEM_0 == 1);
    
    %% get timing of events as % of trimmed data
    clock3_time_c = [clock3_time_c; events(1)/height(data)];
    clock12_time_c = [clock12_time_c; events(2)/height(data)];
    fc_time_c = [fc_time_c; events(3)/height(data)];
    br_time_c = [br_time_c; events(4)/height(data)];
    
    %% time normalizing
    % create original time axis
    time_old = (0:size(data)-1)'/fs; 
    % normalize original time axis
    time_old = time_old/time_old(end);
    % create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);
    
    %% lefty if statments
    if max(data.LwristAV) > max(data.RwristAV)
        data.thor_rot = data.thor_rot*-1;
        data.thor_latflex = data.thor_latflex*-1;
        data.pelv_rot = data.pelv_rot*-1;
        data.pelv_lattilt = data.pelv_lattilt*-1;
        data.LkneeFlex = data.RkneeFlex;
        data.prctRL = data.prctLR;
        data.RelbFlex = data.LelbFlex;
    end
    %% interpolate variables    
        knee_flex_mstr_c = [knee_flex_mstr_c; interp1(time_old,data.LkneeFlex,time,'spline')];
        elb_flex_mstr_c = [elb_flex_mstr_c; interp1(time_old,data.RelbFlex,time,'spline')];
        
        thor_rot_mstr_c = [thor_rot_mstr_c; interp1(time_old,data.thor_rot,time,'spline')];
        thor_flex_mstr_c = [thor_flex_mstr_c; interp1(time_old,data.thor_flex,time,'spline')];
        thor_lflex_mstr_c = [thor_lflex_mstr_c; interp1(time_old,data.thor_latflex,time,'spline')];
        thor_rot_velo_mstr_c = [thor_rot_velo_mstr_c; interp1(time_old,data.thor_rot_velo,time,'spline')];
        
        pelv_rot_mstr_c = [pelv_rot_mstr_c; interp1(time_old,data.pelv_rot,time,'spline')];
        pelv_antpst_mstr_c = [pelv_antpst_mstr_c; interp1(time_old,data.pelv_antpost,time,'spline')];
        pelv_ltilt_mstr_c = [pelv_ltilt_mstr_c; interp1(time_old,data.pelv_lattilt,time,'spline')];
        pelv_rot_velo_mstr_c = [pelv_rot_velo_mstr_c; interp1(time_old,data.pelv_rot_velo,time,'spline')];
        
        com_velo_mstr_c = [com_velo_mstr_c; interp1(time_old,data.COMvelo,time,'spline')];
        prctRL_mstr_c = [prctRL_mstr_c; interp1(time_old,data.prctRL,time,'spline')];
end
fclose('all');
%% compute average (+/-sd) timing of events for plots
clock3mean_c = mean(clock3_time_c);
clock3std_c = std(clock3_time_c);
% patch x-limits for highlighted plot areas
clock3xlims_c = [clock3mean_c-clock3std_c clock3mean_c+clock3std_c clock3mean_c+clock3std_c clock3mean_c-clock3std_c];

clock12mean_c = mean(clock12_time_c);
clock12std_c = std(clock12_time_c);
clock12xlims_c = [clock12mean_c-clock12std_c clock12mean_c+clock12std_c clock12mean_c+clock12std_c clock12mean_c-clock12std_c];

fcmean_c = mean(fc_time_c);
fcstd_c = std(fc_time_c);
fcxlims_c = [fcmean_c-fcstd_c fcmean_c+fcstd_c fcmean_c+fcstd_c fcmean_c-fcstd_c];

brmean_c = mean(br_time_c);
brstd_c = std(br_time_c);
brxlims_c = [brmean_c-brstd_c brmean_c+brstd_c brmean_c+brstd_c brmean_c-brstd_c];

%% average arrays

pStartRow = 1:3:length(thor_lflex_mstr_c)+1;

for i = 1:(length(pStartRow)-1)
    aveknee_flex_mstr_c(i,:) = mean(knee_flex_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    aveelb_flex_mstr_c(i,:) = mean(elb_flex_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avecom_velo_mstr_c(i,:) = mean(com_velo_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    aveprct_mstr_c(i,:) = mean(prctRL_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avethor_rot_mstr_c(i,:) = mean(thor_rot_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avethor_rot_velo_mstr_c(i,:) = mean(thor_rot_velo_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avethor_flex_mstr_c(i,:) = mean(thor_flex_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avethor_lflex_mstr_c(i,:) = mean(thor_lflex_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avepelv_rot_mstr_c(i,:) = mean(pelv_rot_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avepelv_rot_velo_mstr_c(i,:) = mean(pelv_rot_velo_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avepelv_antpst_mstr_c(i,:) = mean(pelv_antpst_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));
    avepelv_ltilt_mstr_c(i,:) = mean(pelv_ltilt_mstr_c(pStartRow(i):pStartRow(i+1)-1,:));  
end