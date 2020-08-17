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
%% 

thor_rot_mstr_y = [];
thor_flex_mstr_y = [];
thor_lflex_mstr_y = [];

pelv_rot_mstr_y = [];
pelv_antpst_mstr_y = [];
pelv_ltilt_mstr_y = [];

clock3_time_y = [];
clock12_time_y = [];
fc_time_y = [];
br_time_y = [];

for i = 1:numFiles
    
    %% extract sampling frequency (fs)
    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});

    %% read in and clean data
    data = extractData(fileNames.fileNames{i},'text',9);
    events = find(data.VEM_0 == 1);
    
    % round to .1 seconds before first and after fourth event
    data = data(events(1)-round(.1*fs):events(4)+round(.1*fs),:);
    % redefine event indicies
    events = find(data.VEM_0 == 1);
    % get timing of events as % of trimmed data
    clock3_time_y = [clock3_time_y; events(1)/height(data)];
    clock12_time_y = [clock12_time_y; events(2)/height(data)];
    fc_time_y = [fc_time_y; events(3)/height(data)];
    br_time_y = [br_time_y; events(4)/height(data)];
    %% time normalizing
    % create original time axis
    time_old = (0:size(data)-1)'/fs; 
    % normalize original time axis
    time_old = time_old/time_old(end);
    % create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);
    %% lefty if loop
    if max(data.LwristAV) > max(data.RwristAV)
        data.thor_rot = data.thor_rot*-1;
        data.thor_latflex = data.thor_latflex*-1;
        data.pelv_rot = data.pelv_rot*-1;
        data.pelv_lattilt = data.pelv_lattilt*-1;
    end
        thor_rot_mstr_y = [thor_rot_mstr_y; interp1(time_old,data.thor_rot,time,'spline')];
        thor_flex_mstr_y = [thor_flex_mstr_y; interp1(time_old,data.thor_flex,time,'spline')];
        thor_lflex_mstr_y = [thor_lflex_mstr_y; interp1(time_old,data.thor_latflex,time,'spline')];
        
        pelv_rot_mstr_y = [pelv_rot_mstr_y; interp1(time_old,data.pelv_rot,time,'spline')];
        pelv_antpst_mstr_y = [pelv_antpst_mstr_y; interp1(time_old,data.pelv_antpost,time,'spline')];
        pelv_ltilt_mstr_y = [pelv_ltilt_mstr_y; interp1(time_old,data.pelv_lattilt,time,'spline')];
end
fclose('all');
%%
%SB_C_KO_02-21-19 has funky pelvis
% pelv_rot_mstr_y(100:102,:) = "Na";
% pelv_antpst_mstr_y(100:102,:) = "Na";
% pelv_ltilt_mstr_y(100:102,:) = "Na";
% 
% thor_rot_mstr_y(100:102,:) = "Na";
% thor_lflex_mstr_y(100:102,:) = "Na";
% thor_flex_mstr_y(100:102,:) = "Na";
%% compute average (+/-sd) timing of events for plots
clock3mean_y = mean(clock3_time_y);
clock3std_y = std(clock3_time_y);
clock3xlims_y = [clock3mean_y-clock3std_y clock3mean_y+clock3std_y clock3mean_y+clock3std_y clock3mean_y-clock3std_y];

clock12mean_y = mean(clock12_time_y);
clock12std_y = std(clock12_time_y);
clock12xlims_y = [clock12mean_y-clock12std_y clock12mean_y+clock12std_y clock12mean_y+clock12std_y clock12mean_y-clock12std_y];

fcmean_y = mean(fc_time_y);
fcstd_y = std(fc_time_y);
fcxlims_y = [fcmean_y-fcstd_y fcmean_y+fcstd_y fcmean_y+fcstd_y fcmean_y-fcstd_y];

brmean_y = mean(br_time_y);
brstd_y = std(br_time_y);
brxlims_y = [brmean_y-brstd_y brmean_y+brstd_y brmean_y+brstd_y brmean_y-brstd_y];


%% plots

plot_distribution(time,thor_flex_mstr_y)
plot_distribution(time,thor_flex_mstr_c)
ax = gca;
color1 = ax.ColorOrder(1,:);
color2 = ax.ColorOrder(2,:);
color3 = ax.ColorOrder(3,:);
ylims = [min(ylim) min(ylim) max(ylim) max(ylim)];
close all
figure('color','w')
box('on')
title(["thorax lateral flexion"; "\it\fontsize{8}{youth = blue\newlinecollege = yellow}"],'fontweight','normal')
xlabel("Time (%)")
ylabel("Joint Angle (degrees)")
xline(mean(clock3_time_y),'linestyle','--','color',color1,'label','3:00','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(clock12_time_y),'linestyle','--','color',color1,'label','12:00','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(fc_time_y),'linestyle','--','color',color1,'label','FC','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(br_time_y),'linestyle','--','color',color1,'label','BR','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(clock3_time_c),'linestyle','--','color',color3,'label','3:00','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(clock12_time_c),'linestyle','--','color',color3,'label','12:00','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(fc_time_c),'linestyle','--','color',color3,'label','FC','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
xline(mean(br_time_c),'linestyle','--','color',color3,'label','BR','labelhorizontalalignment','center',...
'labelorientation','horizontal','fontname','times new roman','fontsize',8);
yline(0);
% patch(clock3xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(clock12xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(fcxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(brxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
plot_distribution(time,thor_flex_mstr_y);
plot_distribution(time,thor_flex_mstr_c);
text((clock3mean_c+clock3mean_y)/2,25,'3:00','fontname','times new roman','fontsize',8)
set(gca,'fontname','times new roman');
xticks([0 .25 .50 .75 1]);
xticklabels(["0" "25" "50" "75" "100"]);

%%
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
%% 

thor_rot_mstr_y = [];
thor_flex_mstr_y = [];
thor_lflex_mstr_y = [];

pelv_rot_mstr_y = [];
pelv_antpst_mstr_y = [];
pelv_ltilt_mstr_y = [];

clock3_time_y = [];
clock12_time_y = [];
fc_time_y = [];
br_time_y = [];

for i = 1:numFiles
    
    %% extract sampling frequency (fs)
    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});

    %% read in and clean data
    data = extractData(fileNames.fileNames{i},'text',9);
    events = find(data.VEM_0 == 1);
    
    % round to .1 seconds before first and after fourth event
    data = data(events(1)-round(.1*fs):events(4)+round(.1*fs),:);
    % redefine event indicies
    events = find(data.VEM_0 == 1);
    % get timing of events as % of trimmed data
    clock3_time_y = [clock3_time_y; events(1)/height(data)];
    clock12_time_y = [clock12_time_y; events(2)/height(data)];
    fc_time_y = [fc_time_y; events(3)/height(data)];
    br_time_y = [br_time_y; events(4)/height(data)];
    %% time normalizing
    % create original time axis
    time_old = (0:size(data)-1)'/fs; 
    % normalize original time axis
    time_old = time_old/time_old(end);
    % create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);
    %% lefty if loop
    if max(data.LwristAV) > max(data.RwristAV)
        data.thor_rot = data.thor_rot*-1;
        data.thor_latflex = data.thor_latflex*-1;
        data.pelv_rot = data.pelv_rot*-1;
        data.pelv_lattilt = data.pelv_lattilt*-1;
    end
        thor_rot_mstr_y = [thor_rot_mstr_y; interp1(time_old,data.thor_rot,time,'spline')];
        thor_flex_mstr_y = [thor_flex_mstr_y; interp1(time_old,data.thor_flex,time,'spline')];
        thor_lflex_mstr_y = [thor_lflex_mstr_y; interp1(time_old,data.thor_latflex,time,'spline')];
        
        pelv_rot_mstr_y = [pelv_rot_mstr_y; interp1(time_old,data.pelv_rot,time,'spline')];
        pelv_antpst_mstr_y = [pelv_antpst_mstr_y; interp1(time_old,data.pelv_antpost,time,'spline')];
        pelv_ltilt_mstr_y = [pelv_ltilt_mstr_y; interp1(time_old,data.pelv_lattilt,time,'spline')];
end
fclose('all');
%%
%SB_C_KO_02-21-19 has funky pelvis
% pelv_rot_mstr_y(100:102,:) = "Na";
% pelv_antpst_mstr_y(100:102,:) = "Na";
% pelv_ltilt_mstr_y(100:102,:) = "Na";
% 
% thor_rot_mstr_y(100:102,:) = "Na";
% thor_lflex_mstr_y(100:102,:) = "Na";
% thor_flex_mstr_y(100:102,:) = "Na";
%% compute average (+/-sd) timing of events for plots
clock3mean_y = mean(clock3_time_y);
clock3std_y = std(clock3_time_y);
clock3xlims_y = [clock3mean_y-clock3std_y clock3mean_y+clock3std_y clock3mean_y+clock3std_y clock3mean_y-clock3std_y];

clock12mean_y = mean(clock12_time_y);
clock12std_y = std(clock12_time_y);
clock12xlims_y = [clock12mean_y-clock12std_y clock12mean_y+clock12std_y clock12mean_y+clock12std_y clock12mean_y-clock12std_y];

fcmean_y = mean(fc_time_y);
fcstd_y = std(fc_time_y);
fcxlims_y = [fcmean_y-fcstd_y fcmean_y+fcstd_y fcmean_y+fcstd_y fcmean_y-fcstd_y];

brmean_y = mean(br_time_y);
brstd_y = std(br_time_y);
brxlims_y = [brmean_y-brstd_y brmean_y+brstd_y brmean_y+brstd_y brmean_y-brstd_y];


%% plots

plot_distribution(time,thor_lflex_mstr_y)
plot_distribution(time,thor_lflex_mstr_c)
ax = gca;
color1 = ax.ColorOrder(1,:);
color2 = ax.ColorOrder(2,:);
color3 = ax.ColorOrder(3,:);
ylims = [min(ylim) min(ylim) max(ylim) max(ylim)];
close all
figure('color','w')
box('on')
title(["thorax lateral flexion"; "\it\fontsize{8}{youth = blue\newlinecollege = yellow\newline.}"],'fontweight','normal')
xlabel("Time (%)")
ylabel(['Joint Angle (' char(176) ')'])
xline(mean(clock3_time_y),'linestyle','--','color',color1)
xline(mean(clock12_time_y),'linestyle','--','color',color1)
xline(mean(fc_time_y),'linestyle','--','color',color1)
xline(mean(br_time_y),'linestyle','--','color',color1)
xline(mean(clock3_time_c),'linestyle','--','color',color3)
xline(mean(clock12_time_c),'linestyle','--','color',color3)
xline(mean(fc_time_c),'linestyle','--','color',color3)
xline(mean(br_time_c),'linestyle','--','color',color3)
yline(0);
% patch(clock3xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(clock12xlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(fcxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
% patch(brxlims_y,ylims,'k','facealpha',.1,'edgecolor','none')
plot_distribution(time,thor_lflex_mstr_y);
plot_distribution(time,thor_lflex_mstr_c);
text((clock3mean_c+clock3mean_y)/2,max(ylim)+1,'3:00','fontname','times new roman',...
    'fontsize',8,'horizontalalignment','center')
text((clock12mean_c+clock12mean_y)/2,max(ylim)+1,'12:00','fontname','times new roman',...
    'fontsize',8,'horizontalalignment','center')
text((fcmean_c+fcmean_y)/2,max(ylim)+1,'FC','fontname','times new roman',...
    'fontsize',8,'horizontalalignment','center')
text((brmean_c+brmean_y)/2,max(ylim)+1,'BR','fontname','times new roman',...
    'fontsize',8,'horizontalalignment','center')
set(gca,'fontname','times new roman');
xticks([0 .25 .50 .75 1]);
xticklabels(["0" "25" "50" "75" "100"]);