%% Energy generation, absorption, & transfer @ the shoulder & elbow
% *************************************************************************
% Script to parse energy flow between the thorax, humerus, and forearm into
% energy transfer, generation, and absorption as outlined in Robinson &
% Winter (1980): Mechanical Energy Generation, Absorption, and Transfer
% Amongst Segments Druing Walking
%
% Writes results to .csv table in selected directory
% 
% Contributors:
%   Anne de Swart
%       Human Movement Sciences
%       Vrije Universiteit 
%       Amsterdam, Netherlands
%   Kyle Wasserberger
%       Sports Medicine and Movement Lab
%       School of Kinesiology
%       Auburn University
%       Alabama, USA
% Last Updated: 2020-06-28
% *************************************************************************
%% Close all previous windows

clear
clc
close all

% *************************************************************************
% sfc = stride foot contact
% mer = maximum external rotation
% br = ball release
% mir = maximum internal rotation
% cp = cocking phase (sfc - mer)
% ap = acceleration phase (mer - br)
% dp = deceleration phase (br - mir)
% jfp = joint force power
% stp = segment torque power
% abs = energy absorption
% gen = energy generation
% trans = energy transfer
% fs = sampling frequency
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

%% Variable initialization

shldr_trans_jfp = nan(100,numFiles);
shldr_trans_stp = nan(100,numFiles);
shldr_gen_stp = nan(100,numFiles);
shldr_trans_net = nan(100,numFiles);

elb_trans_jfp = nan(100,numFiles);
elb_trans_stp = nan(100,numFiles);
elb_gen_stp = nan(100,numFiles);
elb_trans_net = nan(100,numFiles);

shldr_trans_jfp_cp = nan(numFiles,1);
shldr_trans_jfp_ap = nan(numFiles,1);
shldr_trans_jfp_dp = nan(numFiles,1);

shldr_trans_stp_cp = nan(numFiles,1);
shldr_trans_stp_ap = nan(numFiles,1);
shldr_trans_stp_dp = nan(numFiles,1);

shldr_trans_net_cp = nan(numFiles,1);
shldr_trans_net_ap = nan(numFiles,1);
shldr_trans_net_dp = nan(numFiles,1);

shldr_gen_cp = nan(numFiles,1);
shldr_gen_ap = nan(numFiles,1);
shldr_gen_dp = nan(numFiles,1);

shldr_abs_cp = nan(numFiles,1);
shldr_abs_ap = nan(numFiles,1);
shldr_abs_dp = nan(numFiles,1);

elb_trans_jfp_cp = nan(numFiles,1);
elb_trans_jfp_ap = nan(numFiles,1);
elb_trans_jfp_dp = nan(numFiles,1);

elb_trans_stp_cp = nan(numFiles,1);
elb_trans_stp_ap = nan(numFiles,1);
elb_trans_stp_dp = nan(numFiles,1);

elb_trans_net_cp = nan(numFiles,1);
elb_trans_net_ap = nan(numFiles,1);
elb_trans_net_dp = nan(numFiles,1);

elb_gen_cp = nan(numFiles,1);
elb_gen_ap = nan(numFiles,1);
elb_gen_dp = nan(numFiles,1);

elb_abs_cp = nan(numFiles,1);
elb_abs_ap = nan(numFiles,1);
elb_abs_dp = nan(numFiles,1);

%pkh_master = nan(numFiles,1);
%sfc_master = nan(numFiles,1);
mer_master = nan(numFiles,1);
br_master = nan(numFiles,1);
%mir_master = nan(numFiles,1);

%% Shoulder & elbow generation, transfer, & absorption per trial

websave('extractData.m',...
    'https://raw.githubusercontent.com/kww-22/matlab/master/extractData.m');

for i = 1:numFiles
    %% Import i-th file
    data = extractData(fileNames.fileNames{i},'text',9);
    
    %% Get sample frequency from .txt file
    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});
    clear txt

    %% Get event timings and trim trial
    
    % find event indices
    events = find(data.VEM_0 == 1);
    
    % trim data to only include rows between SFC & MIR
    data = data(events(2):events(end),:); 
    % if including PKH in analysis, change "events(2)" to "events(1)"
    % and un-comment necessary lines below
    
    % redefine and label event indices
    events = find(data.VEM_0 == 1);
   %pkh = events(1);    % Peak knee height
    sfc = events(1);    % Stride foot contact
    mer = events(2);    % Max external rotation
    br = events(3);     % Ball release
    mir = events(4);    % Max internal rotation
    
    % event timing as a % between pkh and mir
    % sfc_master(i,1) = (sfc - pkh)/(mir-pkh);
    % mer_master(i,1) = (mer - pkh)/(mir-pkh);
    % br_master(i,1) = (br - pkh)/(mir-pkh);
    
    % event timing as a % between sfc and mir
    mer_master(i,1) = (mer - sfc)/(mir - sfc);
    br_master(i,1) = (br - sfc)/(mir - sfc);   


    %% Calculate energy transfer, generation, & absorption time series at shoulder
    
    pwrtrans_stp = nan(height(data),1);
    pwrgen_stp = nan(height(data),1);
    trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.HumJfpProx';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.ThorStpHum(j).*data.HumStpProx(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.ThorStpHum(j)) abs(data.HumStpProx(j))]);
            
            if abs(data.ThorStpHum(j)) == abs(data.HumStpProx(j))
                pwrtrans_stp(j) = data.HumStpProx(j);
                pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.ThorStpHum(j))
                pwrtrans_stp(j) = -data.ThorStpHum(j);
                pwrgen_stp(j) = data.HumStpProx(j) + data.ThorStpHum(j);
            elseif minSTP == abs(data.HumStpProx(j))
                pwrtrans_stp(j) = data.HumStpProx(j);
                pwrgen_stp(j) = data.HumStpProx(j) + data.ThorStpHum(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            pwrtrans_stp(j) = 0;
            pwrgen_stp(j) = data.HumStpProx(j) + data.ThorStpHum(j);
            
        end
    trans_net(j) = pwrtrans_jfp(j) + pwrtrans_stp(j);
    end
    %% Calculate total trans/gen/abs over phases @ shoulder
        
        shldr_trans_jfp_cp(i,1) = trapz(pwrtrans_jfp(sfc:mer))/fs;
        shldr_trans_jfp_ap(i,1) = trapz(pwrtrans_jfp(mer:br))/fs;
        shldr_trans_jfp_dp(i,1) = trapz(pwrtrans_jfp(br:mir))/fs;
        
        shldr_trans_stp_cp(i,1) = trapz(pwrtrans_stp(sfc:mer))/fs;
        shldr_trans_stp_ap(i,1) = trapz(pwrtrans_stp(mer:br))/fs;
        shldr_trans_stp_dp(i,1) = trapz(pwrtrans_stp(br:mir))/fs;
        
        shldr_trans_net_cp(i,1) = trapz(trans_net(sfc:mer))/fs;
        shldr_trans_net_ap(i,1) = trapz(trans_net(mer:br))/fs;
        shldr_trans_net_dp(i,1) = trapz(trans_net(br:mir))/fs;
        
        shldr_gen_cp(i,1) = (trapz(pwrgen_stp(sfc:mer)) + trapz(abs(pwrgen_stp(sfc:mer)))) / (2*fs);
        shldr_gen_ap(i,1) = (trapz(pwrgen_stp(mer:br)) + trapz(abs(pwrgen_stp(mer:br)))) / (2*fs);
        shldr_gen_dp(i,1) = (trapz(pwrgen_stp(br:mir)) + trapz(abs(pwrgen_stp(br:mir)))) / (2*fs);
        
        shldr_abs_cp(i,1) = (trapz(pwrgen_stp(sfc:mer)) - trapz(abs(pwrgen_stp(sfc:mer)))) / (2*fs);
        shldr_abs_ap(i,1) = (trapz(pwrgen_stp(mer:br)) - trapz(abs(pwrgen_stp(mer:br)))) / (2*fs);
        shldr_abs_dp(i,1) = (trapz(pwrgen_stp(br:mir)) - trapz(abs(pwrgen_stp(br:mir)))) / (2*fs);
        
    %% Inerpolate shoulder power curves
    
    % Create original time axis
    time_old = (0:size(data)-1)'/fs;
    
    % Normalize original time axis
    time_old = time_old/time_old(end);
    
    % Create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);
    
    % Interpolation
    shldr_trans_jfp(:,i) = interp1(time_old,pwrtrans_jfp,time,'spline');
    shldr_trans_stp(:,i) = interp1(time_old,pwrtrans_stp,time,'spline');
    shldr_gen_stp(:,i) = interp1(time_old,pwrgen_stp,time,'spline');
    shldr_trans_net(:,i) = interp1(time_old,trans_net,time,'spline');
    
    clear pwrtans_jfp pwrtrans_stp pwrgen_stp trans_net;
    
    %% Calculate energy transfer, generation, & absorption at elbow
    
    pwrtrans_stp = nan(height(data),1);
    pwrgen_stp = nan(height(data),1);
    trans_net = nan(height(data),1);
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.FaJfpProx';
    
    for k = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.HumStpDist(k).*data.FaStpProx(k));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.HumStpDist(k)) abs(data.FaStpProx(k))]);
            
            if abs(data.HumStpDist(k)) == abs(data.FaStpProx(k))
                pwrtrans_stp(k) = data.FaStpProx(k);
                pwrgen_stp(k) = 0;
            elseif minSTP == abs(data.HumStpDist(k))
                pwrtrans_stp(k) = -data.HumStpDist(k);
                pwrgen_stp(k) = data.FaStpProx(k) + data.HumStpDist(k);
            elseif minSTP == abs(data.FaStpProx(k))
                pwrtrans_stp(k) = data.FaStpProx(k);
                pwrgen_stp(k) = data.FaStpProx(k) + data.HumStpDist(k);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            pwrtrans_stp(k) = 0;
            pwrgen_stp(k) = data.HumStpProx(k) + data.ThorStpHum(k);
            
        end
    trans_net(k) = pwrtrans_jfp(k) + pwrtrans_stp(k);    
    end
    %% Calculate total trans/gen/abs over phases @ elbow
        
        elb_trans_jfp_cp(i,1) = trapz(pwrtrans_jfp(sfc:mer))/fs;
        elb_trans_jfp_ap(i,1) = trapz(pwrtrans_jfp(mer:br))/fs;
        elb_trans_jfp_dp(i,1) = trapz(pwrtrans_jfp(br:mir))/fs;
        
        elb_trans_stp_cp(i,1) = trapz(pwrtrans_stp(sfc:mer))/fs;
        elb_trans_stp_ap(i,1) = trapz(pwrtrans_stp(mer:br))/fs;
        elb_trans_stp_dp(i,1) = trapz(pwrtrans_stp(br:mir))/fs;
        
        elb_trans_net_cp(i,1) = trapz(trans_net(sfc:mer))/fs;
        elb_trans_net_ap(i,1) = trapz(trans_net(mer:br))/fs;
        elb_trans_net_dp(i,1) = trapz(trans_net(br:mir))/fs;
        
        elb_gen_cp(i,1) = (trapz(pwrgen_stp(sfc:mer)) + trapz(abs(pwrgen_stp(sfc:mer)))) / (2*fs);
        elb_gen_ap(i,1) = (trapz(pwrgen_stp(mer:br)) + trapz(abs(pwrgen_stp(mer:br)))) / (2*fs);
        elb_gen_dp(i,1) = (trapz(pwrgen_stp(br:mir)) + trapz(abs(pwrgen_stp(br:mir)))) / (2*fs);
        
        elb_abs_cp(i,1) = (trapz(pwrgen_stp(sfc:mer)) - trapz(abs(pwrgen_stp(sfc:mer)))) / (2*fs);
        elb_abs_ap(i,1) = (trapz(pwrgen_stp(mer:br)) - trapz(abs(pwrgen_stp(mer:br)))) / (2*fs);
        elb_abs_dp(i,1) = (trapz(pwrgen_stp(br:mir)) - trapz(abs(pwrgen_stp(br:mir)))) / (2*fs);
    
    %% Interpolate elbow power curves
    
    elb_trans_jfp(:,i) = interp1(time_old,pwrtrans_jfp,time,'spline');
    elb_trans_stp(:,i) = interp1(time_old,pwrtrans_stp,time,'spline');
    elb_gen_stp(:,i) = interp1(time_old,pwrgen_stp,time,'spline');
    elb_trans_net(:,i) = interp1(time_old,trans_net,time,'spline');
    
    clear pwrtrans_jfp pwrtrans_stp pwrgen_stp trans_net 
end

%% Calculate maxes and mins of interest

peak_shldr_trans_jfp = max(shldr_trans_jfp)';
peak_elb_trans_jfp = max(elb_trans_jfp)';
peak_shldr_trans_stp = max(shldr_trans_stp)';
peak_elb_trans_stp = max(elb_trans_stp)';
peak_shldr_trans_net = max(shldr_trans_net)';
peak_elb_trans_net = max(elb_trans_net)';
peak_shldr_gen = max(shldr_gen_stp)';
peak_elb_gen = max(elb_gen_stp)';
peak_shldr_abs = min(shldr_gen_stp)';
peak_elb_abs = min(elb_gen_stp)';

%% Plot individual trial curves for visual inspection

figure('color','w');
subplot(2,3,1)
plot(shldr_trans_jfp);
title('Shoulder Transfer JFP');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,3,2)
plot(shldr_trans_stp);
title('Shoulder Transfer STP');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,3,3)
plot(shldr_gen_stp);
title('Shoulder Gen/Abs');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,3,4)
plot(elb_trans_jfp);
title('Elbow Transfer JFP');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,3,5)
plot(elb_trans_stp);
title('Elbow Transfer STP');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,3,6)
plot(elb_gen_stp)
title('Elbow Gen/Abs');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);

%% Plot net transfer @ shoulder & elbow

figure('color','w');
subplot(2,1,1)
plot(shldr_trans_net);
title('Net Shoulder Transfer');
%xlabel('Time (%SFC-MIR)');
ylabel('Joint Power (W)');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);
subplot(2,1,2)
plot(elb_trans_net);
title('Net Elbow Transfer');
xlabel('Time (%SFC-MIR)');
ylabel('Joint Power (W)');
xline(mean(mer_master)*100);
xline(mean(br_master)*100);
yline(0);

%% Write results to table
participants = fileNames.fileNames;
genabstrans_master = table();

genabstrans_master = ...
    addvars(genabstrans_master,participants,...
    peak_shldr_abs,...
    peak_shldr_gen,...
    peak_shldr_trans_jfp,...
    peak_shldr_trans_stp,...
    peak_shldr_trans_net,...
    shldr_abs_ap,...
    shldr_abs_cp,...
    shldr_abs_dp,...
    shldr_gen_ap,...
    shldr_gen_cp,...
    shldr_gen_dp,...
    shldr_trans_jfp_ap,...
    shldr_trans_jfp_cp,...
    shldr_trans_jfp_dp,...
    shldr_trans_stp_ap,...
    shldr_trans_stp_cp,...
    shldr_trans_stp_dp,...
    shldr_trans_net_ap,...
    shldr_trans_net_cp,...
    shldr_trans_net_dp,...
    peak_elb_abs,...
    peak_elb_gen,...
    peak_elb_trans_jfp,...
    peak_elb_trans_stp,...
    peak_elb_trans_net,...
    elb_abs_ap,...
    elb_abs_cp,...
    elb_abs_dp,...
    elb_gen_ap,...
    elb_gen_cp,...
    elb_gen_dp,...
    elb_trans_jfp_ap,...
    elb_trans_jfp_cp,...
    elb_trans_jfp_dp,...
    elb_trans_stp_ap,...
    elb_trans_stp_cp,...
    elb_trans_stp_dp,...
    elb_trans_net_ap,...
    elb_trans_net_cp,...
    elb_trans_net_dp);

writetable(genabstrans_master,[path '/genabstrans_master.csv'])
%% Clean up after youself like mama taught you
delete extractData.m;
clear ans br data events fs i j k mer minSTP mir sfc Sign time time_old
clc;
disp('The script has finished');