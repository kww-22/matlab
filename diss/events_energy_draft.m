%% Gather directory information
fclose('all');
clear
clc

% Set directory
% Choose folder path containing exported and necessary matlab files
path = uigetdir();
cd(path);
% cd('C:\ProgramData\Innsport\TMM_xGen\MotionMonitor\User\kww\Export');


% Isolate files with .exp extension
% fileDir = dir('*.exp');
fileDir = dir('*.txt');

% Isolate file names from directory structure
fileNames = {fileDir.name}.';

% Convert cell array to table
% fileNames = cell2table(fileNames);

% Get number of files in directory
numFiles = length(fileNames);

%% find indicies of warmup and condition throws
inds_warm = find(contains(fileNames,"warmup"));

%% Variable initialization

% Get number of warmup files in directory
numFiles = length(inds_warm);

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

fc_master = table('VariableNames',{'fname','idx_raw','idx_interp'},'Size',[numFiles 3],'VariableTypes',{'string','doubleNaN','doubleNaN'});
stride_master = table('VariableNames',{'fname','idx_raw','idx_interp'},'Size',[numFiles 3],'VariableTypes',{'string','doubleNaN','doubleNaN'});
mav_master =  table('VariableNames',{'fname','idx_raw','idx_interp'},'Size',[numFiles 3],'VariableTypes',{'string','doubleNaN','doubleNaN'});
mer_master = table('VariableNames',{'fname','value','idx_raw','idx_interp'},'Size',[numFiles 4],'VariableTypes',{'string','doubleNaN','doubleNaN','doubleNaN'});
br_master = table('VariableNames',{'fname','value','idx_raw','idx_interp'},'Size',[numFiles 4],'VariableTypes',{'string','doubleNaN','doubleNaN','doubleNaN'});
mir_master = table('VariableNames',{'fname','value','idx_raw','idx_interp'},'Size',[numFiles 4],'VariableTypes',{'string','doubleNaN','doubleNaN','doubleNaN'});

events_master = table('VariableNames',{'fname','i_mav','i_fc','i_mer','i_br','i_mir','ckng'},'size',[numFiles 7],'VariableTypes',{'string','doubleNaN','doubleNaN','doubleNaN','doubleNaN','doubleNaN','doubleNaN'});

% ir_master = table('VariableNames',{'fname','ir_max'},'size',[numFiles 2],'VariableTypes',{'string','doubleNaN'});
% var_master = table('VariableNames',{'fname','ir_max'},'size',[numFiles 2],'VariableTypes',{'string','doubleNaN'});
%%
for index = 1:numFiles

%% Get sample frequency from .txt file
    txt = textscan(fopen(fileNames{inds_warm(index)}), '%q','Delimiter','//');
    txt = txt{1};
    fs  = str2double(txt{10});
    clear txt
%%
data = readtable(fileNames{inds_warm(index)},'VariableNamingRule','preserve');

% get max wrist velocity value + index for BR event mark
[max_wrist_velo,i_br] = max(data.WristLinVelo);
i_br = i_br + 2;
% get max shoulder external rotation value + index for MER event mark
[mer,i_mer] = max(data.ShldrER);
% get max shoulder internal rotation value + index for MIR event mark
    % MIR defined as first time after BR that the derivative of shoulder
    % rotation angle crosses zero (index.e. first local minimum following BR) 
i_mir = find(diff(sign(diff(data.ShldrER(i_mer:end))/(1/fs))),1) + i_mer-1;
mir = data.ShldrER(i_mir);
% get foot contact index for FC event mark
    % FC identified using lead ankle linear velocity in the Z (up/down)
    % direction and defined as the first frame after max ankle z velocity
    % that has a value greater than -0.1 m/s
thresh = -0.1;
[max_ank_vel, i_mav] = min(data.AnkVeloZ(round(fs*0.2):i_mer));
i_fc = i_mav + round(fs*0.2)-1;
% i_fc = find(data.AnkVeloZ(i_mav:end)>=thresh,1) + i_mav-1;
% clear i_mav

events_master(index,:) = [fileNames(inds_warm(index)), i_mav, i_fc, i_mer, i_br, i_mir, i_mer-i_fc];
%%
fc_master(index,:) = [fileNames(inds_warm(index)), i_fc, round(i_fc/height(data)*100)];
mav_master(index,:) = [fileNames(inds_warm(index)), i_mav, round(i_mav/height(data)*100)];
mer_master(index,:) = [fileNames(inds_warm(index)), mer, i_mer, round(i_mer/height(data)*100)];
br_master(index,:) = [fileNames(inds_warm(index)), max_wrist_velo, i_br, round(i_br/height(data)*100)];
mir_master(index,:) = [fileNames(inds_warm(index)), mir, i_mir, round(i_mir/height(data)*100)];

if index == 31
%% Override funky foot contacts
events_master.i_fc(31) = 14;
fc_master.idx_raw(31) = 14;
fc_master.idx_interp(31) = round(14/148*100);
end

if index == 35
%% Override funky foot contacts
events_master.i_fc(35) = 23;
fc_master.idx_raw(35) = 23;
fc_master.idx_interp(35) = round(23/164*100);
end

if index == 36
%% Override funky foot contacts
events_master.i_fc(36) = 29;
fc_master.idx_raw(36) = 29;
fc_master.idx_interp(36) = round(29/171*100);
end
%% Calculate energy transfer, generation, & absorption time series at shoulder
    if sum(isnan(data.Var_mom)) >= 1
        disp(['Possible missing data in ' fileNames{inds_warm(index)} '. (loop iteration ' num2str(index) ')'])
        continue
    end
    pwrtrans_stp = nan(height(data),1);
    pwrgen_stp = nan(height(data),1);
    trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.UaJFPp';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.ThorSTPd(j).*data.UaSTPp(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.ThorSTPd(j)) abs(data.UaSTPp(j))]);
            
            if abs(data.ThorSTPd(j)) == abs(data.UaSTPp(j))
                pwrtrans_stp(j) = data.UaSTPp(j);
                pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.ThorSTPd(j))
                pwrtrans_stp(j) = -data.ThorSTPd(j);
                pwrgen_stp(j) = data.UaSTPp(j) + data.ThorSTPd(j);
            elseif minSTP == abs(data.UaSTPp(j))
                pwrtrans_stp(j) = data.UaSTPp(j);
                pwrgen_stp(j) = data.UaSTPp(j) + data.ThorSTPd(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            pwrtrans_stp(j) = 0;
            pwrgen_stp(j) = data.UaSTPp(j) + data.ThorSTPd(j);
            
        end
    trans_net(j) = pwrtrans_jfp(j) + pwrtrans_stp(j);
    end
    
    %% Inerpolate shoulder power curves
    
    % Create original time axis
    time_old = (0:size(data)-1)'/fs;
    
    % Normalize original time axis
    time_old = time_old/time_old(end);
    
    % Create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);
    
    % Calcuate proportional measurement rate and event indexes for integrating interpolated curves
    fs_new = fs*100/height(data);
%     i_fc_interp = i_fc/height(data);
%     i_br_interp = i_br/height(data);
%     i_mer_interp = i_mer/height(data);
%     i_mir_interp = i_mir/height(data);
    
    % Interpolation
    shldr_trans_jfp(:,index) = interp1(time_old,pwrtrans_jfp,time,'spline');
    shldr_trans_stp(:,index) = interp1(time_old,pwrtrans_stp,time,'spline');
    shldr_gen_stp(:,index) = interp1(time_old,pwrgen_stp,time,'spline');
    shldr_trans_net(:,index) = interp1(time_old,trans_net,time,'spline');
    
    clear pwrtans_jfp pwrtrans_stp pwrgen_stp trans_net;
 %% Calculate total trans/gen/abs over phases @ shoulder
        
        shldr_trans_jfp_cp(index,1) = trapz(shldr_trans_jfp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        shldr_trans_jfp_ap(index,1) = trapz(shldr_trans_jfp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        shldr_trans_jfp_dp(index,1) = trapz(shldr_trans_jfp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        shldr_trans_stp_cp(index,1) = trapz(shldr_trans_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        shldr_trans_stp_ap(index,1) = trapz(shldr_trans_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        shldr_trans_stp_dp(index,1) = trapz(shldr_trans_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        shldr_trans_net_cp(index,1) = trapz(shldr_trans_net(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        shldr_trans_net_ap(index,1) = trapz(shldr_trans_net(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        shldr_trans_net_dp(index,1) = trapz(shldr_trans_net(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        shldr_gen_cp(index,1) = (trapz(shldr_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) + trapz(abs(shldr_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)))) / (2*fs_new);
        shldr_gen_ap(index,1) = (trapz(shldr_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) + trapz(abs(shldr_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)))) / (2*fs_new);
        shldr_gen_dp(index,1) = (trapz(shldr_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) + trapz(abs(shldr_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)))) / (2*fs_new);
        
        shldr_abs_cp(index,1) = (trapz(shldr_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) - trapz(abs(shldr_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)))) / (2*fs_new);
        shldr_abs_ap(index,1) = (trapz(shldr_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) - trapz(abs(shldr_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)))) / (2*fs_new);
        shldr_abs_dp(index,1) = (trapz(shldr_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) - trapz(abs(shldr_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)))) / (2*fs_new);   

    %% Calculate energy transfer, generation, & absorption at elbow
    
    pwrtrans_stp = nan(height(data),1);
    pwrgen_stp = nan(height(data),1);
    trans_net = nan(height(data),1);
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.FaJFPp';
    
    for k = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.UaSTPd(k).*data.FaSTPp(k));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.UaSTPd(k)) abs(data.FaSTPp(k))]);
            
            if abs(data.UaSTPd(k)) == abs(data.FaSTPp(k))
                pwrtrans_stp(k) = data.FaSTPp(k);
                pwrgen_stp(k) = 0;
            elseif minSTP == abs(data.UaSTPd(k))
                pwrtrans_stp(k) = -data.UaSTPd(k);
                pwrgen_stp(k) = data.FaSTPp(k) + data.UaSTPd(k);
            elseif minSTP == abs(data.FaSTPp(k))
                pwrtrans_stp(k) = data.FaSTPp(k);
                pwrgen_stp(k) = data.FaSTPp(k) + data.UaSTPd(k);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            pwrtrans_stp(k) = 0;
            pwrgen_stp(k) = data.UaSTPd(k) + data.FaSTPp(k);
            
        end
    trans_net(k) = pwrtrans_jfp(k) + pwrtrans_stp(k);    
    end
      %% Interpolate elbow power curves
    
    elb_trans_jfp(:,index) = interp1(time_old,pwrtrans_jfp,time,'spline');
    elb_trans_stp(:,index) = interp1(time_old,pwrtrans_stp,time,'spline');
    elb_gen_stp(:,index) = interp1(time_old,pwrgen_stp,time,'spline');
    elb_trans_net(:,index) = interp1(time_old,trans_net,time,'spline');
    
    clear pwrtrans_jfp pwrtrans_stp pwrgen_stp trans_net 
    %% Calculate total trans/gen/abs over phases @ elbow
        
        elb_trans_jfp_cp(index,1) = trapz(elb_trans_jfp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        elb_trans_jfp_ap(index,1) = trapz(elb_trans_jfp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        elb_trans_jfp_dp(index,1) = trapz(elb_trans_jfp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        elb_trans_stp_cp(index,1) = trapz(elb_trans_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        elb_trans_stp_ap(index,1) = trapz(elb_trans_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        elb_trans_stp_dp(index,1) = trapz(elb_trans_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        elb_trans_net_cp(index,1) = trapz(elb_trans_net(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) / fs_new;
        elb_trans_net_ap(index,1) = trapz(elb_trans_net(mer_master.idx_interp(index):br_master.idx_interp(index),index)) / fs_new;
        elb_trans_net_dp(index,1) = trapz(elb_trans_net(br_master.idx_interp(index):mir_master.idx_interp(index),index)) / fs_new;
        
        elb_gen_cp(index,1) = (trapz(elb_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) + trapz(abs(elb_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)))) / (2*fs_new);
        elb_gen_ap(index,1) = (trapz(elb_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) + trapz(abs(elb_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)))) / (2*fs_new);
        elb_gen_dp(index,1) = (trapz(elb_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) + trapz(abs(elb_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)))) / (2*fs_new);
        
        elb_abs_cp(index,1) = (trapz(elb_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)) - trapz(abs(elb_gen_stp(fc_master.idx_interp(index):mer_master.idx_interp(index),index)))) / (2*fs_new);
        elb_abs_ap(index,1) = (trapz(elb_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)) - trapz(abs(elb_gen_stp(mer_master.idx_interp(index):br_master.idx_interp(index),index)))) / (2*fs_new);
        elb_abs_dp(index,1) = (trapz(elb_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)) - trapz(abs(elb_gen_stp(br_master.idx_interp(index):mir_master.idx_interp(index),index)))) / (2*fs_new);
    %% Get peak shoulder and elbow moments
%     ir_master(index,:) = [fileNames(inds_warm(index)), max(data.Var_mom(fc_master.idx_raw(index):br_master.idx_raw(index)))];
%     var_master(index,:) = [fileNames(inds_warm(index)), max(data.Ir_mom(fc_master.idx_raw(index):br_master.idx_raw(index)))];
end
fclose('all');
%% Override funky foot contacts
events_master.i_fc(31) = 14;
fc_master.idx_raw(31) = 14;
fc_master.idx_interp(31) = round(14/148*100);
%%
figure('color','w')
plot(time,shldr_trans_net(:,index));
yline(0)
xline([fc_master.idx_interp(index)/100 br_master.idx_interp(index)/100 mer_master.idx_interp(index)/100 mir_master.idx_interp(index)/100],'--',{'fc','br','mer','mir'});
title(fileNames(index))

histogram(fc_master.idx_interp);

%%

[test, itest] = min(fc_master.idx_interp);
[~, itest] = min(test)
fileNames(itest)

%%
data = readtable(fileNames{inds_warm(index)},'VariableNamingRule','preserve');
figure('color','w')
hold on
plot(data.("Sample #"),data.AnkVeloX)
plot(data.("Sample #"),data.AnkVeloZ)
yline(0)
xline([fc_master.idx_raw(index) br_master.idx_raw(index) mer_master.idx_raw(index) mir_master.idx_raw(index)],'--',{'fc','br','mer','mir'});
legend('x','z')
title(fileNames(index))

%%

data = readtable(fileNames{inds_warm(index)},'VariableNamingRule','preserve');
figure('color','w')
plot(data.("Sample #"),data.Var_mom)
yline(0)
xline([fc_master.idx_raw(index) br_master.idx_raw(index) mer_master.idx_raw(index) mir_master.idx_raw(index)],'--',{'fc','br','mer','mir'});
title(fileNames(index))

%% write csv

t = table(fileNames,shldr_trans_net_cp,elb_trans_net_cp,'VariableNames',{'fname','shldr_trans_net_cp','elb_trans_net_cp'});
writetable(t,'energy_transfer.csv');


