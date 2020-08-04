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


%% filter to only 240 Hz files
% files240 = nan(numFiles,1);
% 
% for i = 1:numFiles
%     txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
%     txt = txt{1};
%     fs  = str2double(txt{10});
%     clear txt
%     if fs > 200
%         files240(i) = 1;
%     else
%         files240(i) = 0;
%     end
% end
% fclose('all');
% fileNames240 = fileNames.fileNames(files240 == 1,:);
% numFiles240 = length(fileNames240);

%% preallocate variable arrays

% event timing master arrays
sfc_time = nan(numFiles,1);
mer_time = nan(numFiles,1);
br_time = nan(numFiles,1);
mir_time = nan(numFiles,1);

bodyheight = nan(numFiles,1);
bodymass = nan(numFiles,1);
measRate = nan(numFiles,1);

% total energy gen/abs/trans over phases of interest master arrays
    % back knee
    Bknee_trns_sp = nan(numFiles,1);
    Bknee_trns_cp = nan(numFiles,1);
    Bknee_trns_ap = nan(numFiles,1);

    Bknee_gen_sp = nan(numFiles,1);
    Bknee_gen_cp = nan(numFiles,1);
    Bknee_gen_ap = nan(numFiles,1);

    Bknee_abs_sp = nan(numFiles,1);
    Bknee_abs_cp = nan(numFiles,1);
    Bknee_abs_ap = nan(numFiles,1);

    % front knee
    Fknee_trns_sp = nan(numFiles,1);
    Fknee_trns_cp = nan(numFiles,1);
    Fknee_trns_ap = nan(numFiles,1);

    Fknee_gen_sp = nan(numFiles,1);
    Fknee_gen_cp = nan(numFiles,1);
    Fknee_gen_ap = nan(numFiles,1);

    Fknee_abs_sp = nan(numFiles,1);
    Fknee_abs_cp = nan(numFiles,1);
    Fknee_abs_ap = nan(numFiles,1);
    
    % back hip
    Bhip_trns_sp = nan(numFiles,1);
    Bhip_trns_cp = nan(numFiles,1);
    Bhip_trns_ap = nan(numFiles,1);

    Bhip_gen_sp = nan(numFiles,1);
    Bhip_gen_cp = nan(numFiles,1);
    Bhip_gen_ap = nan(numFiles,1);

    Bhip_abs_sp = nan(numFiles,1);
    Bhip_abs_cp = nan(numFiles,1);
    Bhip_abs_ap = nan(numFiles,1);
    
    % front hip
    Fhip_trns_sp = nan(numFiles,1);
    Fhip_trns_cp = nan(numFiles,1);
    Fhip_trns_ap = nan(numFiles,1);

    Fhip_gen_sp = nan(numFiles,1);
    Fhip_gen_cp = nan(numFiles,1);
    Fhip_gen_ap = nan(numFiles,1);

    Fhip_abs_sp = nan(numFiles,1);
    Fhip_abs_cp = nan(numFiles,1);
    Fhip_abs_ap = nan(numFiles,1);
    
    % thorax-pelvis
    
    PelvThor_trns_sp = nan(numFiles,1);
    PelvThor_trns_cp = nan(numFiles,1);
    PelvThor_trns_ap = nan(numFiles,1);

    PelvThor_gen_sp = nan(numFiles,1);
    PelvThor_gen_cp = nan(numFiles,1);
    PelvThor_gen_ap = nan(numFiles,1);

    PelvThor_abs_sp = nan(numFiles,1);
    PelvThor_abs_cp = nan(numFiles,1);
    PelvThor_abs_ap = nan(numFiles,1);

    % shoulder
    shldr_trns_cp = nan(numFiles,1);
    shldr_trns_ap = nan(numFiles,1);
    shldr_trns_dp = nan(numFiles,1);

    shldr_gen_cp = nan(numFiles,1);
    shldr_gen_ap = nan(numFiles,1);
    shldr_gen_dp = nan(numFiles,1);

    shldr_abs_cp = nan(numFiles,1);
    shldr_abs_ap = nan(numFiles,1);
    shldr_abs_dp = nan(numFiles,1);
    
    % elbow
    elb_trns_cp = nan(numFiles,1);
    elb_trns_ap = nan(numFiles,1);
    elb_trns_dp =nan(numFiles,1);

    elb_gen_cp = nan(numFiles,1);
    elb_gen_ap = nan(numFiles,1);
    elb_gen_dp = nan(numFiles,1);

    elb_abs_cp = nan(numFiles,1);
    elb_abs_ap = nan(numFiles,1);
    elb_abs_dp = nan(numFiles,1);

% joint gen/abs/trans interpolation masters for plotting
Bknee_gen = nan(numFiles,100);
Bknee_trns = nan(numFiles,100);

Bhip_gen = nan(numFiles,100);
Bhip_trns = nan(numFiles,100);

Fknee_gen = nan(numFiles,100);
Fknee_trns = nan(numFiles,100);

Fhip_trns = nan(numFiles,100);
Fhip_gen = nan(numFiles,100);

PelvThor_trns = nan(numFiles,100);
PelvThor_gen = nan(numFiles,100);

Shldr_trns = nan(numFiles,100);
Shldr_gen = nan(numFiles,100);

Elb_trns = nan(numFiles,100);
Elb_gen = nan(numFiles,100);

% peak rates of gen/abs/trans
pBknee_gen = nan(numFiles,1);
pBknee_abs = nan(numFiles,1);
pBknee_trns = nan(numFiles,1);

pBhip_gen = nan(numFiles,1);
pBhip_abs = nan(numFiles,1);
pBhip_trns = nan(numFiles,1);

pFknee_gen = nan(numFiles,1);
pFknee_abs = nan(numFiles,1);
pFknee_trns = nan(numFiles,1);

pFhip_trns = nan(numFiles,1);
pFhip_abs = nan(numFiles,1);
pFhip_gen = nan(numFiles,1);

pPelvThor_trns = nan(numFiles,1);
pPelvThor_abs = nan(numFiles,1);
pPelvThor_gen = nan(numFiles,1);

pShldr_trns = nan(numFiles,1);
pShldr_abs = nan(numFiles,1);
pShldr_gen = nan(numFiles,1);

pElb_trns = nan(numFiles,1);
pElb_abs = nan(numFiles,1);
pElb_gen = nan(numFiles,1);

% segment power interpolated masters for plotting
% back shank (proximal = ankle; distal = knee)
% Bshk_spP_master = nan(numFiles,100);
% Bshk_spD_master =  nan(numFiles,100);
% Bshk_sp_master =  nan(numFiles,100);
% 
% % front shank (proximal = ankle; distal = knee)
% Fshk_spP_master = nan(numFiles,100);
% Fshk_spD_master = nan(numFiles,100);
% Fshk_sp_master = nan(numFiles,100);
% 
% % back thigh (proximal = knee; distal = hip)
% Bthi_spP_master =  nan(numFiles,100);
% Bthi_spD_master = nan(numFiles,100);
% Bthi_sp_master = nan(numFiles,100);
% 
% % front thigh (proximal = knee; distal = hip)
% Fthi_spP_master = nan(numFiles,100);
% Fthi_spD_master = nan(numFiles,100);
% Fthi_sp_master = nan(numFiles,100);
% 
% % pelvis (proximal = back & front hips; distal = thorax)
% pelv_spBhip_master = nan(numFiles,100);
% pelv_spFhip_master = nan(numFiles,100);
% pelv_spThor_master = nan(numFiles,100);
% pelv_sp_master = nan(numFiles,100);

%%
for i = 1:numFiles
    %% get sampling frequency

    txt = textscan(fopen(fileNames.fileNames{i}), '%q','Delimiter','//');
        txt = txt{1};
        fs  = str2double(txt{10});
        clear txt

%% read in data, trim trial, and find events
    
    % load i-th file from directory
    data = extractData(fileNames.fileNames{i},'text',9);
    
    % get height & weight
    bodyheight(i,1) = data.BodyHeight(1);
    bodymass(i,1) = data.BodyMass(1);
    measRate(i,1) = fs;
    
    % find event indices
    events = find(data.VEM_0 == 1);
    
    % trim data to only include rows between pkh & mir
    data = data(events(2):events(5)+round(fs*.1),:);
    
    % redefine and label event indices
    events = find(data.VEM_0 == 1);
%     pkh = events(1);    % Peak knee height
    sfc = events(1);    % Stride foot contact
    mer = events(2);    % Max external rotation
    br = events(3);     % Ball release
    mir = events(4);    % Max internal rotation
    
    % event timing as a % between pkh and mir
%     sfc_time(i,1) = (sfc - pkh)/(height(data)-pkh);
%     mer_time(i,1) = (mer - pkh)/(height(data)-pkh);
%     br_time(i,1) = (br - pkh)/(height(data)-pkh);
%     mir_time(i,1) = (mir - pkh)/(height(data)-pkh);
    
    % event timing as a % between sfc and mir
    mer_time(i,1) = (mer - sfc)/(height(data)-sfc);
    br_time(i,1) = (br - sfc)/(height(data)-sfc);
    mir_time(i,1) = (mir - sfc)/(height(data)-sfc);
%% new times for interpolation
    % Create original time axis
    time_old = (0:size(data)-1)'/fs;
    
    % Normalize original time axis
    time_old = time_old/time_old(end);
    
    % Create new time axis for interpolation
    time = linspace(time_old(1),time_old(end),100);


      %% calculate segment powers
% 
%     % back shank (proximal = ankle; distal = knee)
%     Bshk_spP = data.BshankJfpP + data.BshankStpP;
%     Bshk_spD = data.BshankJfpD + data.BshankStpD;
%     Bshk_sp = Bshk_spP + Bshk_spD;
% 
%     % front shank (proximal = ankle; distal = knee)
%     Fshk_spP = data.FshankJfpP + data.FshankStpP;
%     Fshk_spD = data.FshankJfpD + data.FshankStpD;
%     Fshk_sp = Fshk_spP + Fshk_spD;
% 
%     % back thigh (proximal = knee; distal = hip)
%     Bthi_spP = data.BthighJfpP + data.BthighStpP;
%     Bthi_spD = data.BthighJfpD + data.BthighStpD;
%     Bthi_sp = Bthi_spP + Bthi_spD;
% 
%     % front thigh (proximal = knee; distal = hip)
%     Fthi_spP = data.FthighJfpP + data.FthighStpP;
%     Fthi_spD = data.FthighJfpD + data.FthighStpD;
%     Fthi_sp = Fthi_spP + Fthi_spD;
% 
%     % pelvis (proximal = back & front hips; distal = thorax)
%     pelv_spBhip = data.PelvJfpBhip + data.PelvStpBhip;
%     pelv_spFhip = data.PelvJfpFhip + data.PelvStpFhip;
%     pelv_spThor = data.PelvJfpThor + data.PelvStpThor;
%     pelv_sp = pelv_spBhip + pelv_spFhip + pelv_spThor;
% 
%     % thorax (proximal = pelvis; distal = throwing shoulder)
%     thor_spPelv = data.ThorJfpPelv + data.ThorStpPelv;
%     thor_spThrow = data.ThorJfpThrow + data.ThorStpThrow;
%     thor_spGlove = data.ThorJfpGlove + data.ThorStpGlove;
%     thor_sp = thor_spPelv + thor_spThrow + thor_spGlove;
% 
%     % upper arm (proximal = shoulder; distal = elbow)
%     hum_spP = data.HumJfpP + data.HumStpP;
%     hum_spD = data.HumJfpD + data.HumStpD;
%     hum_sp = hum_spP + hum_spD;
% 
%     % forearm (proximal = elbow; distal = wrist)
%     fa_spP = data.FaJfpP + data.FaStpP;
%     fa_spD = data.FaJfpD + data.FaStpD;
%     fa_sp = fa_spP + fa_spD;

    %% Calculate energy transfer, generation, & absorption time series at back knee
    % (distal back shank; proximal back thigh)
    
    Bknee_pwrtrans_stp = nan(height(data),1);
    Bknee_pwrgen_stp = nan(height(data),1);
    Bknee_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.BthighJfpP';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.BshankStpD(j).*data.BthighStpP(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.BshankStpD(j)) abs(data.BthighStpP(j))]);
            
            if abs(data.BshankStpD(j)) == abs(data.BthighStpP(j))
                Bknee_pwrtrans_stp(j) = data.BthighStpP(j);
                Bknee_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.BshankStpD(j))
                Bknee_pwrtrans_stp(j) = -data.BshankStpD(j);
                Bknee_pwrgen_stp(j) = data.BthighStpP(j) + data.BshankStpD(j);
            elseif minSTP == abs(data.BthighStpP(j))
                Bknee_pwrtrans_stp(j) = data.BthighStpP(j);
                Bknee_pwrgen_stp(j) = data.BthighStpP(j) + data.BshankStpD(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Bknee_pwrtrans_stp(j) = 0;
            Bknee_pwrgen_stp(j) = data.BthighStpP(j) + data.BshankStpD(j);
            
        end
        
        Bknee_trans_net(j) = pwrtrans_jfp(j) + Bknee_pwrtrans_stp(j);
    end
    
    Bknee_trns(i,:) = interp1(time_old, Bknee_trans_net, time, 'spline');
    Bknee_gen(i,:) = interp1(time_old, Bknee_pwrgen_stp, time, 'spline');
    
    %Bknee_trns_sp(i,1) = trapz(Bknee_trans_net(pkh:sfc)) / fs;
    Bknee_trns_cp(i,1) = trapz(Bknee_trans_net(sfc:mer)) / fs;
    Bknee_trns_ap(i,1) = trapz(Bknee_trans_net(mer:br)) / fs;
    
    %Bknee_gen_sp(i,1) = (trapz(Bknee_pwrgen_stp(pkh:sfc)) + trapz(abs(Bknee_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Bknee_gen_cp(i,1) = (trapz(Bknee_pwrgen_stp(sfc:mer)) + trapz(abs(Bknee_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Bknee_gen_ap(i,1) = (trapz(Bknee_pwrgen_stp(mer:br)) + trapz(abs(Bknee_pwrgen_stp(mer:br)))) / (2 * fs);
    
    %Bknee_abs_sp(i,1) = (trapz(Bknee_pwrgen_stp(pkh:sfc)) - trapz(abs(Bknee_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Bknee_abs_cp(i,1) = (trapz(Bknee_pwrgen_stp(sfc:mer)) - trapz(abs(Bknee_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Bknee_abs_ap(i,1) = (trapz(Bknee_pwrgen_stp(mer:br)) - trapz(abs(Bknee_pwrgen_stp(mer:br)))) / (2 * fs);
    
    pBknee_abs(i,1) = min(Bknee_pwrgen_stp);
    pBknee_gen(i,1) = max(Bknee_pwrgen_stp);
    pBknee_trns(i,1) = max(Bknee_trans_net);
    %% Calculate energy transfer, generation, & absorption time series at front knee
    % (distal front shank; proximal front thigh)
    
    Fknee_pwrtrans_stp = nan(height(data),1);
    Fknee_pwrgen_stp = nan(height(data),1);
    Fknee_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.FthighJfpP';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.FshankStpD(j).*data.FthighStpP(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.FshankStpD(j)) abs(data.FthighStpP(j))]);
            
            if abs(data.FshankStpD(j)) == abs(data.FthighStpP(j))
                Fknee_pwrtrans_stp(j) = data.FthighStpP(j);
                Fknee_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.FshankStpD(j))
                Fknee_pwrtrans_stp(j) = -data.FshankStpD(j);
                Fknee_pwrgen_stp(j) = data.FthighStpP(j) + data.FshankStpD(j);
            elseif minSTP == abs(data.FthighStpP(j))
                Fknee_pwrtrans_stp(j) = data.FthighStpP(j);
                Fknee_pwrgen_stp(j) = data.FthighStpP(j) + data.FshankStpD(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Fknee_pwrtrans_stp(j) = 0;
            Fknee_pwrgen_stp(j) = data.FthighStpP(j) + data.FshankStpD(j);
            
        end
        
        Fknee_trans_net(j) = pwrtrans_jfp(j) + Fknee_pwrtrans_stp(j);
    end
    
    Fknee_trns(i,:) = interp1(time_old, Fknee_trans_net, time, 'spline');
    Fknee_gen(i,:) = interp1(time_old, Fknee_pwrgen_stp, time, 'spline');
    
    %Fknee_trns_sp(i,1) = trapz(Fknee_trans_net(pkh:sfc)) / fs;
    Fknee_trns_cp(i,1) = trapz(Fknee_trans_net(sfc:mer)) / fs;
    Fknee_trns_ap(i,1) = trapz(Fknee_trans_net(mer:br)) / fs;
    
    %Fknee_gen_sp(i,1) = (trapz(Fknee_pwrgen_stp(pkh:sfc)) + trapz(abs(Fknee_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Fknee_gen_cp(i,1) = (trapz(Fknee_pwrgen_stp(sfc:mer)) + trapz(abs(Fknee_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Fknee_gen_ap(i,1) = (trapz(Fknee_pwrgen_stp(mer:br)) + trapz(abs(Fknee_pwrgen_stp(mer:br)))) / (2 * fs);
    
    %Fknee_abs_sp(i,1) = (trapz(Fknee_pwrgen_stp(pkh:sfc)) - trapz(abs(Fknee_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Fknee_abs_cp(i,1) = (trapz(Fknee_pwrgen_stp(sfc:mer)) - trapz(abs(Fknee_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Fknee_abs_ap(i,1) = (trapz(Fknee_pwrgen_stp(mer:br)) - trapz(abs(Fknee_pwrgen_stp(mer:br)))) / (2 * fs);
    
    pFknee_abs(i,1) = min(Fknee_pwrgen_stp);
    pFknee_gen(i,1) = max(Fknee_pwrgen_stp);
    pFknee_trns(i,1) = max(Fknee_trans_net);
    %% Calculate energy transfer, generation, & absorption time series at back hip
    % (distal front thigh; proximal pelvis front hip)
    
    Bhip_pwrtrans_stp = nan(height(data),1);
    Bhip_pwrgen_stp = nan(height(data),1);
    Bhip_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.PelvJfpBhip';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.BthighStpD(j).*data.PelvStpBhip(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.BthighStpD(j)) abs(data.PelvStpBhip(j))]);
            
            if abs(data.BthighStpD(j)) == abs(data.PelvStpBhip(j))
                Bhip_pwrtrans_stp(j) = data.PelvStpBhip(j);
                Bhip_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.BthighStpD(j))
                Bhip_pwrtrans_stp(j) = -data.BthighStpD(j);
                Bhip_pwrgen_stp(j) = data.BthighStpD(j) + data.PelvStpBhip(j);
            elseif minSTP == abs(data.PelvStpBhip(j))
                Bhip_pwrtrans_stp(j) = data.PelvStpBhip(j);
                Bhip_pwrgen_stp(j) = data.BthighStpD(j) + data.PelvStpBhip(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Bhip_pwrtrans_stp(j) = 0;
            Bhip_pwrgen_stp(j) = data.BthighStpD(j) + data.PelvStpBhip(j);
            
        end
        
        Bhip_trans_net(j) = pwrtrans_jfp(j) + Bhip_pwrtrans_stp(j);
    end
    
    Bhip_trns(i,:) = interp1(time_old, Bhip_trans_net, time, 'spline');
    Bhip_gen(i,:) = interp1(time_old, Bhip_pwrgen_stp, time, 'spline');
    
    %Bhip_trns_sp(i,1) = trapz(Bhip_trans_net(pkh:sfc)) / fs;
    Bhip_trns_cp(i,1) = trapz(Bhip_trans_net(sfc:mer)) / fs;
    Bhip_trns_ap(i,1) = trapz(Bhip_trans_net(mer:br)) / fs;
    
    %Bhip_gen_sp(i,1) = (trapz(Bhip_pwrgen_stp(pkh:sfc)) + trapz(abs(Bhip_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Bhip_gen_cp(i,1) = (trapz(Bhip_pwrgen_stp(sfc:mer)) + trapz(abs(Bhip_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Bhip_gen_ap(i,1) = (trapz(Bhip_pwrgen_stp(mer:br)) + trapz(abs(Bhip_pwrgen_stp(mer:br)))) / (2 * fs);
    
    %Bhip_abs_sp(i,1) = (trapz(Bhip_pwrgen_stp(pkh:sfc)) - trapz(abs(Bhip_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Bhip_abs_cp(i,1) = (trapz(Bhip_pwrgen_stp(sfc:mer)) - trapz(abs(Bhip_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Bhip_abs_ap(i,1) = (trapz(Bhip_pwrgen_stp(mer:br)) - trapz(abs(Bhip_pwrgen_stp(mer:br)))) / (2 * fs);
    
    pBhip_abs(i,1) = min(Bhip_pwrgen_stp);
    pBhip_gen(i,1) = max(Bhip_pwrgen_stp);
    pBhip_trns(i,1) = max(Bhip_trans_net);
    %% Calculate energy transfer, generation, & absorption time series at front hip
    % (distal front thigh; proximal pelvis front hip)
    
    Fhip_pwrtrans_stp = nan(height(data),1);
    Fhip_pwrgen_stp = nan(height(data),1);
    Fhip_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.PelvJfpFhip';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.FthighStpD(j).*data.PelvStpFhip(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.FthighStpD(j)) abs(data.PelvStpFhip(j))]);
            
            if abs(data.FthighStpD(j)) == abs(data.PelvStpFhip(j))
                Fhip_pwrtrans_stp(j) = data.PelvStpFhip(j);
                Fhip_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.FthighStpD(j))
                Fhip_pwrtrans_stp(j) = -data.FthighStpD(j);
                Fhip_pwrgen_stp(j) = data.FthighStpD(j) + data.PelvStpFhip(j);
            elseif minSTP == abs(data.PelvStpFhip(j))
                Fhip_pwrtrans_stp(j) = data.PelvStpFhip(j);
                Fhip_pwrgen_stp(j) = data.FthighStpD(j) + data.PelvStpFhip(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Fhip_pwrtrans_stp(j) = 0;
            Fhip_pwrgen_stp(j) = data.FthighStpD(j) + data.PelvStpFhip(j);
            
        end
        
        Fhip_trans_net(j) = pwrtrans_jfp(j) + Fhip_pwrtrans_stp(j);
    end
    
    Fhip_trns(i,:) = interp1(time_old, Fhip_trans_net, time, 'spline');
    Fhip_gen(i,:) = interp1(time_old, Fhip_pwrgen_stp, time, 'spline');
    
    %Fhip_trns_sp(i,1) = trapz(Fhip_trans_net(pkh:sfc)) / fs;
    Fhip_trns_cp(i,1) = trapz(Fhip_trans_net(sfc:mer)) / fs;
    Fhip_trns_ap(i,1) = trapz(Fhip_trans_net(mer:br)) / fs;
    
    %Fhip_gen_sp(i,1) = (trapz(Fhip_pwrgen_stp(pkh:sfc)) + trapz(abs(Fhip_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Fhip_gen_cp(i,1) = (trapz(Fhip_pwrgen_stp(sfc:mer)) + trapz(abs(Fhip_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Fhip_gen_ap(i,1) = (trapz(Fhip_pwrgen_stp(mer:br)) + trapz(abs(Fhip_pwrgen_stp(mer:br)))) / (2 * fs);
    
    %Fhip_abs_sp(i,1) = (trapz(Fhip_pwrgen_stp(pkh:sfc)) - trapz(abs(Fhip_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    Fhip_abs_cp(i,1) = (trapz(Fhip_pwrgen_stp(sfc:mer)) - trapz(abs(Fhip_pwrgen_stp(sfc:mer)))) / (2 * fs);
    Fhip_abs_ap(i,1) = (trapz(Fhip_pwrgen_stp(mer:br)) - trapz(abs(Fhip_pwrgen_stp(mer:br)))) / (2 * fs);
    
    pFhip_abs(i,1) = min(Fhip_pwrgen_stp);
    pFhip_gen(i,1) = max(Fhip_pwrgen_stp);
    pFhip_trns(i,1) = max(Fhip_trans_net);
    %% Calculate energy transfer, generation, & absorption time series between pelvis & thorax
    % (distal pelvis; proximal thorax)
    
    PelvThor_pwrtrans_stp = nan(height(data),1);
    PelvThor_pwrgen_stp = nan(height(data),1);
    PelvThor_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.ThorJfpPelv';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.PelvStpThor(j).*data.ThorStpPelv(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.PelvStpThor(j)) abs(data.ThorStpPelv(j))]);
            
            if abs(data.PelvStpThor(j)) == abs(data.ThorStpPelv(j))
                PelvThor_pwrtrans_stp(j) = data.ThorStpPelv(j);
                PelvThor_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.PelvStpThor(j))
                PelvThor_pwrtrans_stp(j) = -data.PelvStpThor(j);
                PelvThor_pwrgen_stp(j) = data.PelvStpThor(j) + data.ThorStpPelv(j);
            elseif minSTP == abs(data.ThorStpPelv(j))
                PelvThor_pwrtrans_stp(j) = data.ThorStpPelv(j);
                PelvThor_pwrgen_stp(j) = data.PelvStpThor(j) + data.ThorStpPelv(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            PelvThor_pwrtrans_stp(j) = 0;
            PelvThor_pwrgen_stp(j) = data.PelvStpThor(j) + data.ThorStpPelv(j);
            
        end
        
        PelvThor_trans_net(j) = pwrtrans_jfp(j) + PelvThor_pwrtrans_stp(j);
    end
    
    PelvThor_trns(i,:) = interp1(time_old, PelvThor_trans_net, time, 'spline');
    PelvThor_gen(i,:) = interp1(time_old, PelvThor_pwrgen_stp, time, 'spline');
    
    %PelvThor_trns_sp(i,1) = trapz(PelvThor_trans_net(pkh:sfc)) / fs;
    PelvThor_trns_cp(i,1) = trapz(PelvThor_trans_net(sfc:mer)) / fs;
    PelvThor_trns_ap(i,1) = trapz(PelvThor_trans_net(mer:br)) / fs;
    
    %PelvThor_gen_sp(i,1) = (trapz(PelvThor_pwrgen_stp(pkh:sfc)) + trapz(abs(PelvThor_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    PelvThor_gen_cp(i,1) = (trapz(PelvThor_pwrgen_stp(sfc:mer)) + trapz(abs(PelvThor_pwrgen_stp(sfc:mer)))) / (2 * fs);
    PelvThor_gen_ap(i,1) = (trapz(PelvThor_pwrgen_stp(mer:br)) + trapz(abs(PelvThor_pwrgen_stp(mer:br)))) / (2 * fs);
    
    %PelvThor_abs_sp(i,1) = (trapz(PelvThor_pwrgen_stp(pkh:sfc)) - trapz(abs(PelvThor_pwrgen_stp(pkh:sfc)))) / (2 * fs);
    PelvThor_abs_cp(i,1) = (trapz(PelvThor_pwrgen_stp(sfc:mer)) - trapz(abs(PelvThor_pwrgen_stp(sfc:mer)))) / (2 * fs);
    PelvThor_abs_ap(i,1) = (trapz(PelvThor_pwrgen_stp(mer:br)) - trapz(abs(PelvThor_pwrgen_stp(mer:br)))) / (2 * fs);
    
    pPelvThor_abs(i,1) = min(PelvThor_pwrgen_stp);
    pPelvThor_gen(i,1) = max(PelvThor_pwrgen_stp);
    pPelvThor_trns(i,1) = max(PelvThor_trans_net);
    %% Calculate energy transfer, generation, & absorption at shoulder
    % (distal thorax; proximal humerus)
    
    Shldr_pwr_trans_stp = nan(height(data),1);
    Shldr_pwrgen_stp = nan(height(data),1);
    Shldr_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.HumJfpP';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.ThorStpThrow(j).*data.HumStpP(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.ThorStpThrow(j)) abs(data.HumStpP(j))]);
            
            if abs(data.ThorStpThrow(j)) == abs(data.HumStpP(j))
                Shldr_pwr_trans_stp(j) = data.HumStpP(j);
                Shldr_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.ThorStpThrow(j))
                Shldr_pwr_trans_stp(j) = -data.ThorStpThrow(j);
                Shldr_pwrgen_stp(j) = data.ThorStpThrow(j) + data.HumStpP(j);
            elseif minSTP == abs(data.HumStpP(j))
                Shldr_pwr_trans_stp(j) = data.HumStpP(j);
                Shldr_pwrgen_stp(j) = data.ThorStpThrow(j) + data.HumStpP(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Shldr_pwr_trans_stp(j) = 0;
            Shldr_pwrgen_stp(j) = data.ThorStpThrow(j) + data.HumStpP(j);    
        end
        Shldr_trans_net(j) = pwrtrans_jfp(j) + Shldr_pwr_trans_stp(j);
    end
    
    Shldr_trns(i,:) = interp1(time_old, Shldr_trans_net, time, 'spline');
    Shldr_gen(i,:) = interp1(time_old, Shldr_pwrgen_stp, time, 'spline');
    
    shldr_trns_cp(i,1) = trapz(Shldr_trans_net(sfc:mer)) / fs;
    shldr_trns_ap(i,1) = trapz(Shldr_trans_net(mer:br)) / fs;
    shldr_trns_dp(i,1) = trapz(Shldr_trans_net(br:mir+round(fs*.1))) / fs;
    
    shldr_gen_cp(i,1) = (trapz(Shldr_pwrgen_stp(sfc:mer)) + trapz(abs(Shldr_pwrgen_stp(sfc:mer)))) / (2 * fs);
    shldr_gen_ap(i,1) = (trapz(Shldr_pwrgen_stp(mer:br)) + trapz(abs(Shldr_pwrgen_stp(mer:br)))) / (2 * fs);
    shldr_gen_dp(i,1) = (trapz(Shldr_pwrgen_stp(br:mir+round(fs*.1))) + trapz(abs(Shldr_pwrgen_stp(br:mir+round(fs*.1))))) / (2 * fs);
    
    shldr_abs_cp(i,1) = (trapz(Shldr_pwrgen_stp(sfc:mer)) - trapz(abs(Shldr_pwrgen_stp(sfc:mer)))) / (2 * fs);
    shldr_abs_ap(i,1) = (trapz(Shldr_pwrgen_stp(mer:br)) - trapz(abs(Shldr_pwrgen_stp(mer:br)))) / (2 * fs);
    shldr_abs_dp(i,1) = (trapz(Shldr_pwrgen_stp(br:mir+round(fs*.1))) - trapz(abs(Shldr_pwrgen_stp(br:mir+round(fs*.1))))) / (2 * fs);
    
    pShldr_abs(i,1) = min(Shldr_pwrgen_stp);
    pShldr_gen(i,1) = max(Shldr_pwrgen_stp);
    pShldr_trns(i,1) = max(Shldr_trans_net);
    %% Calculate energy transfer, generation, & absorption at Elbow
    % (distal humerus; proximal forearm)
    
    Elb_pwr_trans_stp = nan(height(data),1);
    Elb_pwrgen_stp = nan(height(data),1);
    Elb_trans_net = nan(height(data),1);   
    
    % Power transfer by joint force power
    pwrtrans_jfp = data.FaJfpP';
    
    for j = 1:height(data)
        
        % Define if segment torque powers have the same or opposite sign
        Sign = sign(data.HumStpD(j).*data.FaStpP(j));
            
        if Sign == -1   % = opposite sign
            
            % Find the minimum of the two segment torque powers
            minSTP = min([abs(data.HumStpD(j)) abs(data.FaStpP(j))]);
            
            if abs(data.HumStpD(j)) == abs(data.FaStpP(j))
                Elb_pwr_trans_stp(j) = data.FaStpP(j);
                Elb_pwrgen_stp(j) = 0;
            elseif minSTP == abs(data.HumStpD(j))
                Elb_pwr_trans_stp(j) = -data.HumStpD(j);
                Elb_pwrgen_stp(j) = data.HumStpD(j) + data.FaStpP(j);
            elseif minSTP == abs(data.FaStpP(j))
                Elb_pwr_trans_stp(j) = data.FaStpP(j);
                Elb_pwrgen_stp(j) = data.HumStpD(j) + data.FaStpP(j);
            end
            
        elseif Sign == 1 || Sign == 0   % = same sign
            
            Elb_pwr_trans_stp(j) = 0;
            Elb_pwrgen_stp(j) = data.HumStpD(j) + data.FaStpP(j);
            
        end
        
        Elb_trans_net(j) = pwrtrans_jfp(j) + Elb_pwr_trans_stp(j);
    end
    
    Elb_trns(i,:) = interp1(time_old, Elb_trans_net, time, 'spline');
    Elb_gen(i,:) = interp1(time_old, Elb_pwrgen_stp, time, 'spline');
    
    elb_trns_cp(i,1) = trapz(Elb_trans_net(sfc:mer)) / fs;
    elb_trns_ap(i,1) = trapz(Elb_trans_net(mer:br)) / fs;
    elb_trns_dp(i,1) = trapz(Elb_trans_net(br:mir+round(fs*.1))) / fs;
    
    elb_gen_cp(i,1) = (trapz(Elb_pwrgen_stp(sfc:mer)) + trapz(abs(Elb_pwrgen_stp(sfc:mer)))) / (2 * fs);
    elb_gen_ap(i,1) = (trapz(Elb_pwrgen_stp(mer:br)) + trapz(abs(Elb_pwrgen_stp(mer:br)))) / (2 * fs);
    elb_gen_dp(i,1) = (trapz(Elb_pwrgen_stp(br:mir+round(fs*.1))) + trapz(abs(Elb_pwrgen_stp(br:mir+round(fs*.1))))) / (2 * fs);
    
    elb_abs_cp(i,1) = (trapz(Elb_pwrgen_stp(sfc:mer)) - trapz(abs(Elb_pwrgen_stp(sfc:mer)))) / (2 * fs);
    elb_abs_ap(i,1) = (trapz(Elb_pwrgen_stp(mer:br)) - trapz(abs(Elb_pwrgen_stp(mer:br)))) / (2 * fs);
    elb_abs_dp(i,1) = (trapz(Elb_pwrgen_stp(br:mir+round(fs*.1))) - trapz(abs(Elb_pwrgen_stp(br:mir+round(fs*.1))))) / (2 * fs);
    
    pElb_abs(i,1) = min(Elb_pwrgen_stp);
    pElb_gen(i,1) = max(Elb_pwrgen_stp);
    pElb_trns(i,1) = max(Elb_trans_net);
end


% create average interpolated masters
pStartRow = 1:3:length(Bhip_abs_ap)+1;
for i = 1:(length(pStartRow)-1)
    aveBknee_trns(i,:) = mean(Bknee_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveBknee_gen(i,:) = mean(Bknee_gen(pStartRow(i):pStartRow(i+1)-1,:));
    aveFknee_trns(i,:) = mean(Fknee_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveFknee_gen(i,:) = mean(Fknee_gen(pStartRow(i):pStartRow(i+1)-1,:));
    aveBhip_trns(i,:) = mean(Bhip_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveBhip_gen(i,:) = mean(Bhip_gen(pStartRow(i):pStartRow(i+1)-1,:));
    aveFhip_trns(i,:) = mean(Fhip_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveFhip_gen(i,:) = mean(Fhip_gen(pStartRow(i):pStartRow(i+1)-1,:));
    avePelvThor_trns(i,:) = mean(PelvThor_trns(pStartRow(i):pStartRow(i+1)-1,:));
    avePelvThor_gen(i,:) = mean(PelvThor_gen(pStartRow(i):pStartRow(i+1)-1,:));
    aveShldr_trns(i,:) = mean(Shldr_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveShldr_gen(i,:) = mean(Shldr_gen(pStartRow(i):pStartRow(i+1)-1,:));
    aveElb_trns(i,:) = mean(Elb_trns(pStartRow(i):pStartRow(i+1)-1,:));
    aveElb_gen(i,:) = mean(Elb_gen(pStartRow(i):pStartRow(i+1)-1,:));
end

%% plot front knee & hip generation/transfer
close all

fig1 = figure('name', 'Front Leg Energy Generation & Transfer', 'color','w');
    subplot(2, 2, 1)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveFhip_trns);
        title("Front Hip Transfer");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2,2,2)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveFhip_gen);
        title("Front Hip Generation");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 2, 3)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveFknee_trns);
        title("Front Knee Transfer");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 2, 4)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveFknee_gen);
        title("Front Knee Generation");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    ax4 = axes(fig1, 'visible', 'off', 'fontname','times new roman');
    ax4.Title.Visible = 'on';
    ax4.XLabel.Visible = 'on';
    ax4.YLabel.Visible = 'on';
    labelx = xlabel(ax4, 'Time(%)');
    labely = ylabel(ax4, 'Power (W)');
    labely.Position(1) = labely.Position(1)*1.1; % move y label a bit further away from y axis
   
%% plot back knee & hip generation/transfer

fig2 = figure('name', 'Back Leg Energy Generation & Tansfer', 'color','w');
    hold on
    subplot(2, 2, 1)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveBhip_trns);
        title("Back Hip Transfer");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 2, 2)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveBhip_gen);
        title("Back Hip Generation");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 2, 3)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveBknee_trns,'prctile',25);
        title("Back Knee Transfer");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 2, 4)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, aveBknee_gen,'prctile',25);
        title("Back Knee Generation");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1])
        xticklabels(["0" "25" "50" "75" "100"]); 
    
 %% plot thoracopelvic generation/transfer

fig3 = figure('name', 'Thoraco-Pelvic Energy Generation & Transfer', 'color','w');
    hold on
    subplot(2, 1, 1)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, avePelvThor_trns);
        title("Thorax-from-Pelvis Transfer");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 1, 2)
        yline(0);
        xline(mean(sfc_time));
        xline(mean(mer_time));
        xline(mean(br_time));
        plot_distribution_prctile(time, avePelvThor_gen);
        title("Thorax-from-Pelvis Generation");
        set(gca,'fontname','times new roman');
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
        
%% plot shoulder generation/transfer

fig4 = figure('name','Shoulder Energy Generation & Transfer', 'color','w');
    hold on
    subplot(2, 1, 1)
        yline(0);
        xline(mean(mer_time),'linestyle','--','label','MER','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(br_time),'linestyle','--','label','BR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(mir_time),'linestyle','--','label','MIR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        plot_distribution_prctile(time, Shldr_trns);
        title("Net Shoulder Transfer",'fontsize',14);
        set(gca,'fontname','times new roman','fontsize',10);
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 1, 2)
        yline(0);
         xline(mean(mer_time),'linestyle','--','label','MER','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(br_time),'linestyle','--','label','BR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(mir_time),'linestyle','--','label','MIR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        plot_distribution_prctile(time, Shldr_gen);
        title("Shoulder Generation",'fontsize',14);
        set(gca,'fontname','times new roman','fontsize',10);
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    ax4 = axes(fig4, 'visible', 'off', 'fontname','times new roman');
    ax4.Title.Visible = 'on';
    ax4.XLabel.Visible = 'on';
    ax4.YLabel.Visible = 'on';
    labelx = xlabel(ax4, 'Time (%)','fontsize',12);
    labely = ylabel(ax4, 'Power (W)','fontsize',12);
    labely.Position(1) = labely.Position(1)*1.4; % move y label a bit further away from y axis
        
%% plot elbow generation/transfer

fig5 = figure('name', 'Elbow Energy Generation & Transfer', 'color','w');
    hold on
    subplot(2, 1, 1)
        yline(0);
%         xline(mean(sfc_time));
        xline(mean(mer_time),'linestyle','--','label','MER','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(br_time),'linestyle','--','label','BR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        xline(mean(mir_time),'linestyle','--','label','MIR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman','fontsize',8);
        plot_distribution_prctile(time, aveElb_trns);
        title("Net Elbow Transfer",'fontsize',14);
        set(gca,'fontname','times new roman','fontsize',10);
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    subplot(2, 1, 2)
        yline(0);
%         xline(mean(sfc_time));
        xline(mean(mer_time),'linestyle','--','label','MER','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman',...
            'fontsize',8,'labelverticalalignment','bottom');
        xline(mean(br_time),'linestyle','--','label','BR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman',...
            'fontsize',8,'labelverticalalignment','bottom');
        xline(mean(mir_time),'linestyle','--','label','MIR','labelhorizontalalignment','center',...
            'labelorientation','horizontal','fontname','times new roman',...
            'fontsize',8,'labelverticalalignment','bottom');
        plot_distribution_prctile(time, aveElb_gen);
        title("Elbow Generation",'fontsize',14);
        set(gca,'fontname','times new roman','fontsize',10);
        xticks([0 .25 .50 .75 1]);
        xticklabels(["0" "25" "50" "75" "100"]);
    ax5 = axes(fig5, 'visible', 'off', 'fontname','times new roman');
    ax5.Title.Visible = 'on';
    ax5.XLabel.Visible = 'on';
    ax5.YLabel.Visible = 'on';
    labelx = xlabel(ax5, 'Time (%)','fontsize',12);
    labely = ylabel(ax5, 'Power (W)','fontsize',12);
    labely.Position(1) = labely.Position(1)*1.4; % move y label a bit further away from y axis
        
%% write gen/abs/trans table

fullbodyGAT_master = table();
participants = fileNames.fileNames;
pID = (1:numFiles/3)';
pID = repelem(pID,3);

trial = [1 2 3]';
trial = repmat(trial,numFiles/3,1);

fullbodyGAT_master = ...
    addvars(fullbodyGAT_master,...
    participants,...
    pID,...
    trial,...
    measRate,...
    bodyheight,...
    bodymass,...
    measRate,...
    Bknee_trns_sp,...
    Bknee_trns_cp,...
    Bknee_trns_ap,...
    Bknee_gen_sp,...
    Bknee_gen_cp,...
    Bknee_gen_ap,...
    Bknee_abs_sp,...
    Bknee_abs_cp,...
    Bknee_abs_ap,...
    Fknee_trns_sp,...
    Fknee_trns_cp,...
    Fknee_trns_ap,...
    Fknee_gen_sp,...
    Fknee_gen_cp,...
    Fknee_gen_ap,...
    Fknee_abs_sp,...
    Fknee_abs_cp,...
    Fknee_abs_ap,...
    Bhip_trns_sp,...
    Bhip_trns_cp,...
    Bhip_trns_ap,...
    Bhip_gen_sp,...
    Bhip_gen_cp,...
    Bhip_gen_ap,...
    Bhip_abs_sp,...
    Bhip_abs_cp,...
    Bhip_abs_ap,...
    Fhip_trns_sp,...
    Fhip_trns_cp,...
    Fhip_trns_ap,...
    Fhip_gen_sp,...
    Fhip_gen_cp,...
    Fhip_gen_ap,...
    Fhip_abs_sp,...
    Fhip_abs_cp,...
    Fhip_abs_ap,...
    PelvThor_trns_sp,...
    PelvThor_trns_cp,...
    PelvThor_trns_ap,...
    PelvThor_gen_sp,...
    PelvThor_gen_cp,...
    PelvThor_gen_ap,...
    PelvThor_abs_sp,...
    PelvThor_abs_cp,...
    PelvThor_abs_ap,...
    shldr_trns_cp,...
    shldr_trns_ap,...
    shldr_trns_dp,...
    shldr_gen_cp,...
    shldr_gen_ap,...
    shldr_gen_dp,...
    shldr_abs_cp,...
    shldr_abs_ap,...
    shldr_abs_dp,...
    elb_trns_cp,...
    elb_trns_ap,...
    elb_trns_dp,...
    elb_gen_cp,...
    elb_gen_ap,...
    elb_gen_dp,...
    elb_abs_cp,...
    elb_abs_ap,...
    elb_abs_dp,...
    pBknee_gen,...
    pBknee_abs,... 
    pBknee_trns,...
    pBhip_gen,...
    pBhip_abs,...
    pBhip_trns,...
    pFknee_gen,...
    pFknee_abs,...
    pFknee_trns,...
    pFhip_trns,...
    pFhip_abs,...
    pFhip_gen,...
    pPelvThor_trns,...
    pPelvThor_abs,...
    pPelvThor_gen,...
    pShldr_trns,...
    pShldr_abs,...
    pShldr_gen,...
    pElb_trns,...
    pElb_abs,...
    pElb_gen);

fullbodyGAT_master = renamevars(fullbodyGAT_master,...
    ["participants" "bodyheight" "bodymass"],...
    ["fname" "height" "mass"]);
% 
% writetable(fullbodyGAT_master,'GAT_master9-7.csv');
% 
% 
% %% write average table
pStartRow = 1:3:height(fullbodyGAT_master)+1;
aveGAT_master = array2table(nan((numFiles/3),width(fullbodyGAT_master)));

for i = 1:(length(pStartRow)-1)
    aveGAT_master{i,2:end} = mean(fullbodyGAT_master{pStartRow(i):pStartRow(i+1)-1,2:end});
end

% give average table same variable names and file names as full master
aveGAT_master.Properties.VariableNames = fullbodyGAT_master.Properties.VariableNames;

fStartRow = 1:3:length(fileNames.fileNames);
aveGAT_master.fname = fileNames.fileNames(fStartRow);
% 
% writetable(aveGAT_master,'aveGAT_master9-7.csv');
% 
% fclose('all');
