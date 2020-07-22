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

%% preallocate variable arrays

% event timing master arrays
sfc_time = nan(numFiles,1);
mer_time = nan(numFiles,1);
br_time = nan(numFiles,1);

% total energy gen/abs/trans over phases of interest master arrays
Bknee_gen_sp = nan(numFiles,1);
Bknee_trns_sp = nan(numFiles,1);

Bhip_gen_sp = nan(numFiles,1);
Bhip_trans_sp = nan(numFiles,1);

Fknee_gen_cp = nan(numFiles,1);
Fknee_trns_cp = nan(numFiles,1);

Fknee_trns_ap = nan(numFiles,1);
Fknee_gen_ap = nan(numFiles,1);

% joint gen/abs/trans interpolation masters for plotting
Bknee_gen = nan(numFiles,100);
Bknee_trns = nan(numFiles,100);

Bhip_gen = nan(numFiles,100);
Bhip_trans = nan(numFiles,100);

Fknee_gen = nan(numFiles,100);
Fknee_trns = nan(numFiles,100);

Fhip_trns = nan(numFiles,100);
Fhip_gen = nan(numFiles,100);

% segment power interpolated masters for plotting
% back shank (proximal = ankle; distal = knee)
Bshk_spP_master = nan(numFiles,100);
Bshk_spD_master =  nan(numFiles,100);
Bshk_sp_master =  nan(numFiles,100);

% front shank (proximal = ankle; distal = knee)
Fshk_spP_master = nan(numFiles,100);
Fshk_spD_master = nan(numFiles,100);
Fshk_sp_master = nan(numFiles,100);

% back thigh (proximal = knee; distal = hip)
Bthi_spP_master =  nan(numFiles,100);
Bthi_spD_master = nan(numFiles,100);
Bthi_sp_master = nan(numFiles,100);

% front thigh (proximal = knee; distal = hip)
Fthi_spP_master = nan(numFiles,100);
Fthi_spD_master = nan(numFiles,100);
Fthi_sp_master = nan(numFiles,100);

% pelvis (proximal = back & front hips; distal = thorax)
pelv_spBhip_master = nan(numFiles,100);
pelv_spFhip_master = nan(numFiles,100);
pelv_spThor_master = nan(numFiles,100);
pelv_sp_master = nan(numFiles,100);

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
    
    % find event indices
    events = find(data.VEM_0 == 1);
    
    % trim data to only include rows between pkh & mir
    data = data(events(1):events(end),:);
    
    % redefine and label event indices
    events = find(data.VEM_0 == 1);
    pkh = events(1);    % Peak knee height
    sfc = events(2);    % Stride foot contact
    mer = events(3);    % Max external rotation
    br = events(4);     % Ball release
    mir = events(5);    % Max internal rotation
    
    % event timing as a % between pkh and mir
    sfc_time(i,1) = (sfc - pkh)/(mir-pkh);
    mer_time(i,1) = (mer - pkh)/(mir-pkh);
    br_time(i,1) = (br - pkh)/(mir-pkh);
    
    % event timing as a % between sfc and mir
%     mer_time(i,1) = (mer - sfc)/(mir - sfc);
%     br_time(i,1) = (br - sfc)/(mir - sfc);
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
end
    
 %%

figure('color','w');
    hold on
    subplot(1, 2, 1)
    yline(0);
    xline(sfc);
    xline(br);
    p1 = plot(time, Fhip_trns);
    title("Front Hip Transfer");
    subplot(1, 2, 2)
    yline(0);
    xline(sfc);
    xline(br);
    p2 = plot(time, Fhip_gen);
    title("Front Hip Generation");
    set(gca,'fontname','times new roman');
   
    
    