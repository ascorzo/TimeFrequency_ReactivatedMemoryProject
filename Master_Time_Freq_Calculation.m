addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off
%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.05; % try with 0.05
s_fstep = 0.2; % 0.05

cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet'; %'mtmconvol';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  =  0.5:s_fstep:20; % Reduce to 20
cfg.width                   = 12;
cfg_Tf.toi                  = -25:s_tstep:55;
cfg_Tf.keeptrials           = 'yes';


%--------------------------------------------------------------------------
% Parameters for Baseline Correction 
%--------------------------------------------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-15 45];
cfg_Bas.baselinetype        = 'zscore';

%--------------------------------------------------------------------------
% Define clusters of channels 
%--------------------------------------------------------------------------
Clust.left_frontal = {...
    'E15', 'E16', 'E11', 'E18', 'E19', 'E22', 'E23', 'E24', 'E26', ...
    'E27', 'E33', 'E38'};
Clust.right_frontal = {...
    'E15', 'E16', 'E11', 'E10', 'E4', 'E9', 'E3', 'E124', 'E2', ...
    'E123', 'E122', 'E121'};
Clust.frontal = {...
    'E3', 'E4', 'E9', 'E10', 'E11', 'E15', 'E16', 'E18', 'E19', ...
    'E22', 'E23', 'E24', 'E124'};
Clust.left_central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55'};
Clust.right_central = {...
    'E6', 'E55', 'E112', 'E106', 'E105', 'E80', 'E87', 'E79'};
Clust.central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55', 'E79', ...
    'E80', 'E87', 'E105', 'E106', 'E112'};
Clust.left_temporal = {...
    'E46', 'E51', 'E45', 'E50', 'E58', 'E56', 'E63'};
Clust.right_temporal = {...
    'E108', 'E102', 'E101', 'E97', 'E96', 'E99', 'E107'};
Clust.left_parietal = {...
    'E53', 'E61', 'E62', 'E72', 'E67', 'E52', 'E60', 'E59', 'E66', ...
    'E65', 'E64', 'E68'};
Clust.right_patietal = {...
    'E62', 'E72', 'E78', 'E77', 'E86', 'E85', 'E84', 'E92', 'E91', ...
    'E90', 'E95', 'E94'};
Clust.parietal = {...
    'E52', 'E61', 'E62', 'E59', 'E60', 'E67', 'E66', 'E72', 'E78', ...
    'E77', 'E86', 'E85', 'E84', 'E92', 'E91','E53'};
Clust.left_occipital = {...
    'E71', 'E70', 'E75', 'E74', 'E69', 'E73'};
Clust.right_occipital = {...
    'E75', 'E76', 'E82', 'E83', 'E88', 'E89'};
Clust.occipital = {...
    'E71', 'E70', 'E74', 'E69', 'E73', 'E75', 'E76', 'E83', 'E82', ...
    'E89', 'E88'};


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp/OdorD_Night/';

files = dir(strcat(filepath,'*.set'));


for subj = 1:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
    addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    
 
    -----------Selection of data Channels--------------------------------- 
    dataOdor = ft_selectdata(cfg_Sel,dataOdor);
    
    
    -----------Time-Frequency Calculation---------------------------------
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    --------------Select time of interest---------------------------------
    cfg = [];
    cfg.latency = [-15 45];
    Time_Freq_DA{subj}  = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    %-----------baseline Correction ---------------------------------------
    
    Time_Freq_DA_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq_DA{subj} );
    
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_DA_Baseline,  'trialinfo');
    Time_Freq_DA_Mean = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    
    %-----------Separate conditions -----------------------  
    cfg = [];
    cfg.latency = [-15 15];
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_DA_Mean);
    
    cfg = [];
    cfg.latency = [15 45];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_DA_Mean);
    Time_Freq_Vehicle.time = Time_Freq_DA_Mean.time;
    
    %-----------Combine Mean of all frontal channels-----------------------
    [~,~,ind2] = intersect(Clust.frontal,Time_Freq_Cue.label);
    Time_Freq_Cue_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
    
    [~,~,ind2] = intersect(Clust.frontal,Time_Freq_Vehicle.label);
    Time_Freq_Vehicle_Frontal(subj,:,:) = ...
        squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));
    
    %-----------Combine Mean of all central channels-----------------------
    [~,~,ind2] = intersect(Clust.central,Time_Freq_Cue.label);
    Time_Freq_Cue_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
    
    [~,~,ind2] = intersect(Clust.central,Time_Freq_Vehicle.label);
    Time_Freq_Vehicle_Central(subj,:,:) = ...
        squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));

end



%% Calculation of power change in specific bands

SpindleBand     = [12 16];
DeltaBand       = [0.5 4];
time = Time_Freq_DA{1}.time ;
subjects = 1:numel(files);

%------------------------------------------------------
% For Cue Night
%------------------------------------------------------

SpindleIdx = find(Time_Freq_DA{subj}.freq>=SpindleBand(1) &...
    Time_Freq_DA{subj}.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq_DA{subj}.freq>=DeltaBand(1) &...
    Time_Freq_DA{subj}.freq<=DeltaBand(2));

% -----------------For Cue----------------------------------
v_Spindle_Cue_Frontal = squeeze(mean(Time_Freq_Cue_Frontal(subjects,SpindleIdx,:),2));

v_Delta_Cue_Frontal = squeeze(mean(Time_Freq_Cue_Frontal(subjects,DeltaIdx,:),2));


v_Spindle_Cue_Central = squeeze(mean(Time_Freq_Cue_Central(subjects,SpindleIdx,:),2));

v_Delta_Cue_Central = squeeze(mean(Time_Freq_Cue_Central(subjects,DeltaIdx,:),2));
    
    

% -----------------For Vehicle----------------------------------
v_Spindle_Vehicle_Frontal = squeeze(mean(Time_Freq_Vehicle_Frontal(subjects,SpindleIdx,:),2));

v_Delta_Vehicle_Frontal = squeeze(mean(Time_Freq_Vehicle_Frontal(subjects,DeltaIdx,:),2));


v_Spindle_Vehicle_Central = squeeze(mean(Time_Freq_Vehicle_Central(subjects,SpindleIdx,:),2));

v_Delta_Vehicle_Central = squeeze(mean(Time_Freq_Vehicle_Central(subjects,DeltaIdx,:),2));

%% plot CUE Night

%p_plot_TimeFreq_All

 %% plot Wilcoxon test
 
 figure
 
 v_time = Time_Freq_Cue.time ;
 addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
 
% eeglab nogui
 addpath('./Scripts_Wilc')
 
 Title = 'Odor D Night Ref-Interp';
 
 subplot(2,2,1)
 f_WilcTest('Spindle Frontal',...
     'Time(sec)',' ','Cue','Vehicle',...
     v_Spindle_Cue_Frontal,...
     v_Spindle_Vehicle_Frontal,...
     v_time,'-r',1,[])
 
 subplot(2,2,3)
 ValidTimeCue = v_time(~isnan(v_Delta_Cue_Frontal(1,:)));
 ValidTimeVehicle = v_time(~isnan(v_Delta_Vehicle_Frontal(1,:)));
 
 ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
 start_sample = ceil(ValidTime(1)-v_time(1)/s_tstep);
 end_sample = floor(ValidTime(end)-v_time(1)/s_tstep);
 
 f_WilcTest('Delta Frontal',...
     'Time(sec)',' ','Cue','Vehicle',...
     v_Delta_Cue_Frontal,...
     v_Delta_Vehicle_Frontal,...
     v_time,'-r',start_sample,end_sample)
 
 subplot(2,2,2)
 f_WilcTest('Spindle Central',...
     'Time(sec)',' ','Cue','Vehicle',...
     v_Spindle_Cue_Central,...
     v_Spindle_Vehicle_Central,...
     v_time,'-r',1,[])
 
 subplot(2,2,4)
 ValidTimeCue = v_time(~isnan(v_Delta_Cue_Central(1,:)));
 ValidTimeVehicle = v_time(~isnan(v_Delta_Vehicle_Central(1,:)));
 
 ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
 start_sample = ceil(ValidTime(1)-v_time(1)/s_tstep);
 end_sample = floor(ValidTime(end)-v_time(1)/s_tstep);
 
 f_WilcTest('Delta Central',...
     'Time(sec)',' ','Cue','Vehicle',...
     v_Delta_Cue_Central,...
     v_Delta_Vehicle_Central,...
     v_time,'-r',start_sample,end_sample)

 sgtitle(Title)
 
 
 set(gcf,'position',[1,35,1920,926])
 
 saveas(gcf,strcat(Title,'.png'))