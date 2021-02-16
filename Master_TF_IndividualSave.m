addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
load('clusterChans.mat')
ft_warning off

%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.005; % try with 0.005
s_fstep = 0.2; % 0.005
cycles  = 12;

cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  = 0.5:s_fstep:20; 
cfg.width                   = cycles;
cfg_Tf.toi                  = -22:s_tstep:52; 
% toi this is extended before and after to deal with border effect of wavelet
cfg_Tf.keeptrials           = 'yes';



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
    
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-15 45];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear EEGOdor dataOdor
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-15 45];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    clear Time_Freq
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-15 15];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [15 45];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(files(subj).name(1:12),'TF_DN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(files(subj).name(1:12),'TF_DN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Motor Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp/OdorM_Night/';

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
    
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    clear EEGOdor dataOdor
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-15 45];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear Time_Freq_DA_Temp
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-15 45];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-15 15];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [15 45];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(files(subj).name(1:12),'TF_MN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(files(subj).name(1:12),'TF_MN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

