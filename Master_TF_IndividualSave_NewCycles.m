addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')
% addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')

% ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
ft_warning off

%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.1; % try with 0.005
s_fstep = 0.1; % 0.005
cycles  = 12;
minfreq = 0.5;
maxfreq = 20;


cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  = minfreq:s_fstep:maxfreq; 
% cfg_Tf.width                = cycles;
v_timeWindows               = 5:(-4.5/numel(cfg_Tf.foi)):0.5;
cfg_Tf.t_ftimwin            = v_timeWindows;   
cfg_Tf.toi                  = -11:s_tstep:66; 
% toi this is extended before and after to deal with border effect of wavelet
cfg_Tf.keeptrials           = 'yes';



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/DNight/';
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

% filepath = 'D:\Thesis_Publication\Epoched_90SecTrial_MastoidRef-Interp_NewEpochs\DNight\';
% savepath = 'D:\TF_Calculation_90SecTrial\';

files = dir(strcat(filepath,'*.set'));


for subj = 1:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
%     addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    addpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
%     rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    rmpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1\'))
    
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-5 60];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear EEGOdor dataOdor
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-5 60];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    clear Time_Freq
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [25 60];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(savepath,files(subj).name(1:12),'TF_NewCycles_DN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_NewCycles_DN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Motor Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/MNight/';
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

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
    cfg.latency = [-5 60];
    Time_Freq = ft_selectdata(cfg, Time_Freq_DA_Temp);
    
    clear Time_Freq_DA_Temp
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-5 60];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.latency = [25 60];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(savepath,files(subj).name(1:12),'TF_NewCycles_MN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_NewCycles_MN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

