addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
ft_warning off

%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.005; % try with 0.005
s_fstep = 0.005; % 0.005
cycles  = 12;

cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  = 0.5:s_fstep:20; 
cfg_Tf.width                = cycles;
cfg_Tf.toi                  = -5:s_tstep:10; 
% toi this is extended before and after to deal with border effect of wavelet
cfg_Tf.keeptrials           = 'yes';



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/DNight/';
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/DNight/';

files = dir(strcat(filepath,'*.set'));


for subj = 1:numel(files)

    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
    addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))

    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Data_Odor = ft_selectdata(cfg, dataOdor);
    
    cfg = [];
    cfg.latency = [25 60];
    Data_Vehicle = ft_selectdata(cfg, dataOdor);
    Data_Vehicle.time = Data_Odor.time;

    clear EEGOdor dataOdor
    %-----------Select first Seconds with an extension period---------------
    cfg = [];
    cfg.latency = [-5 10];
    Data_Odor = ft_selectdata(cfg, Data_Odor);
    Data_Vehicle = ft_selectdata(cfg, Data_Vehicle);
    
    %-----------Time-Frequency Calculation---------------------------------
    Time_Freq_Odor  = ft_freqanalysis(cfg_Tf, Data_Odor);
    Time_Freq_Vehicle  = ft_freqanalysis(cfg_Tf, Data_Vehicle);

    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-1 5];
    Time_Freq_Odor  = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehicle  = ft_selectdata(cfg, Time_Freq_Vehicle);

    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-1 5];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Odor  = ft_freqbaseline(cfg_Bas, Time_Freq_Odor);
    Time_Freq_Vehicle  = ft_freqbaseline(cfg_Bas, Time_Freq_Vehicle);


    save(strcat(savepath,files(subj).name(1:12),'TF_DN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_DN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Motor Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/andrea/German_Study/Data/PreProcessed/Epoched_90SecTrial_MastoidRef-Interp_NewEpochs/MNight/';
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/MNight/';

files = dir(strcat(filepath,'*.set'));


for subj = 1:numel(files)
    
    
    disp(strcat('Sujeto: ',num2str(subj)))
    disp(strcat('Cue'))
    
    %--------- Load Data --------------------------------------------------
    addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
    EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    rmpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))

    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    cfg = [];
    cfg.latency = [-5 30];
    Data_Odor = ft_selectdata(cfg, dataOdor);
    
    cfg = [];
    cfg.latency = [25 60];
    Data_Vehicle = ft_selectdata(cfg, dataOdor);
    Data_Vehicle.time = Data_Odor.time;
    clear EEGOdor dataOdor
    
    %-----------Select first Seconds with an extension period---------------
    cfg = [];
    cfg.latency = [-5 10];
    Data_Odor = ft_selectdata(cfg, Data_Odor);
    Data_Vehicle = ft_selectdata(cfg, Data_Vehicle);    
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_Odor  = ft_freqanalysis(cfg_Tf, Data_Odor);
    Time_Freq_Vehicle  = ft_freqanalysis(cfg_Tf, Data_Vehicle);

    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-1 5];
    Time_Freq_Odor  = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehicle  = ft_selectdata(cfg, Time_Freq_Vehicle);

    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-1 5];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Odor  = ft_freqbaseline(cfg_Bas, Time_Freq_Odor);
    Time_Freq_Vehicle  = ft_freqbaseline(cfg_Bas, Time_Freq_Vehicle);

    save(strcat(savepath,files(subj).name(1:12),'TF_MN_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(subj).name(1:12),'TF_MN_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end

