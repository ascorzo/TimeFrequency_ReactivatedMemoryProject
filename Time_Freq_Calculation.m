addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')
ft_defaults
addpath(genpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/'))
eeglab nogui

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Cue Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/home/sleep/Documents/DAVID/Datasets/Ori/TrialGroups_DetrendOptim_Off_On/';


filesOdor = dir(strcat(filepath,'*Odor*.set'));
filesSham = dir(strcat(filepath,'*Sham*.set'));

for subj = 1%:numel(filesOdor)
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
    %________________________________________________________
    %Cue Odor
    %________________________________________________________
    
    EEGOdor = pop_loadset(strcat(filepath,filesOdor(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    dataOdor.sampleinfo = zeros(numel(EEGOdor.event),2);
    dataOdor.sampleinfo(:,1) = [EEGOdor.event.latency]-15*dataOdor.fsample;
    dataOdor.sampleinfo(:,2) = [EEGOdor.event.latency]+15*dataOdor.fsample;
    
    cfg                     = [];
    cfg.output              = 'pow';
    cfg.method              = 'mtmconvol';
    cfg.taper               = 'hanning';
    cfg.foi                 = 0.5:0.5:20;
    cfg.t_ftimwin           = ones(length(cfg.foi),1).*0.5;
    cfg.toi                 = -15:0.5:15;
    %cfg.baseline            = [-15 0];
    %cfg.baselinetype        = 'zscore';
    %cfg.channel             = 'E36';
    cfg.keeptrials          = 'yes';
    %Time_Freq_Odor{subj}    = ft_freqanalysis(cfg, dataOdor);
    Time_Freq_Odor    = ft_freqanalysis(cfg, dataOdor);
    
    %------------------------------------------------------
    %Baseline Correction
    %------------------------------------------------------
    begsample = dataOdor.sampleinfo(:,1);
    endsample = dataOdor.sampleinfo(:,2);
    time      = ((begsample+endsample)/2) / dataOdor.fsample;

    freq_continuous           = Time_Freq_Odor;
    freq_continuous.powspctrm = permute(Time_Freq_Odor.powspctrm, [1, 3, 4, 2]);
    freq_continuous.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
    freq_continuous.time      = time;             % add the description of the time dimension
    
    figure
    cfg                = [];
    cfg.baseline       = [-15 0];
    cfg.baselinetype   = 'normchange';
    cfg.zlim           = [-0.2 0.2];
    ft_singleplotTFR(cfg, Time_Freq_Odor);

    %________________________________________________________
    %Sham Odor
    %________________________________________________________
    
%     EEGSham = pop_loadset(strcat(filepath,filesSham(subj).name));
%     dataSham = eeglab2fieldtrip(EEGSham,'raw');
%     
%     cfg                     = [];
%     cfg.output              = 'pow';
%     cfg.method              = 'mtmconvol';
%     cfg.taper               = 'hanning';
%     cfg.foi                 = 0.5:0.5:20;                        
%     cfg.t_ftimwin           = ones(length(cfg.foi),1).*0.5;
%     cfg.toi                 = -15:0.5:15; 
%     %cfg.baseline            = [-15 0];
%     %cfg.baselinetype        = 'zscore';
%     cfg.keeptrials          = 'yes';
%     Time_Freq_Sham_CN{subj}    = ft_freqanalysis(cfg, dataSham);
    
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Placebo Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = '/home/sleep/Documents/DAVID/Datasets/Ori_PlaceboNight/TrialGroups_DetrendOptim_Off_On/';


filesOdor = dir(strcat(filepath,'*Odor*.set'));
filesSham = dir(strcat(filepath,'*Sham*.set'));

for subj = 1:numel(filesOdor)
    
    
    
    %---Placebo Odor---%
    EEGOdor = pop_loadset(strcat(filepath,filesOdor(subj).name));
    dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
    
    cfg                     = [];
    cfg.output              = 'pow';
    cfg.method              = 'mtmconvol';
    cfg.taper               = 'hanning';
    cfg.foi                 = 0.5:0.5:20;
    cfg.t_ftimwin           = ones(length(cfg.foi),1).*0.5;
    cfg.toi                 = -15:0.5:15;
    %cfg.baseline            = [-15 0];
    %cfg.baselinetype        = 'zscore';
    cfg.keeptrials          = 'yes';
    Time_Freq_Placebo{subj}    = ft_freqanalysis(cfg, dataOdor);
    

     %---Sham Odor---%
    EEGSham = pop_loadset(strcat(filepath,filesSham(subj).name));
    dataSham = eeglab2fieldtrip(EEGSham,'raw');
    
    cfg                     = [];
    cfg.output              = 'pow';
    cfg.method              = 'mtmconvol';
    cfg.taper               = 'hanning';
    cfg.foi                 = 0.5:0.5:20;                        
    cfg.t_ftimwin           = ones(length(cfg.foi),1).*0.5;
    cfg.toi                 = -15:0.5:15; 
    %cfg.baseline            = [-15 0];
    %cfg.baselinetype        = 'zscore';
    cfg.keeptrials          = 'yes';
    Time_Freq_Sham_PN{subj}    = ft_freqanalysis(cfg, dataSham);
    
end










