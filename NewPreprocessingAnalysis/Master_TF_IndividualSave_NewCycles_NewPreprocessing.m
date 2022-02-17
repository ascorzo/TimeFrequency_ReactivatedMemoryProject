addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')
% addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
ft_warning off

%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.1; % try with 0.005
s_fstep = 0.1; % 0.005
minfreq = 0.5;
maxfreq = 20;

cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  = minfreq:s_fstep:maxfreq; 
% cfg_Tf.width                = cycles;
v_timeWindows               = 5:(-4.5/numel(cfg_Tf.foi)):0.5;
cfg_Tf.t_ftimwin            = v_timeWindows;   
cfg_Tf.toi                  = -10:s_tstep:20; 
% toi this is extended before and after to deal with border effect of wavelet
cfg_Tf.keeptrials           = 'yes';


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


filepath = '/mnt/disk1/sleep/German_Study/Data/FT_Preprocessing_250/';
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/';

% filepath = 'D:\Thesis_Publication\Epoched_90SecTrial_MastoidRef-Interp_NewEpochs\DNight\';
% savepath = 'D:\TF_Calculation_90SecTrial\';

files = dir(strcat(filepath,'*.mat'));
load('/home/andrea/Documents/Github/rc_preproc/EventsDescription.mat')


for file = 19%:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('File: ',files(file).name))
    
    %--------- Load Data --------------------------------------------------
    
    load(strcat(filepath,files(file).name))
    
    %-----------Time-Frequency Calculation---------------------------------
    
    Time_Freq_Temp  = ft_freqanalysis(cfg_Tf, data_downsamp_250);
    
    %--------------Select time of interest removing borders-----------------
    cfg = [];
    cfg.latency = [-5 15];
    Time_Freq = ft_selectdata(cfg, Time_Freq_Temp);
    
    
    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    cfg_Bas                     = [];
    cfg_Bas.baseline            = [-5 15];
    cfg_Bas.baselinetype        = 'zscore';
    Time_Freq_Baseline          = ft_freqbaseline(cfg_Bas,Time_Freq);
    
    clear Time_Freq
    
    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    % EventsFile = strcmp(AllEvents(1,:),files(file).name(1:6));

    % ThisFileEvents = AllEvents{2,EventsFile};
    
    % ValidTrials  = ([ThisFileEvents.Rejected] == 0);
    % ValidTrials  = ThisFileEvents(ValidTrials);

    % OdorTrials      = strcmp({ValidTrials.stimulation},'ODOR');
    % VehicleTrials   = strcmp({ValidTrials.stimulation},'VEHICLE');

    OdorTrials      = find(data_downsamp_250.trialinfo == 1);
    VehicleTrials   = find(data_downsamp_250.trialinfo == 0);
    
    cfg = [];
    cfg.trials = OdorTrials;
    Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
    
    cfg = [];
    cfg.trials = VehicleTrials;
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
    Time_Freq_Vehicle.time = Time_Freq_Odor.time;

    save(strcat(savepath,files(file).name(1:6),'TF_NewCycles_Odor'),'Time_Freq_Odor','-v7.3')
    save(strcat(savepath,files(file).name(1:6),'TF_NewCycles_Vehicle'),'Time_Freq_Vehicle','-v7.3')

end