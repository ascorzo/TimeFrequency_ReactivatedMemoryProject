% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/fieldtrip/') %stjude server
% addpath('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\fieldtrip\') %stjude computer

% addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')

ft_defaults
% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/fieldtrip/qsub') %stjude server
% addpath('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\fieldtrip\qsub') %stjude computer

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


% filepath = '/mnt/disk1/sleep/German_Study/Data/FT_Preprocessing_250/';
% savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/';

%stjude server
filepath = '/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/FT_Preprocessing_250/';
savepath = '/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/TF_Calculation_NewDatasets_ft/';

%stjude computer
% filepath = 'Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\FT_Preprocessing_250\';
% savepath = 'Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\TF_Calculation_NewDatasets_ft\';

% filepath = 'D:\Thesis_Publication\Epoched_90SecTrial_MastoidRef-Interp_NewEpochs\DNight\';
% savepath = 'D:\TF_Calculation_90SecTrial\';

files = dir(strcat(filepath,'*.mat'));
% load('/home/andrea/Documents/Github/rc_preproc/EventsDescription.mat')
load('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/rc_preproc/EventsDescription.mat')%stjude server
% load('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\rc_preproc\EventsDescription.mat')%stjude computer


%-------------------Run Pairing--------------------------------------------------
% addpath('/home/andrea/Documents/Github/rc_preproc/')
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/rc_preproc/')%stjude server
% addpath('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\rc_preproc\')%stjude computer

p_TrialPairing


for file = 3:numel(files)
    
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
    % pairing
    %----------------------------------------------------------------------

    recording       = find(strcmp(PairingAll(:,1), files(file).name(1:6)));
    pairs           = PairingAll{recording,2};
    idxOdors        = pairs(:,1);
    idxVehicles     = pairs(:,2);
    
    % check pairing
    
    if data_downsamp_250.trialinfo(idxOdors)~= 1
        error('pairing indexes are wrong')
    end
    
    if data_downsamp_250.trialinfo(idxVehicles)~= 0
        error('pairing indexes are wrong')
    end

    %----------------------------------------------------------------------
    % Separate conditions
    %----------------------------------------------------------------------
    
    cfg             = [];
    cfg.trials      = idxOdors;
    Time_Freq_Odor  = ft_selectdata(cfg, Time_Freq);
    
    cfg                 = [];
    cfg.trials          = idxVehicles;
    Time_Freq_Vehicle   = ft_selectdata(cfg, Time_Freq);

    %----------------------------------------------------------------------
    % baseline correction
    %----------------------------------------------------------------------
    %%
    for trial = 1:numel(idxOdors)
        TF_Baseline     = Time_Freq_Odor;
        
        Vehicle_TF      = squeeze(Time_Freq_Vehicle.powspctrm(trial,:,:,:));
        meanVals        = squeeze(mean(Vehicle_TF,3));
        stdVals         = squeeze(std(Vehicle_TF,1,3));
        
        meanVals        = repmat(meanVals,[1 1 size(Vehicle_TF, 3)]);
        stdVals         = repmat(stdVals,[1 1 size(Vehicle_TF, 3)]);
        
        Odor_TF         = squeeze(Time_Freq_Odor.powspctrm(trial,:,:,:));
        TrialTF_Odor    =(Odor_TF-meanVals)./stdVals;
        
        TF_Baseline.powspctrm(trial,:,:,:) = TrialTF_Odor;
    end

    save(strcat(savepath,files(file).name(1:6),'_TF_NewBaseline'),'TF_Baseline','-v7.3')
end