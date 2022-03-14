% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/fieldtrip/') %stjude server
% addpath('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\fieldtrip\') %stjude computer



ft_defaults
% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/fieldtrip/qsub') %stjude server
% addpath('Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\Github\fieldtrip\qsub') %stjude computer


ft_warning off

s_tstep = 0.1; % try with 0.005
s_fstep = 0.15; % 0.005

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Organize files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------

% filepath_DNight = 'Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\TF_Calculation_NewDatasets_ft\DNight\'; %stjude computer
% filepath_MNight = 'Z:\ResearchHome\ClusterHome\asanch24\ReactivatedConnectivity\TF_Calculation_NewDatasets_ft\MNight\'; %stjude computer

filepath_DNight = '/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/TF_Calculation_NewDatasets_ft/DNight/'; %stjude server
filepath_MNight = '/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/TF_Calculation_NewDatasets_ft/MNight/'; %stjude server



%files = dir(strcat(filepath,'*.set'));
files_DNight = dir(strcat(filepath_DNight,'*.mat'));
files_MNight = dir(strcat(filepath_MNight,'*.mat'));
addpath('/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/Github/TimeFrequency_ReactivatedMemoryProject') %stjude server
p_clustersOfInterest
clusters = fieldnames(Clust);

Savepath = '/research/rgs01/home/clusterHome/asanch24/ReactivatedConnectivity/TF_Calculation_NewDatasets_ft/'; %stjude server

load(strcat(Savepath,'Time_Freq_Clust_ToPlot_DNight'))
load(strcat(Savepath,'Time_Freq_MeanAllchans_ToPlot_DNight'))

for subj = 13:numel(files_DNight)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath_DNight,files_DNight(subj).name));
    Time_Freq_DNight = TF_Baseline; clear TF_Baseline
    
    load(strcat(filepath_MNight,files_MNight(subj).name));
    Time_Freq_MNight = TF_Baseline; clear TF_Baseline
   
    if any(isnan(Time_Freq_DNight.powspctrm(:))) || any(isnan(Time_Freq_MNight.powspctrm(:)))
        error(files_DNight(subj).name)
    end
    
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_DNight,  'trialinfo');
    Time_Freq_DNight = ft_selectdata(cfg, Time_Freq_DNight);
    

    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_MNight,  'trialinfo');
    Time_Freq_MNight = ft_selectdata(cfg, Time_Freq_MNight);

    if any(isnan(Time_Freq_DNight.powspctrm(:))) || any(isnan(Time_Freq_MNight.powspctrm(:)))
        error(files_DNight(subj).name)
    end

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                         = [];
        cfg.avgoverchan             = 'yes';
        cfg.channel                 = Clust.(clusters{cluster});


        Time_Freq_DNight_MeanTemp    = ft_selectdata(cfg, Time_Freq_DNight);
        Time_Freq_DNight_clust.(clusters{cluster})(subj,:,:) = Time_Freq_DNight_MeanTemp.powspctrm;

        
        
        Time_Freq_MNight_MeanTemp    = ft_selectdata(cfg, Time_Freq_MNight);
        Time_Freq_MNight_clust.(clusters{cluster})(subj,:,:) = Time_Freq_MNight_MeanTemp.powspctrm;


    end

    cfg = [];
    cfg.avgoverchan = 'yes';
    cfg.channel     = {'all', '-E57', '-E100'};

    Time_Freq_DNight_Mean = ft_selectdata(cfg, Time_Freq_DNight);
    Time_Freq_MNight_Mean= ft_selectdata(cfg, Time_Freq_MNight);

    Time_Freq_DNight_All(subj,:,:) = Time_Freq_DNight_Mean.powspctrm;
    Time_Freq_MNight_All(subj,:,:) = Time_Freq_MNight_Mean.powspctrm;
    
    save(strcat(Savepath,'Time_Freq_Clust_ToPlot_DNight'),'Time_Freq_DNight_clust','Time_Freq_MNight_clust','-v7.3')
    save(strcat(Savepath,'Time_Freq_MeanAllchans_ToPlot_DNight'),'Time_Freq_DNight_All','Time_Freq_MNight_All','-v7.3')


    save(strcat(Savepath,'Time_Freq_Parms'),'Time_Freq_DNight')
        
end

