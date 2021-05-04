%% Analyze Time Frequency Calculation
% This notebook takes the Time-frequency calculation for the Reactivated connectivity 
% project, loading each file individually and calculating time series and statistics
%% Importing Fieldtrip libraries
% In this section we add the the Fieldtrip folder in the path and run ft_defaults 
% to have access to all functions of fieldtrip

addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/')

ft_defaults
addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/qsub')

ft_warning off

%% Initial parameters
%addpath('/mnt/disk1/andrea/German_Study/');
p_clustersOfInterest


cfg_Sel = [];
cfg_Sel.channel = Clust.occipital;

%--------------------------------------------------------------------------
% Start MVPA-Light
%--------------------------------------------------------------------------

addpath('/gpfs01/born/group/Andrea/Github/MVPA-Light/startup/')
startup_MVPA_Light

%--------------------------------------------------------------------------
% Parameters for Baseline Correction
%--------------------------------------------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-15 45];
cfg_Bas.baselinetype        = 'zscore';

%--------------------------------------------------------------------------
% Files to work with
%--------------------------------------------------------------------------

% For Odor D Night
filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/';

% For Odor M Night
%filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_OdorD%_Night/';

filesOdor = dir(strcat(filepath,'*Odor.mat'));
filesVehicle = dir(strcat(filepath,'*Vehicle.mat'));
%% MVPA TF analysis
% 
subjects = 1:numel(filesOdor);

% load('Time-Freq_DNight_MastoidRef_Interp_12cycles_90secTrial.mat');

% Time_Freq = Time_Freq_DA;

%clear Time_Freq_DA

for subj = subjects
    % display current subject name
    disp(strcat('Sujeto: ',num2str(subj)))
    
    % Load Data
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));

%     %----------------------------------------------------------------------
%     % baseline correction
%     %----------------------------------------------------------------------
%     cfg_Bas                     = [];
%     cfg_Bas.baseline            = [-15 45];
%     cfg_Bas.baselinetype        = 'zscore';
%     Time_Freq_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq);
    
%     %----------------------------------------------------------------------
%     % Separate conditions
%     %----------------------------------------------------------------------
%     cfg = [];
%     cfg.latency = [-15 15];
%     Time_Freq_Odor = ft_selectdata(cfg, Time_Freq_Baseline);
%     
%     cfg = [];
%     cfg.latency = [15 45];
%     Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Baseline);
%     Time_Freq_Vehicle.time = Time_Freq_Odor.time;
    


    
    clabel = [ones(1,size(Time_Freq_Odor.powspctrm,1)),...
    ones(1,size(Time_Freq_Vehicle.powspctrm,1))+1];


%     -- For Clust.left_frontal
    cfg = []; cfg.channel = Clust.left_frontal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.left_frontal{subj,1} = result_freq;


    %-- For Clust.right_frontal 
    cfg = []; cfg.channel = Clust.right_frontal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.right_frontal{subj,1} = result_freq;


    %-- For Clust.frontal 
    cfg = []; cfg.channel = Clust.frontal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.frontal{subj,1} = result_freq;


    %-- For Clust.left_central
    cfg = []; cfg.channel = Clust.left_central; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.left_central{subj,1} = result_freq;


    %-- For Clust.right_central
    cfg = []; cfg.channel = Clust.right_central; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.right_central{subj,1} = result_freq;


    %-- For Clust.central
    cfg = []; cfg.channel = Clust.central; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.central{subj,1} = result_freq;


    %-- For Clust.left_temporal
    cfg = []; cfg.channel = Clust.left_temporal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.left_temporal{subj,1} = result_freq;


    %-- For Clust.right_temporal
    cfg = []; cfg.channel = Clust.right_temporal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.right_temporal{subj,1} = result_freq;


    %-- For Clust.left_parietal 
    cfg = []; cfg.channel = Clust.left_parietal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.left_parietal{subj,1} = result_freq;


    %-- For Clust.right_parietal
    cfg = []; cfg.channel = Clust.right_parietal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.right_parietal{subj,1} = result_freq;


    %-- For Clust.parietal
    cfg = []; cfg.channel = Clust.parietal; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.parietal{subj,1} = result_freq;


    %-- For Clust.left_occipital
    cfg = []; cfg.channel = Clust.left_occipital; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.left_occipital{subj,1} = result_freq;


    %-- For Clust.right_occipital 
    cfg = []; cfg.channel = Clust.right_occipital; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.right_occipital{subj,1} = result_freq;


    %-- For Clust.occipital
    cfg = []; cfg.channel = Clust.occipital; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.occipital{subj,1} = result_freq;

    %-- For all chans
    cfg = []; cfg.channel = 'all'; cfg.avgoverchan = 'yes';
    Time_Freq_Odor_Temp = ft_selectdata(cfg,Time_Freq_Odor);
    Time_Freq_Vehicle_Temp = ft_selectdata(cfg,Time_Freq_Vehicle);

    p_MVPA_TFclassifier
    resultsfreq.Allchans{subj,1} = result_freq;
 
    %Time_Freq{subj} = [];
    save(strcat('/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/','resultsfreq_ONOFF.mat'),'resultsfreq','-v7.3')
    
end
<<<<<<< Updated upstream

=======
>>>>>>> Stashed changes
%%
%-- For Clust.left_frontal
result_freq = resultsfreq.left_frontal;
p_MVPA_TFClusterStats

stat2.left_frontal = stat_level2;
resultAverage.left_frontal = result_average;

%-- For Clust.right_frontal 
result_freq = resultsfreq.right_frontal;
p_MVPA_TFClusterStats

stat2.right_frontal = stat_level2;
resultAverage.right_frontal = result_average;

%-- For Clust.frontal 
result_freq = resultsfreq.frontal;
p_MVPA_TFClusterStats

stat2.frontal = stat_level2;
resultAverage.frontal = result_average;


%-- For Clust.left_central
result_freq = resultsfreq.left_central;
p_MVPA_TFClusterStats

stat2.left_central = stat_level2;
resultAverage.left_central = result_average;

%-- For Clust.right_central
result_freq = resultsfreq.right_central;
p_MVPA_TFClusterStats

stat2.right_central = stat_level2;
resultAverage.right_central = result_average;

%-- For Clust.central
result_freq = resultsfreq.central;
p_MVPA_TFClusterStats

stat2.central = stat_level2;
resultAverage.central = result_average;

%-- For Clust.left_temporal
result_freq = resultsfreq.left_temporal;
p_MVPA_TFClusterStats

stat2.left_temporal = stat_level2;
resultAverage.left_temporal = result_average;

%-- For Clust.right_temporal
result_freq = resultsfreq.right_temporal;
p_MVPA_TFClusterStats

stat2.right_temporal = stat_level2;
resultAverage.right_temporal = result_average;

%-- For Clust.left_parietal 
result_freq = resultsfreq.left_parietal;
p_MVPA_TFClusterStats

stat2.left_parietal = stat_level2;
resultAverage.left_parietal = result_average;

%-- For Clust.right_parietal
result_freq = resultsfreq.right_parietal;
p_MVPA_TFClusterStats

stat2.right_parietal = stat_level2;
resultAverage.right_parietal = result_average;

%-- For Clust.parietal
result_freq = resultsfreq.parietal;
p_MVPA_TFClusterStats

stat2.parietal = stat_level2;
resultAverage.parietal = result_average;

%-- For Clust.left_occipital
result_freq = resultsfreq.left_occipital;
p_MVPA_TFClusterStats

stat2.left_occipital = stat_level2;
resultAverage.left_occipital = result_average;

%-- For Clust.right_occipital 
result_freq = resultsfreq.right_occipital;
p_MVPA_TFClusterStats

stat2.right_occipital = stat_level2;
resultAverage.right_occipital = result_average;

%-- For Clust.occipital
result_freq = resultsfreq.occipital;
p_MVPA_TFClusterStats

stat2.occipital = stat_level2;
resultAverage.occipital = result_average;

%-- For AllChans
result_freq = resultsfreq.Allchans;
p_MVPA_TFClusterStats

stat2.Allchans = stat_level2;
resultAverage.Allchans = result_average;

save(strcat('/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/','stat2_TF_ONOFF'),'stat2','-v7.3')
save(strcat('/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/','resultAverage_ONOFF'),'resultAverage','-v7.3')

%%
clusters = fieldnames(resultAverage);

for clust = 1:numel(clusters)
    
    cluster = clusters{clust};
    mv_plot_result(resultAverage.(cluster),(-5:0.1:30),(0.5:0.1:20),'mask', stat2.(cluster).mask)
    
    l = title(strcat(cluster,'  Average'));
    set(l,'interpreter','none')
    xlabel('time(s)')
    ylabel('Frequency(Hz)')
    
    
%     saveas(gcf,strcat(cluster,'_Mask_OdorD_Night.png'))
    
    mv_plot_result(resultAverage.(cluster),(-5:0.1:30),(0.5:0.1:20))
    
    l = title(strcat(cluster,'  Average'));
    set(l,'interpreter','none')
    xlabel('time(s)')
    ylabel('Frequency(Hz)')
    
%     saveas(gcf,strcat(cluster,'_NoMask_OdorD_Night.png'))
end
% close all

