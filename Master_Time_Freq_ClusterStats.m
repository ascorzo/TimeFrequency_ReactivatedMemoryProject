addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

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

%filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Filtered-David/CueNight/';

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    

    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                 = [];
        cfg.avgoverchan     = 'yes';
        cfg.channel         = Clust.(clusters{cluster});
        Time_Freq_OdorD_clust.(clusters{cluster}){subj} = ...
        ft_selectdata(cfg,Time_Freq_Cue);


        cfg                 = [];
        cfg.avgoverchan     = 'yes';
        cfg.channel         = Clust.(clusters{cluster});
        Time_Freq_VehicleD_clust.(clusters{cluster}){subj} = ...
        ft_selectdata(cfg,Time_Freq_Vehicle);
    end
     
end
savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/';
save(strcat(savepath,'FT_Time_Freq_Clust_ToPlot_DNight_7sec'),'Time_Freq_OdorD_clust','Time_Freq_VehicleD_clust')

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/MNight/';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor)

    disp(strcat('Sujeto: ',num2str(subj)))


    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    

    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                 = [];
        cfg.avgoverchan     = 'yes';
        cfg.channel         = Clust.(clusters{cluster});
        Time_Freq_OdorM_clust.(clusters{cluster}){subj} = ...
        ft_selectdata(cfg,Time_Freq_Cue);


        cfg                 = [];
        cfg.avgoverchan     = 'yes';
        cfg.channel         = Clust.(clusters{cluster});
        Time_Freq_VehicleM_clust.(clusters{cluster}){subj} = ...
        ft_selectdata(cfg,Time_Freq_Vehicle);
    end
     
end

savepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_7SecTrial/';
save(strcat(savepath,'FT_Time_Freq_Clust_ToPlot_MNight_7sec'),'Time_Freq_OdorM_clust','Time_Freq_VehicleM_clust')


%%

for cluster = 1:numel(clusters)
    
    cfg = [];
    cfg.method            = 'montecarlo';
    cfg.statistic         = 'indepsamplesT';
    cfg.correctm          = 'cluster';
    cfg.tail              = 1;
    cfg.alpha             = 0.05;
    cfg.clustertail       = 1;
    cfg.clusteralpha      = 0.05;  
    cfg.clusterstatistic  = 'maxsum';
    cfg.correcttail       = 'prob'; 
    cfg.numrandomization  = 1000; 
    cfg.neighbours        = [];
    cfg.design            = [ones(1,numel(Time_Freq_OdorD_clust.(clusters{cluster}))),...
        ones(1,numel(Time_Freq_OdorM_clust.(clusters{cluster})))+1];
    % design matrix, note the transpose
    cfg.ivar              = 1;
    
    statpar_fdr =           ft_freqstatistics(cfg, Time_Freq_OdorD_clust.(clusters{cluster}){:},...
        Time_Freq_OdorM_clust.(clusters{cluster}){:});
    
    
    cfg = [];
    [grandavg_odor] = ft_freqgrandaverage(cfg,  Time_Freq_OdorD_clust.(clusters{cluster}){:});
    grandavg_odor.mask = statpar_fdr.mask;
    
    sum(statpar_fdr.mask(:))

    figure
    cfg = [];
    cfg.maskparameter = 'mask';
    ft_singleplotTFR(cfg, grandavg_odor);
    
end