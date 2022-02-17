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

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/DNight/';


%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
addpath('/home/andrea/Documents/Github/TimeFrequency_ReactivatedMemoryProject')
p_clustersOfInterest
clusters = fieldnames(Clust);


for subj = 1:numel(filesOdor)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
   
    if (any(isnan(Time_Freq_Odor.powspctrm)) | any(isnan(Time_Freq_Vehicle.powspctrm)))
        error(filesOdor(subj).name)
    end

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Odor);
    

    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);

    if (any(isnan(Time_Freq_Cue.powspctrm)) | any(isnan(Time_Freq_Vehicle.powspctrm)))
        error(filesOdor(subj).name)
    end

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                         = [];
        cfg.avgoverchan             = 'yes';
        cfg.channel                 = Clust.(clusters{cluster});


        Time_Freq_OdorD_MeanTemp    = ft_selectdata(cfg, Time_Freq_Cue);
        Time_Freq_OdorD_clust.(clusters{cluster})(subj,:,:) = Time_Freq_OdorD_MeanTemp.powspctrm;

        
        
        Time_Freq_VehicleD_MeanTemp    = ft_selectdata(cfg, Time_Freq_Vehicle);
        Time_Freq_VehicleD_clust.(clusters{cluster})(subj,:,:) = Time_Freq_VehicleD_MeanTemp.powspctrm;


    end

    cfg = [];
    cfg.avgoverchan = 'yes';
    cfg.channel     = {'all', '-E57', '-E100'};

    Time_Freq_OdorD_Mean = ft_selectdata(cfg, Time_Freq_Cue);
    Time_Freq_VehicleD_Mean= ft_selectdata(cfg, Time_Freq_Vehicle);

    Time_Freq_OdorD_All(subj,:,:) = Time_Freq_OdorD_Mean.powspctrm;
    Time_Freq_VehicleD_All(subj,:,:) = Time_Freq_VehicleD_Mean.powspctrm;
        
     
end

save(strcat(filepath,'Time_Freq_Clust_ToPlot_DNight'),'Time_Freq_OdorD_clust','Time_Freq_VehicleD_clust')
save(strcat(filepath,'Time_Freq_MeanAllchans_ToPlot_DNight'),'Time_Freq_OdorD_All','Time_Freq_VehicleD_All')

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_NewDatasets_ft/MNight/';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

for subj = 1:numel(filesOdor)

    disp(strcat('Sujeto: ',num2str(subj)))


    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));

    if (any(isnan(Time_Freq_Odor.powspctrm)) | any(isnan(Time_Freq_Vehicle.powspctrm)))
        error(filesOdor(subj).name)
    end
   

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_Odor);
    

    % Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Vehicle,  'trialinfo');
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_Vehicle);



    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)
        cfg                         = [];
        cfg.avgoverchan             = 'yes';
        cfg.channel                 = Clust.(clusters{cluster});


        Time_Freq_OdorM_MeanTemp    = ft_selectdata(cfg, Time_Freq_Cue);
        Time_Freq_OdorM_clust.(clusters{cluster})(subj,:,:) = Time_Freq_OdorM_MeanTemp.powspctrm;
        
        
        Time_Freq_VehicleM_MeanTemp    = ft_selectdata(cfg, Time_Freq_Vehicle);
        Time_Freq_VehicleM_clust.(clusters{cluster})(subj,:,:) = Time_Freq_VehicleM_MeanTemp.powspctrm;

    end

    cfg = [];
    cfg.avgoverchan = 'yes';

    Time_Freq_OdorM_Mean = ft_selectdata(cfg, Time_Freq_Cue);
    Time_Freq_VehicleM_Mean= ft_selectdata(cfg, Time_Freq_Vehicle);

    Time_Freq_OdorM_All(subj,:,:) = Time_Freq_OdorM_Mean.powspctrm;
    Time_Freq_VehicleM_All(subj,:,:) = Time_Freq_VehicleM_Mean.powspctrm;
     
end

save(strcat(filepath,'Time_Freq_Clust_ToPlot_MNight'),'Time_Freq_OdorM_clust','Time_Freq_VehicleM_clust')
save(strcat(filepath,'Time_Freq_MeanAllchans_ToPlot_MNight'),'Time_Freq_OdorM_All','Time_Freq_VehicleM_All')


save(strcat(filepath,'Time_Freq_Parms'),'Time_Freq_Cue')