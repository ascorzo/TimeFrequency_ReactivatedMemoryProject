% p_TFOverTrials

addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

p_clustersOfInterest
clusters = fieldnames(Clust);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % D Night
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));


DNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_DNight');

for subj = 1:numel(filesOdor) 
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
    v_cycles =  DNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
    [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);

    s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
    s_Lastepoch = s_MaxNumStimsIdx;
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));


    cfg                         = [];
    cfg.trials                  = s_Firstepoch;
    TF_OdorFirstTrial{subj}     = ft_selectdata(cfg, Time_Freq_Odor);
    TF_VehicleFirstTrial{subj}  = ft_selectdata(cfg, Time_Freq_Vehicle);
    

    cfg                         = [];
    cfg.trials                  = s_Lastepoch;
    TF_OdorLastTrial{subj}      = ft_selectdata(cfg, Time_Freq_Odor);
    TF_VehicleLastTrial{subj}   = ft_selectdata(cfg, Time_Freq_Vehicle);

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                         = [];
        cfg.channel                 = Clust.(clusters{cluster});
        cfg.avgoverchan             = 'yes';

        TF_OdorFirstTrial_chan           = ft_selectdata(cfg, TF_OdorFirstTrial);
        TF_VehicleFirstTrial_chan        = ft_selectdata(cfg, TF_VehicleFirstTrial);
        TF_OdorLastTrial_chan            = ft_selectdata(cfg, TF_OdorLastTrial);
        TF_VehicleLastTrial_chan         = ft_selectdata(cfg, TF_VehicleLastTrial);  

        TF_OdorFirstTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_OdorFirstTrial_chan.powspctrm);
        TF_VehicleFirstTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_VehicleFirstTrial_chan.powspctrm);
        TF_OdorLastTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_OdorLastTrial_chan.powspctrm);
        TF_VehicleLastTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_VehicleLastTrial_chan.powspctrm);

    end


    cfg             = [];
    cfg.channel     = 'all';
    cfg.avgoverchan = 'yes';

    TF_OdorFirstTrial_chan           = ft_selectdata(cfg, TF_OdorFirstTrial);
    TF_VehicleFirstTrial_chan        = ft_selectdata(cfg, TF_VehicleFirstTrial);
    TF_OdorLastTrial_chan            = ft_selectdata(cfg, TF_OdorLastTrial);
    TF_VehicleLastTrial_chan         = ft_selectdata(cfg, TF_VehicleLastTrial); 

    TF_OdorFirstTrial_All.All(subj,:,:) = ...
        squeeze(TF_OdorFirstTrial_chan.powspctrm);
    TF_VehicleFirstTrial_All.All(subj,:,:) = ...
        squeeze(TF_VehicleFirstTrial_chan.powspctrm);
    TF_OdorLastTrial_All.All(subj,:,:) = ...
        squeeze(TF_OdorLastTrial_chan.powspctrm);
    TF_VehicleLastTrial_All.All(subj,:,:) = ...
        squeeze(TF_VehicleLastTrial_chan.powspctrm);


end

savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

save(strcat(savepath,'TF_DNight_FirstLastTrial_clusters'),...
    'TF_OdorFirstTrial_All','TF_VehicleFirstTrial_All','TF_OdorLastTrial_All','TF_VehicleLastTrial_All')

save(strcat(savepath,'TF_DNight_FirstLastTrial'),...
    'TF_OdorFirstTrial','TF_VehicleFirstTrial','TF_OdorLastTrial','TF_VehicleLastTrial')


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % M Night
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));


MNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_MNight');

for subj = 1:numel(filesOdor) 
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
    v_cycles =  MNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
    [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);

    s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
    s_Lastepoch = s_MaxNumStimsIdx;
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));


    cfg                         = [];
    cfg.trials                  = s_Firstepoch;
    TF_OdorFirstTrial{subj}     = ft_selectdata(cfg, Time_Freq_Odor);
    TF_VehicleFirstTrial{subj}  = ft_selectdata(cfg, Time_Freq_Vehicle);
    

    cfg                         = [];
    cfg.trials                  = s_Lastepoch;
    TF_OdorLastTrial{subj}      = ft_selectdata(cfg, Time_Freq_Odor);
    TF_VehicleLastTrial{subj}   = ft_selectdata(cfg, Time_Freq_Vehicle);

    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    
    for cluster = 1:numel(clusters)

        cfg                         = [];
        cfg.channel                 = Clust.(clusters{cluster});
        cfg.avgoverchan             = 'yes';

        TF_OdorFirstTrial_chan           = ft_selectdata(cfg, TF_OdorFirstTrial);
        TF_VehicleFirstTrial_chan        = ft_selectdata(cfg, TF_VehicleFirstTrial);
        TF_OdorLastTrial_chan            = ft_selectdata(cfg, TF_OdorLastTrial);
        TF_VehicleLastTrial_chan         = ft_selectdata(cfg, TF_VehicleLastTrial);  

        TF_OdorFirstTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_OdorFirstTrial_chan.powspctrm);
        TF_VehicleFirstTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_VehicleFirstTrial_chan.powspctrm);
        TF_OdorLastTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_OdorLastTrial_chan.powspctrm);
        TF_VehicleLastTrial_All.(clusters{cluster})(subj,:,:) = ...
            squeeze(TF_VehicleLastTrial_chan.powspctrm);

    end


    cfg             = [];
    cfg.channel     = 'all';
    cfg.avgoverchan = 'yes';

    TF_OdorFirstTrial_chan           = ft_selectdata(cfg, TF_OdorFirstTrial);
    TF_VehicleFirstTrial_chan        = ft_selectdata(cfg, TF_VehicleFirstTrial);
    TF_OdorLastTrial_chan            = ft_selectdata(cfg, TF_OdorLastTrial);
    TF_VehicleLastTrial_chan         = ft_selectdata(cfg, TF_VehicleLastTrial); 

    TF_OdorFirstTrial_All.All(subj,:,:) = ...
        squeeze(TF_OdorFirstTrial_chan.powspctrm);
    TF_VehicleFirstTrial_All.All(subj,:,:) = ...
        squeeze(TF_VehicleFirstTrial_chan.powspctrm);
    TF_OdorLastTrial_All.All(subj,:,:) = ...
        squeeze(TF_OdorLastTrial_chan.powspctrm);
    TF_VehicleLastTrial_All.All(subj,:,:) = ...
        squeeze(TF_VehicleLastTrial_chan.powspctrm);

end

savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

save(strcat(savepath,'TF_MNight_FirstLastTrial_clusters'),...
    'TF_OdorFirstTrial_All','TF_VehicleFirstTrial_All','TF_OdorLastTrial_All','TF_VehicleLastTrial_All')

save(strcat(savepath,'TF_MNight_FirstLastTrial'),...
    'TF_OdorFirstTrial','TF_VehicleFirstTrial','TF_OdorLastTrial','TF_VehicleLastTrial')