addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off


filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

frequency_bands = [0.5 2; 1 4; 4 8; 8 13; 13 16];

for subj = 1:numel(filesOdor)

    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));

    cfg                         = [];
    cfg.latency                 = [-1 1];
    cfg.avgovertime             = 'yes';

    Time_Freq_Odor_aroundCero   = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_aroundCero   = ft_selectdata(cfg, Time_Freq_Vehicle);

    for cluster = 1:numel(clusters) 
        cfg                 = [];
        cfg.avgoverchan     = 'yes';
        cfg.channel         = Clust.(clusters{cluster});

        Time_Freq_OdorM_clust.(clusters{cluster}){subj} = ...
            ft_selectdata(cfg,Time_Freq_Odor_aroundCero);
        
        Time_Freq_VehicleM_clust.(clusters{cluster}){subj} = ...
            ft_selectdata(cfg,Time_Freq_Vehi_aroundCero);

    end

end

%%
mintrials = 500;

for subj = 1:numel(filesOdor)

    subj_trials = size(Time_Freq_OdorM_clust.central{subj}.powspctrm,1);
    mintrials = min(mintrials,subj_trials);

end


for subj = 1:numel(filesOdor)
    
    for cluster = 1:numel(clusters)
        for freq = 1:length(frequency_bands)
            
            cfg                 = [];
            cfg.frequency       = frequency_bands(freq,:);
            cfg.avgoverfreq     = 'yes';
            
            Time_Freq_OdorM_freq.(clusters{cluster}){subj} = ...
                ft_selectdata(cfg,Time_Freq_OdorM_clust.(clusters{cluster}){subj});
            
            Time_Freq_VehicleM_freq.(clusters{cluster}){subj}= ...
                ft_selectdata(cfg,Time_Freq_VehicleM_clust.(clusters{cluster}){subj});
            
            Time_Freq_OdorM_evol.(clusters{cluster})(freq,:,subj) = ...
                Time_Freq_OdorM_freq.(clusters{cluster}){subj}.powspctrm(1:mintrials);
            
            Time_Freq_VehicleM_evol.(clusters{cluster})(freq,:,subj) = ...
                Time_Freq_VehicleM_freq.(clusters{cluster}){subj}.powspctrm(1:mintrials);
            
        end
    end
    
end
