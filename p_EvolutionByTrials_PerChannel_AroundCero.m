addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Odor Night
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

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

    Time_Freq_Odor_aroundCero{subj}   = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_aroundCero{subj}   = ft_selectdata(cfg, Time_Freq_Vehicle);
end

%%
mintrials = 500;

for subj = 1:numel(filesOdor)

    subj_trials = size(Time_Freq_Odor_aroundCero{subj}.powspctrm,1);
    mintrials = min(mintrials,subj_trials);

end

Time_Freq_OdorD_evol = [];
Time_Freq_VehicleD_evol = [];

for subj = 1:numel(filesOdor)

    for freq = 1:length(frequency_bands)
        
        cfg                 = [];
        cfg.frequency       = frequency_bands(freq,:);
        cfg.avgoverfreq     = 'yes';
        
        Time_Freq_OdorD_freq{subj} = ...
            ft_selectdata(cfg,Time_Freq_Odor_aroundCero{subj});
        
        Time_Freq_VehicleD_freq{subj} = ...
            ft_selectdata(cfg,Time_Freq_Vehi_aroundCero{subj});
        
        Time_Freq_OdorD_evol(:,:,freq,subj) = ...
            Time_Freq_OdorD_freq{subj}.powspctrm(1:mintrials,:);
        
        Time_Freq_VehicleD_evol(:,:,freq,subj) = ...
            Time_Freq_VehicleD_freq{subj}.powspctrm(1:mintrials,:);
        
    end
    
end

savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

save(strcat(savepath,'firstSec_TFcount_DNight_chans'),...
    'Time_Freq_OdorD_evol','Time_Freq_VehicleD_evol','frequency_bands')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

    Time_Freq_Odor_aroundCero{subj}   = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_aroundCero{subj}   = ft_selectdata(cfg, Time_Freq_Vehicle);

end

%%
mintrials = 500;

for subj = 1:numel(filesOdor)

    subj_trials = size(Time_Freq_Odor_aroundCero{subj}.powspctrm,1);
    mintrials = min(mintrials,subj_trials);

end

Time_Freq_OdorM_evol = [];
Time_Freq_VehicleM_evol = [];

for subj = 1:numel(filesOdor)
    
    for freq = 1:length(frequency_bands)
        
        cfg                 = [];
        cfg.frequency       = frequency_bands(freq,:);
        cfg.avgoverfreq     = 'yes';
        
        Time_Freq_OdorM_freq{subj} = ...
            ft_selectdata(cfg,Time_Freq_Odor_aroundCero{subj});
        
        Time_Freq_VehicleM_freq{subj} = ...
            ft_selectdata(cfg,Time_Freq_Vehi_aroundCero{subj});
        
        Time_Freq_OdorM_evol(:,:,freq,subj)  = ...
            Time_Freq_OdorM_freq{subj}.powspctrm(1:mintrials,:);
        
        Time_Freq_VehicleM_evol(:,:,freq,subj)  = ...
            Time_Freq_VehicleM_freq{subj}.powspctrm(1:mintrials,:);
        
    end
    
end

savepath ='/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/';

save(strcat(savepath,'firstSec_TFcount_MNight_chans'),...
    'Time_Freq_OdorM_evol','Time_Freq_VehicleM_evol','frequency_bands')
