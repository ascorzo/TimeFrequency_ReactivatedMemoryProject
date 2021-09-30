addpath('C:\Users\lanan\Documents\Github\TimeFrequency_ReactivatedMemoryProject\')

p_clustersOfInterest
clusters = fieldnames(Clust);
clusters{length(clusters)+1}  = 'all';

% filepathD = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/'; %server
% filepathM = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/'; %server

filepathD = 'D:\TF_Calculation_90SecTrial\DNight\'; %Windows
filepathM = 'D:\TF_Calculation_90SecTrial\MNight\'; %Windows

filesOdorD = dir(strcat(filepathD,'*_Odor.mat'));
filesVehicleD = dir(strcat(filepathD,'*_Vehicle.mat'));

filesOdorM = dir(strcat(filepathM,'*_Odor.mat'));
filesVehicleM = dir(strcat(filepathM,'*_Vehicle.mat'));


%%

bands = {'SW';'Delta';'Theta';'Spindle'};

SpindleBand     = [12 16];
DeltaBand       = [1 4];
ThetaBand       = [4 8];
SWBand          = [0.5 2];

load('D:\Thesis_Publication\Using\Time_Freq_Parms')
Time_Series_OdorD =[];
%--------------------------------------------------------------------------
% Analysis by bands
%--------------------------------------------------------------------------
for band = 1:numel(bands)
    
    cfg = [];
    
    % select frequency band
    if strcmp(bands(band),'Spindle')
        cfg.frequency = SpindleBand;
    end
    
    if strcmp(bands(band),'Theta')
        cfg.frequency = ThetaBand;
    end
    
    if strcmp(bands(band),'Delta')
        cfg.frequency = DeltaBand;
    end
    
    if strcmp(bands(band),'SW')
        cfg.frequency = SWBand;
    end
    
    cfg.avgoverrpt      = 'yes';
    cfg.avgoverfreq     = 'yes';
    cfg.channel         = {'all','-E49', '-E48', '-E17', '-E128', '-E32', ...
        '-E1', '-E125', '-E119', '-E113', '-E56', ...
        '-E63', '-E68', '-E73', '-E81', '-E88', ...
        '-E94', '-E99', '-E107'};
    
    Time_Series_OdorD.(bands{band}) = Time_Freq_Parms;
    
    for subj = 1:numel(filesOdorD)
        
        disp(strcat('Sujeto: ',num2str(subj)))
        
        % For D Night
        DNight.Odor = load(strcat(filepathD,filesOdorD(subj).name));
        DNight.Vehicle = load(strcat(filepathD,filesVehicleD(subj).name));
        
        Time_Freq_OdorD = rmfield(DNight.Odor.Time_Freq_Odor,  'trialinfo');
        Time_Freq_VehicleD = rmfield(DNight.Vehicle.Time_Freq_Vehicle  ,  'trialinfo');
        
        Sel_Freq_OdorD = ft_selectdata(cfg, Time_Freq_OdorD);
        Sel_Freq_VehicleD = ft_selectdata(cfg, Time_Freq_VehicleD);
        
        Time_Series_OdorD.(bands{band}).trial{subj} = squeeze(Sel_Freq_OdorD.powspctrm)';
        Time_Series_VehicleD.(bands{band}).trial{subj} = squeeze(Sel_Freq_VehicleD.powspctrm)';
        
        % For M Night
        MNight.Odor = load(strcat(filepathM,filesOdorM(subj).name));
        MNight.Vehicle = load(strcat(filepathM,filesVehicleM(subj).name));
        
        Time_Freq_OdorM = rmfield(MNight.Odor.Time_Freq_Odor,  'trialinfo');
        Time_Freq_VehicleM = rmfield(MNight.Vehicle.Time_Freq_Vehicle  ,  'trialinfo');
        
        Sel_Freq_OdorM = ft_selectdata(cfg, Time_Freq_OdorM);
        Sel_Freq_VehicleM = ft_selectdata(cfg, Time_Freq_VehicleM);
        
        Time_Series_OdorM.(bands{band}).trial{subj} = squeeze(Sel_Freq_OdorM.powspctrm)';
        Time_Series_VehicleM.(bands{band}).trial{subj} = squeeze(Sel_Freq_VehicleM.powspctrm)';
        
    end
    
end

cfg = [];
cfg.channel = 'E2';
figure; ft_singleplotER(cfg,Time_Series_OdorM.Delta);

