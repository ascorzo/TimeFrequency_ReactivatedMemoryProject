addpath('/home/andrea/Documents/Github/TimeFrequency_ReactivatedMemoryProject/')

% addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\') %Windows
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828') %server

ft_defaults


p_clustersOfInterest
clusters = fieldnames(Clust);
clusters{length(clusters)+1}  = 'all';

filepathD = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/'; %server
filepathM = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/'; %server

% filepathD = 'D:\TF_Calculation_90SecTrial\DNight\'; %Windows
% filepathM = 'D:\TF_Calculation_90SecTrial\MNight\'; %Windows

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

% load('D:\Thesis_Publication\Using\Time_Freq_Parms')
Time_Series_OdorD =[];
%--------------------------------------------------------------------------
% Analysis by bands
%--------------------------------------------------------------------------
for band = 1%:numel(bands)
    
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
    
    % Time_Series_OdorD.(bands{band}) = Time_Freq_Parms;
    
    for subj = 1:4%numel(filesOdorD)
        
        disp(strcat('Sujeto: ',num2str(subj)))
        
        % For D Night
        DNight.Odor = load(strcat(filepathD,filesOdorD(subj).name));
        DNight.Vehicle = load(strcat(filepathD,filesVehicleD(subj).name));
        
        Sel_Freq_OdorD = ft_selectdata(cfg, DNight.Odor.Time_Freq_Odor);
        Sel_Freq_VehicleD = ft_selectdata(cfg, DNight.Vehicle.Time_Freq_Vehicle);
        
        OdorD_trial{1,subj} = squeeze(Sel_Freq_OdorD.powspctrm);
        VehicleD_trial{1,subj} = squeeze(Sel_Freq_VehicleD.powspctrm);

        
        % For M Night
        MNight.Odor = load(strcat(filepathM,filesOdorM(subj).name));
        MNight.Vehicle = load(strcat(filepathM,filesVehicleM(subj).name));
        
        Sel_Freq_OdorM = ft_selectdata(cfg, MNight.Odor.Time_Freq_Odor);
        Sel_Freq_VehicleM = ft_selectdata(cfg, MNight.Vehicle.Time_Freq_Vehicle);
        
        OdorM_trial{1,subj} = squeeze(Sel_Freq_OdorM.powspctrm);
        VehicleM_trial{1,subj} = squeeze(Sel_Freq_VehicleM.powspctrm);
        
        time_all{1,subj} = Sel_Freq_OdorM.time;
    end
    
    Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD,'freq');
%     Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD, 'trialinfo');
    Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD,'powspctrm');
    
    Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD,'freq');
%     Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD, 'trialinfo');
    Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD,'powspctrm');
    
    Sel_Freq_OdorD.dimord = 'chan_time';
    Sel_Freq_VehicleD.dimord = 'chan_time';
    
    Sel_Freq_OdorD.time      = time_all;
    Sel_Freq_VehicleD.time   = time_all;
    
    Time_Series_OdorD.(bands{band}) = Sel_Freq_OdorD;
    Time_Series_VehicleD.(bands{band}) = Sel_Freq_VehicleD;
    
    Time_Series_OdorD.(bands{band}).trial = OdorD_trial;
    Time_Series_VehicleD.(bands{band}).trial = VehicleD_trial;

    %--- M Night---
    
    Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM,'freq');
%     Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM, 'trialinfo');
    Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM,'powspctrm');
    
    Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM,'freq');
%     Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM, 'trialinfo');
    Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM,'powspctrm');
    
    Sel_Freq_OdorM.dimord = 'chan_time';
    Sel_Freq_VehicleM.dimord = 'chan_time';
    
    Sel_Freq_OdorM.time      = time_all;
    Sel_Freq_VehicleM.time   = time_all;
    
    Time_Series_OdorM.(bands{band}) = Sel_Freq_OdorM;
    Time_Series_VehicleM.(bands{band}) = Sel_Freq_VehicleM;
    
    Time_Series_OdorM.(bands{band}).trial = OdorM_trial;
    Time_Series_VehicleM.(bands{band}).trial = VehicleM_trial;
     
end

%filename     = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Using\dummyfile.set'; %windows
filename     = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/dummyfile.set'; %server
sensors = ft_read_sens(filename);

Nsub = subj;

cfg                  = [];
cfg.method           = 'analytic';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'bonferroni'; % correction methodcfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2;
cfg_neighb.elec      = sensors;
cfg_neighb.method    = 'distance';
cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;


Nsubj  = 10;
design = zeros(2, Nsubj*2);
design(1,:) = [1:Nsubj 1:Nsubj];
design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];

cfg.design = [ones(1,Nsubj) ones(1,Nsubj)*2];;

cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;

stat_t = ft_timelockstatistics(cfg,Time_Series_OdorD.SW,...
    Time_Series_OdorM.SW);

% cfg = [];
% cfg.channel = 'E2';
% figure; ft_singleplotER(cfg,Time_Series_OdorD.Spindle);
% figure; ft_singleplotER(cfg,Time_Series_OdorM.Spindle);
