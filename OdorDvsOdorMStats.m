addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load and Organize files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------

%filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Filtered-David/CueNight/';

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));

for subj = 1:numel(filesOdor)
    
  
    disp(strcat('Sujeto: ',num2str(subj)))

    load(strcat(filepath,filesOdor(subj).name));   

    %p_plotSubject_ConMartin
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_OdorD{subj} = ft_selectdata(cfg, Time_Freq_Cue_baseline2);

end

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));

for subj = 1:numel(filesOdor)

    disp(strcat('Sujeto: ',num2str(subj)))


    load(strcat(filepath,filesOdor(subj).name));   
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    %cfg.latency = [-5 25];
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_Odor,  'trialinfo');
    Time_Freq_OdorM{subj} = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
end

%%

filename     = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/dummyfile.set';
sensors = ft_read_sens(filename);

%%

cfg                     = [];
cfg.latency             = 'all';
cfg.frequency           = 'all';
cfg.channel             = 'all';
cfg.correctm            = 'cluster';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';     % use actvsblT for activation against baseline
cfg.clusterstatistic    = 'maxsum';            % statistic used to decide cluster significance (sum of t-values within a cluster)
cfg.minnbchan           = 2;

cfg_neighb.method       = 'distance';
cfg_neighb.channel      = 'all';
cfg_neighb.elec         = sensors; 
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb);

cfg.tail                = 0;
cfg.clustertail         = cfg.tail;
cfg.alpha               = 0.025;
cfg.numrandomization    = 500;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

design                  = [];
design(1,:)             = [1:length(Time_Freq_OdorD) 1:length(Time_Freq_OdorM)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(Time_Freq_OdorD)) ones(1,length(Time_Freq_OdorM))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

%%
stats1                  = ft_freqstatistics(cfg, Time_Freq_OdorD{:}, Time_Freq_OdorM{:});

%%
% pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(2,:,:)));
%shading interp;
colorbar;

cfg = [];

[grandavg_odor] = ft_freqgrandaverage(cfg,  Time_Freq_OdorD{:});
[grandavg_sham] = ft_freqgrandaverage(cfg,  Time_Freq_OdorM{:});

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)';
tfr_difference = ft_math(cfg, grandavg_odor, grandavg_sham);
%tfr_difference.powspctrm = tfr_difference.powspctrm.*double(stats1.negclusterslabelmat==1);

cfg = [];
cfg.frequency = [stats1.freq(1) stats1.freq(end)];
cfg.latency   = [stats1.time(1) stats1.time(end)];

% probs = tfr_difference;
% probs.powspctrm = stats1.posclusterslabelmat;
% probs.time = stats1.time;
% 
% 


this_tfr = ft_selectdata(cfg, tfr_difference);
this_tfr.mask = stats1.mask;

cfg      = []; 
cfg.elec = Time_Freq_OdorD{1}.elec;
cfg.layout    = ft_prepare_layout(cfg);
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';

cfg.zlim = [-0.5 0.5];
ft_multiplotTFR(cfg, this_tfr);
c = colorbar('location', 'southoutside');
c.Label.String = '(Power ratio odor over sham)';


% cfg.zlim = [-0.2 0.2];
% ft_singleplotTFR(cfg, this_tfr);
% c = colorbar('location', 'southoutside');
% c.Label.String = '(Power ratio odor over sham)';