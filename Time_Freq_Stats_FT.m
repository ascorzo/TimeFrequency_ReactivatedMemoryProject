% load('Time_Freq_All_JensData.mat')


addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/')

ft_defaults
addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/qsub')

ft_warning off

%--------------------------------------------------------------------------
% Parameters for Baseline Correction 
%--------------------------------------------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-15 45];
cfg_Bas.baselinetype        = 'zscore';

%load('Time-Freq_DA.mat')

%% Load dummy file

%commonChans = Time_Freq_Cue_ChanSel{1, 1}.label;

filename     = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/dummyfile.set';
sensors = ft_read_sens(filename);

%% Baseline correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative associated odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for subj = 1:length(Time_Freq_DA)
%     
%     Time_Freq_baseline = ...
%         ft_freqbaseline(cfg_Bas, Time_Freq_DA{subj});
%     
%     % Average over trials
%     cfg = []; 
%     cfg.avgoverrpt = 'yes';
%     %cfg.avgoverchan = 'yes';
%     cfg.channel = cluster;
%     
%     Time_Freq_baseline2 = rmfield(Time_Freq_baseline,  'trialinfo');
%     avg_Time_Freq = ft_selectdata(cfg, Time_Freq_baseline2);
%     
%     % Separate Conditions
%     cfg = [];
%     cfg.latency = [-15 15];
%     Time_Freq_Cue{subj} = ft_selectdata(cfg, avg_Time_Freq);
%     
%     cfg = [];
%     cfg.latency = [15 45];
%     Time_Freq_Vehicle{subj} = ft_selectdata(cfg, avg_Time_Freq);
%     
%     Time_Freq_Vehicle{subj}.time = Time_Freq_Cue{subj}.time;
%      
% end

%clear TimeFreq_DA 
filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_ONOFF_OdorD_Night/';
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));

subjects = 1:numel(filesOdor);
for subj = subjects
    temp = load(strcat(filepath,filesOdor(subj).name));
    
    %Average over trials
    cfg = []; 
    cfg.avgoverrpt = 'yes';
    Time_Freq_Cue{subj} = ft_selectdata(cfg, temp.Time_Freq_Odor);
    
    temp = load(strcat(filepath,filesVehicle(subj).name));
    Time_Freq_Vehicle{subj} =  ft_selectdata(cfg, temp.Time_Freq_Vehicle);
    
end

%% statistics


%--------------------
% Parameters
%--------------------

cfg                     = [];
cfg.latency             = [0 30];
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
cfg.numrandomization    = 1000;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

design                  = [];
design(1,:)             = [1:length(Time_Freq_Cue) 1:length(Time_Freq_Vehicle)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(Time_Freq_Cue)) ones(1,length(Time_Freq_Vehicle))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

%%
stats1                  = ft_freqstatistics(cfg, Time_Freq_Cue{:}, Time_Freq_Vehicle{:});

%%

% pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(2,:,:)));
%shading interp;
colorbar;

cfg = [];

[grandavg_odor] = ft_freqgrandaverage(cfg,  Time_Freq_Cue{:});
[grandavg_sham] = ft_freqgrandaverage(cfg,  Time_Freq_Vehicle{:});

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
cfg.elec = Time_Freq_Cue{1}.elec;
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

