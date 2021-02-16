% load('Time_Freq_All_JensData.mat')


addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

load('clusterChans.mat')
ft_warning off

central_channels_small = {'E13';'E6';'E112';'E30';'E7';'E106';'E105';'E37';'E31';...
    'E80';'E87';'E54';'E55';'E79'};

% Clusters of interest
Clust.left_frontal = {...
    'E15', 'E16', 'E11', 'E18', 'E19', 'E22', 'E23', 'E24', 'E26', ...
    'E27', 'E33', 'E38'};
Clust.right_frontal = {...
    'E15', 'E16', 'E11', 'E10', 'E4', 'E9', 'E3', 'E124', 'E2', ...
    'E123', 'E122', 'E121'};
Clust.frontal = {...
    'E3', 'E4', 'E9', 'E10', 'E11', 'E15', 'E16', 'E18', 'E19', ...
    'E22', 'E23', 'E24', 'E124'};
Clust.left_central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55'};
Clust.right_central = {...
    'E6', 'E55', 'E112', 'E106', 'E105', 'E80', 'E87', 'E79'};
Clust.central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55', 'E79', ...
    'E80', 'E87', 'E105', 'E106', 'E112'};
Clust.left_temporal = {...
    'E46', 'E51', 'E45', 'E50', 'E58', 'E56', 'E63'};
Clust.right_temporal = {...
    'E108', 'E102', 'E101', 'E97', 'E96', 'E99', 'E107'};
Clust.left_parietal = {...
    'E53', 'E61', 'E62', 'E72', 'E67', 'E52', 'E60', 'E59', 'E66', ...
    'E65', 'E64', 'E68'};
Clust.right_patietal = {...
    'E62', 'E72', 'E78', 'E77', 'E86', 'E85', 'E84', 'E92', 'E91', ...
    'E90', 'E95', 'E94'};
Clust.parietal = {...
    'E52', 'E61', 'E62', 'E59', 'E60', 'E67', 'E66', 'E72', 'E78', ...
    'E77', 'E86', 'E85', 'E84', 'E92', 'E91','E53'};
Clust.left_occipital = {...
    'E71', 'E70', 'E75', 'E74', 'E69', 'E73'};
Clust.right_occipital = {...
    'E75', 'E76', 'E82', 'E83', 'E88', 'E89'};
Clust.occipital = {...
    'E71', 'E70', 'E74', 'E69', 'E73', 'E75', 'E76', 'E83', 'E82', ...
    'E89', 'E88'};

cluster = Clust.occipital;
%--------------------------------------------------------------------------
% Parameters for Baseline Correction 
%--------------------------------------------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-15 45];
cfg_Bas.baselinetype        = 'zscore';

%load('Time-Freq_DA.mat')

%% Load dummy file

%commonChans = Time_Freq_Cue_ChanSel{1, 1}.label;

filename     = 'dummyfile.set';
sensors = ft_read_sens(filename);

%% Baseline correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative associated odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for subj = 1:length(Time_Freq_DA)
    
    Time_Freq_baseline = ...
        ft_freqbaseline(cfg_Bas, Time_Freq_DA{subj});
    
    % Average over trials
    cfg = []; 
    cfg.avgoverrpt = 'yes';
    %cfg.avgoverchan = 'yes';
    cfg.channel = cluster;
    
    Time_Freq_baseline2 = rmfield(Time_Freq_baseline,  'trialinfo');
    avg_Time_Freq = ft_selectdata(cfg, Time_Freq_baseline2);
    
    % Separate Conditions
    cfg = [];
    cfg.latency = [-15 15];
    Time_Freq_Cue{subj} = ft_selectdata(cfg, avg_Time_Freq);
    
    cfg = [];
    cfg.latency = [15 45];
    Time_Freq_Vehicle{subj} = ft_selectdata(cfg, avg_Time_Freq);
    
    Time_Freq_Vehicle{subj}.time = Time_Freq_Cue{subj}.time;
     
end

%clear TimeFreq_DA 


%% statistics


%--------------------
% Parameters
%--------------------

cfg                     = [];
cfg.latency             = [-15 15];
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
cfg.numrandomization    = 10000;             
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

