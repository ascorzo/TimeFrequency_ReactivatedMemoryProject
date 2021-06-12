% p_FirstLastTrial statistics

% load('Time_Freq_All_JensData.mat')


addpath('C:\Users\lanan\Documents\MATLAB\fieldtrip-20190828\')

ft_defaults
addpath('C:\Users\lanan\Documents\MATLAB\fieldtrip-20190828\qsub')

ft_warning off

%% Load dummy file

%commonChans = Time_Freq_Cue_ChanSel{1, 1}.label;

filename     = 'G:\Mi unidad\2021\AnalysisTemp\dummyfile.set';
sensors = ft_read_sens(filename);

%% statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative associated odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%--------------------
% Parameters
%--------------------

cfg                     = [];
cfg.latency             = [-5 5];
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
cfg.alpha               = 0.05;
cfg.numrandomization    = 500;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

design                  = [];
design(1,:)             = [1:length(TF_OdorLastTrial) 1:length(TF_VehicleLastTrial)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(TF_OdorLastTrial)) ones(1,length(TF_VehicleLastTrial))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

%%
stats1                  = ft_freqstatistics(cfg, TF_OdorLastTrial{:}, TF_VehicleLastTrial{:});

%%

% pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(2,:,:)));
%shading interp;
colorbar;

cfg = [];

[grandavg_odor] = ft_freqgrandaverage(cfg,  TF_OdorLastTrial{:});
[grandavg_sham] = ft_freqgrandaverage(cfg,  TF_VehicleLastTrial{:});

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
cfg.elec = TF_OdorLastTrial{1}.elec;
cfg.layout    = ft_prepare_layout(cfg);
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';

cfg.zlim = [-0.2 0.2];
ft_multiplotTFR(cfg, this_tfr);
c = colorbar('location', 'southoutside');
c.Label.String = '(Power ratio odor over sham)';


% cfg.zlim = [-0.2 0.2];
% ft_singleplotTFR(cfg, this_tfr);
% c = colorbar('location', 'southoutside');
% c.Label.String = '(Power ratio odor over sham)';
