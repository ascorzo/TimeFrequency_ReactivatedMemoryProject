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
cfg.numrandomization    = 100;
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

design                  = [];
design(1,:)             = [1:size(Time_Freq_Odor.powspctrm,1) 1:size(Time_Freq_Vehicle.powspctrm,1)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,size(Time_Freq_Odor.powspctrm,1)) ones(1,size(Time_Freq_Vehicle.powspctrm,1))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

%%
stats1                  = ft_freqstatistics(cfg, Time_Freq_Odor, Time_Freq_Vehicle);

%% Plot Statistics

%Average over trials
cfg = []; 
cfg.avgoverrpt = 'yes';
Time_Freq_Odor_avg = ft_selectdata(cfg, Time_Freq_Odor);
Time_Freq_Vehicle_avg =  ft_selectdata(cfg, Time_Freq_Vehicle);

cfg = [];
cfg.parameter    = 'powspctrm';
cfg.operation    = '(x1-x2)';
tfr_difference = ft_math(cfg, Time_Freq_Odor_avg, Time_Freq_Vehicle_avg);
%tfr_difference.powspctrm = tfr_difference.powspctrm.*double(stats1.negclusterslabelmat==1);

cfg = [];
cfg.frequency = [stats1.freq(1) stats1.freq(end)];
cfg.latency   = [stats1.time(1) stats1.time(end)];
this_tfr      = ft_selectdata(cfg, tfr_difference);
this_tfr.mask = stats1.mask;

cfg                 = []; 
cfg.elec            = Time_Freq_Odor.elec;
cfg.layout          = ft_prepare_layout(cfg);
cfg.parameter       = 'powspctrm';
cfg.maskparameter   = 'mask';
cfg.maskstyle       = 'outline';

%cfg.zlim = [-0.5 0.5];
ft_multiplotTFR(cfg, this_tfr);
c = colorbar('location', 'southoutside');
c.Label.String = '(Power ratio Odor over Vehicle)';