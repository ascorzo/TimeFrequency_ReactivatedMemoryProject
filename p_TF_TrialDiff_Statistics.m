%% Statistics


for cluster = 1:numel(clusters)
    
    cfg                     = [];
    cfg.method              = 'montecarlo';
    cfg.statistic           = 'depsamplesT'; 
    cfg.numrandomization    = 1000;
    cfg.uvar                = 1;                 % condition (uvar would be the subjects)
    cfg.ivar                = 2;
    
    design                  = [];
    design(1,:)             = [1:length(Time_Freq_Diff_D_clust.(clusters{cluster})) 1:length(Time_Freq_Diff_M_clust.(clusters{cluster}))];        % conditions, eg:   1 1 1 1 2 2 2 2
    design(2,:)             = [ones(1,length(Time_Freq_Diff_D_clust.(clusters{cluster}))) ones(1,length(Time_Freq_Diff_M_clust.(clusters{cluster})))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
    cfg.design              = design;

    %%
    stats1                  = ft_freqstatistics(cfg, Time_Freq_Diff_D_clust.(clusters{cluster}){:}, Time_Freq_Diff_M_clust.(clusters{cluster}){:});


    %% Plot
    cfg = [];

    [grandavg_odor] = ft_freqgrandaverage(cfg,  Time_Freq_Diff_D_clust.(clusters{cluster}){:});
    [grandavg_sham] = ft_freqgrandaverage(cfg,  Time_Freq_Diff_M_clust.(clusters{cluster}){:});

    cfg                 = [];
    cfg.parameter       = 'powspctrm';
    cfg.operation       = '(x1-x2)';
    tfr_difference      = ft_math(cfg, grandavg_odor, grandavg_sham);
    tfr_difference.mask = stats1.mask;

    cfg                = [];
    cfg.parameter      = 'powspctrm';
    cfg.maskparameter  = 'mask';
    ft_singleplotTFR(cfg,tfr_difference)

end

%%

filename     = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/dummyfile.set';

filename     = 'dummyfile.set';
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
cfg.alpha               = 0.05;
cfg.numrandomization    = 500;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

%% Subject by Subject

for subj = 1:length(Time_Freq_Diff_D_All)
    design                  = [];
    design(1,:)             = [1:length(Time_Freq_Diff_D_All{subj}) 1:length(Time_Freq_Diff_M_All{subj})];        % conditions, eg:   1 1 1 1 2 2 2 2
    design(2,:)             = [ones(1,length(Time_Freq_Diff_D_All{subj})) ones(1,length(Time_Freq_Diff_M_All{subj}))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
    cfg.design              = design;
    
    
    stats1{subj}            = ft_freqstatistics(cfg, Time_Freq_Diff_D_All{subj}, Time_Freq_Diff_M_All{subj});
    
end