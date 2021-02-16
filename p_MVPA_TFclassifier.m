cfg                     = [];
cfg.classifier          = 'lda';
cfg.metric              = {'auc'};
cfg.repeat              = 2;
cfg.sample_dimension    = 1;
cfg.feature_dimension   = 2;
cfg.dimension_names     = {'samples','channels','frequencies','time points'};

[cf_freq{subj,1}, result_freq{subj,1}] = mv_classify(cfg, ...
    cat(1,Time_Freq_Odor.powspctrm,...
    Time_Freq_Vehicle.powspctrm), clabel);