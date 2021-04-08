cfg                             = [];
cfg.classifier                  = 'svm';
cfg.metric                      = {'auc'};
cfg.repeat                      = 2;
cfg.k                           = 5;
cfg.sample_dimension            = 1;
cfg.feature_dimension           = [2,4];
cfg.generalization_dimension    = 3;
cfg.dimension_names     = {'samples','channels','frequencies','time points'};

[cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
    cat(1,Time_Freq_Odor.powspctrm,...
    Time_Freq_Vehicle.powspctrm), clabel);