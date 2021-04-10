cfg                             = [];
cfg.classifier                  = 'libsvm'; % svm for two classes, and libsvm for multi-class
cfg.metric                      = {'accuracy'}; % The following metrics work only for 2 classes: auc dval tval
cfg.repeat                      = 2;
cfg.k                           = 5;
cfg.sample_dimension            = 1;
cfg.feature_dimension           = [2,4];
cfg.generalization_dimension    = 3;
cfg.dimension_names             = {'samples','channels','frequencies','time points'};


[cf_freqxfreq ,result_freq] = mv_classify(cfg, ...
    cat(1,D_Night_Time_Freq_Odor.powspctrm(1:s_minTrials,:,:,:),...
    D_Night_Time_Freq_Vehicle.powspctrm(1:s_minTrials,:,:,:),...
    M_Night_Time_Freq_Odor.powspctrm(1:s_minTrials,:,:,:),...
    M_Night_Time_Freq_Vehicle.powspctrm(1:s_minTrials,:,:,:)), clabel);