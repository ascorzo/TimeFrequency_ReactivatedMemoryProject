


p_clustersOfInterest
clusters = fieldnames(Clust);

cfg                         = [];
cfg.latency                 = [-1 1];
cfg.avgovertime             = 'yes';

Time_Freq_Odor_aroundCero   = ft_selectdata(cfg, Time_Freq_Odor);
Time_Freq_Vehi_aroundCero   = ft_selectdata(cfg, Time_Freq_Vehicle);

for cluster = 1:numel(clusters) 
    cfg                 = [];
    cfg.avgoverchan     = 'yes';
    cfg.channel         = Clust.(clusters{cluster});
    cfg.frequency       = [4 8];
    cfg.avgoverfreq     = 'yes';
    
    Time_Freq_OdorD_clust.(clusters{cluster}) = ...
        ft_selectdata(cfg,Time_Freq_Odor_aroundCero);
    
    Time_Freq_VehicleD_clust.(clusters{cluster}){subj} = ...
        ft_selectdata(cfg,Time_Freq_Vehi_aroundCero);
end



