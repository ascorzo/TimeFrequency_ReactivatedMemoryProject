% load('Time_Freq_All_JensData.mat')
% 
% for subj = 1:length(Time_Freq_Cue)
%     
%     cfg                     = [];
%     cfg.baseline            = [-2 0];
%     cfg.baselinetype        = 'relative';
%     Time_Freq_Odor_baseline{subj} = ...
%         ft_freqbaseline(cfg, Time_Freq_Cue{subj});
%     Time_Freq_Sham_baseline{subj} = ...
%         ft_freqbaseline(cfg, Time_Freq_Sham_CN{subj});
%     
%     cfg = [];
%     cfg.channel = commonChans;
%     
%     Time_Freq_Cue_ChanSel{subj}   = ft_selectdata(cfg, Time_Freq_Odor_baseline{subj});
%     Time_Freq_Sham_CN_ChanSel{subj}  = ft_selectdata(cfg, Time_Freq_Sham_baseline{subj});
% 
% end
%%

commonChans = Time_Freq_Cue_ChanSel{1, 1}.label;

cfg = []; 
filename     = 'D:\GermanData\PSD\dummyfile.set';
sensors = ft_read_sens(filename);

%%


for n = 1:length(Time_Freq_Sham_CN_ChanSel) 
    cfg = []; 
    cfg.avgoverrpt = 'yes';
    avg_sham_data{n} = ft_selectdata(cfg, Time_Freq_Sham_CN_ChanSel{n});

    cfg = []; 
    cfg.avgoverrpt = 'yes';
    avg_cue_data{n} = ft_selectdata(cfg, Time_Freq_Cue_ChanSel{n});

end

%%


%--------------------
% Parameters
%--------------------

cfg                     = [];
cfg.latency             = [-2 15];
cfg.frequency           = 'all';
cfg.channel             = commonChans;
cfg.correctm            = 'cluster';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';     % use actvsblT for activation against baseline
cfg.clusterstatistic    = 'maxsum';            % statistic used to decide cluster significance (sum of t-values within a cluster)
cfg.minnbchan           = 0;
cfg_neighb.method       = 'triangulation';
cfg_neighb.channel      = commonChans;
cfg_neighb.elec         = sensors; 
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb);
cfg.tail                = 1;
cfg.clustertail         = cfg.tail;

cfg.alpha               = 0.05;
cfg.numrandomization    = 5000;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast

design                  = [];
design(1,:)             = [1:length(avg_sham_data) 1:length(avg_cue_data)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(avg_sham_data)) ones(1,length(avg_cue_data))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;

%%
stats1                  = ft_freqstatistics(cfg, avg_cue_data{:}, avg_sham_data{:});

%%

pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(5,:,:)));
shading interp;
colorbar;





