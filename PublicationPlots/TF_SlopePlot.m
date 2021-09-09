addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')

ft_defaults
addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\qsub')

ft_warning off

p_clustersOfInterest
clusters = fieldnames(Clust);


%%

load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat')

%--------------------------------------------------------------------------
% D Night: Odor vs Vehicle
%--------------------------------------------------------------------------

for cluster = 1:numel(clusters)
    
    for subj = 1:size(Time_Freq_Odor_Slope,2)
        SlopeTemp_Odor = Time_Freq_Odor_Slope{subj};
        SlopeTemp_Vehi = Time_Freq_Vehi_Slope{subj};
        
        cfg             = [];
        cfg.channel     = intersect(SlopeTemp_Odor.label,Clust.(clusters{cluster}));
        cfg.avgoverchan = 'yes';
        SlopeTemp_Odor       = ft_selectdata(cfg,SlopeTemp_Odor);
        SlopeTemp_Vehi       = ft_selectdata(cfg,SlopeTemp_Vehi);
        
        SlopeClust_Odor(subj,:,:)= squeeze(SlopeTemp_Odor.powspctrm);
        SlopeClust_Vehi(subj,:,:)= squeeze(SlopeTemp_Vehi.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_Odor.freq)
        for t = 1:length(SlopeTemp_Odor.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_Odor(:,f,t),SlopeClust_Vehi(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_Odor,1))-squeeze(mean(SlopeClust_Vehi,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_Odor.time,SlopeTemp_Odor.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\D_Night\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_DNight.png'))
    close all
end

%%
%--------------------------------------------------------------------------
% M Night: Odor vs Vehicle
%--------------------------------------------------------------------------
load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat')

for cluster = 1:numel(clusters)
    
    for subj = 1:size(Time_Freq_Odor_Slope,2)
        SlopeTemp_Odor = Time_Freq_Odor_Slope{subj};
        SlopeTemp_Vehi = Time_Freq_Vehi_Slope{subj};
        
        cfg             = [];
        cfg.channel     = intersect(SlopeTemp_Odor.label,Clust.(clusters{cluster}));
        cfg.avgoverchan = 'yes';
        SlopeTemp_Odor       = ft_selectdata(cfg,SlopeTemp_Odor);
        SlopeTemp_Vehi       = ft_selectdata(cfg,SlopeTemp_Vehi);
        
        SlopeClust_Odor(subj,:,:)= squeeze(SlopeTemp_Odor.powspctrm);
        SlopeClust_Vehi(subj,:,:)= squeeze(SlopeTemp_Vehi.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_Odor.freq)
        for t = 1:length(SlopeTemp_Odor.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_Odor(:,f,t),SlopeClust_Vehi(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_Odor,1))-squeeze(mean(SlopeClust_Vehi,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_Odor.time,SlopeTemp_Odor.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\M_Night\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_MNight.png'))
    close all
end


%%
% clear Time_Freq_Odor_Slope Time_Freq_Vehi_Slope
%--------------------------------------------------------------------------
% Odor D vs Odor M
%--------------------------------------------------------------------------
D_Night = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat');
M_Night = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat');

for cluster = 1:numel(clusters)
    
    for subj = 1:size(D_Night.Time_Freq_Odor_Slope,2)
        SlopeTemp_OdorD = D_Night.Time_Freq_Odor_Slope{subj};
        SlopeTemp_OdorM = M_Night.Time_Freq_Odor_Slope{subj};
        
        cfg             = [];
        cfg.channel     = intersect(SlopeTemp_OdorD.label,Clust.(clusters{cluster}));
        cfg.avgoverchan = 'yes';
        SlopeTemp_OdorD       = ft_selectdata(cfg,SlopeTemp_OdorD);
        SlopeTemp_OdorM       = ft_selectdata(cfg,SlopeTemp_OdorM);
        
        SlopeClust_OdorD(subj,:,:)= squeeze(SlopeTemp_OdorD.powspctrm);
        SlopeClust_OdorM(subj,:,:)= squeeze(SlopeTemp_OdorM.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_OdorD.freq)
        for t = 1:length(SlopeTemp_OdorD.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_OdorD(:,f,t),SlopeClust_OdorM(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_OdorD,1))-squeeze(mean(SlopeClust_OdorM,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_OdorD.time,SlopeTemp_OdorD.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\BothOdors\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_DvsM.png'))
    close all
end


%% Fieldtrip cluster permutation statistics
% 
filename     = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Using\dummyfile.set';
sensors = ft_read_sens(filename);
%%
%--------------------
% Parameters
%--------------------

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
design(1,:)             = [1:length(D_Night.Time_Freq_Odor_Slope) 1:length(M_Night.Time_Freq_Odor_Slope)];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,length(D_Night.Time_Freq_Odor_Slope)) ones(1,length(M_Night.Time_Freq_Odor_Slope))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;



%% Simple statistical analysis
% 
% cfg = [];
% cfg.method = 'analytic'; %% do a mass-univariatetest
% cfg.alpha = 0.05; %% critical value around ± 2.09
% cfg.statistic = 'depsamplesT'; % use a dependent samples t-test
% 
% design                  = [];
% design(1,:)             = [1:length(D_Night.Time_Freq_Odor_Slope) 1:length(M_Night.Time_Freq_Odor_Slope)];        % conditions, eg:   1 1 1 1 2 2 2 2
% design(2,:)             = [ones(1,length(D_Night.Time_Freq_Odor_Slope)) ones(1,length(M_Night.Time_Freq_Odor_Slope))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
% cfg.design              = design;
% cfg.uvar                = 1;                 % condition (uvar would be the subjects)
% cfg.ivar                = 2;


%%
stats1                  = ft_freqstatistics(cfg, D_Night.Time_Freq_Odor_Slope{:}, M_Night.Time_Freq_Odor_Slope{:});

%%

% pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(2,:,:)));
%shading interp;
colorbar;

cfg = [];

[grandavg_odor] = ft_freqgrandaverage(cfg,  D_Night.Time_Freq_Odor_Slope{:});
[grandavg_sham] = ft_freqgrandaverage(cfg,  M_Night.Time_Freq_Odor_Slope{:});

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
cfg.elec = D_Night.Time_Freq_Odor_Slope{1}.elec;
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



%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Just around 5 sec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat')

%--------------------------------------------------------------------------
% D Night: Odor vs Vehicle
%--------------------------------------------------------------------------

for cluster = 1:numel(clusters)
    
    for subj = 1:size(Time_Freq_Odor_Slope,2)
        SlopeTemp_Odor = Time_Freq_Odor_Slope{subj};
        SlopeTemp_Vehi = Time_Freq_Vehi_Slope{subj};
        
        cfg                     = [];
        cfg.channel             = intersect(SlopeTemp_Odor.label,Clust.(clusters{cluster}));
        cfg.avgoverchan         = 'yes';
        cfg.latency             = [-5 5];
        SlopeTemp_Odor          = ft_selectdata(cfg,SlopeTemp_Odor);
        SlopeTemp_Vehi          = ft_selectdata(cfg,SlopeTemp_Vehi);
        
        SlopeClust_Odor(subj,:,:)= squeeze(SlopeTemp_Odor.powspctrm);
        SlopeClust_Vehi(subj,:,:)= squeeze(SlopeTemp_Vehi.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_Odor.freq)
        for t = 1:length(SlopeTemp_Odor.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_Odor(:,f,t),SlopeClust_Vehi(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_Odor,1))-squeeze(mean(SlopeClust_Vehi,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_Odor.time,SlopeTemp_Odor.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\D_Night\5sec\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_DNight.png'))
    close all
end

%%
%--------------------------------------------------------------------------
% M Night: Odor vs Vehicle
%--------------------------------------------------------------------------
load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat')

for cluster = 1:numel(clusters)
    
    for subj = 1:size(Time_Freq_Odor_Slope,2)
        SlopeTemp_Odor = Time_Freq_Odor_Slope{subj};
        SlopeTemp_Vehi = Time_Freq_Vehi_Slope{subj};
        
        cfg             = [];
        cfg.channel     = intersect(SlopeTemp_Odor.label,Clust.(clusters{cluster}));
        cfg.avgoverchan = 'yes';
        cfg.latency             = [-5 5];
        SlopeTemp_Odor       = ft_selectdata(cfg,SlopeTemp_Odor);
        SlopeTemp_Vehi       = ft_selectdata(cfg,SlopeTemp_Vehi);
        
        SlopeClust_Odor(subj,:,:)= squeeze(SlopeTemp_Odor.powspctrm);
        SlopeClust_Vehi(subj,:,:)= squeeze(SlopeTemp_Vehi.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_Odor.freq)
        for t = 1:length(SlopeTemp_Odor.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_Odor(:,f,t),SlopeClust_Vehi(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_Odor,1))-squeeze(mean(SlopeClust_Vehi,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_Odor.time,SlopeTemp_Odor.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\M_Night\5sec\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_MNight.png'))
    close all
end


%%

%--------------------------------------------------------------------------
% Odor D vs Odor M
%--------------------------------------------------------------------------
D_Night = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat');
M_Night = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat');

for cluster = 1:numel(clusters)
    
    for subj = 1:size(D_Night.Time_Freq_Odor_Slope,2)
        SlopeTemp_OdorD = D_Night.Time_Freq_Odor_Slope{subj};
        SlopeTemp_OdorM = M_Night.Time_Freq_Odor_Slope{subj};
        
        cfg             = [];
        cfg.channel     = intersect(SlopeTemp_OdorD.label,Clust.(clusters{cluster}));
        cfg.avgoverchan = 'yes';
        cfg.latency           = [-5 5];
        SlopeTemp_OdorD       = ft_selectdata(cfg,SlopeTemp_OdorD);
        SlopeTemp_OdorM       = ft_selectdata(cfg,SlopeTemp_OdorM);
        
        SlopeClust_OdorD(subj,:,:)= squeeze(SlopeTemp_OdorD.powspctrm);
        SlopeClust_OdorM(subj,:,:)= squeeze(SlopeTemp_OdorM.powspctrm);
    end
    
    for f = 1:length(SlopeTemp_OdorD.freq)
        for t = 1:length(SlopeTemp_OdorD.time)
            [~,TF_ttest(f,t)] = ttest(SlopeClust_OdorD(:,f,t),SlopeClust_OdorM(:,f,t));
        end
    end
     
%     zmin = 0.93;
%     zmax = 1;

    Diff = squeeze(mean(SlopeClust_OdorD,1))-squeeze(mean(SlopeClust_OdorM,1));
    lim = max(abs(Diff(:)));
    zmin = -lim;
    zmax = lim;
    
    FigureHandle = figure;
    f_ImageMatrix(Diff,SlopeTemp_OdorD.time,SlopeTemp_OdorD.freq,[zmin zmax])
    titl = title({'Odor - Vehicle';clusters{cluster}},'interpreter','none');
    
    colorbar
    
    colormap(FigureHandle, parula)
    %=======================
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\BothOdors\5sec\';
    saveas(gcf,strcat(savepath,clusters{cluster},'_slope_Diff_DvsM.png'))
    close all
end
%%
for subj = 1:size(D_Night.Time_Freq_Odor_Slope,2)
    
    SlopeTemp_OdorD = D_Night.Time_Freq_Odor_Slope{subj};
    SlopeTemp_OdorM = M_Night.Time_Freq_Odor_Slope{subj};
    
    cfg                     = [];
    cfg.latency             = [-5 5];
    
    Slope_OdorD{subj}       = ft_selectdata(cfg,SlopeTemp_OdorD);
    Slope_OdorM{subj}       = ft_selectdata(cfg,SlopeTemp_OdorM);
end

%%
%--------------------
% Parameters
%--------------------

cfg                     = [];
cfg.latency             = 'all';
cfg.frequency           = 'all';
cfg.channel             = 'all';
cfg.correctm            = 'cluster';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';        % use actvsblT for activation against baseline
cfg.clusterstatistic    = 'maxsum';            % statistic used to decide cluster significance (sum of t-values within a cluster)
cfg.minnbchan           = 2;

cfg_neighb.method       = 'distance';
cfg_neighb.channel      = 'all';
cfg_neighb.elec         = sensors; 
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb);

cfg.tail                = 0;
cfg.clustertail         = cfg.tail;
cfg.alpha               = 0.05;
cfg.numrandomization    = 100;             
cfg.clusteralpha        = 0.05;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.uvar                = 1;                 % condition (uvar would be the subjects)
cfg.ivar                = 2;
% Design the statistical contrast
subjects                = length(Slope_OdorD);

design                  = [];
design(1,:)             = [1:subjects 1:subjects];        % conditions, eg:   1 1 1 1 2 2 2 2
design(2,:)             = [ones(1,subjects) ones(1,subjects)*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;
%%
stats1                  = ft_freqstatistics(cfg, Slope_OdorD{:}, Slope_OdorM{:});

%%

% pcolor(stats1.time, stats1.freq, squeeze(stats1.stat(2,:,:)));
%shading interp;
colorbar;

cfg = [];

[grandavg_odor] = ft_freqgrandaverage(cfg,  D_Night.Time_Freq_Odor_Slope{:});
[grandavg_sham] = ft_freqgrandaverage(cfg,  M_Night.Time_Freq_Odor_Slope{:});

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
cfg.elec = D_Night.Time_Freq_Odor_Slope{1}.elec;
cfg.layout    = ft_prepare_layout(cfg);
cfg.parameter = 'powspctrm';
cfg.maskparameter = 'mask';
cfg.maskstyle = 'outline';


ft_multiplotTFR(cfg, this_tfr);
c = colorbar('location', 'southoutside');
c.Label.String = '(Power ratio odor over sham)';


% cfg.zlim = [-0.2 0.2];
% ft_singleplotTFR(cfg, this_tfr);
% c = colorbar('location', 'southoutside');
% c.Label.String = '(Power ratio odor over sham)';
