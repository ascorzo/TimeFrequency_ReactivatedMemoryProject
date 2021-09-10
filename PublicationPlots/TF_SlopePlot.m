addpath(genpath('C:\Users\asanch24\Documents\MATLAB\eeglab2019_1'))
addpath(genpath('C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject'))

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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Same TF analysis by bands but with slope
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------
%           Combine Mean of channels by clusters
%----------------------------------------------------------------------
for subj = 1:numel(D_Night.Time_Freq_Odor_Slope)
    for cluster = 1:numel(clusters)
        
        SlopeSubjOdorD = D_Night.Time_Freq_Odor_Slope{subj};
        SlopeSubjVehiD = D_Night.Time_Freq_Vehi_Slope{subj};
        SlopeSubjOdorM = M_Night.Time_Freq_Odor_Slope{subj};
        SlopeSubjVehiM = M_Night.Time_Freq_Vehi_Slope{subj};
        
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),...
            SlopeSubjOdorD.label);
        
        Slope_OdorD_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(SlopeSubjOdorD.powspctrm(ind2,:,:),1));
        
        Slope_VehiD_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(SlopeSubjVehiD.powspctrm(ind2,:,:),1));
        
        Slope_OdorM_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(SlopeSubjOdorM.powspctrm(ind2,:,:),1));
        
        Slope_VehiM_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(SlopeSubjVehiM.powspctrm(ind2,:,:),1));
        
    end
    
    Slope_OdorD_clust.all(subj,:,:) = ...
        squeeze(mean(SlopeSubjOdorD.powspctrm,1));
    
    Slope_VehiD_clust.all(subj,:,:) = ...
        squeeze(mean(SlopeSubjVehiD.powspctrm,1));
    
    Slope_OdorM_clust.all(subj,:,:) = ...
        squeeze(mean(SlopeSubjOdorM.powspctrm,1));
    
    Slope_VehiM_clust.all(subj,:,:) = ...
        squeeze(mean(SlopeSubjVehiM.powspctrm,1));
end

SpindleBand     = [12 16];
DeltaBand       = [1 4];
Around_2_5      = [2 8];
ThetaBand       = [4 8];
SWBand          = [0.5 2];
%BetaBand       = [18 30];
subjects        = 1:23; %numel(filesOdor);

frequencies = D_Night.Time_Freq_Odor_Slope{1, 1}.freq;

SpindleIdx = find(frequencies>=SpindleBand(1) &...
    frequencies<=SpindleBand(2));
DeltaIdx = find(frequencies>=DeltaBand(1) &...
    frequencies<=DeltaBand(2));
ThetaIdx = find(frequencies>=ThetaBand(1) &...
    frequencies<=ThetaBand(2));
Around_2_5Idx = find(frequencies>=Around_2_5(1) &...
    frequencies<=Around_2_5(2));
SWIdx = find(frequencies>=SWBand(1) &...
    frequencies<=SWBand(2));

clusters = fieldnames(Clust);
clusters{length(clusters)+1}  = 'all';

for cluster = 1:numel(clusters)
    % -----------------For Odor D ----------------------------------
    v_Spindle_OdorD.(clusters{cluster}) = squeeze(mean(Slope_OdorD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorD.(clusters{cluster}) = squeeze(mean(Slope_OdorD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorD.(clusters{cluster}) = squeeze(mean(Slope_OdorD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_OdorD.(clusters{cluster}) = squeeze(mean(Slope_OdorD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorD.(clusters{cluster}) = squeeze(mean(Slope_OdorD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle D----------------------------------
    v_Spindle_VehicleD.(clusters{cluster}) = squeeze(mean(Slope_VehiD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleD.(clusters{cluster}) = squeeze(mean(Slope_VehiD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleD.(clusters{cluster}) = squeeze(mean(Slope_VehiD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleD.(clusters{cluster}) = squeeze(mean(Slope_VehiD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_VehicleD.(clusters{cluster}) = squeeze(mean(Slope_VehiD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));
    
    % -----------------For Odor M ----------------------------------
    v_Spindle_OdorM.(clusters{cluster}) = squeeze(mean(Slope_OdorM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorM.(clusters{cluster}) = squeeze(mean(Slope_OdorM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorM.(clusters{cluster}) = squeeze(mean(Slope_OdorM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_OdorM.(clusters{cluster}) = squeeze(mean(Slope_OdorM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorM.(clusters{cluster}) = squeeze(mean(Slope_OdorM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle M----------------------------------
    v_Spindle_VehicleM.(clusters{cluster}) = squeeze(mean(Slope_VehiM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleM.(clusters{cluster}) = squeeze(mean(Slope_VehiM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleM.(clusters{cluster}) = squeeze(mean(Slope_VehiM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleM.(clusters{cluster}) = squeeze(mean(Slope_VehiM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_VehicleM.(clusters{cluster}) = squeeze(mean(Slope_VehiM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

for cluster = 1:numel(clusters)
    
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = D_Night.Time_Freq_Odor_Slope{1}.time;
    total_subplots_row      = 6;
    total_subplots_column   = 3;
    count                   = 1;
    frequencies             = D_Night.Time_Freq_Odor_Slope{1}.freq;
    v_xlim                  = [-5 20];
    v_xticks                = [-5 0 5 10 15 20 25];
    
    
    TF_OdorD = squeeze(mean(Slope_OdorD_clust.(clusters{cluster}),1));
    TF_VehicleD = squeeze(mean(Slope_VehiD_clust.(clusters{cluster}),1));
    
    TF_OdorM = squeeze(mean(Slope_OdorM_clust.(clusters{cluster}),1));
    TF_VehicleM = squeeze(mean(Slope_VehiM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD(:)),min(TF_VehicleD(:)),...
        min(TF_OdorM(:)),min(TF_VehicleM(:))])*0.6;
    top_TF     = max([max(TF_OdorD(:)),max(TF_VehicleD(:)),...
        max(TF_OdorM(:)),max(TF_VehicleM(:))])*0.6;
    
    max_TF = max(abs(bottom_TF),abs(top_TF));
    
    % -- ylims Spindles
    bottom_Spindles  = min([min(v_Spindle_OdorD.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        min(v_Spindle_OdorM.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Spindles     = max([max(v_Spindle_OdorD.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        max(v_Spindle_OdorM.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Theta
    bottom_Theta  = min([min(v_Theta_OdorD.(clusters{cluster})(:)),...
        min(v_Theta_VehicleD.(clusters{cluster})(:)),...
        min(v_Theta_OdorM.(clusters{cluster})(:)),...
        min(v_Theta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Theta     = max([max(v_Theta_OdorD.(clusters{cluster})(:)),...
        max(v_Theta_VehicleD.(clusters{cluster})(:)),...
        max(v_Theta_OdorM.(clusters{cluster})(:)),...
        max(v_Theta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Delta
    bottom_Delta  = min([min(v_Delta_OdorD.(clusters{cluster})(:)),...
        min(v_Delta_VehicleD.(clusters{cluster})(:)),...
        min(v_Delta_OdorM.(clusters{cluster})(:)),...
        min(v_Delta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Delta     = max([max(v_Delta_OdorD.(clusters{cluster})(:)),...
        max(v_Delta_VehicleD.(clusters{cluster})(:)),...
        max(v_Delta_OdorM.(clusters{cluster})(:)),...
        max(v_Delta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims SW
    bottom_SW  = min([min(v_SW_OdorD.(clusters{cluster})(:)),...
        min(v_SW_VehicleD.(clusters{cluster})(:)),...
        min(v_SW_OdorM.(clusters{cluster})(:)),...
        min(v_SW_VehicleM.(clusters{cluster})(:))])*0.2;
    top_SW     = max([max(v_SW_OdorD.(clusters{cluster})(:)),...
        max(v_SW_VehicleD.(clusters{cluster})(:)),...
        max(v_SW_OdorM.(clusters{cluster})(:)),...
        max(v_SW_VehicleM.(clusters{cluster})(:))])*0.12;
    
    
    start_sample = 1;
    end_sample = length(time_parcial);
    

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
   
    
    % -- M Night ----  
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,2,...
        'Position',[0.5233798195242,0.855206042570128,0.305988515176375,0.089345535800919]);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Odor M')
    xlabel('')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,5,...
        'Position',[0.525840853158327,0.72179636669222,0.305988515176375,0.089345535800919]);
    f_ImageMatrix(TF_VehicleM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Vehicle')
    xlabel('Time (s)')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
        % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,1,...
        'Position',[0.171222313371618,0.855206042570114,0.306218211648892,0.089345535800919]);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    ylabel('Frequency (Hz)','Position',[-6.884988259404048,-3.446376583590895,1])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Odor D')
    xlabel('')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,4,...
        'Position',[0.171452009844135,0.719076938180593,0.305373256767845,0.0923576200811]);
    
    f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylim([0 20])
    colormap(parula)
    caxis manual
    caxis([-max_TF max_TF]);
    title('Vehicle')
    xlabel('Time (s)')
    hold on 
    plot([0,15],[0.25,0.25],'k','LineWidth',3);
    hold off
    
%     % -- Both Odors ----  
%     
%     tf(count) = subplot(total_subplots_row,total_subplots_column,3);
%     f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
%     xlim(v_xlim)
%     colormap(parula)
%     caxis manual
%     caxis([bottom_TF top_TF]);
%     title('Odor D')
%     xlabel('')
%     
%     
%     count = count+1;
%     tf(count) = subplot(total_subplots_row,total_subplots_column,6);
%     f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
%     xlim(v_xlim)
%     colormap(parula)
%     caxis manual
%     caxis([bottom_TF top_TF]);
%     title('Odor M')
%     xlabel('')
    a = colorbar('Position',[0.846735996404601,0.713822894168464,0.011321300281366,0.235421162817687]);
    a.Label.String = 'Power (Z)';
    
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,7);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,8);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Spindle_OdorM.(clusters{cluster}),...
        v_Spindle_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,9);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),{' '},'-',{' '},num2str(SpindleBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    hold on 
    plot([0,15],[bottom_Spindles,bottom_Spindles],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot Theta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,10);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    ylabel('Power (Z)','position',[-9.497365001570415,-0.116481838204675,-1])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,11);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Theta_OdorM.(clusters{cluster}),...
        v_Theta_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,12);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),{' '},'-',{' '},num2str(ThetaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Theta,bottom_Theta],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,13);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,14);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_Delta_OdorM.(clusters{cluster}),...
        v_Delta_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
   
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,15);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),{' '},'-',{' '},num2str(DeltaBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_Delta,bottom_Delta],'k','LineWidth',3);
    hold off
    
    %----------------------------------------------------------------------
    % Plot SW Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,16);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Vehicle',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_VehicleD.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],'-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.265927658827671,0.038916064881278,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xticks([-5 0 5 10 15 20 25]) 
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,17);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor M','Vehicle',...
        v_SW_OdorM.(clusters{cluster}),...
        v_SW_VehicleM.(clusters{cluster}),...
        time_parcial,[ 0.27059 ,0.45882,0.70588],'-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.545970189311835,0.041075892258916,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    
    
    % -- Both Odors ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,18);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),{' '},'-',{' '},num2str(SWBand(2)),' Hz)'),...
        'Time (s)',' ','Odor D','Odor M',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_OdorM.(clusters{cluster}),...
        time_parcial,[0.84314  0.18824  0.15294],[ 0.27059 ,0.45882,0.70588],start_sample,end_sample);
    legend('Position',...
        [0.825707679057529,0.039995978652004,0.085110747295877,0.040631750283952],...
        'AutoUpdate','off');
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    xticks([-5 0 5 10 15 20 25])
    hold on 
    plot([0,15],[bottom_SW,bottom_SW],'k','LineWidth',3);
    

    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    set(gcf,'position',[345.8,41.8,975.2,740.8])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    savepath = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Figures\TF_Publication\Slope\ByBands\';
    saveas(gcf,strcat(savepath,clusters{cluster},'.png'))
    close all
end
