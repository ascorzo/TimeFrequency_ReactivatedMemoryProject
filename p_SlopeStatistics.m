addpath('C:\Users\lanan\Documents\MATLAB\fieldtrip-20190828\')
ft_defaults

Subj = Time_Freq_Odor;
Subj = rmfield(Subj,'freq');
Subj = rmfield(Subj,'time');
Subj = rmfield(Subj,'trialinfo');
Subj = rmfield(Subj,'powspctrm');
Subj.time = 1;

for subj = 1:size(slope_DNight_Odor,2)
    Subj.theta_slope = squeeze(slope_DNight_Odor(:,subj));
    AllThetaSlopes_OdorD{subj} =  Subj;
end

for subj = 1:size(slope_MNight_Odor,2)
    Subj.theta_slope = squeeze(slope_MNight_Odor(:,subj));
    AllThetaSlopes_OdorM{subj} =  Subj;
end

for subj = 1:size(slope_DNight_Vehicle,2)
    Subj.theta_slope = squeeze(slope_DNight_Vehicle(:,subj));
    AllThetaSlopes_VehicleD{subj} =  Subj;
end

for subj = 1:size(slope_MNight_Vehicle,2)
    Subj.theta_slope = squeeze(slope_MNight_Vehicle(:,subj));
    AllThetaSlopes_VehicleM{subj} =  Subj;
end

for subj = 1:size(slope_DNight_Vehicle,2)
    Subj.theta_slope = squeeze(slope_DNight_Diff(:,subj));
    AllThetaSlopes_DiffD{subj} =  Subj;
end

for subj = 1:size(slope_MNight_Vehicle,2)
    Subj.theta_slope = squeeze(slope_MNight_Diff(:,subj));
    AllThetaSlopes_DiffM{subj} =  Subj;
end

filename     = 'G:\Mi unidad\2021\AnalysisTemp\dummyfile.set';
sensors = ft_read_sens(filename);

% define the parameters for the statistical comparison
cfg = [];
cfg.parameter           = 'theta_slope';
cfg.method              = 'montecarlo';
cfg.statistic           = 'depsamplesT';
cfg.alpha               = 0.025;
cfg.correctm            = 'bonferroni';
%cfg.correcttail         = 0;
cfg.clusteralpha        = 0.05;
cfg.numrandomization    = 10000;
cfg.feedback            = 'yes';   
cfg_neighb.method       = 'distance';
cfg.minnbchan           = 2;
cfg_neighb.elec         = sensors; 
cfg.neighbours          = ft_prepare_neighbours(cfg_neighb);

Nsub = 23;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat_OdorDvsM = ft_timelockstatistics(cfg, AllThetaSlopes_OdorD{:}, AllThetaSlopes_OdorM{:});   % don't forget the {:}!
stat_VehicleDvsM = ft_timelockstatistics(cfg, AllThetaSlopes_VehicleD{:}, AllThetaSlopes_VehicleM{:});   % don't forget the {:}!
stat_DNight_OdorvsVehicle = ft_timelockstatistics(cfg, AllThetaSlopes_OdorD{:}, AllThetaSlopes_VehicleD{:});   % don't forget the {:}!
stat_MNight_OdorvsVehicle = ft_timelockstatistics(cfg, AllThetaSlopes_OdorM{:}, AllThetaSlopes_VehicleM{:});   % don't forget the {:}!
stat_DiffDvsM = ft_timelockstatistics(cfg, AllThetaSlopes_DiffD{:}, AllThetaSlopes_DiffM{:});   % don't forget the {:}!

%%

load('G:\Mi unidad\2021\AnalysisTemp\Using\reducedChanlocs.mat');
eeglab nogui

figure
subplot(2,2,1);
confidence = 1-stat_OdorDvsM.prob;
topoplot(confidence,...
    reducedchanlocs,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula);
a = colorbar;
a.Label.String = 'Confidence';

caxis([0.9 1])
hold on
no_results = zeros(numel(reducedchanlocs), 1);
idx_clusters = find((stat_OdorDvsM.prob <= 0.05));
topoplot(no_results, reducedchanlocs, ...
    'style', 'blank', ...
    'electrodes', 'pts', ...
    'shading', 'interp', ...
    'headcolor', [0, 0, 0], ...
    'plotchans', idx_clusters, ...
    'hcolor','none',...
    'whitebk','on',...
    'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
title({'Odor D vs Odor M';'significant slope differences'})
hold off

subplot(2,2,3);
topoplot(squeeze(mean(slope_DNight_Odor,2)-mean(slope_MNight_Odor,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_DNight_Odor(:)-slope_MNight_Odor(:)),...
    max(slope_DNight_Odor(:)-slope_MNight_Odor(:))]*0.1,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula); 
a = colorbar;
a.Label.String = 'Slope difference';

title('Mean OdorD - Odor M Slope')

subplot(2,2,2);
confidence = 1-stat_VehicleDvsM.prob;
topoplot(confidence,...
    reducedchanlocs,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula);
a = colorbar;
a.Label.String = 'Confidence';

caxis([0.9 1])
hold on
no_results = zeros(numel(reducedchanlocs), 1);
idx_clusters = find((stat_VehicleDvsM.prob <= 0.05));
topoplot(no_results, reducedchanlocs, ...
    'style', 'blank', ...
    'electrodes', 'pts', ...
    'shading', 'interp', ...
    'headcolor', [0, 0, 0], ...
    'plotchans', idx_clusters, ...
    'hcolor','none',...
    'whitebk','on',...
    'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
title({'Vehicle D vs Vehicle M';'significant slope differences'})
hold off
hold off

subplot(2,2,4);
topoplot(squeeze(mean(slope_DNight_Vehicle,2)-mean(slope_MNight_Vehicle,2))',...
    reducedchanlocs,... 'pmask',(p_VehicleSlope<0.05),...
    'maplimits',[min(slope_DNight_Vehicle(:)-slope_MNight_Vehicle(:)),...
    max(slope_DNight_Vehicle(:)- slope_MNight_Vehicle(:))]*0.1,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula); 
a = colorbar;
a.Label.String = 'Slope difference';
title({'Mean Vehicle Slope';'D Night - M Night'})


%% figure
figure
subplot(2,2,1);
confidence = 1-stat_DNight_OdorvsVehicle.prob;
topoplot(confidence,...
    reducedchanlocs,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula);
a = colorbar;
a.Label.String = 'Confidence';

caxis([0.9 1])
hold on
no_results = zeros(numel(reducedchanlocs), 1);
idx_clusters = find((stat_DNight_OdorvsVehicle.prob <= 0.05));
topoplot(confidence, reducedchanlocs, ...
    'style', 'blank', ...
    'electrodes', 'on', ...
    'shading', 'interp', ...
    'headcolor', [0, 0, 0], ...
    'plotchans', idx_clusters, ...    
    'whitebk','on',...
    'hcolor','none',...
    'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
title({'Odor vs Vehicle';'significant slope differences';'D night'})
hold off


subplot(2,2,3);
topoplot(squeeze(mean(slope_DNight_Odor,2)-mean(slope_DNight_Vehicle,2))',...
    reducedchanlocs,...
    'maplimits',[min(slope_DNight_Odor(:)-slope_DNight_Vehicle(:)),...
    max(slope_DNight_Odor(:)-slope_DNight_Vehicle(:))]*0.1,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula); 
a = colorbar;
a.Label.String = 'Slope difference';
title({'Mean Odor-Vehicle Slope'; 'D Night'})

subplot(2,2,2);
confidence = 1-stat_MNight_OdorvsVehicle.prob;
topoplot(confidence,...
    reducedchanlocs,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula);
a = colorbar;
a.Label.String = 'Confidence';

caxis([0.9 1])
hold on
no_results = zeros(numel(reducedchanlocs), 1);
idx_clusters = find((stat_MNight_OdorvsVehicle.prob <= 0.05));
topoplot(no_results, reducedchanlocs, ...
    'style', 'blank', ...
    'electrodes', 'pts', ...
    'shading', 'interp', ...
    'headcolor', [0, 0, 0], ...
    'plotchans', idx_clusters, ...
    'hcolor','none',...
    'whitebk','on',...
    'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
title({'Odor vs Vehicle';'significant slope differences'; 'M night'})
hold off

subplot(2,2,4);
topoplot(squeeze(mean(slope_MNight_Odor,2)-mean(slope_MNight_Vehicle,2))',...
    reducedchanlocs,...
    'maplimits',[min(slope_MNight_Odor(:)-slope_MNight_Vehicle(:)),...
    max(slope_MNight_Odor(:)-slope_MNight_Vehicle(:))]*0.1,...
    'conv', 'on', ...
    'whitebk','on',...
    'electrodes','on','colormap',parula);
    %'plotrad', 0.519); 
a = colorbar;
a.Label.String = 'Slope difference';
title({'Mean Odor-Vehicle Slope';'M Night'})

