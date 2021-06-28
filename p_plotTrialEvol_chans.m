

load('G:\Mi unidad\2021\AnalysisTemp\firstSec_TFcount_DNight_chans.mat');
load('G:\Mi unidad\2021\AnalysisTemp\firstSec_TFcount_MNight_chans.mat');

p_clustersOfInterest
clusters = fieldnames(Clust);

for cluster = 1:numel(clusters)
    Time_Freq_OdorD_evol = Time_Freq_OdorD_evol(1:7,:,:,:);
    Time_Freq_OdorM_evol = Time_Freq_OdorM_evol(1:7,:,:,:);
    Time_Freq_VehicleD_evol = Time_Freq_VehicleD_evol(1:7,:,:,:);
    Time_Freq_VehicleM_evol= Time_Freq_VehicleM_evol(1:7,:,:,:);
end


%% Implementar esto por sujeto para ver la regresión del cambio de theta

for chan = 1:size(Time_Freq_OdorD_evol,2)
    
    for subj = 1:23
        Xt = 1:size(Time_Freq_OdorD_evol,1);
        X = Xt;
        
        OdorEvol    = squeeze(Time_Freq_OdorD_evol(:,chan,3,subj));
        VehicleEvol = squeeze(Time_Freq_VehicleD_evol(:,chan,3,subj));
        DiffEvol    = OdorEvol-VehicleEvol;
        
        mdl_Odor = fitlm(X,OdorEvol, 'linear');
        mdl_Vehicle  = fitlm(X,VehicleEvol, 'linear');
        mdl_Diff  = fitlm(X,DiffEvol, 'linear');
        
        % Get slope of regression
        slope_DNight_Odor(chan,subj) = (mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_DNight_Vehicle(chan,subj) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_DNight_Diff(chan,subj) = (mdl_Diff.Fitted(end) - mdl_Diff.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
    end
    
    
    %% Implementar esto por sujeto para ver la regresión del cambio de theta
    
    for subj = 1:23
        
        Xt = 1:size(Time_Freq_OdorM_evol,1);
        X = Xt;
        
        OdorEvol    = squeeze(Time_Freq_OdorM_evol(:,chan,3,subj));
        VehicleEvol = squeeze(Time_Freq_VehicleM_evol(:,chan,3,subj));
        DiffEvol    = OdorEvol-VehicleEvol;
        
        X = Xt;
        
        mdl_Odor = fitlm(X,OdorEvol, 'linear');
        mdl_Vehicle  = fitlm(X,VehicleEvol, 'linear');
        mdl_Diff  = fitlm(X,DiffEvol, 'linear');
        
        
        % Get slope of regression
        slope_MNight_Odor(chan,subj) = ( mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Vehicle(chan,subj) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Diff(chan,subj) = ( mdl_Diff.Fitted(end) - mdl_Diff.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
    end
    
    %%

    [h,p_OdorSlope(chan)] = ttest(slope_DNight_Odor(chan,:),slope_MNight_Odor(chan,:));
    [h,p_VehicleSlope(chan)] = ttest(slope_DNight_Vehicle(chan,:),slope_MNight_Vehicle(chan,:));
    [h,p_DiffSlope(chan)] = ttest(slope_DNight_Diff(chan,:),slope_MNight_Diff(chan,:));

    % Slope odor vs Vehicle

    [h,p_OvsVSlope_D(chan)] = ttest(slope_DNight_Odor(chan,:),slope_DNight_Vehicle(chan,:));
    [h,p_OvsVSlope_M(chan)] = ttest(slope_MNight_Odor(chan,:),slope_MNight_Vehicle(chan,:));

end

%%

load('G:\Mi unidad\2021\AnalysisTemp\reducedChanlocs.mat');
eeglab nogui


figure
subplot(3,3,1);
topoplot((1-p_OdorSlope),reducedchanlocs,...'pmask',(p_OdorSlope<0.05),...
    'maplimits',...
    [0 1],'electrodes','on',...
    'plotrad', 0.519); colorbar
title('Odor Slope 1-pval D night vs M Night')

subplot(3,3,4);
topoplot(squeeze(mean(slope_DNight_Odor,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_DNight_Odor(:)),max(slope_DNight_Odor(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar

title('Mean Odor Slope D night')

subplot(3,3,7);

topoplot(squeeze(mean(slope_MNight_Odor,2))',...%.*(p_OdorSlope<0.05),...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_MNight_Odor(:)),max(slope_MNight_Odor(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Odor Slope M Night')

subplot(3,3,2);
topoplot((1-p_VehicleSlope),reducedchanlocs,... 'pmask',(p_VehicleSlope<0.05),...
    'maplimits',...
    [0 1],'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Vehicle Slope 1-pval D night vs M Night')

subplot(3,3,5);
topoplot(squeeze(mean(slope_DNight_Vehicle,2))',...
    reducedchanlocs,... 'pmask',(p_VehicleSlope<0.05),...
    'maplimits',[min(slope_DNight_Vehicle(:)),max(slope_DNight_Vehicle(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Vehicle Slope D Night')

subplot(3,3,8);

topoplot(squeeze(mean(slope_MNight_Vehicle,2))',...
    reducedchanlocs,...'pmask',(p_VehicleSlope<0.05),...
    'maplimits',[min(slope_MNight_Vehicle(:)),max(slope_MNight_Vehicle(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Vehicle Slope M Night')


subplot(3,3,3);
topoplot((1-p_DiffSlope),reducedchanlocs,...'pmask',(p_DiffSlope<0.05),...
    'maplimits',...
    [0 1],'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Odor-Vehicle Slope 1-pval D night vs M Night')

subplot(3,3,6);
topoplot(squeeze(mean(slope_DNight_Diff,2))',...
    reducedchanlocs,... 'pmask',(p_DiffSlope<0.05),...
    'maplimits',[min(slope_DNight_Diff(:)),max(slope_DNight_Diff(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Diff Slope D night')

subplot(3,3,9);

topoplot(squeeze(mean(slope_MNight_Diff,2))',...
    reducedchanlocs,...'pmask',(p_DiffSlope<0.05),...
    'maplimits',[min(slope_MNight_Diff(:)),max(slope_MNight_Diff(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Diff Slope M Night')



%% figure
figure
subplot(3,2,1);
topoplot((1-p_OvsVSlope_D),reducedchanlocs,...'pmask',(p_OvsVSlope_D < 0.05),...
    'maplimits',...
    [0 1],'electrodes','on',...
    'plotrad', 0.519); colorbar
title('Odor vs Vehicle Slope 1-pval D night')

subplot(3,2,3);
topoplot(squeeze(mean(slope_DNight_Odor,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_DNight_Odor(:)),max(slope_DNight_Odor(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Odor Slope D Night')

subplot(3,2,5);
topoplot(squeeze(mean(slope_DNight_Vehicle,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_DNight_Vehicle(:)),max(slope_DNight_Vehicle(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Vehicle Slope D Night')

subplot(3,2,2);
topoplot((1-p_OvsVSlope_M),reducedchanlocs,...'pmask',(p_OvsVSlope_M<0.05),...
    'maplimits',...
    [0 1],'electrodes','on',...
    'plotrad', 0.519); colorbar
title('Odor vs Vehicle Slope 1-pval M night')

subplot(3,2,4);
topoplot(squeeze(mean(slope_MNight_Odor,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_MNight_Odor(:)),max(slope_MNight_Odor(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Odor Slope M Night')


subplot(3,2,6);
topoplot(squeeze(mean(slope_MNight_Vehicle,2))',...
    reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
    'maplimits',[min(slope_MNight_Vehicle(:)),max(slope_MNight_Vehicle(:))]*0.1,...
    'electrodes','on','colormap',parula,...
    'plotrad', 0.519); colorbar
title('Mean Vehicle Slope M Night')



%%

% 'plotchans', idx_clusters, ...



