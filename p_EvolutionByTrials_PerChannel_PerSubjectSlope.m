addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/')

ft_defaults
addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828/qsub')

ft_warning off

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Odor Night
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

frequency_bands = [0.5 2; 1 4; 4 8; 8 13; 13 16];

DNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_DNight');
MNight_Cycles = load('/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/OdorCycles_MNight');

for subj = 1:numel(filesOdor)
    clear Time_Freq_OdorD_evol Time_Freq_VehicleD_evol  
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
    v_cycles =  DNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
    [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);

    s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
    s_Lastepoch = s_MaxNumStimsIdx;
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
    
    cfg                         = [];
    cfg.trials                  = s_Firstepoch:s_Lastepoch;
    cfg.latency                 = [-1 1];
    cfg.avgovertime             = 'yes';
    
    Time_Freq_Odor_aroundCero   = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_aroundCero   = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    for freq = 1:length(frequency_bands)
        
        cfg                 = [];
        cfg.frequency       = frequency_bands(freq,:);
        cfg.avgoverfreq     = 'yes';
        
        Time_Freq_OdorD_freq = ...
            ft_selectdata(cfg,Time_Freq_Odor_aroundCero);
        
        Time_Freq_VehicleD_freq = ...
            ft_selectdata(cfg,Time_Freq_Vehi_aroundCero);
        
        Time_Freq_OdorD_evol(freq,:,:) = Time_Freq_OdorD_freq.powspctrm;

        Time_Freq_VehicleD_evol(freq,:,:) = Time_Freq_VehicleD_freq.powspctrm;
    end
    
    for chan = 1:size(Time_Freq_OdorD_evol,3)
        % Evolution in Theta
        
        Xt = 1:size(Time_Freq_OdorD_freq.powspctrm,1);
        X = Xt;
        
        OdorEvol    = squeeze(Time_Freq_OdorD_evol(3,:,chan));
        VehicleEvol = squeeze(Time_Freq_VehicleD_evol(3,:,chan));
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
    
    
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Motor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/MNight/';

%files = dir(strcat(filepath,'*.set'));
filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));
p_clustersOfInterest
clusters = fieldnames(Clust);

frequency_bands = [0.5 2; 1 4; 4 8; 8 13; 13 16];

for subj = 1:numel(filesOdor)
    clear Time_Freq_OdorM_evol Time_Freq_VehicleM_evol 
    
    disp(strcat('Sujeto: ',num2str(subj)))

    v_cycles =  MNight_Cycles.OdorCyles.(filesOdor(subj).name(1:6));
    
    [s_MaxNumStims, s_MaxNumStimsIdx] = max(v_cycles);

    s_Firstepoch  = s_MaxNumStimsIdx-s_MaxNumStims+1;
    s_Lastepoch = s_MaxNumStimsIdx;
    
    load(strcat(filepath,filesOdor(subj).name));
    load(strcat(filepath,filesVehicle(subj).name));
    
    cfg                         = [];
    cfg.trials                  = s_Firstepoch:s_Lastepoch;
    cfg.latency                 = [-1 1];
    cfg.avgovertime             = 'yes';
    
    Time_Freq_Odor_aroundCero   = ft_selectdata(cfg, Time_Freq_Odor);
    Time_Freq_Vehi_aroundCero   = ft_selectdata(cfg, Time_Freq_Vehicle);
    
    for freq = 1:length(frequency_bands)
        
        cfg                 = [];
        cfg.frequency       = frequency_bands(freq,:);
        cfg.avgoverfreq     = 'yes';
        
        Time_Freq_OdorM_freq = ...
            ft_selectdata(cfg,Time_Freq_Odor_aroundCero);
        
        Time_Freq_VehicleM_freq = ...
            ft_selectdata(cfg,Time_Freq_Vehi_aroundCero);

        Time_Freq_OdorM_evol(freq,:,:) = Time_Freq_OdorM_freq.powspctrm;
        Time_Freq_VehicleM_evol(freq,:,:) = Time_Freq_VehicleM_freq.powspctrm;
    end
    
    % Evolution in Theta
    
    for chan = 1:size(Time_Freq_OdorM_evol,3)
        Xt = 1:size(Time_Freq_OdorM_freq.powspctrm,1);
        X = Xt;
        
        OdorEvol    = squeeze(Time_Freq_OdorM_evol(3,:,chan));
        VehicleEvol = squeeze(Time_Freq_VehicleM_evol(3,:,chan));
        DiffEvol    = OdorEvol-VehicleEvol;
        
        mdl_Odor = fitlm(X,OdorEvol, 'linear');
        mdl_Vehicle  = fitlm(X,VehicleEvol, 'linear');
        mdl_Diff  = fitlm(X,DiffEvol, 'linear');
        
        % Get slope of regression
        slope_MNight_Odor(chan,subj) = (mdl_Odor.Fitted(end) - mdl_Odor.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Vehicle(chan,subj) = (mdl_Vehicle.Fitted(end) - mdl_Vehicle.Fitted(1) ) / ...
            (numel(X) - 1) ;
        
        slope_MNight_Diff(chan,subj) = (mdl_Diff.Fitted(end) - mdl_Diff.Fitted(1) ) / ...
            (numel(X) - 1) ;
    end
end

for chan = 1:size(slope_MNight_Odor,1)
    
    [h,p_OdorSlope(chan)] = ttest(slope_DNight_Odor(chan,:),slope_MNight_Odor(chan,:));
    [h,p_VehicleSlope(chan)] = ttest(slope_DNight_Vehicle(chan,:),slope_MNight_Vehicle(chan,:));
    [h,p_DiffSlope(chan)] = ttest(slope_DNight_Diff(chan,:),slope_MNight_Diff(chan,:));

    % Slope odor vs Vehicle

    [h,p_OvsVSlope_D(chan)] = ttest(slope_DNight_Odor(chan,:),slope_DNight_Vehicle(chan,:));
    [h,p_OvsVSlope_M(chan)] = ttest(slope_MNight_Odor(chan,:),slope_MNight_Vehicle(chan,:));
end

save('slopeResults.mat','slope_DNight_Odor','slope_MNight_Odor','slope_DNight_Vehicle','slope_MNight_Vehicle',...
    'slope_DNight_Diff','slope_MNight_Diff','p_OdorSlope','p_VehicleSlope','p_DiffSlope','p_OvsVSlope_D','p_OvsVSlope_M')

% %% Plot
% load('reducedChanlocs.mat');
% addpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1')
% eeglab nogui


% figure
% subplot(3,3,1);
% topoplot((1-p_OdorSlope),reducedchanlocs,...
%     'pmask',(p_OdorSlope<0.05),...
%     'maplimits',...
%     [0 1],'electrodes','on'); colorbar
% title('Odor Slope p-val D night vs M Night')

% subplot(3,3,4);
% topoplot(squeeze(mean(slope_DNight_Odor,2))',...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_DNight_Odor(:)),max(slope_DNight_Odor(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar

% title('Mean Odor Slope D night')

% subplot(3,3,7);

% topoplot(squeeze(mean(slope_MNight_Odor,2))',...%.*(p_OdorSlope<0.05),...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_MNight_Odor(:)),max(slope_MNight_Odor(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Odor Slope M Night')

% subplot(3,3,2);
% topoplot((1-p_VehicleSlope),reducedchanlocs,...
%     'pmask',(p_VehicleSlope<0.05),...
%     'maplimits',...
%     [0 1],'electrodes','on','colormap',parula); colorbar
% title('Vehicle Slope p-val D night vs M Night')

% subplot(3,3,5);
% topoplot(squeeze(mean(slope_DNight_Vehicle,2))',...
%     reducedchanlocs,... 'pmask',(p_VehicleSlope<0.05),...
%     'maplimits',[min(slope_DNight_Vehicle(:)),max(slope_DNight_Vehicle(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Vehicle Slope D Night')

% subplot(3,3,8);

% topoplot(squeeze(mean(slope_MNight_Vehicle,2))',...
%     reducedchanlocs,...'pmask',(p_VehicleSlope<0.05),...
%     'maplimits',[min(slope_MNight_Vehicle(:)),max(slope_MNight_Vehicle(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Vehicle Slope M Night')


% subplot(3,3,3);
% topoplot((1-p_DiffSlope),reducedchanlocs,...
%     'pmask',(p_DiffSlope<0.05),...
%     'maplimits',...
%     [0 1],'electrodes','on','colormap',parula); colorbar
% title('Odor-Vehicle Slope p-val D night vs M Night')

% subplot(3,3,6);
% topoplot(squeeze(mean(slope_DNight_Diff,2))',...
%     reducedchanlocs,... 'pmask',(p_DiffSlope<0.05),...
%     'maplimits',[min(slope_DNight_Diff(:)),max(slope_DNight_Diff(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Diff Slope D night')

% subplot(3,3,9);

% topoplot(squeeze(mean(slope_MNight_Diff,2))',...
%     reducedchanlocs,...'pmask',(p_DiffSlope<0.05),...
%     'maplimits',[min(slope_MNight_Diff(:)),max(slope_MNight_Diff(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Diff Slope M Night')



% %% figure
% figure
% subplot(3,2,1);
% topoplot((1-p_OvsVSlope_D),reducedchanlocs,...
%     'pmask',(p_OvsVSlope_D < 0.05),...
%     'maplimits',...
%     [0 1],'electrodes','on'); colorbar
% title('Odor vs Vehicle Slope p-val D night')

% subplot(3,2,3);
% topoplot(squeeze(mean(slope_DNight_Odor,2))',...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_DNight_Odor(:)),max(slope_DNight_Odor(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Odor Slope D Night')

% subplot(3,2,5);
% topoplot(squeeze(mean(slope_DNight_Vehicle,2))',...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_DNight_Vehicle(:)),max(slope_DNight_Vehicle(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Vehicle Slope D Night')

% subplot(3,2,2);
% topoplot((1-p_OvsVSlope_M),reducedchanlocs,...
%     'pmask',(p_OvsVSlope_M<0.05),...
%     'maplimits',...
%     [0 1],'electrodes','on'); colorbar
% title('Odor vs Vehicle Slope p-val M night')

% subplot(3,2,4);
% topoplot(squeeze(mean(slope_MNight_Odor,2))',...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_MNight_Odor(:)),max(slope_MNight_Odor(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Odor Slope M Night')


% subplot(3,2,6);
% topoplot(squeeze(mean(slope_MNight_Vehicle,2))',...
%     reducedchanlocs,...%'pmask',(p_OdorSlope<0.05),...
%     'maplimits',[min(slope_MNight_Vehicle(:)),max(slope_MNight_Vehicle(:))]*0.1,...
%     'electrodes','on','colormap',parula); colorbar
% title('Mean Vehicle Slope M Night')
