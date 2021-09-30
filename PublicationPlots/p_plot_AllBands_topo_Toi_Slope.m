% plot topo of toi all bands for slope
addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')
ft_defaults

DNight = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat');
MNight = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat');
%=========================================================================
% In order to be more organized with the results, and not have multiple 
% TF figures for each cluster of channels, I would like to check the
% topographical map of the power for the defined frequency bands in the
% critical points where I can see some significance
%--------------------------------------------------------------------------

% The moments of interest are:
% SO: [-1 1],[9 12],[9.5 10.5],[11 13],[15 20]
% Delta: [-0.5 0.5],[3.5 4.5],[7 8],[9 11],[15 20]
% Theta: [-0.5 0.5],[-4.5 -3.5],[7 8],[10 12]
% Spindles: [2 3]

% Toi.SW          = [-5 -3;-3 -1.5;-1.5 0.5;0.5 1;1 2];
% Toi.Delta       = [-5 -3;-3 -1.5;-1.5 0.5;0.5 1;1 2];
% Toi.Theta       = [-5 -3;-3 -1.5;-1.5 0.5;0.5 1;1 2];

Toi.SW          = [-1 0;0 1;1 2;2 3];
Toi.Delta       = [-3 -1; -2 -1;-1 0;0 1;1 2;2 3];
Toi.Theta       = [-5 -4;-3 -2;-1.5 -0.5;-1 0;0 1;1 2;2 3;3 4];
Toi.Spindle     = [-4.5 -3.5;-3.5 -2.5;-2 -1;0 1];

bands = fieldnames(Toi);

SpindleBand     = [12 20];
DeltaBand       = [1 4];
ThetaBand       = [4 8];
SWBand          = [0.5 2];

filename     = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Using\dummyfile.set';
sensors = ft_read_sens(filename);

load('D:\Thesis_Publication\Using\reducedChanlocs_2.mat');

eeglab nogui

%%
for band = 1:numel(bands)
    cfg = [];
    % select frequency band
    if strcmp(bands(band),'Spindle')
        cfg.frequency = SpindleBand;
    end
    
    if strcmp(bands(band),'Theta')
        cfg.frequency = ThetaBand;
    end
    
    if strcmp(bands(band),'Delta')
        cfg.frequency = DeltaBand;
    end
    
    if strcmp(bands(band),'SW')
        cfg.frequency = SWBand;
    end
    
    if strcmp(bands(band),'Spindle')
        cfg.frequency = SpindleBand;
    end
    
    % select Time of interest
    
    for toi = 1:size(Toi.(bands{band}),1)
        
        cfg.avgovertime     = 'yes';
        %cfg.avgoverfreq     = 'yes';
        cfg.channel         = {'all','-E49', '-E48', '-E17', '-E128', '-E32', ...
            '-E1', '-E125', '-E119', '-E113', '-E56', ...
            '-E63', '-E68', '-E73', '-E81', '-E88', ...
            '-E94', '-E99', '-E107'};
        cfg.latency = Toi.(bands{band})(toi,:);
        latency = cfg.latency;
        
        
        for subj = 1:numel(DNight.Time_Freq_Odor_Slope)
                       
            disp(strcat('Sujeto: ',num2str(subj)))
            
            % For D Night
            SlopeOdor = DNight.Time_Freq_Odor_Slope{subj};
            SlopeVehicle = DNight.Time_Freq_Vehi_Slope{subj};
            
            Sel_Time_Freq_OdorD{subj} = ft_selectdata(cfg, SlopeOdor);
            Sel_Time_Freq_VehicleD{subj} = ft_selectdata(cfg, SlopeVehicle);
            
            Sel_Time_Freq_OdorD{subj}.powspctrm = mean(Sel_Time_Freq_OdorD{subj}.powspctrm,2);
            Sel_Time_Freq_VehicleD{subj}.powspctrm = mean(Sel_Time_Freq_VehicleD{subj}.powspctrm,2);
            
            Sel_Time_Freq_OdorD{subj} = rmfield(Sel_Time_Freq_OdorD{subj},'freq');
            Sel_Time_Freq_OdorD{subj} = rmfield(Sel_Time_Freq_OdorD{subj},'time');
            
            Sel_Time_Freq_VehicleD{subj} = rmfield(Sel_Time_Freq_VehicleD{subj},'freq');
            Sel_Time_Freq_VehicleD{subj} = rmfield(Sel_Time_Freq_VehicleD{subj},'time');
            
            Sel_Time_Freq_OdorD{subj}.time = 1;
            Sel_Time_Freq_VehicleD{subj}.time = 1;
            
            
            % For M Night
            SlopeOdor = MNight.Time_Freq_Odor_Slope{subj};
            SlopeVehicle = MNight.Time_Freq_Vehi_Slope{subj};
            
            Sel_Time_Freq_OdorM{subj} = ft_selectdata(cfg, SlopeOdor);
            Sel_Time_Freq_VehicleM{subj} = ft_selectdata(cfg, SlopeVehicle);
            
            Sel_Time_Freq_OdorM{subj}.powspctrm = mean(Sel_Time_Freq_OdorM{subj}.powspctrm,2);
            Sel_Time_Freq_VehicleM{subj}.powspctrm = mean(Sel_Time_Freq_VehicleM{subj}.powspctrm,2);
            
            Sel_Time_Freq_OdorM{subj} = rmfield(Sel_Time_Freq_OdorM{subj},'freq');
            Sel_Time_Freq_OdorM{subj} = rmfield(Sel_Time_Freq_OdorM{subj},'time');
            
            Sel_Time_Freq_VehicleM{subj} = rmfield(Sel_Time_Freq_VehicleM{subj},'freq');
            Sel_Time_Freq_VehicleM{subj} = rmfield(Sel_Time_Freq_VehicleM{subj},'time');
            
            Sel_Time_Freq_OdorM{subj}.time = 1;
            Sel_Time_Freq_VehicleM{subj}.time = 1;
        end
        
        % define the parameters for the statistical comparison
        cfg = [];
        cfg.parameter           = 'powspctrm';
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
        
        Nsub = numel(Sel_Time_Freq_OdorD);
        cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
        cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
        cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
        cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
        
        stat_OdorDvsVehicle      = ft_timelockstatistics(cfg, ...
            Sel_Time_Freq_OdorD{:}, Sel_Time_Freq_VehicleD{:});   % don't forget the {:}!
        
        stat_OdorMvsVehicle      = ft_timelockstatistics(cfg, ...
            Sel_Time_Freq_OdorM{:}, Sel_Time_Freq_VehicleM{:});   % don't forget the {:}!
        
        stat_OdorDvsOdorM      = ft_timelockstatistics(cfg, ...
            Sel_Time_Freq_OdorD{:}, Sel_Time_Freq_OdorM{:});   % don't forget the {:}!
        
        
        % Plot odor D vs Vehicle
        for subj = 1:numel(Sel_Time_Freq_OdorD)
            Difference(subj,:) = Sel_Time_Freq_OdorD{subj}.powspctrm -...
                Sel_Time_Freq_VehicleD{subj}.powspctrm;
        end
        
        Diff = mean(Difference,1);
        
        figure
        confidence = 1-stat_OdorDvsVehicle.prob;
        
        lim = max(abs(Diff));
        topoplot(Diff,...
            reducedchanlocs,...
            'conv', 'on', ...
            'whitebk','on',...
            'electrodes','on',...
            'colormap',parula,...
            'maplimits',[-lim lim]);%if plotting difference
        a = colorbar;
        a.Label.String = 'OdorD - Vehicle';
      
        
        %caxis([0.9 1]) %if plotting confidence

        
        
        hold on
        no_results = zeros(numel(reducedchanlocs), 1);
        idx_clusters = find((stat_OdorDvsVehicle.prob <= 0.05));
        if ~isempty(idx_clusters)
            topoplot(no_results, reducedchanlocs, ...
                'style', 'blank', ...
                'electrodes', 'pts', ...
                'shading', 'interp', ...
                'headcolor', [0, 0, 0], ...
                'plotchans', idx_clusters, ...
                'hcolor','none',...
                'whitebk','on',...
                'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
        end
        
        Text = strcat(bands(band),'[',num2str(latency),']','s');
        title({'Odor D vs Vehicle';Text{1}})
        
        filename_save = strcat('Slope_DNight_',bands(band),'Toi',num2str(toi));
        saveas(gcf,strcat(filename_save{1},'.png'))
        
        hold off
        close all
        
        % Plot odor M vs Vehicle
        figure
        for subj = 1:numel(Sel_Time_Freq_OdorM)
            Difference(subj,:) = Sel_Time_Freq_OdorM{subj}.powspctrm -...
                Sel_Time_Freq_VehicleM{subj}.powspctrm;
        end
        
        Diff = mean(Difference,1);
        
        confidence = 1-stat_OdorMvsVehicle.prob;
        lim = max(abs(Diff));
        topoplot(Diff,...
            reducedchanlocs,...
            'conv', 'on', ...
            'whitebk','on',...
            'electrodes','on',...
            'colormap',parula,...
            'maplimits',[-lim lim]);%if plotting difference
        a = colorbar;
        
        a.Label.String = 'Odor M - Vehicle';
        
        %caxis([0.9 1]) %if plotting confidence

        
        hold on
        no_results = zeros(numel(reducedchanlocs), 1);
        idx_clusters = find((stat_OdorMvsVehicle.prob <= 0.05));
        if ~isempty(idx_clusters)
            topoplot(no_results, reducedchanlocs, ...
                'style', 'blank', ...
                'electrodes', 'pts', ...
                'shading', 'interp', ...
                'headcolor', [0, 0, 0], ...
                'plotchans', idx_clusters, ...
                'hcolor','none',...
                'whitebk','on',...
                'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
        end
        Text = strcat(bands(band),'[',num2str(latency),']','s');
        title({'Odor M vs Vehicle';Text{1}})
        
        filename_save = strcat('Slope_MNight_',bands(band),'Toi',num2str(toi));
        saveas(gcf,strcat(filename_save{1},'.png'))
        hold off
        
        close all
        % Plot odor D vs Odor M
        figure
        
        for subj = 1:numel(Sel_Time_Freq_OdorD)
            Difference(subj,:) = Sel_Time_Freq_OdorD{subj}.powspctrm -...
                Sel_Time_Freq_OdorM{subj}.powspctrm;
        end
        
        Diff = mean(Difference,1);
        
        confidence = 1-stat_OdorDvsOdorM.prob;
        lim = max(abs(Diff));
        topoplot(Diff,...
            reducedchanlocs,...
            'conv', 'on', ...
            'whitebk','on',...
            'electrodes','on',...
            'colormap',parula,...
            'maplimits',[-lim lim]);%if plotting difference
        a = colorbar;
        
        a.Label.String = 'Odor D - Odor M';
        
        %caxis([0.9 1]) %if plotting confidence

        
        hold on
        no_results = zeros(numel(reducedchanlocs), 1);
        idx_clusters = find((stat_OdorDvsOdorM.prob <= 0.05));
        if ~isempty(idx_clusters)
            topoplot(no_results, reducedchanlocs, ...
                'style', 'blank', ...
                'electrodes', 'pts', ...
                'shading', 'interp', ...
                'headcolor', [0, 0, 0], ...
                'plotchans', idx_clusters, ...
                'hcolor','none',...
                'whitebk','on',...
                'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
        end
        
        Text = strcat(bands(band),'[',num2str(latency),']','s');
        title({'Odor D vs Odor M';Text{1}})
        
        filename_save = strcat('Slope_DvsM_',bands(band),'Toi',num2str(toi));
        saveas(gcf,strcat(filename_save{1},'.png'))
        hold off
        
        clear stat_OdorDvsVehicle stat_OdorMvsVehicle stat_OdorDvsOdorM 
        clear Sel_Time_Freq_OdorM Sel_Time_Freq_OdorD 
        clear Sel_Time_Freq_VehicleM Sel_Time_Freq_VehicleD
        close all
    end
    
end
