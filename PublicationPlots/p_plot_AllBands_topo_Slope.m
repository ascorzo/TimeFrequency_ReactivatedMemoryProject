% plot topo of toi all bands
addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\')
ft_defaults


%filepath = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/TF_Calculation_90SecTrial/DNight/';
filepath = 'D:\TF_Calculation_90SecTrial\DNight\';

filesOdor = dir(strcat(filepath,'*_Odor.mat'));
filesVehicle = dir(strcat(filepath,'*_Vehicle.mat'));



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

Toi.SW          = [-1 1;9 12;9.5 10.5;11 13;15 20];
Toi.Delta       = [-0.5 0.5;3.5 4.5;7 8;9 11;15 20];
Toi.Theta       = [-0.5 0.5;-4.5 -3.5;7 8;10 12];
Toi.Spindle     = [2 3];

bands = fieldnames(Toi);

SpindleBand     = [12 16];
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
    cfg.avgoverrpt      = 'yes';
    cfg.avgovertime     = 'yes';
    cfg.avgoverfreq     = 'yes'; 
    cfg.channel         = {'all','-E49', '-E48', '-E17', '-E128', '-E32', ...
                             '-E1', '-E125', '-E119', '-E113', '-E56', ...
                             '-E63', '-E68', '-E73', '-E81', '-E88', ...
                             '-E94', '-E99', '-E107'};
    
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
    
    % select Time of interest
    
    for toi = 1:size(Toi.(bands{band}),1)
        cfg.latency = Toi.(bands{band})(toi,:);
        latency = cfg.latency;
        
        
        for subj = 1:numel(filesOdor)
                       
            disp(strcat('Sujeto: ',num2str(subj)))
            
            load(strcat(filepath,filesOdor(subj).name));
            load(strcat(filepath,filesVehicle(subj).name));
            
            Time_Freq_Odor = rmfield(Time_Freq_Odor,  'trialinfo');
            Time_Freq_Vehicle = rmfield(Time_Freq_Vehicle,  'trialinfo');
            
            Sel_Time_Freq_Odor{subj} = ft_selectdata(cfg, Time_Freq_Odor);
            Sel_Time_Freq_Vehi{subj} = ft_selectdata(cfg, Time_Freq_Vehicle);
            
            Sel_Time_Freq_Odor{subj}.powspctrm = squeeze(Sel_Time_Freq_Odor{subj}.powspctrm);
            Sel_Time_Freq_Vehi{subj}.powspctrm = squeeze(Sel_Time_Freq_Vehi{subj}.powspctrm);
            
            Sel_Time_Freq_Odor{subj} = rmfield(Sel_Time_Freq_Odor{subj},'freq');
            Sel_Time_Freq_Odor{subj} = rmfield(Sel_Time_Freq_Odor{subj},'time');
            
            Sel_Time_Freq_Vehi{subj} = rmfield(Sel_Time_Freq_Vehi{subj},'freq');
            Sel_Time_Freq_Vehi{subj} = rmfield(Sel_Time_Freq_Vehi{subj},'time');
            
            Sel_Time_Freq_Odor{subj}.time = 1;
            Sel_Time_Freq_Vehi{subj}.time = 1;
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
        
        Nsub = numel(Sel_Time_Freq_Odor);
        cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
        cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
        cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
        cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
        
        stat_OdorvsVehicle      = ft_timelockstatistics(cfg, ...
            Sel_Time_Freq_Odor{:}, Sel_Time_Freq_Vehi{:});   % don't forget the {:}!
        
        
        figure
        confidence = 1-stat_OdorvsVehicle.prob;
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
        idx_clusters = find((stat_OdorvsVehicle.prob <= 0.05));
        topoplot(no_results, reducedchanlocs, ...
            'style', 'blank', ...
            'electrodes', 'pts', ...
            'shading', 'interp', ...
            'headcolor', [0, 0, 0], ...
            'plotchans', idx_clusters, ...
            'hcolor','none',...
            'whitebk','on',...
            'emarker', {'.', [0.1, 0.1, 0.1], 20, 1});
        
        Text = strcat(bands(band),'[',num2str(latency),']','s');
        title({'Odor D vs Vehicle';Text{1}})
        
        filename_save = strcat('DNight_',bands(band),'Toi',num2str(toi));
        saveas(gcf,strcat(filename_save{1},'.png'))
        
        hold off
        
        close all
        
    end
    
end
