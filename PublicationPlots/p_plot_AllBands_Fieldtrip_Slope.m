% addpath('/home/andrea/Documents/Github/TimeFrequency_ReactivatedMemoryProject/')

addpath('C:\Users\asanch24\Documents\MATLAB\fieldtrip-20190828\') %Windows
% addpath('/home/andrea/Documents/MatlabFunctions/fieldtrip-20200828') %server

ft_defaults

DNight = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_DNight_TFSlope.mat');
MNight = load('C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\TF_Slope\TF_MNight_TFSlope.mat');

%%

bands = {'SW';'Delta';'Theta';'Spindle'};

SpindleBand     = [12 16];
DeltaBand       = [1 4];
ThetaBand       = [4 8];
SWBand          = [0.5 2];

% load('D:\Thesis_Publication\Using\Time_Freq_Parms')
Time_Series_OdorD =[];

%% Set electrode layout

cfg = [];
cfg.layout = 'C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject\GSN-HydroCel-128.mat';
layout = ft_prepare_layout(cfg);


%% customized colormap




%%
%--------------------------------------------------------------------------
% Analysis by bands
%--------------------------------------------------------------------------
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
    
%     cfg.avgoverrpt      = 'yes';
    cfg.avgoverfreq     = 'yes';
    cfg.channel         = {'all','-E49', '-E48', '-E17', '-E128', '-E32', ...
        '-E1', '-E125', '-E119', '-E113', '-E56', ...
        '-E63', '-E68', '-E73', '-E81', '-E88', ...
        '-E94', '-E99', '-E107'};
    
    % Time_Series_OdorD.(bands{band}) = Time_Freq_Parms;
    
    for subj = 1:numel(DNight.Time_Freq_Odor_Slope)
        
        disp(strcat('Sujeto: ',num2str(subj)))
        
        % For D Night
        SlopeOdor = DNight.Time_Freq_Odor_Slope{subj};
        SlopeVehicle = DNight.Time_Freq_Vehi_Slope{subj};
        
        Sel_Freq_OdorD = ft_selectdata(cfg, SlopeOdor);
        Sel_Freq_VehicleD = ft_selectdata(cfg, SlopeVehicle);
        
        OdorD_trial{1,subj} = squeeze(Sel_Freq_OdorD.powspctrm);
        VehicleD_trial{1,subj} = squeeze(Sel_Freq_VehicleD.powspctrm);

        
        % For M Night
        SlopeOdor = MNight.Time_Freq_Odor_Slope{subj};
        SlopeVehicle = MNight.Time_Freq_Vehi_Slope{subj};

        
        Sel_Freq_OdorM = ft_selectdata(cfg, SlopeOdor);
        Sel_Freq_VehicleM = ft_selectdata(cfg, SlopeVehicle);
        
        OdorM_trial{1,subj} = squeeze(Sel_Freq_OdorM.powspctrm);
        VehicleM_trial{1,subj} = squeeze(Sel_Freq_VehicleM.powspctrm);
        
        time_all{1,subj} = Sel_Freq_OdorM.time;
    end
    
    Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD,'freq');
%     Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD, 'trialinfo');
    Sel_Freq_OdorD = rmfield(Sel_Freq_OdorD,'powspctrm');
    
    Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD,'freq');
%     Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD, 'trialinfo');
    Sel_Freq_VehicleD = rmfield(Sel_Freq_VehicleD,'powspctrm');
    
    Sel_Freq_OdorD.dimord = 'chan_time';
    Sel_Freq_VehicleD.dimord = 'chan_time';
    
    Sel_Freq_OdorD.time      = time_all;
    Sel_Freq_VehicleD.time   = time_all;
    
    Time_Series_OdorD.(bands{band}) = Sel_Freq_OdorD;
    Time_Series_VehicleD.(bands{band}) = Sel_Freq_VehicleD;
    
    Time_Series_OdorD.(bands{band}).trial = OdorD_trial;
    Time_Series_VehicleD.(bands{band}).trial = VehicleD_trial;

    %--- M Night---
    
    Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM,'freq');
%     Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM, 'trialinfo');
    Sel_Freq_OdorM = rmfield(Sel_Freq_OdorM,'powspctrm');
    
    Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM,'freq');
%     Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM, 'trialinfo');
    Sel_Freq_VehicleM = rmfield(Sel_Freq_VehicleM,'powspctrm');
    
    Sel_Freq_OdorM.dimord = 'chan_time';
    Sel_Freq_VehicleM.dimord = 'chan_time';
    
    Sel_Freq_OdorM.time      = time_all;
    Sel_Freq_VehicleM.time   = time_all;
    
    Time_Series_OdorM.(bands{band}) = Sel_Freq_OdorM;
    Time_Series_VehicleM.(bands{band}) = Sel_Freq_VehicleM;
    
    Time_Series_OdorM.(bands{band}).trial = OdorM_trial;
    Time_Series_VehicleM.(bands{band}).trial = VehicleM_trial;
     
end
%%
filename     = 'C:\Users\asanch24\OneDrive - St. Jude Children''s Research Hospital\Thesis_Publication\Using\dummyfile.set'; %windows
%filename     = '/mnt/disk1/andrea/German_Study/Time_Frequency_FT/dummyfile.set'; %server
sensors        = ft_read_sens(filename);

Nsub = 23;



%%
for band = 1:numel(bands)
    
    cfg                  = [];
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.correctm         = 'cluster'; % correction
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
    cfg.minnbchan        = 2;
    cfg_neighb.elec      = 'GSN-HydroCel-128.sfp';%sensors;
    cfg_neighb.method    = 'distance';
    cfg.neighbours       = ft_prepare_neighbours(cfg_neighb);
    cfg.tail             = 0;
    cfg.clustertail      = 0;
    cfg.alpha            = 0.025;
    cfg.numrandomization = 10000;
    
    
    Nsubj  = Nsub;
    design = zeros(2, Nsubj*2);
    design(1,:) = [1:Nsubj 1:Nsubj];
    design(2,:) = [ones(1,Nsubj) ones(1,Nsubj)*2];
    
    cfg.design = [ones(1,Nsubj) ones(1,Nsubj)*2];
    
    cfg.design = design;
    cfg.uvar   = 1;
    cfg.ivar   = 2;
    
    stat_OdorD_VehicleD = ft_timelockstatistics(cfg,Time_Series_OdorD.(bands{band}),...
        Time_Series_VehicleD.(bands{band}));
    
    stat_OdorM_VehicleM = ft_timelockstatistics(cfg,Time_Series_OdorM.(bands{band}),...
        Time_Series_VehicleM.(bands{band}));
    
    stat_OdorD_OdorM = ft_timelockstatistics(cfg,Time_Series_OdorD.(bands{band}),...
        Time_Series_OdorM.(bands{band}));
    
    band_stat_OdorD_VehicleD.(bands{band})=stat_OdorD_VehicleD;
    band_stat_OdorM_VehicleM.(bands{band})=stat_OdorM_VehicleM;
    band_stat_OdorD_OdorM.(bands{band})=stat_OdorD_OdorM;
    
    % Take the difference using ft_math subject per subject (data is organized
    % such that each subjects is a trial in the data)
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'trial';
    Diff_OdorDvsVehicleD.(bands{band}) =...
        ft_math(cfg, Time_Series_OdorD.(bands{band}), Time_Series_VehicleD.(bands{band}));
    Diff_OdorMvsVehicleM.(bands{band}) =...
        ft_math(cfg, Time_Series_OdorM.(bands{band}), Time_Series_VehicleM.(bands{band}));
    Diff_OdorDvsOdorM.(bands{band}) =...
        ft_math(cfg, Time_Series_OdorD.(bands{band}), Time_Series_OdorM.(bands{band}));
    
    % Do average for all subjects subject (each subject is organized as a
    % trial)
    cfg                     = [];
    cfg.avgoverrpt          = 'yes';
    Diff_OdorDvsVehicleD.(bands{band}) = ...
        ft_selectdata(cfg, Diff_OdorDvsVehicleD.(bands{band}));
    
    Diff_OdorMvsVehicleM.(bands{band}) = ...
        ft_selectdata(cfg, Diff_OdorMvsVehicleM.(bands{band}));
    
    Diff_OdorDvsOdorM.(bands{band}) = ...
        ft_selectdata(cfg, Diff_OdorDvsOdorM.(bands{band}));
    
    
    
    % Make a vector of all p-values associated with the clusters from ft_timelockstatistics.
    pos_cluster_pvals_OdorD_VehicleD = [stat_OdorD_VehicleD.posclusters(:).prob];
    pos_cluster_pvals_OdorM_VehicleM = [stat_OdorM_VehicleM.posclusters(:).prob];
    pos_cluster_pvals_OdorD_OdorM    = [stat_OdorD_OdorM.posclusters(:).prob];
    
    % Then, find which clusters are deemed interesting to visualize, here we use a cutoff criterion based on the
    % cluster-associated p-value, and take a 5% two-sided cutoff (i.e. 0.025 for the positive and negative clusters,
    % respectively
    pos_clust_OdorD_VehicleD = find(pos_cluster_pvals_OdorD_VehicleD < 0.025);
    pos_OdorD_VehicleD       = ismember(stat_OdorD_VehicleD.posclusterslabelmat, pos_clust_OdorD_VehicleD);
    
    pos_clust_OdorM_VehicleM = find(pos_cluster_pvals_OdorM_VehicleM < 0.025);
    pos_OdorM_VehicleM       = ismember(stat_OdorM_VehicleM.posclusterslabelmat, pos_clust_OdorM_VehicleM);
    
    pos_clust_OdorD_OdorM = find(pos_cluster_pvals_OdorD_OdorM < 0.025);
    pos_OdorD_OdorM       = ismember(stat_OdorD_OdorM.posclusterslabelmat, pos_clust_OdorD_OdorM);
    
    % and now for the negative clusters...
    neg_cluster_pvals_OdorD_VehicleD = [stat_OdorD_VehicleD.negclusters(:).prob];
    neg_clust_OdorD_VehicleD         = find(neg_cluster_pvals_OdorD_VehicleD < 0.025);
    neg_OdorD_VehicleD               = ismember(stat_OdorD_VehicleD.negclusterslabelmat, neg_clust_OdorD_VehicleD);
    
    neg_cluster_pvals_OdorM_VehicleM = [stat_OdorM_VehicleM.negclusters(:).prob];
    neg_clust_OdorM_VehicleM         = find(neg_cluster_pvals_OdorM_VehicleM < 0.025);
    neg_OdorM_VehicleM               = ismember(stat_OdorM_VehicleM.negclusterslabelmat, neg_clust_OdorM_VehicleM);
    
    neg_cluster_pvals_OdorD_OdorM  = [stat_OdorD_OdorM.negclusters(:).prob];
    neg_clust_OdorD_OdorM          = find(neg_cluster_pvals_OdorD_OdorM < 0.025);
    neg_OdorD_OdorM                = ismember(stat_OdorD_OdorM.negclusterslabelmat, neg_clust_OdorD_OdorM );
    
    % If we only want to see the extent of the first (i.e. most significant)
    % positive and negative clusters:
    
    pos_OdorD_VehicleD = stat_OdorD_VehicleD.posclusterslabelmat == 1; % or == 2, or 3, etc.
    neg_OdorD_VehicleD = stat_OdorD_VehicleD.negclusterslabelmat == 1;
    
    pos_OdorM_VehicleM = stat_OdorM_VehicleM.posclusterslabelmat == 1; % or == 2, or 3, etc.
    neg_OdorM_VehicleM = stat_OdorM_VehicleM.negclusterslabelmat == 1;
    
    pos_OdorD_OdorM = stat_OdorD_OdorM.posclusterslabelmat == 1; % or == 2, or 3, etc.
    neg_OdorD_OdorM = stat_OdorD_OdorM.negclusterslabelmat == 1;
    
    
    %----------------------------------
    % plot for OdorD_VehicleD
    %----------------------------------
    
    timestep      = 0.5; % timestep between time windows for each subplot (in seconds)
    sampling_rate = 10; % Data has a temporal resolution of 300 Hz
    sample_count  = length(stat_OdorD_VehicleD.time);
    % number of temporal samples in the statistics object
    j = -5:timestep:5; % Temporal endpoints (in seconds) of the ERP average computed in each subplot
    m = 1:timestep*sampling_rate:sample_count; % temporal endpoints in M/EEG samples


    
    % First ensure the channels to have the same order in the average and in the statistical output.
    % This might not be the case, because ft_math might shuffle the order
    [i1,i2] = match_str(Time_Series_OdorD.(bands{band}).label, stat_OdorD_VehicleD.label);
    
    for k = 1:numel(j)-1
        %subplot(4,8,k);
        subplot(3,7,k);
%         subplot(6,9,k);
        
        cfg = [];
        cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
%         cfg.zlim = [-2.5e-13 2.5e-13];
        % If a channel is in a to-be-plotted cluster, then
        % the element of pos_int with an index equal to that channel
        % number will be set to 1 (otherwise 0).
        
        % Next, check which channels are in the clusters over the
        % entire time interval of interest.
        pos_int = zeros(numel(Time_Series_OdorD.(bands{band}).label),1);
        neg_int = zeros(numel(Time_Series_OdorD.(bands{band}).label),1);
        pos_int(i1) = all(pos_OdorD_VehicleD(i2, m(k):m(k+1)), 2);
        neg_int(i1) = all(neg_OdorD_VehicleD(i2, m(k):m(k+1)), 2);
        
        if find(pos_int==1)
            cfg.highlight   = 'on';
            % Get the index of the to-be-highlighted channel
            cfg.highlightchannel = find(pos_int | neg_int);
        end
        
        data = Diff_OdorDvsVehicleD.(bands{band});
        data.avg = data.trial;
        data.mask               = stat_OdorD_VehicleD.mask;
        
        cfg.zlim            = [-max(max(abs(data.avg{1}))),...
            max(max(abs(data.avg{1})))];
        cfg.parameter       = 'avg';
        cfg.comment         = 'xlim';
        cfg.commentpos      = 'title';
        cfg.layout          = layout;
        cfg.colormap        = 'parula';
%         cfg.maskparameter   = 'mask';
%         cfg.colorbar        = 'yes';
        cfg.style           = 'straight';
        

        ft_topoplotER(cfg, data);
    end
    colorbar
    set(gcf,'position',[1,41,1920,963])
    filename_save = strcat('C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject\PublicationPlots\DNight_',bands{band});
    saveas(gcf,strcat(filename_save,'.png'))
    
    
    %----------------------------------
    % plot for OdorM_VehicleM
    %----------------------------------
    
    sampling_rate = 10; % Data has a temporal resolution of 300 Hz
    sample_count  = length(stat_OdorM_VehicleM.time);
    % number of temporal samples in the statistics object
    m = 1:timestep*sampling_rate:sample_count; % temporal endpoints in M/EEG samples


    
    % First ensure the channels to have the same order in the average and in the statistical output.
    % This might not be the case, because ft_math might shuffle the order
    [i1,i2] = match_str(Time_Series_OdorM.(bands{band}).label, stat_OdorM_VehicleM.label);
    
    for k = 1:numel(j)-1
        %subplot(4,8,k);
        subplot(3,7,k);
%         subplot(6,9,k);
        
        cfg = [];
        cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
%         cfg.zlim = [-2.5e-13 2.5e-13];
        % If a channel is in a to-be-plotted cluster, then
        % the element of pos_int with an index equal to that channel
        % number will be set to 1 (otherwise 0).
        
        % Next, check which channels are in the clusters over the
        % entire time interval of interest.
        pos_int = zeros(numel(Time_Series_OdorM.(bands{band}).label),1);
        neg_int = zeros(numel(Time_Series_OdorM.(bands{band}).label),1);
        pos_int(i1) = all(pos_OdorM_VehicleM(i2, m(k):m(k+1)), 2);
        neg_int(i1) = all(neg_OdorM_VehicleM(i2, m(k):m(k+1)), 2);
        
        if find(pos_int==1)
            cfg.highlight   = 'on';
            % Get the index of the to-be-highlighted channel
            cfg.highlightchannel = find(pos_int | neg_int);
        end
        
        data = Diff_OdorMvsVehicleM.(bands{band});
        data.avg = data.trial;
        data.mask               = stat_OdorM_VehicleM.mask;
        
         cfg.zlim            = [-max(max(abs(data.avg{1}))),...
            max(max(abs(data.avg{1})))];
        cfg.parameter       = 'avg';
        cfg.comment         = 'xlim';
        cfg.commentpos      = 'title';
        cfg.layout          = layout;
        cfg.colormap        = 'parula';
%         cfg.maskparameter   = 'mask';
%         cfg.colorbar        = 'yes';
        cfg.style           = 'straight';
        

        
        ft_topoplotER(cfg, data);
    end
    colorbar
    set(gcf,'position',[1,41,1920,963])
    filename_save = strcat('C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject\PublicationPlots\MNight_',bands{band});
    saveas(gcf,strcat(filename_save,'.png'))
    
    %----------------------------------
    % plot for OdorD_OdorM
    %----------------------------------
    
    sampling_rate = 10; % Data has a temporal resolution of 300 Hz
    sample_count  = length(stat_OdorD_OdorM.time);
    % number of temporal samples in the statistics object
    m = 1:timestep*sampling_rate:sample_count; % temporal endpoints in M/EEG samples


    
    % First ensure the channels to have the same order in the average and in the statistical output.
    % This might not be the case, because ft_math might shuffle the order
    [i1,i2] = match_str(Time_Series_OdorM.(bands{band}).label, stat_OdorD_OdorM.label);
    
    for k = 1:numel(j)-1
        %subplot(4,8,k);
        subplot(3,7,k);
%         subplot(6,9,k);
        
        cfg = [];
        cfg.xlim = [j(k) j(k+1)];   % time interval of the subplot
        % If a channel is in a to-be-plotted cluster, then
        % the element of pos_int with an index equal to that channel
        % number will be set to 1 (otherwise 0).
        
        % Next, check which channels are in the clusters over the
        % entire time interval of interest.
        pos_int = zeros(numel(Time_Series_OdorM.(bands{band}).label),1);
        neg_int = zeros(numel(Time_Series_OdorM.(bands{band}).label),1);
        pos_int(i1) = all(pos_OdorD_OdorM(i2, m(k):m(k+1)), 2);
        neg_int(i1) = all(neg_OdorD_OdorM(i2, m(k):m(k+1)), 2);
        
        if sum((pos_int==1) + (neg_int==1))>0
            cfg.highlight   = 'on';
            % Get the index of the to-be-highlighted channel
            cfg.highlightchannel = find(pos_int | neg_int);
        end
        
        data = Diff_OdorDvsOdorM.(bands{band});
        data.avg                = data.trial;
        data.mask               = stat_OdorD_OdorM.mask;
        
         cfg.zlim            = [-max(max(abs(data.avg{1}))),...
            max(max(abs(data.avg{1})))];
        cfg.parameter       = 'avg';
        cfg.comment         = 'xlim';
        cfg.commentpos      = 'title';
        cfg.layout          = layout;
        cfg.colormap        = 'parula';
%         cfg.maskparameter   = 'mask';
%         cfg.colorbar        = 'yes';
%         cfg.style           = 'straight';
        
        ft_topoplotER(cfg, data);
      
    end
    colorbar
    set(gcf,'position',[1,41,1920,963])
    filename_save = strcat('C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject\PublicationPlots\DvsM_',bands{band});
    saveas(gcf,strcat(filename_save,'.png'))
end


%%

dirlist  = dir('C:\Users\asanch24\Documents\Github\TimeFrequency_ReactivatedMemoryProject\GSN-HydroCel-128.mat');
filename = {dirlist(~[dirlist.isdir]).name}';

for i=1:length(filename)
  cfg = [];
  cfg.layout = filename{i};
  layout = ft_prepare_layout(cfg);

  figure
  ft_plot_layout(layout);
  title(filename{i}, 'Interpreter', 'none');

  [p, f, x] = fileparts(filename{i});
  print([lower(f) '.png'], '-dpng');
end