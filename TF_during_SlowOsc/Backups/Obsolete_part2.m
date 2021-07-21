% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   P A R T 2   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% Here, we average channels' TFs for pre-defined clusters and compute
% permutation-based analysis on each of those and extract clusters that are
% found.


%% TF matrix statistics and wave forms for pre-defined clusters
%  ------------------------------------------------------------------------

% Clusters of interest
% --------------------
PM.Clust.left_frontal = {...
    'E15', 'E16', 'E11', 'E18', 'E19', 'E22', 'E23', 'E24', 'E26', ...
    'E27', 'E33', 'E38'};
PM.Clust.right_frontal = {...
    'E15', 'E16', 'E11', 'E10', 'E4', 'E9', 'E3', 'E124', 'E2', ...
    'E123', 'E122', 'E121'};
PM.Clust.frontal = {...
    'E3', 'E4', 'E9', 'E10', 'E11', 'E15', 'E16', 'E18', 'E19', ...
    'E22', 'E23', 'E24', 'E124'};
PM.Clust.left_central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55'};
PM.Clust.right_central = {...
    'E6', 'E55', 'E112', 'E106', 'E105', 'E80', 'E87', 'E79'};
PM.Clust.central = {...
    'E6', 'E7', 'E13', 'E30', 'E31', 'E37', 'E54', 'E55', 'E79', ...
    'E80', 'E87', 'E105', 'E106', 'E112'};
PM.Clust.left_temporal = {...
    'E46', 'E51', 'E45', 'E50', 'E58', 'E56', 'E63'};
PM.Clust.right_temporal = {...
    'E108', 'E102', 'E101', 'E97', 'E96', 'E99', 'E107'};
PM.Clust.left_parietal = {...
    'E53', 'E61', 'E62', 'E72', 'E67', 'E52', 'E60', 'E59', 'E66', ...
    'E65', 'E64', 'E68'};
PM.Clust.right_parietal = {...
    'E62', 'E72', 'E78', 'E77', 'E86', 'E85', 'E84', 'E92', 'E91', ...
    'E90', 'E95', 'E94'};
PM.Clust.parietal = {...
    'E52', 'E61', 'E62', 'E59', 'E60', 'E67', 'E66', 'E72', 'E78', ...
    'E77', 'E86', 'E85', 'E84', 'E92', 'E91', 'E53'};
PM.Clust.left_occipital = {...
    'E71', 'E70', 'E75', 'E74', 'E69', 'E73'};
PM.Clust.right_occipital = {...
    'E75', 'E76', 'E82', 'E83', 'E88', 'E89'};
PM.Clust.occipital = {...
    'E71', 'E70', 'E74', 'E69', 'E73', 'E75', 'E76', 'E83', 'E82', ...
    'E89', 'E88'};

Clusters = fieldnames(PM.Clust);

for i_defined_clust = 1:numel(Clusters)
    
    
    % Prepare data structures and plot parameters
    % -------------------------------------------
    str_cluster     = char(Clusters(i_defined_clust));
    idx_defchans    = ...
        find(ismember(dummyFile.PM.Clust.all, PM.Clust.(str_cluster)));
    
    clust_sensors   = sensors;
    clust_sensors.chanpos = ...
        clust_sensors.chanpos(idx_defchans, :);
    clust_sensors.chantype = ...
        clust_sensors.chantype(idx_defchans);
    clust_sensors.chanunit = ...
        clust_sensors.chanunit(idx_defchans);
    clust_sensors.elecpos  = ...
        clust_sensors.elecpos(idx_defchans, :);
    clust_sensors.label    = ...
        clust_sensors.label(idx_defchans);
    
    % We generate to options for later: Keep channel dimension by only
    % taking into account channels of cluster or average channels in order
    % to remove the channel dimension
    % ---------------------------------------------------------------------
    TF_Chans.ShamOn              = TF_subject.ShamOn;
    TF_Chans.ShamOn.label        = PM.Clust.(str_cluster);
    TF_Chans.ShamOn.dimord       = 'rpt_chan_freq_time';
    TF_Chans.ShamOn.powspctrm    = ...
        TF_subject.ShamOn.powspctrm(:, idx_defchans, :, :);
    TF_Chans.ShamOn.elec.chanpos = ...
        TF_Chans.ShamOn.elec.chanpos(idx_defchans, :);
    TF_Chans.ShamOn.elec.chantype = ...
        TF_Chans.ShamOn.elec.chantype(idx_defchans);
    TF_Chans.ShamOn.elec.chanunit = ...
        TF_Chans.ShamOn.elec.chanunit(idx_defchans);
    TF_Chans.ShamOn.elec.elecpos  = ...
        TF_Chans.ShamOn.elec.elecpos(idx_defchans, :);
    TF_Chans.ShamOn.elec.label    = ...
        TF_Chans.ShamOn.elec.label(idx_defchans);
    
    TF_Chans.OdorOn              = TF_subject.OdorOn;
    TF_Chans.OdorOn.label        = PM.Clust.(str_cluster);
    TF_Chans.OdorOn.dimord       = 'rpt_chan_freq_time';
    TF_Chans.OdorOn.powspctrm    = ...
        TF_subject.OdorOn.powspctrm(:, idx_defchans, :, :);
    TF_Chans.OdorOn.elec.chanpos = ...
        TF_Chans.OdorOn.elec.chanpos(idx_defchans, :);
    TF_Chans.OdorOn.elec.chantype = ...
        TF_Chans.OdorOn.elec.chantype(idx_defchans);
    TF_Chans.OdorOn.elec.chanunit = ...
        TF_Chans.OdorOn.elec.chanunit(idx_defchans);
    TF_Chans.OdorOn.elec.elecpos  = ...
        TF_Chans.OdorOn.elec.elecpos(idx_defchans, :);
    TF_Chans.OdorOn.elec.label    = ...
        TF_Chans.OdorOn.elec.label(idx_defchans);
    
    
    TF_meanChans.ShamOn           = TF_subject.ShamOn;
    TF_meanChans.ShamOn           = rmfield(TF_meanChans.ShamOn, 'elec');
    TF_meanChans.ShamOn.label     = {str_cluster};
    TF_meanChans.ShamOn.dimord    = 'rpt_freq_time';
    TF_meanChans.ShamOn.powspctrm = squeeze(mean(...
        TF_subject.ShamOn.powspctrm(:, idx_defchans, :, :), 2));
    TF_meanChans.OdorOn           = TF_subject.OdorOn;
    TF_meanChans.OdorOn           = rmfield(TF_meanChans.OdorOn, 'elec');
    TF_meanChans.OdorOn.label     = {str_cluster};
    TF_meanChans.OdorOn.dimord    = 'rpt_freq_time';
    TF_meanChans.OdorOn.powspctrm = squeeze(mean(...
        TF_subject.OdorOn.powspctrm(:, idx_defchans, :, :), 2));
    
    WF_meanChans.ShamOn          = WF_subject.ShamOn(:, idx_defchans, :);
    WF_meanChans.OdorOn          = WF_subject.OdorOn(:, idx_defchans, :);
    
    
    
    %% Permutation-based statistics for cluster
    %  --------------------------------------------------------------------
    
    addpath(fieldtrippath)
    ft_defaults
    
    cfg_stats                   = [];
    cfg_stats.frequency         = 'all';
    cfg_stats.channel           = 'all';
    cfg_stats.correctm          = 'cluster';
    cfg_stats.method            = 'montecarlo';
    cfg_stats.statistic         = 'depsamplesT';
    % use actvsblT for activation against baseline
    cfg_stats.clusterstatistic  = 'maxsum';
    cfg_stats.minnbchan         = 1;
    cfg_stats.tail              = 0; % -/+1 or 0 for one-sided or two-sided
    cfg_stats.clustertail       = cfg_stats.tail;
    % cfg_stats.alpha                 = 0.1;
    cfg_stats.alpha             = 0.05; % = False detection rate ?
    % Commented out in order to get p vals of channels that would otherwise
    % be set to p = 1 in case the p value is higher than the alpha here.
    % Cluster perm stats --> Which clusters are "significant"
    cfg_stats.correcttail       = 'alpha'; % {'prob', 'alpha'}
    % Adapts alpha value corresponding to one-tailed or two-tailed test
    cfg_stats.numrandomization  = 10000;
    % cfg_stats.clusteralpha          = 0.1;
    cfg_stats.clusteralpha      = 0.05;  % "which points for clustering"
    % aka include more or less channels, freqs and times
    % threshold over which a triplet is chosen
    % Design needs to be adapted to cfg.uvar and cfg.ivar
    cfg_stats.uvar              = 1;
    % condition (uvar would be the subjects)
    % Iactivate for cfg.statistic = 'indepsamplesT'
    cfg_stats.ivar              = 2;
    cfg_stats.design            = [1:numel(c_kept_subj), ...
                                    1:numel(c_kept_subj); ...
                                    repmat(1, 1, numel(c_kept_subj)), ...
                                    repmat(2, 1, numel(c_kept_subj))];
    
    
    % Depending on options: channel dimension kept or averaged
%     cfg_stats.neighbours        = [];
%     stats  = ft_freqstatistics(cfg_stats, ...
%         TF_meanChans.OdorOn, TF_meanChans.ShamOn);

    cfg_neighb.method               = 'distance'; %'triangulation';
    cfg_neighb.neighbourdist        = 3.5; % Default is "smart guess".
    % 3.5 gets rid of edge clusters while maintaining high statistical
    % significance in central clusters.
    cfg_neighb.channel              = 'all';
    cfg_neighb.elec                 = clust_sensors;
    cfg_stats.neighbours            = ft_prepare_neighbours(cfg_neighb);
    stats  = ft_freqstatistics(cfg_stats, ...
        TF_Chans.OdorOn, TF_Chans.ShamOn);
    
    allstats.(str_cluster) = stats;
    
    
    
    %% Retrieve statistical information
    %  --------------------------------------------------------------------
    
    dimensions      = size(stats.prob);
    cluster_labels  = zeros(dimensions);
    cluster_pvals   = [];
    
    % ---------------------------------------------------------------------
    %                           /!\ /!\ /!\
    % Here, better circumvent the use of stats.mask since we will be
    % limited for cluster labels according to set alpha during permutation
    % analysis. Only cluster with p values lower than set alpha would be
    % represented in mask. We get more freedom with pos/negclusterlabelmat.
    % ---------------------------------------------------------------------
    
    % Positive cluster
    % cluster_labels(stats.mask)  = stats.posclusterslabelmat(stats.mask);
    cluster_labels(stats.posclusterslabelmat ~= 0) = ...
        stats.posclusterslabelmat(stats.posclusterslabelmat ~= 0);
    
    % Extract p values of positive clusters: Labels indicate position in
    % stats.posclusters.prob array
    idx_labels_pos                      = unique(cluster_labels(:));
    idx_labels_pos(idx_labels_pos == 0) = []; % zeros are "no cluster"
    cluster_pvals                       = [cluster_pvals, ...
        stats.posclusters(idx_labels_pos).prob];
    
    % Negative cluster
    s_raise_label               = max(cluster_labels(:));
    matrix_raise                = repmat(s_raise_label, dimensions);
    stats.negclusterslabelmat   = stats.negclusterslabelmat + matrix_raise;
    % cluster_labels(stats.mask)  = stats.negclusterslabelmat(stats.mask);
    cluster_labels(stats.negclusterslabelmat ~= matrix_raise) = ...
        stats.negclusterslabelmat(stats.negclusterslabelmat ~= matrix_raise);
    
    % Extract p values of negative clusters
    idx_labels_neg              = ...
        unique(cluster_labels(cluster_labels > s_raise_label));
    
    idx_labels_neg              = idx_labels_neg - s_raise_label;
    % Point to correct position in neg cluster structure
    
    cluster_pvals               = [cluster_pvals, ...
        stats.negclusters(idx_labels_neg).prob];
    
    
    
    %% Plot time-frequency clusters of channels belonging to them
    %  --------------------------------------------------------------------
    
    significant_clusters = unique(cluster_labels(:))';
    significant_clusters(significant_clusters == 0) = [];
    
    if numel(cluster_pvals) ~= numel(significant_clusters)
        % Checkpoint
        error('Something went wrong')
    end
    
    for i_clust = significant_clusters
        
        
        if cluster_pvals(significant_clusters(i_clust)) >= s_thrld
            % Not "significant enough"
            continue
        end
        
        
        true_labels = NaN(size(cluster_labels));
        true_labels(cluster_labels == i_clust) = 1;
        
        % Find channels belonging to cluster
        % ----------------------------------
        % Channels belonging to clusters should be converted to 1, the rest
        % is NaN
        idx_clustchans = nanmean(nanmean(true_labels, 3), 2);
        idx_clustchans = find(idx_clustchans == 1);
        
        
        TF_Odor = TF_subject.OdorOn.powspctrm;
        TF_Odor = mean(TF_Odor, 1);
        TF_Odor = mean(TF_Odor, 2);
        TF_Odor = squeeze(TF_Odor);
        
        TF_Sham = TF_subject.ShamOn.powspctrm;
        TF_Sham = mean(TF_Sham, 1);
        TF_Sham = mean(TF_Sham, 2);
        TF_Sham = squeeze(TF_Sham);
        
        minVal = min([TF_Odor(:); TF_Sham(:)]);
        maxVal = max([TF_Odor(:); TF_Sham(:)]);
        
        limitsTF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        average_WF.ShamOn = squeeze(mean(mean(WF_meanChans.ShamOn, 1), 2));
        average_WF.OdorOn = squeeze(mean(mean(WF_meanChans.OdorOn, 1), 2));
        
        minVal      = min([average_WF.ShamOn(:); average_WF.OdorOn(:)]);
        maxVal      = max([average_WF.ShamOn(:); average_WF.OdorOn(:)]);
        
        limitsWF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        
        figure('units','normalized','outerposition', [0 0 1 0.3]);
        
        % Plot all channels as small dots and overlay bigger size
        % electrodes for channels belonging to cluster
        % -------------------------------------------------------
        
        subplot(1, 5, 1)
        no_results = zeros(numel(chanlocs), 1);
        topoplot(no_results, chanlocs, ...
            'style', 'blank', ...
            'electrodes', 'pts', ...
            'shading', 'interp', ...
            'headcolor', [0, 0, 0], ...
            'plotchans', idx_defchans, ...
            'conv', 'on', ...
            'emarker', {'.', [0, 0, 0], 5, 1}, ...
            'plotrad', max(max([chanlocs.radius]),0.5));
        hold on
        topoplot(no_results, chanlocs, ...
            'style', 'blank', ...
            'electrodes', 'pts', ...
            'shading', 'interp', ...
            'headcolor', [0, 0, 0], ...
            'conv', 'off', ...
            'plotchans', idx_defchans(idx_clustchans), ...
            'emarker', {'.', [0.1, 0.1, 0.1], 20, 1}, ...
            'plotrad', max(max([chanlocs.radius]), 0.5));
        hold off
        
        
        colororder({'k', 'r'})
        i_plot = 2;
        for condition = PM.Conditions
            
            subplot(1, 5, i_plot)
            
            % Time frequency plot
            % -------------------
            
            % There must be a more efficient way to determine v_frequencies ...
            v_frequencies = 1:length(TF_series.OdorOn.freq);
            v_frequencies = (max(v_frequencies) - v_frequencies) / ...
                (max(v_frequencies) - min(v_frequencies));
            v_frequencies = v_frequencies * PM.FrRange;
            v_frequencies = v_frequencies - (PM.FrRange / 2);
            v_frequencies = flip(v_frequencies);
            
            v_times       = ...
                PM.cfg_seldat.latency(1) * ...
                PM.Info.TrialParameters.s_fs: 1/PM.s_tstep/2 :...
                PM.cfg_seldat.latency(2) * ...
                PM.Info.TrialParameters.s_fs;
            v_times       = v_times ./ PM.Info.TrialParameters.s_fs;
            
            % Chose channel dimension kept or averaged
%             TF_meanSubj.(char(condition)) = ...
%                 squeeze(mean(TF_meanChans.(char(condition)).powspctrm, 1));
            TF_meanSubj.(char(condition)) = ...
                squeeze(mean(TF_Chans.(char(condition)).powspctrm, 1));
            
            yyaxis left
            pcolor(v_times, v_frequencies, ...
                squeeze(mean(TF_meanSubj.(...
                char(condition))(idx_clustchans, :, :), 1)));
            shading interp
            colorbar;
            set(gca, 'clim', limitsTF)
            ylabel('Distance from spindle peak (Hz)')
            xlabel('Time (s)')
            title(char(condition))
            
            hold on
            
            % Wave form plot
            % --------------
            yyaxis right
            v_times = PM.cfg_seldat.latency(1) * ...
                PM.Info.TrialParameters.s_fs  + 1 : 1 : ...
                PM.cfg_seldat.latency(2) * ...
                PM.Info.TrialParameters.s_fs;
            v_times = v_times ./ PM.Info.TrialParameters.s_fs;
            plot(v_times, average_WF.(char(condition)), ...
                'Color',        [0, 0, 0], ...
                'LineWidth',    2, ...
                'Color',        'r')
            ylim(limitsWF)
            
            i_plot = i_plot + 1;
        end
        
        % Overlay difference in TF
        TF_difference = TF_meanSubj.OdorOn - TF_meanSubj.ShamOn;
        
        TF_cluster = squeeze(mean(TF_difference(idx_clustchans, :, :), 1));
        
        minVal = min(min(TF_cluster));
        maxVal = max(max(TF_cluster));
        
        limitsTF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        subplot(1, 5, i_plot)
        v_times       = ...
            PM.cfg_seldat.latency(1) * ...
            PM.Info.TrialParameters.s_fs: 1/PM.s_tstep/2 :...
            PM.cfg_seldat.latency(2) * ...
            PM.Info.TrialParameters.s_fs;
        v_times       = v_times ./ PM.Info.TrialParameters.s_fs;
        pcolor(v_times, v_frequencies, TF_cluster);
        shading interp
        colorbar;
        set(gca, 'clim', limitsTF)
        ylabel('Distance from spindle peak (Hz)')
        xlabel('Time (s)')
        title('Odor - Sham')
        colormap(parula)
        
        
        subplot(1, 5, 5)
        v_times = PM.cfg_seldat.latency(1) * ...
            PM.Info.TrialParameters.s_fs : 1/PM.s_tstep/2 : ...
            PM.cfg_seldat.latency(2) * ...
            PM.Info.TrialParameters.s_fs;
        v_times = v_times ./ PM.Info.TrialParameters.s_fs;
        
        TF_selective    = NaN(dimensions);
        %        TF_selective(cluster_pvals == p_clusters(i_clust)) = ...
        %          TF_difference(cluster_pvals == p_clusters(i_clust));
        TF_selective(true_labels == 1) = ...
            TF_difference(true_labels == 1);
        
        figSign = pcolor(v_times, v_frequencies, ...
            squeeze(nanmean(TF_selective(idx_clustchans, :, :), 1)));
        shading interp
        colorbar;
        set(gca, 'clim', limitsTF)
        ylabel('Distance from spindle peak (Hz)')
        xlabel('Time (s)')
        
        % figSign = imagesc(flipud(...
        %     squeeze(nanmean(...
        %     true_labels(idx_clustchans, :, :), 1))), [0, 1]);
        
        
        % Put p value of cluster as title
        % -------------------------------
        
        pVal = cluster_pvals(i_clust);
        pVal = round(pVal, 4);
        if pVal < 0.001
            str_p = 'p < 0.001';
        elseif isempty(pVal)
            str_p = 'No cluster found';
        else
            str_convert = num2str(pVal(1));
            if numel(str_convert) < 6
                idx_end = numel(str_convert);
            else
                idx_end = 6;
            end
            str_p = char(strcat('p =', {' '}, str_convert(1:idx_end)));
        end
        title(str_p)
        
        
        
        saveas(gcf, [str_folderOut, filesep, ...
            char(PM_stats.Grp_subj), 'Subjects_', ...
            str_cluster, '_', num2str(i_clust), '_Cluster.png']);
        
        close
    end
    
    
    PM.Stats        = PM_stats;
    PM.Cfgs.neighb  = cfg_neighb;
    PM.Cfgs.stats   = cfg_stats;
    save([str_folderOut, filesep, ...
        'PermutationStats_', str_cluster, '_', ...
        char(PM_stats.Grp_subj), 'Subjects.mat'], ...
        'PM', 'stats', '-v7')
    
end