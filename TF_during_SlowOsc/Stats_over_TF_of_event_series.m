% Compute cluster-based statistical analysis over subjects time-frequencies

% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   N O T E !   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% As of now, the script is non-functional and still in work!


% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% 1. Files by subject containing time series by channel and by condition 
%    prepared by Compute_TF_of_event_series.m as well as the average wave
%    form of the time window
%       - PM        = parameters used during time-frequency computation
%       - TF_series = time-freq matrices
%       - SO_wave   = raw amplitude of event time window
%
% 2. Parameters
%       - File paths
seriespath          = ['D:\germanStudyData\datasetsSETS\', ...
                        'Ori_CueNight\preProcessing\', ...
                        'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
                        '05-Mar-2021_Cue_NEW\SO_timeSeries\TF_matrices\'];
%       - Channel cluster
PM.ClustOI          = 'all';
% Cell array of channels of interest. See Cluster list below. Can also be
% 'all' to include all channels
%       - subject group if interest
grp_subj            = {'BL_DvsM'}; % {'All', 'GL_DvsM', 'BL_DvsM'}
%       - Conditions to analyze
PM.Conditions       = {'ShamOn', 'OdorOn'};
%       - channel locations in fieldtrip and EEGLAB structures
chanlocspath        = ['D:\Gits\Spindle_analysis\subjectsChanlocs\', ...
                        'fieldtrip_chanlocs.mat'];
chanlocspathEEGLAB  = ['D:\Gits\Spindle_analysis\subjectsChanlocs\', ...
                        'chanlocs.mat'];



%% Set up userland
%  ------------------------------------------------------------------------

files = dir(strcat(seriespath,'*.mat'));

% Dummy file is important because most parameters will be read from this
% file
if ~exist('dummy_loaded', 'var') || dummy_loaded == 0
    dummyFile = load([seriespath, files(1).name]);
    dummy_loaded = 1;
end

% Define channels to go through
Cluster = dummyFile.PM.Clust.(PM.ClustOI)';

% TF_dimensions = trials by channels by freqs by times
v_TF_dimensions = ...
    size(dummyFile.TF_series.(char(PM.Conditions(1))).powspctrm);
% First dimension varies between subjects, but it does not matter as this
% will never be taken into account during this script. We just average over
% all trials (events).

% Prepare subject time-freq matrix
for condition = PM.Conditions
    TF_subject.(char(condition)).powspctrm = NaN(...
        length(files), ...
        numel(Cluster), ...
        v_TF_dimensions(3), ...
        v_TF_dimensions(4));
end

% Prepare subject wave form series
s_windowSize = ...
    diff(dummyFile.PM.cfg_seldat.latency) * ...
    dummyFile.PM.Info.TrialParameters.s_fs;
for condition = PM.Conditions
    WF_subject.(char(condition)) = ...
        NaN(length(files), numel(Cluster), s_windowSize);
end



%% Fruitloops
%  ------------------------------------------------------------------------

% Here, we:
% 1. Average by subject the TF matrices over channels belonging to cluster
%    of interest by condition.
% 2. Catenate all subjects mean matrices to ft_freqstatistics-interpretable
%    data structure.
% 3. Extract the average wave form over cluster channels which should
%    ressemble a theoretical slow oscillation.

for i_subj = 1:length(files)
    
    if ~exist('file_loaded', 'var') || file_loaded == 0
        load([seriespath, files(i_subj).name])
        file_loaded = 1; %#ok<NASGU>
    end
    
    disp(strcat('Subject:',     {' '},    files(i_subj).name))
    
    
    for condition = PM.Conditions
        
          
        for i_chan = 1:numel(Cluster)
            
            channel = Cluster(i_chan);
            
            idx_chan = find(strcmp(...
                dummyFile.PM.Info.ROIs.str_chans, channel));
            
            
            
            %% Average time series' TF over trials
            %  ------------------------------------------------------------
            
            TF_chan = nanmean(TF_series.(...
                char(condition)).powspctrm(:, idx_chan, :, :), 1);
            TF_subject.(char(condition)).powspctrm(i_subj, i_chan, :, :) = ...
                TF_chan;
            
            
            
            %% Average wave form of cluster
            %  ------------------------------------------------------------
        
            WF_subject.(char(condition))(i_subj, i_chan, :) = ...
                SO_wave.(char(condition))(idx_chan, :);
            
        end
        
    end
    
    file_loaded = 0;
end



%% Verify TF matrices and wave forms
%  ------------------------------------------------------------------------

% % Find minima and maxima of time-freq matrices
% % --------------------------------------------
% TF_Odor = TF_subject.OdorOn.powspctrm;
% TF_Odor = mean(TF_Odor, 1);
% TF_Odor = mean(TF_Odor, 2);
% TF_Odor = squeeze(TF_Odor);
% 
% TF_Sham = TF_subject.ShamOn.powspctrm;
% TF_Sham = mean(TF_Sham, 1);
% TF_Sham = mean(TF_Sham, 2);
% TF_Sham = squeeze(TF_Sham);
% 
% minVal = min([TF_Odor(:); TF_Sham(:)]);
% maxVal = max([TF_Odor(:); TF_Sham(:)]);
% 
% limits = [- max(abs(minVal), abs(maxVal)), ...
%     max(abs(minVal), abs(maxVal))];
% 
% 
% figure('units','normalized','outerposition', [0 0 1 0.5]);
% colororder({'b', 'k'})
% i_plot = 1;
% for condition = PM.Conditions
%     
%     subplot(1, 3, i_plot)
%     
%     % Time frequency plot
%     % -------------------
%     TF_meanChans  = mean(TF_subject.(char(condition)).powspctrm, 2);
%     
%     % There must be a more efficient way to determine v_frequencies ...
%     v_frequencies = 1:length(TF_series.OdorOn.freq);
%     v_frequencies = (max(v_frequencies) - v_frequencies) / ...
%         (max(v_frequencies) - min(v_frequencies));
%     v_frequencies = v_frequencies * PM.FrRange;
%     v_frequencies = v_frequencies - (PM.FrRange / 2);
%     v_frequencies = flip(v_frequencies);
%     
%     v_times       = 1:size(TF_meanChans, 4);
%     v_times       = v_times - (size(TF_meanChans, 4) / 2);
%     
%     TF_meanSubj.(char(condition))   = squeeze(mean(TF_meanChans, 1));
%     yyaxis left
%     pcolor(v_times, v_frequencies, TF_meanSubj.(char(condition)));
%     shading interp
%     colorbar;
%     set(gca, 'clim', limits)
%     ylabel('Distance from spindle peak (Hz)')
%     xlabel('Time (s)')
%     title(char(condition))
%     
%     hold on
%     
%     % Wave form plot
%     % --------------
%     yyaxis right
%     v_times = PM.cfg_seldat.latency(1) * ...
%         PM.Info.TrialParameters.s_fs  + 1 : 1 : ...
%         PM.cfg_seldat.latency(2) * ...
%         PM.Info.TrialParameters.s_fs;
%     WF_meanChans = mean(WF_subject.(char(condition)), 2);
%     plot(v_times, squeeze(mean(WF_meanChans, 1)), ...
%         'Color',        [0, 0, 0], ...
%         'LineWidth',    2)
%     
%     i_plot = i_plot + 1;
% end
% 
% % Overlay difference in TF
% subplot(1, 3, 3)
% v_times       = 1:size(TF_meanChans, 4);
% v_times       = v_times - (size(TF_meanChans, 4) / 2);
% pcolor(v_times, v_frequencies, TF_meanSubj.OdorOn - TF_meanSubj.ShamOn);
% shading interp
% colorbar;
% % set(gca, 'clim', [minVal, maxVal])
% ylabel('Distance from spindle peak (Hz)')
% xlabel('Time (s)')
% title(char(condition))



%% Prepare structures for stats
%  ------------------------------------------------------------------------

load(chanlocspath) % contains 'sensors'

c_fields = fieldnames(TF_series.OdorOn);
for i_field = 1:numel(c_fields)
    
    if any(strcmp(c_fields(i_field), {'powspctrm', 'cfg', 'elec'}))
        % This field is already generated and contains averaged data
        continue
    end
    
    TF_subject.OdorOn.(char(c_fields(i_field))) = ...
        TF_series.OdorOn.(char(c_fields(i_field)));
    TF_subject.ShamOn.(char(c_fields(i_field))) = ...
        TF_series.ShamOn.(char(c_fields(i_field)));
    
end
TF_subject.OdorOn.elec = sensors;
TF_subject.ShamOn.elec = sensors;



%% Select subject group
%  ------------------------------------------------------------------------

% Generate pseudo data in order to make function work
v_pseudo = NaN(1, numel(files));
[~, c_kept_subj, v_kept_subj] = ...
    f_skip_subject(v_pseudo, ...
    extractBefore(PM.Info.Subjects, '_sleep'), grp_subj);

TF_subject.OdorOn.powspctrm = ...
    TF_subject.OdorOn.powspctrm(v_kept_subj, :, :, :);
TF_subject.ShamOn.powspctrm = ...
    TF_subject.ShamOn.powspctrm(v_kept_subj, :, :, :);



%% Parameters for stats
%  ------------------------------------------------------------------------

ft_defaults

cfg_stats                       = [];
% cfg_stats.latency               = [1]; % Setting this to one rejects the
%time vector. Probably, this indicates the time vector, but it doesn't seem
%necessary to specify since fieldtrip is taking it from structure.time
cfg_stats.frequency             = 'all';
cfg_stats.channel               = 'all';
cfg_stats.correctm              = 'cluster';
cfg_stats.method                = 'montecarlo';
cfg_stats.statistic             = 'depsamplesT';
% use actvsblT for activation against baseline
cfg_stats.clusterstatistic      = 'maxsum';
cfg_stats.minnbchan             = 2;

cfg_neighb.method               = 'distance'; %'triangulation';
cfg_neighb.neighbourdist        = 3.5; % Default is "smart guess".
% 3.5 gets rid of edge clusters while maintaining high statistical 
% significance in central clusters.
cfg_neighb.channel              = 'all';
cfg_neighb.elec                 = sensors;
cfg_stats.neighbours            = ft_prepare_neighbours(cfg_neighb);

cfg_stats.tail                  = 0; % -1, +1 or 0 for one-sided or two-sided
cfg_stats.clustertail           = cfg_stats.tail;
cfg_stats.alpha                 = 0.05; % = False detection rate ?
% Commented out in order to get p values of channels that would otherwise
% be set to p = 1 in case the p value is higher than the alpha here.
% Cluster perm stats --> Which clusters are "significant"
cfg_stats.correcttail           = 'alpha'; % {'prob', 'alpha'}
% Adapts alpha value corresponding to one-tailed or two-tailed test
cfg_stats.numrandomization      = 100;
cfg_stats.clusteralpha          = 0.05;  % "which chans to use for clustering" aka
% include more or less channels
% threshold over which a triplet is chosen
% Design needs to be adapted to cfg.uvar and cfg.ivar
cfg_stats.uvar                  = 1;    % condition (uvar would be the subjects)
                                % Iactivate for cfg.statistic = 'indepsamplesT'
cfg_stats.ivar                  = 2;
cfg_stats.design                = [1:numel(c_kept_subj), ...
                                    1:numel(c_kept_subj); ...
                                    repmat(1, 1, numel(c_kept_subj)), ...
                                    repmat(2, 1, numel(c_kept_subj))];

cfg_stats.feedback = 'on';
stats  = ft_freqstatistics(cfg_stats, ...
    TF_subject.OdorOn, TF_subject.ShamOn);



%% Extract statistical information
%  ------------------------------------------------------------------------

% Here, we:
% 1. Extract information about channel positions for EEGLAB's topoplot that
%    will draw cluster-belonging electrodes over the scalp
% 2. Look if any of the clusters has a significant difference (p <
%    threshold)
% 3. Retrieve the cluster labels (1, 2, ...) and the corresponding p values
% 4. Go through each cluster found by fieldtrip and
%           - Plot the channels that belong to cluster
%           - Plot the time-frequency matrix averaged over these channels
%           - Plot the difference between the time-frequency matrices of
%             the conditions
%           - Selectively plot only the part of the "difference matrix"
%             that belongs to the current cluster



%% Look up the locations of channels that were analyzed
%  ------------------------------------------------------------------------

load(chanlocspathEEGLAB); % Creates chanlocs variable that is
% used below

idx_ROI = NaN(1, numel(Cluster));
for i_chan = 1 :numel(Cluster)
    idx_ROI(i_chan) = find(strcmp({chanlocs.labels}, ...
        Cluster(i_chan)));
end
chanlocs = chanlocs(idx_ROI);


% Clusters and probabilities are given by channel by freq by time.
% We will extract the channels that show significant differences and plot
% their average TF. Then, we will overlay the time-frequency cluster part
% that shows significant differences.

% Look for clusters with probability of <= 0.2 (permissive threshold)
s_thrld = 0.2;
p_all = stats.prob(:);

if any(p_all < s_thrld)
    
    
    
    %% Retrieve information about clusters
    %  --------------------------------------------------------------------
    
    dimensions      = size(stats.prob);
    
    % First, we move everything over to one labelmat array since we are
    % equally interested in pos and neg clusters
    cluster_labels  = zeros(dimensions);
    cluster_pvals   = zeros(dimensions);
    
    for i_chan = 1:dimensions(1)
        for i_freq = 1:dimensions(2)
            for i_time = 1:dimensions(3)
                pos_label = ...
                    stats.posclusterslabelmat(i_chan, i_freq, i_time);
                % We fill our cluster matrix only with clusters that show
                % significant p below threshold value.
                if pos_label > 0 && ...
                        stats.posclusters(pos_label).prob < s_thrld
                    cluster_labels(i_chan, i_freq, i_time) = ...
                        pos_label;
                    cluster_pvals(i_chan, i_freq, i_time) = ...
                        stats.posclusters(pos_label).prob;
                end
            end
        end
    end
    
    % Overlay negative clusters with new labels
    s_raise_label = max(cluster_labels(:));
    for i_chan = 1:dimensions(1)
        for i_freq = 1:dimensions(2)
            for i_time = 1:dimensions(3)
                neg_label = ...
                    stats.negclusterslabelmat(i_chan, i_freq, i_time);
                if neg_label > 0 && ...
                        stats.negclusters(neg_label).prob < s_thrld
                    cluster_labels(i_chan, i_freq, i_time) = ...
                        neg_label + s_raise_label;
                    cluster_pvals(i_chan, i_freq, i_time) = ...
                        stats.negclusters(neg_label).prob;
                end
            end
        end
    end
    
    
    
    %% Plot time-frequency clusters of channels belonging to them
    %  --------------------------------------------------------------------
    
    significant_clusters = unique(cluster_labels(:))';
    significant_clusters(significant_clusters == 0) = [];
    for i_clust = significant_clusters
    
        
        % Find channels, frequencies and times that are involved in cluster
        % -----------------------------------------------------------------
        
        idx_clustchans = [];
        idx_clusttimes = [];
        idx_clustfreqs = [];
        for i_chan = 1:dimensions(1)
            for i_freq = 1:dimensions(2)
                for i_time = 1:dimensions(3)
                    if cluster_labels(i_chan, i_freq, i_time) == i_clust
                        idx_clustchans = [idx_clustchans, i_chan];
                        idx_clusttimes = [idx_clusttimes, i_time];
                        idx_clustfreqs = [idx_clustfreqs, i_freq];
                    end
                end
            end
        end
        idx_clustchans = unique(idx_clustchans);
        idx_clusttimes = unique(idx_clusttimes);
        idx_clustfreqs = unique(idx_clustfreqs);
        
        figure('units','normalized','outerposition', [0 0 1 0.5]);
        
        
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
            'plotchans', [], ...
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
            'plotchans', idx_clustchans, ...
            'emarker', {'.', [0.1, 0.1, 0.1], 20, 1}, ...
            'plotrad', max(max([chanlocs.radius]), 0.5));
        hold off
        
        
        % Find minima and maxima of time-freq matrices
        % --------------------------------------------
        
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
        
        limits = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        i_plot = 2;
        for condition = PM.Conditions
            
            TF_cluster = TF_subject.(...
                char(condition)).powspctrm(:, idx_clustchans, :, :);
            
            subplot(1, 5, i_plot)
            
            
            % Time frequency plot
            % -------------------
            
            TF_meanChans  = mean(TF_cluster, 2);
            
            v_frequencies = 1:dimensions(2);
            v_frequencies = (max(v_frequencies) - v_frequencies) / ...
                (max(v_frequencies) - min(v_frequencies));
            v_frequencies = v_frequencies * PM.FrRange;
            v_frequencies = v_frequencies - (PM.FrRange / 2);
            v_frequencies = flip(v_frequencies);
            
            v_times       = 1:dimensions(3);
            v_times       = v_times - (size(TF_meanChans, 4) / 2);
            
            TF_meanSubj.(char(condition)) = squeeze(mean(TF_meanChans, 1));
            yyaxis left
            pcolor(v_times, v_frequencies, TF_meanSubj.(char(condition)));
            shading interp
            colorbar;
            colormap parula
            set(gca, 'clim', limits)
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
            WF_meanClust.(char(condition)) = mean(WF_subject.(...
                char(condition))(:, idx_clustchans, :), 2);
            plot(v_times, ...
                squeeze(mean(WF_meanClust.(char(condition)), 1)), ...
                'Color',        [0, 0, 0], ...
                'LineWidth',    2)
            
            i_plot = i_plot + 1;
        end
        
        
        % Plot difference between the two conditions' TF
        % ----------------------------------------------
        
        TF_difference   = TF_meanSubj.OdorOn - TF_meanSubj.ShamOn;
        
        minVal = min(TF_difference(:));
        maxVal = max(TF_difference(:));
        
        limits = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        
        subplot(1, 5, 4)
        v_times         = 1:size(TF_meanChans, 4);
        v_times         = v_times - (size(TF_meanChans, 4) / 2);
        pcolor(v_times, v_frequencies, TF_difference);
        shading interp
        colorbar;
        set(gca, 'clim', limits)
        ylabel('Distance from spindle peak (Hz)')
        xlabel('Time (s)')
        title('Odor - Sham')
        
        
        % Plot time-frequencies of significant difference
        % -----------------------------------------------
        
        subplot(1, 5, 5)
        v_times         = 1:size(TF_meanChans, 4);
        v_times         = v_times - (size(TF_meanChans, 4) / 2);
        TF_selective    = NaN(dimensions(2), dimensions(3));
        TF_selective(idx_clustfreqs, idx_clusttimes) = ...
            TF_difference(idx_clustfreqs, idx_clusttimes);
        pcolor(v_times, v_frequencies, TF_selective);
        shading interp
        colorbar;
        set(gca, 'clim', limits)
        ylabel('Distance from spindle peak (Hz)')
        xlabel('Time (s)')
        
        
        % Put p value of cluster as title
        % -------------------------------
        
        all_labels = cluster_labels(:);
        all_labels = find(all_labels == i_clust);
        all_pvals = cluster_pvals(:);
        all_pvals = all_pvals(all_labels);
        
        if numel(unique(all_pvals)) ~= 1
            error('P value does not match cluster label!')
        end
        
        pVal = all_pvals(1);
        pVal = round(pVal, 4);
        if pVal < 0.001
            str_p = 'p < 0.001';
        else
            str_convert = num2str(all_pvals(1));
            if numel(str_convert) < 6
                idx_end = numel(str_convert);
            else
                idx_end = 6;
            end
            str_p = char(strcat('p =', {' '}, str_convert(1:idx_end)));
        end
        title(str_p)
        
        
    end
    
end




