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
                        '05-Mar-2021_Cue\SO_timeSeries\TF_matrices\'];
% seriespath          = ['/mnt/disk1/sleep/Datasets/', ...
%                         'CueD_SO_TimeSeires/TF_matrices/'];
%       - Channel cluster
PM_stats.ClustOI    = 'all';
% Cell array of channels of interest. See Cluster list below. Can also be
% 'all' to include all channels
%       - subject group if interest
PM_stats.Grp_subj   = {'All'}; % {'All', 'GL_DvsM', 'BL_DvsM'}
%       - Conditions to analyze
PM_stats.Conditions = {'ShamOn', 'OdorOn'};
%       - channel locations in fieldtrip and EEGLAB structures
chanlocspath        = ['D:\Gits\Spindle_analysis\subjectsChanlocs\', ...
                        'fieldtrip_chanlocs.mat'];
% chanlocspath        = ['/home/sleep/Documents/DAVID/Gits/', ...
%                         'Spindle_analysis/subjectsChanlocs/', ...
%                         'fieldtrip_chanlocs.mat'];
chanlocspathEEGLAB  = ['D:\Gits\Spindle_analysis\subjectsChanlocs\', ...
                        'chanlocs.mat'];
% chanlocspathEEGLAB  = ['/home/sleep/Documents/DAVID/Gits/', ...
%                         'Spindle_analysis/subjectsChanlocs/', ...
%                         'chanlocs.mat'];
%       - paths to toolboxes
fieldtrippath       = 'D:\MATLAB\fieldtrip-20200831';
% fieldtrippath       = '/home/sleep/Documents/MATLAB/fieldtrip-20200831';
eeglabpath          = 'D:\MATLAB\eeglab2019_1';
% eeglabpath          = '/home/sleep/Documents/MATLAB/eeglab2019_1';



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
Cluster = dummyFile.PM.Clust.(PM_stats.ClustOI)';

% TF_dimensions = trials by channels by freqs by times
v_TF_dimensions = ...
    size(dummyFile.TF_series.(char(PM_stats.Conditions(1))).powspctrm);
% First dimension varies between subjects, but it does not matter as this
% will never be taken into account during this script. We just average over
% all trials (events).

% Prepare subject time-freq matrix
for condition = PM_stats.Conditions
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
for condition = PM_stats.Conditions
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
        load([seriespath, files(i_subj).name]);
        file_loaded = 1; %#ok<NASGU>
    end
    
    disp(strcat('Subject:', {' '}, files(i_subj).name))
    
    
    for condition = PM_stats.Conditions
        
          
        for i_chan = 1:numel(Cluster)
            
            channel = Cluster(i_chan);
            
            idx_chan = find(strcmp(...
                dummyFile.PM.Info.ROIs.str_chans, channel));
            
            
            
            %% Average time series' TF over trials
            %  ------------------------------------------------------------
            
            TF_chan = nanmean(TF_series.(...
                char(condition)).powspctrm(:, idx_chan, :, :), 1);
            TF_subject.(...
                char(condition)).powspctrm(i_subj, i_chan, :, :) = ...
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

limitsTF = [- max(abs(minVal), abs(maxVal)), ...
    max(abs(minVal), abs(maxVal))];

average_WF.ShamOn = squeeze(mean(mean(WF_subject.ShamOn, 1), 2));
average_WF.OdorOn = squeeze(mean(mean(WF_subject.OdorOn, 1), 2));

minVal      = min([average_WF.ShamOn(:); average_WF.OdorOn(:)]);
maxVal      = max([average_WF.ShamOn(:); average_WF.OdorOn(:)]);

limitsWF = [- max(abs(minVal), abs(maxVal)), ...
    max(abs(minVal), abs(maxVal))];



figure('units','normalized','outerposition', [0 0 1 0.5]);
colororder({'b', 'k'})
i_plot = 1;
for condition = PM.Conditions
    
    subplot(1, 3, i_plot)
    
    % Time frequency plot
    % -------------------
    TF_meanChans  = mean(TF_subject.(char(condition)).powspctrm, 2);
    
    % There must be a more efficient way to determine v_frequencies ...
    v_frequencies = 1:length(TF_series.OdorOn.freq);
    v_frequencies = (max(v_frequencies) - v_frequencies) / ...
        (max(v_frequencies) - min(v_frequencies));
    v_frequencies = v_frequencies * PM.FrRange;
    v_frequencies = v_frequencies - (PM.FrRange / 2);
    v_frequencies = flip(v_frequencies);
    
    v_times       = 1:size(TF_meanChans, 4);
    v_times       = v_times - (size(TF_meanChans, 4) / 2);
    
    TF_meanSubj.(char(condition))   = squeeze(mean(TF_meanChans, 1));
    yyaxis left
    pcolor(v_times, v_frequencies, TF_meanSubj.(char(condition)));
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
    plot(v_times, average_WF.(char(condition)), ...
        'Color',        [0, 0, 0], ...
        'LineWidth',    2)
    ylim(limitsWF)
    
    i_plot = i_plot + 1;
end

% Overlay difference in TF
subplot(1, 3, 3)
v_times       = 1:size(TF_meanChans, 4);
v_times       = v_times - (size(TF_meanChans, 4) / 2);
pcolor(v_times, v_frequencies, TF_meanSubj.OdorOn - TF_meanSubj.ShamOn);
shading interp
colorbar;
% set(gca, 'clim', [minVal, maxVal])
ylabel('Distance from spindle peak (Hz)')
xlabel('Time (s)')
title(char(condition))

close



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
    extractBefore(PM.Info.Subjects, '_sleep'), PM_stats.Grp_subj);

TF_subject.OdorOn.powspctrm = ...
    TF_subject.OdorOn.powspctrm(v_kept_subj, :, :, :);
TF_subject.ShamOn.powspctrm = ...
    TF_subject.ShamOn.powspctrm(v_kept_subj, :, :, :);



% % Prepare Andrea-like structure
% % -----------------------------
% 
% % Place each rpt in separate cell. This does not change the outcome as
% % compared to having all rpt inside same powspctrm field along the first
% % dimension.
% 
% retained_subj = find(v_kept_subj)';
% for i_rpt = 1:numel(retained_subj)
%    
%     struct_tmp = TF_subject.OdorOn;
%     struct_tmp = rmfield(struct_tmp, 'powspctrm');
%     struct_tmp.powspctrm = ...
%         TF_subject.OdorOn.powspctrm(retained_subj(i_rpt), :, :, :);
%     TF_subject_cell.OdorOn{i_rpt} = struct_tmp;
%     
%     struct_tmp = TF_subject.ShamOn;
%     struct_tmp = rmfield(struct_tmp, 'powspctrm');
%     struct_tmp.powspctrm = ...
%         TF_subject.ShamOn.powspctrm(retained_subj(i_rpt), :, :, :);
%     TF_subject_cell.ShamOn{i_rpt} = struct_tmp;
%     
% end



%% Parameters for stats
%  ------------------------------------------------------------------------

addpath(fieldtrippath)
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
% cfg_stats.alpha                 = 0.1;
cfg_stats.alpha                 = 0.05; % = False detection rate ?
% Commented out in order to get p values of channels that would otherwise
% be set to p = 1 in case the p value is higher than the alpha here.
% Cluster perm stats --> Which clusters are "significant"
cfg_stats.correcttail           = 'alpha'; % {'prob', 'alpha'}
% Adapts alpha value corresponding to one-tailed or two-tailed test
cfg_stats.numrandomization      = 125;
% cfg_stats.clusteralpha          = 0.1;
cfg_stats.clusteralpha          = 0.05;  % "which points for clustering" 
% aka include more or less channels, freqs and times
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

% stats  = ft_freqstatistics(cfg_stats, ...
%     TF_subject_cell.OdorOn{:}, TF_subject_cell.ShamOn{:});



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


addpath(eeglabpath)
eeglab nogui % Point to eeglab functions (topoplot, ...)



%% Look for clusters with significant differences between conditions
%  ------------------------------------------------------------------------

% Clusters and probabilities are given by channel by freq by time.
% We will extract the channels that show significant differences and plot
% their average TF. Then, we will overlay the time-frequency cluster part
% that shows significant differences.

% Look for clusters with probability of <= 0.2 (permissive threshold)
s_thrld = 0.2;



%% Retrieve information about clusters
%  --------------------------------------------------------------------

% We go through the label matrices for cluster matrices and search for
% labels higher than 0. Each cluster found by fieldtrip is assigned an
% individual integer > 0. At the coordinates of the cluster matrices
% belonging to such labels, we extract the corresponding labels as well
% as p values and put them into cluster_labels and cluster_pvals.
% We move everything over to one labelmat array since we are
% equally interested in pos and neg clusters.

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
    
    
    
    %% Plot
    %  ----------------------------------------------------------------
    
    % We will generate a 5 panel plot:
    % 1. Channels belonging to current cluster
    % 2. The TF averaged over cluster channels for vehicle condition
    %    with average wave form of time series
    % 3. The TF averaged over cluster channels for odor cue condition
    %    with average wave form of time series
    % 4. The difference between the TF of conditions
    % 5. Only the section of the TF that belongs to the cluster
    
    
    figure('units','normalized','outerposition', [0 0 1 0.3]);
    colororder({'k', 'r'})
    
    
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
    
    limitsTF = [- max(abs(minVal), abs(maxVal)), ...
        max(abs(minVal), abs(maxVal))];
    
    average_WFCluster.ShamOn = squeeze(...
        mean(mean(WF_subject.ShamOn(:, idx_clustchans, :), 2), 1));
    average_WFCluster.OdorOn = squeeze(...
        mean(mean(WF_subject.OdorOn(:, idx_clustchans, :), 2), 1));
    
    minVal = min([average_WFCluster.ShamOn(:); ...
        average_WFCluster.OdorOn(:)]);
    maxVal = max([average_WFCluster.ShamOn(:); ...
        average_WFCluster.OdorOn(:)]);
    
    limitsWF = [- max(abs(minVal), abs(maxVal)), ...
        max(abs(minVal), abs(maxVal))];
    
    
    i_plot = 2;
    for condition = PM_stats.Conditions
        
        TF_meanSubj.(char(condition))  = ...
            squeeze(mean(TF_subject.(char(condition)).powspctrm, 1));
        
        subplot(1, 5, i_plot)
        
        
        % Time frequency plot
        % -------------------
        
        TF_meanCluster.(char(condition)) = squeeze(mean(TF_meanSubj.(...
            char(condition))(idx_clustchans, :, :), 1));
        
        v_frequencies = 1:dimensions(2);
        v_frequencies = (max(v_frequencies) - v_frequencies) / ...
            (max(v_frequencies) - min(v_frequencies));
        v_frequencies = v_frequencies * PM.FrRange;
        v_frequencies = v_frequencies - (PM.FrRange / 2);
        v_frequencies = flip(v_frequencies);
        
        v_times       = 1:dimensions(3);
        v_times       = v_times - (dimensions(3) / 2);
        v_times       = v_times ./ PM.Info.TrialParameters.s_fs;
        
        yyaxis left
        pcolor(v_times, v_frequencies, TF_meanCluster.(char(condition)));
        shading interp
        colorbar;
        colormap parula
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
        WF_meanClust.(char(condition)) = mean(WF_subject.(...
            char(condition))(:, idx_clustchans, :), 2);
        plot(v_times, average_WFCluster.(char(condition)), ...
            'Color',        'r', ...
            'LineWidth',    2)
        ylim(limitsWF)
        
        i_plot = i_plot + 1;
    end
    
    
    TF_difference   = TF_meanSubj.OdorOn - TF_meanSubj.ShamOn;
    
    
    % Plot difference between the two conditions' TF
    % ----------------------------------------------
    
    TF_clustDiff   = TF_meanCluster.OdorOn - TF_meanCluster.ShamOn;
    
    minVal = min(TF_clustDiff(:));
    maxVal = max(TF_clustDiff(:));
    
    limitsTF = [- max(abs(minVal), abs(maxVal)), ...
        max(abs(minVal), abs(maxVal))];
    
    
    subplot(1, 5, 4)
    v_times         = 1:dimensions(3);
    v_times         = v_times - (dimensions(3) / 2);
    v_times         = v_times ./ PM.Info.TrialParameters.s_fs;
    figDiff = pcolor(v_times, v_frequencies, TF_clustDiff);
    shading interp
    colorbar;
    set(gca, 'clim', limitsTF)
    ylabel('Distance from spindle peak (Hz)')
    xlabel('Time (s)')
    title('Odor - Sham')
    
    
    % Plot time-frequencies of significant difference (cluster)
    % ---------------------------------------------------------
    
    subplot(1, 5, 5)
    v_times         = 1:dimensions(3);
    v_times         = v_times - (dimensions(3) / 2);
    v_times         = v_times ./ PM.Info.TrialParameters.s_fs;
    TF_selective    = NaN(dimensions);
    
    TF_selective(true_labels == 1) = TF_difference(true_labels == 1);
    TF_selective = squeeze(...
        nanmean(TF_selective(idx_clustchans, :, :), 1));
    
    figSign = pcolor(v_times, v_frequencies, TF_selective);
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
    
    saveas(gcf, [cd, filesep, 'PermutationStats_', ...
        char(PM_stats.Grp_subj), 'Subjects_', 'Cluster', ...
        num2str(i_clust), '.png']);
    
    
end


%% Save statistical output
%  ------------------------------------------------------------------------

PM.Stats        = PM_stats;
PM.Cfgs.neighb  = cfg_neighb;
PM.Cfgs.stats   = cfg_stats;
save([cd, filesep, 'PermutationStats_', ...
    char(PM_stats.Grp_subj), 'Subjects.mat'], ...
    'PM', 'stats', '-v7')

