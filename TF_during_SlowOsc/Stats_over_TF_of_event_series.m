% Compute cluster-based statistical analysis over subjects time-frequencies

% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% 1. Files by subject containing time series by channel and by condition 
%    prepared by Compute_TF_of_event_series.m as well as the average wave
%    form of the time window
%       - PM            = parameters used during time-frequency computation
%       - TF_series     = time-freq matrices
%       - SO_wave       = raw amplitude of event time window
%
% 2. Files by subject containing the slow ocillation time series. This file
%    should contain spindle latencies as follows:
%       - SS_latencies  = spindle latencies where each field and subfield 
%                         is channel and condition, respectively
%
% 3. Parameters
%       - File paths
SOseriespath        = ['D:\germanStudyData\datasetsSETS\', ...
                         'Ori_CueNight\preProcessing\', ...
                         'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
                         '12-Jun-2021_Cue\SO_timeSeries_Upstate_15s\'];
seriespath          = [SOseriespath, '\TF_matrices\'];
% SOseriespath contains datasets with spindle latencies
str_folderOut       = [cd, filesep, 'Clusters'];
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
%       - path to channel clustering function: All this function does is
%         create cell arrays of electrode labels belonging to a group. The
%         repository containing the function can be found at
%         https://github.com/davidmarcelbaum/EEG_channels
chanclusterpath     = 'D:\Gits\EEG_channels';



%% Set up userland
%  ------------------------------------------------------------------------

files           = dir(strcat(seriespath,'*.mat'));
filesSOseries   = dir(strcat(SOseriespath,'*.mat'));

% Dummy file is important because most parameters will be read from this
% file
if ~exist('dummy_loaded', 'var') || dummy_loaded == 0
    dummyFile = load([seriespath, files(1).name]);
    dummy_loaded = 1;
end

% Define channels to go through
if ~strcmp(PM_stats.ClustOI, 'all')
    warning('Make sure only desired channels are part of the selected cluster')
end
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

% Look up the locations of channels that were analyzed
load(chanlocspathEEGLAB); % Creates chanlocs variable that is used below
idx_ROI = NaN(1, numel(Cluster));
for i_chan = 1 :numel(Cluster)
    idx_ROI(i_chan) = find(strcmp({chanlocs.labels}, ...
        Cluster(i_chan)));
end
chanlocs = chanlocs(idx_ROI);

load(chanlocspath) % contains 'sensors'

addpath(eeglabpath)
eeglab nogui % Point to eeglab functions (topoplot, ...)

% Prepare outputs
mkdir(str_folderOut)


%% Fruitloops
%  ------------------------------------------------------------------------

% Here, we:
% 1. Average by subject the TF matrices over channels belonging to cluster
%    of interest by condition.
% 2. Catenate all subjects mean matrices to ft_freqstatistics-interpretable
%    data structure.
% 3. Extract the average wave form over cluster channels which should
%    ressemble a theoretical slow oscillation.
% 4. Extract the spindle latencies for each channel and catenate them. This
%    serves for identifying outliers in erms of their spindle - SO coupling

for i_subj = 1:length(files)
    
    
    disp(strcat('Subject:', {' '}, files(i_subj).name))
    
    if ~exist('file_loaded', 'var') || file_loaded == 0
        load([seriespath, files(i_subj).name]);
        
        % Holds spindle latencies
        SOseries = load([SOseriespath, filesSOseries(i_subj).name]);
        file_loaded = 1; %#ok<NASGU>
    end
    
    
    % We will catenate spindle latencies of all channels since the
    % number of detected spindles that are coupled to slow oscillations
    % is really small (between 5 and 10 most of time)
    latencies_allChans  = [];
    
    
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
            
            
            
            %% Extract spindle latencies around slow oscillations
            %  ------------------------------------------------------------
            
            latencies_allChans = [latencies_allChans, ...
                SOseries.SS_latencies.(char(channel)).(char(condition))];
            
        end
        
        SS_latencies.(char(condition)){i_subj}     = latencies_allChans;
        
    end
    
    file_loaded = 0;
end



%% Prepare structures for stats
%  ------------------------------------------------------------------------

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



%% Intermediate storage of outputs
%  ------------------------------------------------------------------------

Night.SS_latencies = SS_latencies;
Night.WF_subject   = WF_subject;
Night.TF_subject   = TF_subject;

save([cd, filesep, 'Night.mat'], 'Night', '-v7')



%% Select subject group and removal of outliers
%  ------------------------------------------------------------------------

% =========================================================================
% warning('Next lines of code are there to combine nights ...')
% clearvars TF_subject SS_latencies WF_subject
% 
% NightD = load('NightD.mat');
% NightM = load('NightM.mat');
% 
% SS_latencies.OdorOn = NightD.Night.SS_latencies.OdorOn;
% SS_latencies.ShamOn = NightM.Night.SS_latencies.OdorOn;
% WF_subject.OdorOn   = NightD.Night.WF_subject.OdorOn;
% WF_subject.ShamOn   = NightM.Night.WF_subject.OdorOn;
% TF_subject.OdorOn   = NightD.Night.TF_subject.OdorOn;
% TF_subject.ShamOn   = NightM.Night.TF_subject.OdorOn;
% =========================================================================


% Generate pseudo data in order to make function work
v_pseudo = NaN(1, numel(files));
[~, c_kept_subj, v_kept_subj] = ...
    f_skip_subject(v_pseudo, ...
    extractBefore(PM.Info.Subjects, '_sleep'), PM_stats.Grp_subj);

% Spindle latencies are normalized to the chosen center point of slow
% oscillations. We will have a look at each mean latency of each subject
% and take out the ones that are "off" compared to the group
meanLatencies.OdorOn    = NaN(numel(c_kept_subj), 1);
meanLatencies.ShamOn    = NaN(numel(c_kept_subj), 1);
for i_subj = 1:numel(c_kept_subj)
   
    meanLatencies.OdorOn(i_subj) = median(SS_latencies.OdorOn{i_subj});
    meanLatencies.ShamOn(i_subj) = mean(SS_latencies.ShamOn{i_subj});
    
end
[~, v_outlier.Odor]     = rmoutliers(meanLatencies.OdorOn, 'mean');
[~, v_outlier.Sham]     = rmoutliers(meanLatencies.ShamOn, 'mean');

% If at least one condition detects and outlier, we remove the
% corresponding subject entirely
v_outlier.ConditionIndep = logical(v_outlier.Odor + v_outlier.Sham);

v_kept_subj(v_outlier.ConditionIndep) = false;
c_kept_subj(v_outlier.ConditionIndep) = [];

% Adapt data arrays to outliers detected
TF_subject.OdorOn.powspctrm = ...
    TF_subject.OdorOn.powspctrm(v_kept_subj, :, :, :);
TF_subject.ShamOn.powspctrm = ...
    TF_subject.ShamOn.powspctrm(v_kept_subj, :, :, :);

SS_latencies.OdorOn = SS_latencies.OdorOn(v_kept_subj);
SS_latencies.ShamOn = SS_latencies.ShamOn(v_kept_subj);

WF_subject.ShamOn   = WF_subject.ShamOn(v_kept_subj, :, :);
WF_subject.OdorOn   = WF_subject.OdorOn(v_kept_subj, :, :);



% Here, we have data input defined as repetitions x channels x
% times x frequencies and let the permutation-based analysis find channel
% clusters itself. We vary channels according to the selected cluster of
% interest.

% Clusters of interest
% --------------------
addpath(chanclusterpath)
PM.Clust        = f_chan_clusters;
PM.Clust.all    = dummyFile.PM.Clust.all;
PM.Clust.pseudo = [];

Clusters = fieldnames(PM.Clust);

for i_defined_clust = 1:numel(Clusters)
    
    clearvars TF_Chans
    
    
    
    %% General parameters for stats
    %  --------------------------------------------------------------------
    
    addpath(fieldtrippath)
    ft_defaults
    
    cfg_stats                   = [];
    % cfg_stats.latency           = [1]; % Setting this to one rejects 
    % the time vector. Probably, this indicates the time vector, but it
    % doesn't seemnecessary to specify since fieldtrip is taking it from 
    % structure.time
    cfg_stats.frequency         = 'all';
    cfg_stats.channel           = 'all';
    cfg_stats.correctm          = 'cluster';
    cfg_stats.method            = 'montecarlo';
    cfg_stats.statistic         = 'depsamplesT';
    % use actvsblT for activation against baseline
    cfg_stats.clusterstatistic  = 'maxsum';
    cfg_stats.tail              = 0; % -/+1 = one-sided 0 = two-sided
    cfg_stats.clustertail       = cfg_stats.tail;
    % cfg_stats.alpha             = 0.1;
    cfg_stats.alpha             = 0.05; % = False detection rate ?
    % Commented out in order to get p values of channels that would otherwise
    % be set to p = 1 in case the p value is higher than the alpha here.
    % Cluster perm stats --> Which clusters are "significant"
    cfg_stats.correcttail       = 'alpha'; % {'prob', 'alpha'}
    % Adapts alpha value corresponding to one-tailed or two-tailed test
    cfg_stats.clusteralpha      = 0.05; % "which points for clustering"
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
    
    
    
    %% Prepare data structures
    %  --------------------------------------------------------------------
    
    str_cluster     = char(Clusters(i_defined_clust));
    
    % Since our change in electrode rejection, some clusters such as
    % temporal are listing channels in the PM.Clust structure that are
    % not available anymore. This generates issues when powspctr for
    % example is only including 5 channels in case of left_temporal but
    % assigning the 7 channels listed in PM.Clust to the
    % TF_Chans.(...).label structure. We correct this here:
    Clust2Correct = PM.Clust.(str_cluster);
    if ~isempty(Clust2Correct)
        Clust2Correct = Clust2Correct(...
            ismember(Clust2Correct, dummyFile.PM.Clust.all));
        PM.Clust.(str_cluster) = Clust2Correct;
    end
    
    
    if ~strcmp(str_cluster, 'pseudo')
        
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
        
        
        % Additional paramaeters for stats
        % --------------------------------
        
        cfg_stats.minnbchan             = 2;
        cfg_neighb                      = [];
        cfg_neighb.method               = 'distance'; %'triangulation';
        cfg_neighb.neighbourdist        = 3.5; % Default is "smart guess".
        % 3.5 gets rid of edge clusters while maintaining high statistical
        % significance in central clusters.
        cfg_neighb.channel              = 'all';
        cfg_neighb.elec                 = clust_sensors;
        
        cluster_chanlocs                = chanlocs;
        
    else
        
        % Here, we build a pseudo channel group where each channel is an
        % electrode in the center (hand-picked) of the pre-defined channel
        % groups (left_frontal, ...). The position of these channels is not
        % that important as it will only serve for satisfying the
        % neighbouring during permutation analysis, but minnbchan will be
        % set to 0 so that each individual channel is auto-sufficient in 
        % order to determine ignificant changes inside each group.
        
        % Visualize channel locations
        % for curr_cluster = Clusters'
        %     
        %     if any(strcmp(curr_cluster, {'all', 'pseudo'}))
        %         continue
        %     end
        %     
        %     idx_defchans    = find(ismember(...
        %         dummyFile.PM.Clust.all, PM.Clust.(char(curr_cluster))));
        %     
        %     no_results = zeros(numel(chanlocs), 1);
        %     figure
        %     topoplot(no_results, chanlocs, ...
        %         'style', 'blank', ...
        %         'electrodes', 'labels', ...
        %         'shading', 'interp', ...
        %         'headcolor', [0, 0, 0], ...
        %         'plotchans', idx_defchans, ...
        %         'conv', 'on', ...
        %         'emarker', {'.', [0, 0, 0], 5, 1}, ...
        %         'plotrad', max(max([chanlocs.radius]),0.5));
        %     title(curr_cluster)
        % end
        % close all
        pseudo.occipital        = 'E75';
        pseudo.right_occipital  = 'E89';
        pseudo.left_occipital   = 'E69';
        pseudo.parietal         = 'E62';
        pseudo.right_parietal   = 'E91';
        pseudo.left_parietal    = 'E59';
        pseudo.right_temporal   = 'E108';
        pseudo.left_temporal    = 'E45';
        pseudo.central          = 'E6';
        pseudo.right_central    = 'E105';
        pseudo.left_central     = 'E30';
        pseudo.frontal          = 'E16';
        pseudo.right_frontal    = 'E2';
        pseudo.left_frontal     = 'E26';
        
        % Build data structures according to channel (groups): Mean cluster
        % channels and extract chan position
        c_pseudos   = fieldnames(pseudo);
        i_clust     = 1;
        cat_pseudos = NaN(1, numel(c_pseudos));
        for curr_cluster = c_pseudos'
        
            idx_defchans    = find(ismember(dummyFile.PM.Clust.all, ...
                                PM.Clust.(char(curr_cluster))));
            pseudo_chan     = pseudo.(char(curr_cluster));
            idx_pseudo      = find(strcmp(dummyFile.PM.Clust.all, ...
                                pseudo_chan));
            cat_pseudos(i_clust) = idx_pseudo;
            
            TF_Chans.ShamOn.powspctrm(:, i_clust, :, :) = mean(...
                TF_subject.ShamOn.powspctrm(:, idx_defchans, :, :), 2);
            TF_Chans.ShamOn.label(i_clust, 1)           = curr_cluster;
            TF_Chans.ShamOn.elec.chanpos(i_clust, :)    = ...
                TF_subject.ShamOn.elec.chanpos(idx_pseudo, :);
            TF_Chans.ShamOn.elec.chantype(i_clust, 1)   = ...
                TF_subject.ShamOn.elec.chantype(idx_pseudo);
            TF_Chans.ShamOn.elec.chanunit(i_clust, 1)   = ...
                TF_subject.ShamOn.elec.chanunit(idx_pseudo, 1);
            TF_Chans.ShamOn.elec.elecpos(i_clust, :)    = ...
                TF_subject.ShamOn.elec.elecpos(idx_pseudo, :);
            TF_Chans.ShamOn.elec.label(1, i_clust)      = curr_cluster;
            
            TF_Chans.OdorOn.powspctrm(:, i_clust, :, :) = mean(...
                TF_subject.OdorOn.powspctrm(:, idx_defchans, :, :), 2);
            TF_Chans.OdorOn.label(i_clust, 1)           = curr_cluster;
            TF_Chans.OdorOn.elec.chanpos(i_clust, :)    = ...
                TF_subject.OdorOn.elec.chanpos(idx_pseudo, :);
            TF_Chans.OdorOn.elec.chantype(i_clust, 1)   = ...
                TF_subject.OdorOn.elec.chantype(idx_pseudo);
            TF_Chans.OdorOn.elec.chanunit(i_clust, 1)   = ...
                TF_subject.OdorOn.elec.chanunit(idx_pseudo, 1);
            TF_Chans.OdorOn.elec.elecpos(i_clust, :)    = ...
                TF_subject.OdorOn.elec.elecpos(idx_pseudo, :);
            TF_Chans.OdorOn.elec.label(1, i_clust)      = curr_cluster;
            
            
            i_clust = i_clust + 1;
        end
        
        TF_Chans.ShamOn.freq        = TF_subject.ShamOn.freq;
        TF_Chans.ShamOn.time        = TF_subject.ShamOn.time;
        TF_Chans.ShamOn.elec.type   = TF_subject.ShamOn.elec.type;
        TF_Chans.ShamOn.elec.unit   = TF_subject.ShamOn.elec.unit;
        TF_Chans.ShamOn.dimord      = 'rpt_chan_freq_time';
        
        TF_Chans.OdorOn.freq        = TF_subject.OdorOn.freq;
        TF_Chans.OdorOn.time        = TF_subject.OdorOn.time;
        TF_Chans.OdorOn.elec.type   = TF_subject.OdorOn.elec.type;
        TF_Chans.OdorOn.elec.unit   = TF_subject.OdorOn.elec.unit;
        TF_Chans.OdorOn.dimord      = 'rpt_chan_freq_time';
        
        clust_sensors               = TF_Chans.ShamOn.elec;
        
        idx_defchans                = sort(cat_pseudos); % reset for later
        
        
        % Additional paramaeters for stats
        % --------------------------------
        
        cfg_stats.minnbchan         = 0;
        cfg_neighb                  = [];
        cfg_neighb.method           = 'distance'; %'triangulation';
        cfg_neighb.neighbourdist    = 12;
        cfg_neighb.channel          = 'all';
        cfg_neighb.elec             = clust_sensors;
                
    end
    
    
    
    %% Run stats
    %  --------------------------------------------------------------------
    
    if strcmp(str_cluster, 'all')
        cfg_stats.numrandomization  = 1000;
    else
        cfg_stats.numrandomization  = 10000;
    end
    cfg_stats.neighbours            = ft_prepare_neighbours(cfg_neighb);
    
    
    % ---------------------------------------------------------------------
    %                             /!\/!\/!\ 
    % Main function to perform statistics: If you get an error saying could
    % not determine critical cluster value or such, always check first that
    % dimensions and numbers of channels in labels etc correspond to the
    % rest of the structure for TF_Chans. Most of the time, if everything
    % is set correctly in the cfg_stats, the error comes from
    % inconsistencies in the TF_Chans strucure...
    stats                           = ft_freqstatistics(cfg_stats, ...
                                        TF_Chans.OdorOn, TF_Chans.ShamOn);
    % ---------------------------------------------------------------------
    
    
    
    %% Extract statistical information
    %  --------------------------------------------------------------------
    
    % Here, we:
    % 1. Extract information about channel positions for EEGLAB's topoplot 
    %    that will draw cluster-belonging electrodes over the scalp
    % 2. Look if any of the clusters has a significant difference (p <
    %    threshold)
    % 3. Retrieve the cluster labels (1, 2, ...) and the corresponding p
    %    values
    % 4. Go through each cluster found by fieldtrip and
    %           - Plot the channels that belong to cluster
    %           - Plot the time-frequency matrix averaged over these
    %             channels
    %           - Plot the difference between the time-frequency matrices 
    %             of the conditions
    %           - Selectively plot only the part of the "difference matrix"
    %             that belongs to the current cluster
    
    
    
    %% Look for clusters with significant differences between conditions
    %  --------------------------------------------------------------------
    
    % Clusters and probabilities are given by channel by freq by time.
    % We will extract the channels that show significant differences and 
    % plot their average TF. Then, we will overlay the time-frequency 
    % cluster part that shows significant differences.
    
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
    %                             /!\ /!\ /!\
    % Here, better circumvent the use of stats.mask since we will be
    % limited for cluster labels according to set alpha during permutation
    % analysis. Only cluster with p values lower than set alpha would be
    % represented in mask. We get more freedom with pos/negclusterlabelmat.
    % ---------------------------------------------------------------------
    
    if isfield(stats, 'posclusterslabelmat') && isstruct(stats.posclusters)
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
    end
    
    if isfield(stats, 'negclusterslabelmat') && isstruct(stats.negclusters)
        % Negative cluster
        s_raise_label                 = max(cluster_labels(:));
        matrix_raise                  = repmat(s_raise_label, dimensions);
        stats.negclusterslabelmat_mod = ...
            stats.negclusterslabelmat + matrix_raise;
        % cluster_labels(stats.mask) = stats.negclusterslabelmat(stats.mask);
        cluster_labels(stats.negclusterslabelmat_mod ~= matrix_raise) = ...
            stats.negclusterslabelmat_mod(stats.negclusterslabelmat_mod ~= ...
            matrix_raise);
        
        % Extract p values of negative clusters
        idx_labels_neg            = ...
            unique(cluster_labels(cluster_labels > s_raise_label));
        
        idx_labels_neg            = idx_labels_neg - s_raise_label;
        % Point to correct position in neg cluster structure
        
        cluster_pvals             = [cluster_pvals, ...
            stats.negclusters(idx_labels_neg).prob];
    end
    
    
    
    %% Plot multipanel figure for each cluster found inside channel group
    %  --------------------------------------------------------------------
    
    significant_clusters = unique(cluster_labels(:))';
    significant_clusters(significant_clusters == 0) = [];
    
    if numel(cluster_pvals) ~= numel(significant_clusters)
        % Checkpoint
        error('Something went wrong')
    end
    
    
    for i_clust = significant_clusters
        
        
        if cluster_pvals(i_clust) >= s_thrld
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
        true_labels    = true_labels(idx_clustchans, :, :);
        
        
        
        %% Plot topographic distribution of channels and current cluster
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
        no_results = zeros(numel(cluster_chanlocs), 1);
        topoplot(no_results, cluster_chanlocs, ...
            'style', 'blank', ...
            'electrodes', 'pts', ...
            'shading', 'interp', ...
            'headcolor', [0, 0, 0], ...
            'plotchans', idx_defchans, ...
            'conv', 'on', ...
            'emarker', {'.', [0, 0, 0], 5, 1}, ...
            'plotrad', max(max([cluster_chanlocs.radius]),0.5));
        hold on
        topoplot(no_results, cluster_chanlocs, ...
            'style', 'blank', ...
            'electrodes', 'pts', ...
            'shading', 'interp', ...
            'headcolor', [0, 0, 0], ...
            'conv', 'off', ...
            'plotchans', idx_defchans(idx_clustchans), ...
            'emarker', {'.', [0.1, 0.1, 0.1], 20, 1}, ...
            'plotrad', max(max([cluster_chanlocs.radius]), 0.5));
        hold off
        
        title(strrep(str_cluster, '_', ' '))
        
        
        
        %% Plot parameters for time-frequency and wave form
        %  ----------------------------------------------------------------
        TF_Cluster.OdorOn = ...
            TF_Chans.OdorOn.powspctrm(:, idx_clustchans, :, :);
        
        TF_Cluster.ShamOn = ...
            TF_Chans.ShamOn.powspctrm(:, idx_clustchans, :, :);
        
        TF_Cluster.OdorOn = mean(TF_Cluster.OdorOn, 1);
        TF_Cluster.OdorOn = mean(TF_Cluster.OdorOn, 2);
        TF_Cluster.OdorOn = squeeze(TF_Cluster.OdorOn);
        
        TF_Cluster.ShamOn = mean(TF_Cluster.ShamOn, 1);
        TF_Cluster.ShamOn = mean(TF_Cluster.ShamOn, 2);
        TF_Cluster.ShamOn = squeeze(TF_Cluster.ShamOn);
        
        minVal = min([TF_Cluster.OdorOn(:); TF_Cluster.ShamOn(:)]);
        maxVal = max([TF_Cluster.OdorOn(:); TF_Cluster.ShamOn(:)]);
        
        limitsTF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        average_WFCluster.ShamOn = squeeze(mean(mean(...
            WF_subject.ShamOn(:, idx_defchans(idx_clustchans), :), 2), 1));
        average_WFCluster.OdorOn = squeeze(mean(mean(...
            WF_subject.OdorOn(:, idx_defchans(idx_clustchans), :), 2), 1));
        
        minVal = min([average_WFCluster.ShamOn(:); ...
            average_WFCluster.OdorOn(:)]);
        maxVal = max([average_WFCluster.ShamOn(:); ...
            average_WFCluster.OdorOn(:)]);
        
        limitsWF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        % There must be a more efficient way to build the freq vector...
        v_frequencies = 1:dimensions(2);
        v_frequencies = (max(v_frequencies) - v_frequencies) / ...
            (max(v_frequencies) - min(v_frequencies));
        v_frequencies = v_frequencies * PM.FrRange;
        v_frequencies = v_frequencies - (PM.FrRange / 2);
        v_frequencies = flip(v_frequencies);
        
        v_times_TF = PM.cfg_seldat.latency(1) : PM.s_tstep : ...
            PM.cfg_seldat.latency(2);
        
        v_times_WF = PM.cfg_seldat.latency(1) * ...
                PM.Info.TrialParameters.s_fs  + 1 : 1 : ...
                PM.cfg_seldat.latency(2) * ...
                PM.Info.TrialParameters.s_fs;
        v_times_WF = v_times_WF ./ PM.Info.TrialParameters.s_fs;
        
        
        i_plot = 2;
        for condition = PM_stats.Conditions
                        
            subplot(1, 5, i_plot)
            
            
            % Time frequency plot
            % -------------------
                        
            yyaxis left
            pcolor(v_times_TF, v_frequencies, TF_Cluster.(char(condition)));
            shading interp
            colorbar;
            colormap parula
            set(gca, 'clim', limitsTF)
            ylabel('Distance from spindle peak (Hz)')
            xlabel('Time (s)')
            
            str_title = char(extractBefore(condition, 'On'));
            str_title = strrep(str_title, 'Sham', 'Vehicle');
            title(str_title)
            
            hold on
            
            
            % Wave form plot
            % --------------
            
            yyaxis right
            
            WF_meanClust.(char(condition)) = mean(WF_subject.(...
                char(condition))(:, idx_clustchans, :), 2);
            plot(v_times_WF, average_WFCluster.(char(condition)), ...
                'Color',        'r', ...
                'LineWidth',    2)
            ylim(limitsWF)
            
            i_plot = i_plot + 1;
        end
        
        
        
        %% Plot difference between the two conditions' TF
        %  ----------------------------------------------------------------
        
        TF_difference   = ...
            TF_Chans.OdorOn.powspctrm(:, idx_clustchans, :, :) - ...
            TF_Chans.ShamOn.powspctrm(:, idx_clustchans, :, :);
        TF_difference   = squeeze(mean(TF_difference, 1)); % Avg over subj
        
        TF_clustDiff    = squeeze(mean(TF_difference, 1)); % Avg over cluster
        
        minVal = min(TF_clustDiff(:));
        maxVal = max(TF_clustDiff(:));
        
        limitsTF = [- max(abs(minVal), abs(maxVal)), ...
            max(abs(minVal), abs(maxVal))];
        
        
        subplot(1, 5, 4)
        
        figDiff = pcolor(v_times_TF, v_frequencies, TF_clustDiff);
        shading interp
        colorbar;
        set(gca, 'clim', limitsTF)
        ylabel('Distance from spindle peak (Hz)')
        xlabel('Time (s)')
        title('Odor - Vehicle')
        
        
        
        %% Plot time-frequencies of significant difference (cluster)
        %  ----------------------------------------------------------------
        
        subplot(1, 5, 5)

        TF_selective    = NaN(numel(idx_clustchans), ...
            dimensions(2), dimensions(3));
        TF_selective(true_labels == 1) = TF_difference(true_labels == 1);
        TF_selective = squeeze(nanmean(TF_selective, 1));
        
        figSign = pcolor(v_times_TF, v_frequencies, TF_selective);
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
        
        saveas(gcf, [str_folderOut, filesep, ...
            char(PM_stats.Grp_subj), 'Subjects_', ...
            str_cluster, '_', num2str(i_clust), '_Cluster.png']);
        
        close
        
    end
    
    
    %% Save statistical output
    %  --------------------------------------------------------------------
    
    PM.Stats        = PM_stats;
    PM.Cfgs.neighb  = cfg_neighb;
    PM.Cfgs.stats   = cfg_stats;
    save([str_folderOut, filesep, ...
        'PermutationStats_', str_cluster, '_', ...
        char(PM_stats.Grp_subj), 'Subjects.mat'], ...
        'PM', 'stats', '-v7')
    
end
