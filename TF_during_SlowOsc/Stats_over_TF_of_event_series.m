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
%       - Conditions to analyze
PM.Conditions       = {'ShamOn', 'OdorOn'};
%       - channel locations in fieldtrip structure
chanlocspath        = ['D:\Gits\Spindle_analysis\subjectsChanlocs\', ...
                        'fieldtrip_chanlocs.mat'];



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
    WF_subject.(char(condition)) = NaN(length(files), s_windowSize);
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
        
            % /!\ The time series of the wave forms has the borders still
            % attached that were used to circumvent border effects in TF
            % computation. We hard-code the time selection here. This has 
            % been corrected in the TF computation script and therefore, 
            % v_times here can be erased in future.
            
            warning(['Andrea, if you have run the Compute_TF_of_event', ...
                '_series script again, you can remove this time', ...
                'vector selection here and just take SO_wave(idx_chan, :)'])
            pause(100)
            
            v_times = 201:1:1000;
            WF_subject.(char(condition))(i_subj, :) = ...
                SO_wave.(char(condition))(idx_chan, v_times);
            
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


figure
colororder({'b', 'k'})
i_plot = 1;
for condition = PM.Conditions
    
    subplot(1, 2, i_plot)
    
    % Time frequency plot
    % -------------------
    TF_meanChans  = mean(TF_subject.(char(condition)).powspctrm, 2);
    
    v_frequencies = 1:length(TF_series.OdorOn.freq);
    v_frequencies = (max(v_frequencies) - v_frequencies) / ...
        (max(v_frequencies) - min(v_frequencies));
    v_frequencies = v_frequencies * PM.FrRange;
    v_frequencies = v_frequencies - (PM.FrRange / 2);
    v_frequencies = flip(v_frequencies);
    
    v_times       = 1:size(TF_meanChans, 4);
    v_times       = v_times - (size(TF_meanChans, 4) / 2);
    
    yyaxis left
    pcolor(v_times, v_frequencies, squeeze(mean(TF_meanChans, 1)));
    shading interp
    colorbar;
    % set(gca, 'clim', [minVal, maxVal])
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
    plot(v_times, mean(WF_subject.(char(condition)), 1), ...
        'Color',        [0, 0, 0], ...
        'LineWidth',    2)
    
    i_plot = i_plot + 1;
end



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



%% Parameters for stats
%  ------------------------------------------------------------------------

ft_defaults

cfg_stats                       = [];
%cfg_stats.latency               = [1]; % Setting this to one rejects the
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
% cfg.alpha               = 0.05; % = False detection rate ?
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
cfg_stats.design                = [1:numel(files), ...
                                    1:numel(files); ...
                                    repmat(1, 1, numel(files)), ...
                                    repmat(2, 1, numel(files))];

cfg_stats.feedback = 'off';
stats  = ft_freqstatistics(cfg_stats, ...
    TF_subject.OdorOn, TF_subject.ShamOn);



%% Extract statistical information /!\ under construction /!\
%  ------------------------------------------------------------------------

% Clusters and probabilities are given by channel by freq by time.
% We will extract the channels that show significant differences and plot
% their average TF. Then, we will overlay the time-frequency cluster part
% that shows significant differences.

% Look for clusters with probability of <= 0.2 (permissive threshold)
p_all = stats.prob(:);

