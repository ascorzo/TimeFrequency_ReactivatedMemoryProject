% Compute time-frequency around subject's maximum spindle band peak

% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   N O T E !   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% Events' TF are averaged (post-TF computation) for each channel. If you
% want to maintain individual event TF matrices for each channel, comment
% out line "data_TF_norm.powspctrm  = mean(data_TF_norm.powspctrm, 1);"


% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% 1. File containing information about maximum spindle peaks
%
% 2. Slow oscillation time series as prepared by
%    DataPrep_extract_event_series.m
%
% 3. Parameters
%       - General
PM.Conditions           = {'ShamOn', 'OdorOn'};
% Off periods ignored since baseline done during TF
PM.ClustOI              = 'all';
% PM.ClustOI              = 'frontal';
% See Clust list below. Can be 'all' and should be left 'all' since
% channels of interest can be redefined later in other scripts
PM.Spindleband          = 'fast';
PM.FrRange              = 10; % Hz, range of frequencies to take into
                              % account around spindle band peak 
PM.FrCenter             = 'mean'; % Whether to center around 
                                            % subject's "individual" or 
                                            % "mean" peak of all subjects
%       - Time frequency parameters
PM.s_tstep              = 0.05; % FASTER!
% Probably a good idea to have highest time resolution possible since the
% shift in prefered angles of spindles during SO is rather small.
PM.s_fstep              = 0.05; % 0.005
PM.cycles               = 12;
PM.cfg_Tf.method        = 'wavelet';
PM.cfg_Tf.output        = 'pow';
PM.cfg_Tf.width         = PM.cycles;
PM.cfg_Tf.keeptrials    = 'yes';
%       - Data window selection
PM.cfg_seldat.latency   = [-2.25 1.75]; % [-2 2], [-2.25 1.75]
%       - Baseline parameters
PM.cfg_Bas.baseline      = PM.cfg_seldat.latency;
% Purposely choosing the whole window as baseline as we must assume that
% most SOs are preceeded by another SO (S4) and therefore, baseline 
% calculation with a time window before the SO might neutralize desired
% effects.
PM.cfg_Bas.baselinetype  = 'zscore';
%       - File paths
filepath = ['D:\germanStudyData\datasetsANDREA\DetectionNative\Cue\', ...
    '18-Mar-2022_Jens_Cue\SO_timeSeries\'];
savepath = strcat(filepath, 'TF_matrices');
peakpath = ['D:\Gits\SO_Spindle_Detection_Coupling\', ...
           'SubjectSpecific\Max_spindlebands_byEye.mat'];
%       - paths to toolboxes
fieldtrippath       = 'D:\MATLAB\fieldtrip-20200831';
%       - path to f_chan_clusters function: All this function does is
%         create cell arrays of electrode labels belonging to a group. The
%         repository containing the function can be found at
%         https://github.com/davidmarcelbaum/EEG_channels
chanclusterpath     = 'D:\Gits\EEG_channels';



% -------------------------------------------------------------------------
% --------------------------                 ------------------------------
% --------------------------   O U T P U T   ------------------------------
% --------------------------                 ------------------------------
% -------------------------------------------------------------------------

% 1. TF_series
%    Contains the time-frequency matrices and channel information
%
% 2. SO_wave
%    These are the average wave forms by channel
%
% 3. PM
%    Used parameters



%% Set up userland
%  ------------------------------------------------------------------------

files = dir(strcat(filepath,'*.mat'));

load(peakpath);

if exist(savepath, 'dir') ~= 7
    mkdir(savepath)
end

% Clusters of interest
% --------------------
addpath(chanclusterpath)
PM.Clust = f_chan_clusters;

if ~strcmp(PM.ClustOI, 'all')
    Cluster = PM.Clust.(PM.ClustOI);
    % Otherwise, will be defined later when a subject file has been loaded
    % containing channel information.
end

addpath(fieldtrippath)
ft_defaults
ft_warning off



%% Fruitloops
%  ------------------------------------------------------------------------

for i_subj = 1:numel(files)
    
    % Load Data
    if ~exist('file_loaded', 'var') || file_loaded == 0
        SOseries = load([filepath, files(i_subj).name]);
        % Has to be loaded into new variable, otherwise PM will get
        % overwritten.
        PM.filename = [filepath, files(i_subj).name];
        file_loaded = 1; %#ok<NASGU>
    end
    
    
    % If ClusterOI set to all, we define this here
    if strcmp(PM.ClustOI, 'all')
        Cluster = fieldnames(SOseries.SO_timeSeries)';
    end
    
    
    % Determine frequency range to analyze for subject
    % Get subjects spindle peak to center freq range of TF around
    subject_short    = extractBefore(files(i_subj).name, '_sleep');
    idx_spindle_peak = find(strcmp({spindle_max.subjects}, subject_short));
    
    if strcmp(PM.FrCenter, 'individual')
        PM.peakFr        = spindle_max(idx_spindle_peak).(PM.Spindleband);
    elseif strcmp(PM.FrCenter, 'mean')
        PM.peakFr        = mean([spindle_max.(PM.Spindleband)]);
                                % Wavelet-based TF has different
                                % resolutions for different frequency bands
                                % and we might biais the outcome by
                                % centering the TF matrices around
                                % subjects' individual spindle fr peaks. We
                                % circumvent this here by centering around 
                                % the grand average of spindle peaks
    end
    
    v_freqs          = ...
        PM.peakFr-PM.FrRange/2:PM.s_fstep:PM.peakFr+PM.FrRange/2;
    
    
    for condition = PM.Conditions
        
        
        % We will store the average wave by channel of cluster
        SO_wave_chans       = NaN(numel(Cluster), ...
                              diff(PM.cfg_seldat.latency * ...
                              SOseries.PM.Info.TrialParameters.s_fs));
        TF_condition        = struct(); % By far the heaviest variable here
        TF_condition_restr  = struct();
        
        % Used for determining the size of power spectrum matrix when
        % catenating TF matrices
        numTrials           = NaN(1, numel(Cluster));
        
        
        for i_chan = 1:numel(Cluster)
            % There was probably a better way to prepare data, but the way
            % it has been done (*) forces us to do by channel separately.
            % Reason for this was that each channel would have different
            % numbers of trials. Does this interfere during TF?
            % (* every condition and every channel in separate field)
            
            channel             = char(Cluster(i_chan));
            
            disp(strcat('Subject:',     {' '},    files(i_subj).name))
            disp(strcat('Condition:',   {' '},    condition))
            disp(strcat('Channel:',     {' '},    channel))
            
            
            % Select data channels according to condition and channel
            data_raw = SOseries.SO_timeSeries.(channel).(char(condition));
            
            
            
            %% Computing and processing time-frequency spectrum
            %  ------------------------------------------------------------
            
            % Time-Frequency Calculation
            PM.cfg_Tf.foi       = v_freqs;
%             PM.cfg_Tf.toi       = 'all';
            PM.cfg_Tf.toi       = -SOseries.PM.s_timeWindow / 2:...
                                    PM.s_tstep:...
                                    SOseries.PM.s_timeWindow / 2;
                                 % 'all', Extended before and after to deal 
                                 % with border effect of wavelet
            PM.cfg_Tf.feedback  = 'off'; % Speed up
            data_TF             = ft_freqanalysis(PM.cfg_Tf, data_raw);
            
            % Select time of interest removing borders
            data_TF_sel         = ft_selectdata(PM.cfg_seldat, data_TF);
            
            % Baseline correction
            data_TF_norm        = ft_freqbaseline(PM.cfg_Bas, data_TF_sel);
            
            
            % Store results in channel structure
            % -------------------------------------------------------------
            % You have the option to command the line out in order to
            % preserve TF matrices of all trials.
            % Without meaning, script gave OOM errors after some time since
            % TF_condition will be around 14GB! We have no choice but to 
            % mean the events' TF matrices at this stage already.
%             data_TF_norm.powspctrm  = mean(data_TF_norm.powspctrm, 1);
            % -------------------------------------------------------------
            
            TF_condition.(channel)  = data_TF_norm;
            
            
            % Used for creating catenated data. If meaned in line 229, all 
            % channels will be 1 in trials.
            numTrials(i_chan) = size(data_TF_norm.powspctrm, 1);
            
            
            
            %% Extract average wave form of slow osc. of channel
            %  ------------------------------------------------------------
            
            SO_wave_trl =   NaN(numel(data_raw.trial), ...
                            SOseries.PM.s_timeWindow * ...
                            SOseries.PM.Info.TrialParameters.s_fs);
            
            for i_trl = 1:numel(data_raw.trial)
                SO_wave_trl(i_trl, :) = data_raw.trial{i_trl};
            end
            
            % -------------------------------------------------------------
            % The window of the wave form was the length of the raw data
            % which had extended times at both orders in order to
            % circumvent TF border effects. We need to strip these borders
            % from the wave form time series here.
            SO_wave_wholeWindow = mean(SO_wave_trl, 1);
            
            v_wave_times_whole = ...
                - length(SO_wave_wholeWindow) / 2 + 1 : 1 : ...
                length(SO_wave_wholeWindow) / 2;
            
            v_wave_times = ...
                PM.cfg_seldat.latency(1) * ...
                SOseries.PM.Info.TrialParameters.s_fs + 1 : 1 : ...
                PM.cfg_seldat.latency(2) * ...
                SOseries.PM.Info.TrialParameters.s_fs;
            % -------------------------------------------------------------
            
            SO_wave_chans(i_chan, :) = SO_wave_wholeWindow(...
                ismember(v_wave_times_whole, v_wave_times));
            
        end
        
        
        %% Restructure data so that channels are inside same array
        %  ----------------------------------------------------------------
        
        % This allows for easier data preparation when calling
        % permutation-based stats.
        
        % TF_condition_restr
        %   cfg             = various information about data processing
        %   label           = Cell(channel labels x 1)
        %   freq            = 1 x frequency vector output
        %   time            = 1 x time vector
        %   dimord          = 'rpt_chan_freq_time'
        %   elec            IGNORED AS LOADED IN SUBSEQUENT SCRIPT
        %       chanpos     = coordinates of channel [channels by x, y, z]
        %       chantype    = repmat({'eeg'}, numel(Cluster), 1)
        %       chanunit    = repmat({'V'}, numel(Cluster), 1)
        %       elecpos     = TF_condition_restr.elec.chanpos
        %       label       = TF_condition_restr.label' (permutated)
        %       unit        = 'cm'
        %   trialinfo       = Event triggers (vary between channels.
        %                     Therefore, take this field out)
        %   powspctrm       = rpt_chan_freq_time matrix, this is where we
        %                     canetate trials in dimension 2
        
        
        TF_condition_restr.label  = Cluster';
        TF_condition_restr.cfg    = TF_condition.(char(Cluster(1))).cfg;
        TF_condition_restr.freq   = TF_condition.(char(Cluster(1))).freq;
        TF_condition_restr.time   = TF_condition.(char(Cluster(1))).time;
        TF_condition_restr.dimord = TF_condition.(char(Cluster(1))).dimord;
        
        TF_condition_restr.powspctrm =  NaN(...
                                        max(numTrials), ...
                                        numel(Cluster), ...
                                        numel(TF_condition_restr.freq), ...
                                        numel(TF_condition_restr.time));
%         TF_condition_restr.elec = [];
        for i_chan = 1:numel(Cluster)
            TF_condition_restr.powspctrm(...
                1:numTrials(i_chan), i_chan, :, :) = ...
                TF_condition.(char(Cluster(i_chan))).powspctrm;
            
%             TF_condition_restr.elec = [TF_condition_restr.elec; ...
%                 TF_condition.(char(Cluster(i_chan))).elec];
        end

        % -----------------------------------------------------------------
        % ------------------   I M P O R T A N T   ------------------------
        %
        % Should we randomly select trials and balance the trial number 
        % over channels? Right now, every channel has its own number of 
        % trials based on how many slow osc. were detected.
        % Does not seem necessary, since what we feed to ft_freqstatistics
        % are the subjects' means over their respective trials for each 
        % channel. Therefore, the dimension rpt in the input array for 
        % ft_statistics will be subjects (x channels x freqs x times) where
        % the initially different amount of trials per subject does not
        % matter since we will mean over them to create a subject rpt.
        % -----------------------------------------------------------------
        
        SO_wave.(char(condition))   = SO_wave_chans;
        TF_series.(char(condition)) = TF_condition_restr;
    end

    
    
    %% Parsing last parameters to output and save results for subjects
    %  --------------------------------------------------------------------
    
    PM.datapreparation  = SOseries.PM;
    PM.Clust.all        = SOseries.PM.Info.ROIs.str_chans;
    PM.Info             = SOseries.PM.Info;
    save([savepath, filesep, subject_short, '_TF.mat'], ...
        'TF_series', 'SO_wave', 'PM', '-v7.3')
    % Has to be saved as Version 7.3 because bigger than 2GB. This means
    % generated files will not be Octave-compatible.
    
    file_loaded = 0; % Load next subject file
end

% Had a weird error once atthe very end of the script. Make sure here that
% everything went through!
fprintf('\nScript is done!\n')
