% ========================   C A U T I O N   ==============================

% THIS SCRIPT IS NOT FUNCTIONAL!
error('Do not use this script')

% =========================================================================

%% Point to files ---------------------------------------------------------

pathSpindlePeaks = ['D:\Gits\SO_Spindle_Detection_Coupling\', ...
    'SubjectSpecific\Max_spindlebands_ByEye.mat'];


% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=   C U E   =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
pathSOInfo = ['D:\germanStudyData\datasetsSETS\Ori_CueNight\', ...
    'preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
    '05-Mar-2021_Cue_NEW\'];

fileSOInfo = '05-Mar-2021_22-23-41_AllData.mat';

pathTimeSeries = ['D:\germanStudyData\datasetsSETS\', ...
    'Ori_CueNight\preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW'];

pathAverageWave = ['D:\germanStudyData\datasetsSETS\Ori_CueNight\', ...
    'preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
    'EventValidation\_EventCharact'];


% =-=-=-=-=-=-=-=-=-=-=-=-=   P L A C E B O   =-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% pathSOInfo = ['D:\germanStudyData\datasetsSETS\Ori_PlaceboNight\', ...
%     'preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
%     '06-Mar-2021_Placebo\'];
% 
% fileSOInfo = '06-Mar-2021_18-08-41_AllData.mat';
% 
% pathTimeSeries = ['D:\germanStudyData\datasetsSETS\', ...
%     'Ori_PlaceboNight\preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW'];
% 
% pathAverageWave = ['D:\germanStudyData\datasetsSETS\Ori_PlaceboNight\', ...
%     'preProcessing\EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
%     'EventValidation\_EventCharact'];



%% Parameters -------------------------------------------------------------

% Clusters of interest
PM.ClustOI     = 'all'; % See "Clust" structure in userland

PM.centerPnt   = 7;
% Data location of latency of the time stamp around which time series will
% be centered. Refer to ReadMe.txt of SO_Spindle_Detection_Coupling code
% 2 = startTime: begin of event
% 3 = midTime:  center of event
% 4 = endTime:  end of event
% 5 = duration: duration from start to end
% 6 = maxTime:  time stamp of maximum peak
% 7 = minTime:  time stamp of mininum peak

PM.timeSize    = 2; % seconds of time series to extract and to take TF of

PM.Conditions  = {'OdorOn', 'ShamOn', 'OdorOff', 'ShamOff'};

PM.s_series_length = 15; % seconds of length of each trial for each 
                         % condition. Careful, loaded trial series are 30s 
                         % trials since they are [Off On].

PM.Spindleband = 'fast'; % Chose whole, slow or fast spindle band to 
                         % extract peak information from

PM.FrRange     = 10; % Hz range for plotting TF. This is not the band edge 
                     % but the range around the spindle peak frequency

PM.s_freqstep  = 0.5;

PM.FrResol     = PM.FrRange / PM.s_freqstep + 1;

compute_TF     = 1; % [0, 1], Select to compute TF (1) or not (0). If not, 
                    % you'll have to load a frequency matrix by hand



%% Prepare userland -------------------------------------------------------

Datasets = dir(pathTimeSeries);

if ~exist('loaded', 'var') || loaded == 0
    load([pathSOInfo, fileSOInfo]);
    load(pathSpindlePeaks);
    AvgData = load(pathAverageWave);
    loaded = 1;
end

% Clusters of interest
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
PM.Clust.right_patietal = {...
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

PM.Clust.all = Info.ROIs.str_chans';

% Parameter-derived variables. Auto-adjusted
Cluster     = PM.Clust.(PM.ClustOI);
Subjects    = Info.Subjects;
PM.s_time   = Info.TrialParameters.s_fs * PM.timeSize;



%% Start ------------------------------------------------------------------
if compute_TF == 1
% We go through:
% 1.   Subjects in order to load trial series and extract information about 
%      spindle bandfrequency peak
% 2.   Conditions in order to extract TF for each condition separately
% 3.   Channels belonging to selected cluster of interest in order to
%      extract and average TF of these --> Cluster TF
% 4.   Trials in order to average later the TF of all trials
% 5.   Slow oscillations in each trial in order to average their TF

% Output is average of channels of average of trials of average of slow
% oscillations stored in TF_cond. TF_cond will be saved so that the last
% part of script (TF normalization) can be executed alone without having to
% go through the time consuming loops.
% Output will also be stored in permutation_analysis-compatible structure


TF_cond     = struct();
for condition = PM.Conditions
    TF_cond.(char(condition)) = [];
    
    % Used for permutation analysis
    TF_channels_subjects.(char(condition)) = {};%...
        %NaN(numel(Subjects), numel(Cluster), PM.FrResol, PM.s_time);

end


for i_subj = 1:numel(Subjects)
    
    subject = Subjects(i_subj);
    
    
    % Get subjects spindle peak to center freq range of TF around
    subject_short    = extractBefore(subject, '_sleep');
    idx_spindle_peak = find(strcmp({spindle_max.subjects}, subject_short));
    
    peakFr           = spindle_max(idx_spindle_peak).(PM.Spindleband);
    
    
    for condition = PM.Conditions
        
        if contains(condition, 'Off')
            conditionShort = extractBefore(condition, 'Off');
        elseif contains(condition, 'On')
            conditionShort = extractBefore(condition, 'On');
        end
        
        % Load dataset: Sham and Odor
        idx_subj = find(contains({Datasets.name}, subject));
        idx_cond = find(contains({Datasets.name}, char(conditionShort)));
        
        idx_data = intersect(idx_subj, idx_cond);
        
        data_cond = load([pathTimeSeries, filesep, Datasets(idx_data).name]);
        
        
        % Initialize channel * trial time series
        TF_channels = NaN(numel(Cluster), PM.FrResol, PM.s_time);
        
        for i_chan = 1:numel(Cluster)
            
            curr_chan   = char(Cluster(i_chan));
            idx_chan    = find(strcmp(Info.ROIs.str_chans, curr_chan));
            
            chanSlowOsc = OverallSlowOsc.(curr_chan).(char(condition))(:, i_subj);
            
            chanTimeSer = data_cond.data;
            
            if isempty(chanSlowOsc{size(chanTimeSer, 3)}) || ...
                    ( size(chanTimeSer, 3) ~= length(chanSlowOsc) && ...
                    ~isempty(chanSlowOsc{size(chanTimeSer, 3)+1}) )
                % Intensive check of data compatibility: time series and 
                % slow osc data
                error('incompatible data')
            else
                s_trials = size(chanTimeSer, 3);
            end
            
            
            
            TF_trials = NaN(s_trials, PM.FrResol, PM.s_time);
            
            for i_trl = 1:s_trials
                
                
                % Time series of trial
                v_series    = chanTimeSer(idx_chan, :, i_trl); 
                
                
                % We extract the array that holds information about
                % detected slow oscillations: SO x 13 matrix, 13 columns
                % containing various information about Slow osc.
                v_SO        = chanSlowOsc{i_trl};
                
                if all(isnan(v_SO))
                    % Leave NaNs in TF_trials for this trial
                    continue
                end
                
                
                % Go through each detected slow oscillation in trial
                eventSpect = NaN(size(v_SO, 1), PM.FrResol, PM.s_time);
                for i_SO = 1:size(v_SO, 1)
                    
                    v_SOseries = NaN(1, PM.s_time);
                    s_stamp = v_SO(i_SO, PM.centerPnt);
                    
                    
                    % Extracted trial series is [Off On] as one trial for
                    % both conditions. Therefore, we adjust the latency of
                    % slow oscillations accordingly
                    if contains(condition, 'On')
                        s_stamp = s_stamp + ...
                            PM.s_series_length * Info.TrialParameters.s_fs;
                    end
                    
                    s_start = s_stamp-PM.s_time/2+1;
                    s_stop  = s_stamp+PM.s_time/2;
                    
                    
                    % =====================================================
                    % Avoid out of bounds by only assigning positive
                    % indices to output vector in case the slow osccilation
                    % window is reaching outside of trial series
                    s_startAdj = 1;
                    s_stopAdj  = PM.s_time;
                    if s_start < 1
                        s_startAdj = 2 - s_start;
                        s_start    = 1;
                    end
                    if s_stop > numel(v_series)
                        s_stopAdj  = PM.s_time - (s_stop - numel(v_series));
                        s_stop     = numel(v_series);
                    end
                    
                    v_SOseries(s_startAdj:s_stopAdj) = ...
                        v_series(s_start:s_stop);
                    
                    v_SOseriesFinite = v_SOseries(~isnan(v_SOseries));
                    % =====================================================

                    
                    
                    % =====================================================
                    %                     OPTIONS 1:
%                     v_freqs = [peakFr-PM.FrRange/2, peakFr+PM.FrRange/2];
%                     
%                     [TF_out, f, ~] = ...
%                         cwt(double(v_SOseriesFinite), ...
%                         'bump', Info.TrialParameters.s_fs, ...
%                         'FrequencyLimits', v_freqs);
                    % =====================================================

                    error('AAAAAHHH FIELDTRIP')
                    % ANDREA
                    Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
                    cfg_Tf                      = [];
                    cfg_Tf.method               = 'wavelet'; %'mtmconvol';
                    cfg_Tf.output               = 'pow';
                    cfg_Tf.foi                  =  0.5:s_fstep:20; % Reduce to 20
                    cfg.width                   = 12;
                    cfg_Tf.toi                  = -25:s_tstep:55;
                    cfg_Tf.keeptrials           = 'yes';

                    
                    % =====================================================
                    %                    OPTION 2:
                    v_freqs = ...
                        peakFr-PM.FrRange/2:PM.s_freqstep:peakFr+PM.FrRange/2;
                    
                    [TF_out, f, ~] = spectrogram(v_SOseriesFinite, ...
                        10, 9, v_freqs, Info.TrialParameters.s_fs);
                    
                    s_stopAdj = s_startAdj+size(TF_out, 2)-1;
                    % =====================================================
                    
                    
                    
                    eventSpect(i_SO, :, s_startAdj:s_stopAdj) = TF_out;
                    
                    t = 1:size(eventSpect, 3);
                    t = PM.s_time/2*(t - PM.s_time/2)/Info.TrialParameters.s_fs;
                    
                end
                
                spectOut = squeeze(nanmean(real(eventSpect), 1));
                TF_trials(i_trl, :, :) = spectOut;
            end
            
            TF_channels(i_chan, :, :) = nanmean(TF_trials, 1);
            
            TF_channels_subjects(i_chan).(char(condition))(i_subj, :, :) = ...
                {squeeze(nanmean(TF_trials, 1))};
        end
        
        TF_cond(i_subj).(char(condition)) = ...
            squeeze(nanmean(TF_channels, 1));
    end
end



%% Save TF extractions ----------------------------------------------------

PM.frequenciesOut = f; % Used to draw frequency values on y axis
save([pathSOInfo, 'TF_matrices'], 'TF_cond', 'TF_channels_subjects', ...
    'Info', 'PM', '-V7');

else
    warning('Think of loading TF matrix before going further')    
end



%% Post-TF parameters -----------------------------------------------------

subjgroups  = {'All'};
Subjects    = Info.Subjects;
t           = 1:PM.s_time;
t           = PM.s_time/2*(t - PM.s_time/2)/Info.TrialParameters.s_fs;
f_plot      = 1:numel(PM.frequenciesOut);
f_plot      = f_plot - numel(PM.frequenciesOut)/2;


%% Select subjects --------------------------------------------------------

% Make it easy and simulate data array since we are only interested in
% v_subjKept
dat_array = ones(1, numel(Subjects));
[dat_array, c_listSubj, v_subjKept] = f_skip_subject(dat_array, ...
    extractBefore(Subjects, '_sleep'), subjgroups);

TF_cond_sel = TF_cond(v_subjKept);



%% Normalize TF -----------------------------------------------------------

TF_subj_sham        = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_sham_off    = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_sham_on     = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_odor_off    = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_odor_on     = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_odor        = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));
TF_subj_normalized  = NaN(PM.FrResol, PM.s_time, length(TF_cond_sel));


figure
for i_subj = 1:length(TF_cond_sel)
   
    % Pseudo-baseline conditions relative to Off period
    TF_cond_sel(i_subj).Odor = ...
        TF_cond_sel(i_subj).OdorOn - TF_cond_sel(i_subj).OdorOff;
    TF_cond_sel(i_subj).Sham = ...
        TF_cond_sel(i_subj).ShamOn - TF_cond_sel(i_subj).ShamOff;
    
    % Substract control conditon from odor
    TF_cond_sel(i_subj).Normalized = ...
        TF_cond_sel(i_subj).Odor - TF_cond_sel(i_subj).Sham;
    
    TF_subj_odor(:, :, i_subj)       = TF_cond_sel(i_subj).Odor;
    TF_subj_sham(:, :, i_subj)       = TF_cond_sel(i_subj).Sham;
    TF_subj_odor_off(:, :, i_subj)   = TF_cond_sel(i_subj).OdorOff;
    TF_subj_sham_off(:, :, i_subj)   = TF_cond_sel(i_subj).ShamOff;
    TF_subj_odor_on(:, :, i_subj)    = TF_cond_sel(i_subj).OdorOn;
    TF_subj_sham_on(:, :, i_subj)    = TF_cond_sel(i_subj).ShamOn;
    TF_subj_normalized(:, :, i_subj) = TF_cond_sel(i_subj).Normalized;
    
    
    
    maxVal = max(max(TF_subj_normalized(:, :, i_subj)));
    minVal = min(min(TF_subj_normalized(:, :, i_subj)));
    
    if abs(minVal) > maxVal
        limits = [minVal, -minVal];
    elseif abs(minVal) < maxVal
        limits = [-maxVal, maxVal];
    else
        limits = [minVal, maxVal];
    end
    
    subplot(4, 6, i_subj)
    pcolor(t, f_plot, TF_subj_normalized(:, :, i_subj));
    shading interp
    colorbar
    set(gca, 'clim', limits([1,end]))
    
end


figure

TF_mean_odor        = squeeze(nanmean(TF_subj_odor, 3));
TF_mean_sham        = squeeze(nanmean(TF_subj_sham, 3));
TF_mean_normalized  = squeeze(nanmean(TF_subj_normalized, 3));


% Center colorbar around zero and normalize scales for plotting conditions 
% next to each other
maxVal = max([TF_mean_sham(:); TF_mean_odor(:)]);
minVal = min([TF_mean_sham(:); TF_mean_odor(:)]);

if abs(minVal) > maxVal
    limits = [minVal, -minVal];
elseif abs(minVal) < maxVal
    limits = [-maxVal, maxVal];
else
    limits = [minVal, maxVal];
end


subplot(1, 3, 1)
pcolor(t, f_plot, TF_mean_sham);
shading interp
colorbar;
set(gca, 'clim', limits([1,end]))
ylabel('Max frequency +/- 2.5 (Hz)')
% yticklabels({'', '-2', '', '-1', '', '0', '', '+1', '', '+2', ''})
xlabel('Time (s)')
title('Vehicle')

subplot(1, 3, 2)
pcolor(t, f_plot, TF_mean_odor);
shading interp
colorbar;
set(gca, 'clim', limits([1,end]))
ylabel('Max frequency +/- 2.5 (Hz)')
% yticklabels({'', '-2', '', '-1', '', '0', '', '+1', '', '+2', ''})
xlabel('Time (s)')
title('Odor')


% Center colorbar around zero for plotting TF difference between conditions
maxVal = max(TF_mean_normalized(:));
minVal = min(TF_mean_normalized(:));

if abs(minVal) > maxVal
    limits = [minVal, -minVal];
elseif abs(minVal) < maxVal
    limits = [-maxVal, maxVal];
else
    limits = [minVal, maxVal];
end

subplot(1, 3, 3)
pcolor(t, f_plot, TF_mean_normalized);
shading interp
colorbar;
set(gca, 'clim', limits([1,end]))
ylabel('Max frequency +/- 2.5 (Hz)')
% yticklabels({'', '-2', '', '-1', '', '0', '', '+1', '', '+2', ''})
xlabel('Time (s)')
title('Odor - Vehicle')



%% Permutation-based statistical analysis ---------------------------------

ST.FrRange = 10; % Frequency range for statistical analysis centered around
                 % spindle peak
                 
FreqMid = ST.FrRange/2;

idx_fr4stat = intersect(find(f_plot > -FreqMid), find(f_plot < FreqMid));     

% Reject unwanted frequencies and subjects and normalize conditions
for i_chan = 1:numel(Cluster)
    
    chan_odoroff        = TF_channels_subjects(i_chan).OdorOff(v_subjKept);
    chan_shamoff        = TF_channels_subjects(i_chan).ShamOff(v_subjKept);
    chan_odoron         = TF_channels_subjects(i_chan).OdorOn(v_subjKept);
    chan_shamon         = TF_channels_subjects(i_chan).ShamOn(v_subjKept);
    
    numSubj             = numel(v_subjKept(v_subjKept == 1));
    c_subj_odor         = cell(1, numSubj);
    c_subj_sham         = cell(1, numSubj);
    c_subj_odor         = zeros(1, numSubj, numel(idx_fr4stat), PM.s_time);
    c_subj_sham         = zeros(1, numSubj, numel(idx_fr4stat), PM.s_time);
    
    for i_subj = 1:numSubj
        
        subj_odoroff    = chan_odoroff{i_subj};
        subj_shamoff    = chan_shamoff{i_subj};
        subj_odoron     = chan_odoron{i_subj};
        subj_shamon     = chan_shamon{i_subj};
        
        subj_odor = ...
            subj_odoron(idx_fr4stat, :) - subj_odoroff(idx_fr4stat, :);
        subj_sham = ...
            subj_shamon(idx_fr4stat, :) - subj_shamoff(idx_fr4stat, :);
        
        c_subj_odor(1, i_subj, :, :) = subj_odor;
        c_subj_sham(1, i_subj, :, :) = subj_sham;
        
    end
    
    TF_channels_subjects(i_chan).Odor = c_subj_odor;
    TF_channels_subjects(i_chan).Sham = c_subj_sham;
end


% Load chanlocs file
load('D:\Gits\Spindle_analysis\subjectsChanlocs\chanlocs.mat'); 

c_channels = Info.ROIs.str_chans;
idx_ROI = [];
for i_chan = 1 :numel(c_channels)
    idx_ROI = ...
        [idx_ROI, find(strcmp({chanlocs.labels}, c_channels(i_chan)))];
end
chanlocs = chanlocs(idx_ROI);


ST.freqdim = 1:numel(idx_fr4stat);
ST.timedim = 1:PM.s_time;

% Adapt subject array to selected subject groups
Info.Subjects = Subjects(v_subjKept);


[stats] = f_permutation_scalp(...
    [{TF_channels_subjects.Odor}', {TF_channels_subjects.Sham}'], ...
    chanlocs, Info.TrialParameters.Conditions, Info, ST);

