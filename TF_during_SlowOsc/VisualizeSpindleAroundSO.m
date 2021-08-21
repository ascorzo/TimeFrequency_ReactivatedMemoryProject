% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

%       - File paths
SOseriespath            = ['D:\germanStudyData\datasetsSETS\', ...
                            'Ori_CueNight\preProcessing\', ...
                            'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
                            '12-Jun-2021_Cue\SO_timeSeries_Upstate_15s\'];
% SOseriespath            = ['D:\germanStudyData\datasetsSETS\', ...
%                             'Ori_PlaceboNight\preProcessing\', ...
%                             'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
%                             '13-Jun-2021_Placebo\SO_timeSeries_Upstate_15s\'];
TFseriespath            = [SOseriespath, '\TF_matrices\'];
Channelpath             = 'D:\Gits\EEG_channels';
%       - General
PM.Conditions           = {'ShamOn', 'OdorOn'};
% Off periods ignored since baseline done during TF
PM.ClustOI              = 'central';
%       - Center frequency to filter spindle band
centerFr                = 13.7947; % Avg of all subjets
%       - Data window selection
PM.cfg_seldat.latency   = [-2 2];



% -------------------------------------------------------------------------
% --------------------------                 ------------------------------
% --------------------------   O U T P U T   ------------------------------
% --------------------------                 ------------------------------
% -------------------------------------------------------------------------

% Structure "out" with fields (subfields are conditions):
%   1.  SO_wave: Contains wave form of Slow osc. averaged over channels of
%       interest. First dimension is subjects
%   2.  SO_tf: Time-frequency matrices averaged over channels. First
%       dimension is subjects
%   3.  SS_lats: Cells of latencies of spindles that occured between start 
%       and end of slow osc.. Latencies are catenated over all channels of
%       interest. Each cell is a subject.



%% Prepare userland
%  ------------------------------------------------------------------------

filesSOseries = dir(strcat(SOseriespath,'*.mat'));
filesTFseries = dir(strcat(TFseriespath,'*.mat'));
addpath(Channelpath)

% Clusters of interest
% --------------------
PM.Clust = f_chan_clusters;


%% Fruitloops
%  ------------------------------------------------------------------------

% We will for each subject and each condition:
% 1. Extract the average wave form over channels of selected cluster
% 2. Extract the spindle latencies for the channels of interest
% 3. Extract the average time-frequency matrix over channels of interest

for i_subj = 1:numel(filesSOseries)
    
    disp(strcat('Subject:', {' '}, filesSOseries(i_subj).name))
    
    
    
    %% Load Data
    %  --------------------------------------------------------------------
    
    if ~exist('file_loaded', 'var') || file_loaded == 0
        % Check for compatible datasets
        if ~strcmp(extractBefore(filesTFseries(i_subj).name, '_TF'), ...
                extractBefore(filesSOseries(i_subj).name, '_sleep'))
            error('/!\ You were about to load incompatible datasets')
        end
        
        % Have to be loaded into new variables, otherwise PM will get
        % overwritten.
        SOseries = load([SOseriespath, filesSOseries(i_subj).name]);
        TFseries = load([TFseriespath, filesTFseries(i_subj).name]);
        file_loaded = 1;
    end
    
    
    
    %% Define channels to go through
    %  --------------------------------------------------------------------
    
    if strcmp(PM.ClustOI, 'all')
        Cluster = fieldnames(SOseries.SO_timeSeries)';
    else
        Cluster = PM.Clust.(PM.ClustOI);
    end
    
    
    
    for condition = PM.Conditions
        
        %% Prepare outputs
        %  ----------------------------------------------------------------
        
        % We will catenate spindle latencies of all channels since the
        % number of detected spindles that are coupled to slow oscillations
        % is really small (between 5 and 10 most of time)
        latencies_allChans  = [];
        
        TF_chans            = NaN(numel(Cluster), ...
            numel(TFseries.TF_series.(char(condition)).freq), ...
            numel(TFseries.TF_series.(char(condition)).time));
        
        % We will store the average wave by channel of cluster
        SO_wave_chans       = NaN(numel(Cluster), ...
            diff(PM.cfg_seldat.latency * ...
            SOseries.PM.Info.TrialParameters.s_fs));
        
        
        
        for i_chan = 1:numel(Cluster)
            
            
            channel     = char(Cluster(i_chan));
            
            idx_chan    = find(strcmp(...
                TFseries.PM.Info.ROIs.str_chans, channel));
            
            
            
            %% Information about data series
            %  ------------------------------------------------------------
            
            % Select data channels according to condition and channel
            data_raw = SOseries.SO_timeSeries.(channel).(char(condition));
            
            
            
            %% Extract average wave form of slow osc. of channel
            %  ------------------------------------------------------------
            
            SO_wave_trl = NaN(numel(data_raw.trial), ...
                SOseries.PM.s_timeWindow * ...
                SOseries.PM.Info.TrialParameters.s_fs);
            
            for i_trl = 1:numel(data_raw.trial)
                
%                 signal      = data_raw.trial{i_trl};
%                 min_freq    = centerFr - 1.5;
%                 max_freq    = centerFr + 1.5;
%                 srate       = SOseries.PM.Info.TrialParameters.s_fs;
%                 
%                 [FilteredEEGEnvelope, FilteredEEG, ~] = ...
%                     f_fir2(signal, srate, min_freq, max_freq);
%                 SO_wave_trl(i_trl, :) = [FilteredEEG', NaN];
                

                signal      = double(data_raw.trial{i_trl});
                filt_order  = SOseries.PM.Info.SOparameters.butter_order;
                min_freq    = SOseries.PM.Info.SOparameters.HighPassFr;
                max_freq    = SOseries.PM.Info.SOparameters.LowPassFr;
                srate       = SOseries.PM.Info.TrialParameters.s_fs;
                
                FilteredEEG = f_filter_deltaband(signal, filt_order, ...
                    srate, max_freq, min_freq);
                SO_wave_trl(i_trl, :) = FilteredEEG;
            end
            
            % -------------------------------------------------------------
            % The window of the wave form was the length of the raw data
            % which had extended times at both orders in order to
            % circumvent TF border effects. We need to strip these borders
            % from the wave form time series here.
            SO_wave_wholeWindow = nanmean(SO_wave_trl, 1);
            
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
            
            
            
            %% Extract spindle latencies around slow oscillations
            %  ------------------------------------------------------------
            
            latencies_allChans = [latencies_allChans, ...
                SOseries.SS_latencies.(channel).(char(condition))];
            
            
            
            %% Extract time-frequency matrices
            %  ------------------------------------------------------------
            
            TF_chans(i_chan, :, :) = squeeze(...
                nanmean(TFseries.TF_series.(...
                char(condition)).powspctrm(:, idx_chan, :, :), 1));
            
        end
        
        out.SO_wave.(char(condition))(i_subj, :)  = ...
            squeeze(mean(SO_wave_chans, 1));
        out.SO_tf.(char(condition))(i_subj, :, :) = squeeze(...
            mean(TF_chans, 1));
        out.SS_lats.(char(condition)){i_subj}     = latencies_allChans;
        
    end
    
    file_loaded = 0; % Load next subject file
    
end



%% Plots
%  ------------------------------------------------------------------------

% We will plot a two-panel figure with left being Sham and right being Odor
% We will plot average TF matrix with average wave form as overlay. Below,
% we will plot subjects individual spindle latencies as scatter and
% boxplot.


% Prepare time-frequency matrices and wave forms
% ----------------------------------------------
average_TF.Sham = squeeze(mean(out.SO_tf.ShamOn, 1));
average_TF.Odor = squeeze(mean(out.SO_tf.OdorOn, 1));
average_WF.Sham = squeeze(mean(out.SO_wave.ShamOn, 1));
average_WF.Odor = squeeze(mean(out.SO_wave.OdorOn, 1));


% Prepare spindle time stamps
% ---------------------------
% Boxplots work by column. We transform latencies into a matrix where each
% column represents a subject.
sizes_lats      = cell2mat(cellfun(@size, out.SS_lats.ShamOn, ...
                    'UniformOutput', false));

all_Lats.Sham   = NaN(max(sizes_lats), numel(out.SS_lats.ShamOn));

for i_subj = 1:numel(out.SS_lats.ShamOn)
    num_lats = size(out.SS_lats.ShamOn{i_subj}, 2);
    all_Lats.Sham(1:num_lats, i_subj) = out.SS_lats.ShamOn{i_subj}';
end

sizes_lats      = cell2mat(cellfun(@size, out.SS_lats.OdorOn, ...
                    'UniformOutput', false));

all_Lats.Odor   = NaN(max(sizes_lats), numel(out.SS_lats.OdorOn));

for i_subj = 1:numel(out.SS_lats.OdorOn)
    num_lats = size(out.SS_lats.OdorOn{i_subj}, 2);
    all_Lats.Odor(1:num_lats, i_subj) = out.SS_lats.OdorOn{i_subj}';
end

all_Lats.Supersubject.Sham = all_Lats.Sham(:);
all_Lats.Supersubject.Odor = all_Lats.Odor(:);


% Combine conditions
% ------------------
conditions_TF(:, :, 1) = average_TF.Sham;
conditions_TF(:, :, 2) = average_TF.Odor;

all_Lats.Condition     = [all_Lats.Sham; all_Lats.Odor];

all_Lats.Supersubject.Cond = [all_Lats.Sham(:); all_Lats.Odor(:)];
average_TF.Cond            = mean(conditions_TF, 3);

average_WF.Cond        = mean([average_WF.Sham; average_WF.Odor], 1);


% Prepare plot parameters
% -----------------------
num_freqs     = size(average_TF.Sham, 1);
num_times     = size(average_TF.Sham, 2);

v_frequencies = 1:num_freqs;
v_frequencies = (max(v_frequencies) - v_frequencies) / ...
                    (max(v_frequencies) - min(v_frequencies));
v_frequencies = v_frequencies * TFseries.PM.FrRange;
v_frequencies = v_frequencies - (TFseries.PM.FrRange / 2);
v_frequencies = flip(v_frequencies);

v_times_TF    = 1:num_times;
v_times_TF    = v_times_TF - (num_times / 2);

v_times_WF    = -2 * TFseries.PM.Info.TrialParameters.s_fs  + 1 : 1 : ...
                  2 * TFseries.PM.Info.TrialParameters.s_fs;
              
upscale_times = numel(v_times_WF) / num_times;

v_times_TF    = v_times_TF * upscale_times;

minTF         = min([average_TF.Sham(:); average_TF.Odor(:)]);
maxTF         = max([average_TF.Sham(:); average_TF.Odor(:)]);
limitsTF      = [-max([abs(minTF), abs(maxTF)]), ...
                    max([abs(minTF), abs(maxTF)])];

% minWF         = min([average_WF.Sham(:); average_WF.Odor(:)]);
% maxWF         = max([average_WF.Sham(:); average_WF.Odor(:)]);
minWF         = min(average_WF.Cond(:));
maxWF         = max(average_WF.Cond(:));
limitsWF      = [-max([abs(minWF), abs(maxWF)])-5, ...
                    max([abs(minWF), abs(maxWF)])+5];
limitsFR      = [min(v_frequencies), max(v_frequencies)];
                
% Scatter plot: x are latencies, y are "jittered" values around 0
% Yodor = rand(numel(all_Lats.Supersubject.Odor), 1) * ...
%     2*limitsWF(2) - limitsWF(2);
% Ysham = rand(numel(all_Lats.Supersubject.Sham), 1) * ...
%     2*limitsWF(2) - limitsWF(2);
% Ycond = rand(numel(all_Lats.Supersubject.Cond), 1) * ...
%     2*limitsWF(2) - limitsWF(2);
Ycond = rand(numel(all_Lats.Supersubject.Cond), 1) * ...
    2*limitsFR(2) - limitsFR(2);



figure('Units', 'Normalized', 'Position', [0 0 0.5 0.5])


% Time-frequency Conditions
% -------------------------
subplot(2, 1, 1)

yyaxis left
pcolor(v_times_TF-50, v_frequencies, average_TF.Cond);
shading interp
colorbar('Location', 'WestOutside')
set(gca, 'clim', limitsTF)
ylabel('Distance from spindle peak (Hz)')
xlabel('Time (samples)')
xlim([-350, 250])

hold on

yyaxis right
plot(v_times_WF, average_WF.Cond, ...
    'Color',        [1, 0, 0], ...
    'LineWidth',    2)
ylim(limitsWF)

xlim([-350, 250])
xticks([])
ylabel('Amplitude (uv)')


% Overlay spindle latencies from supersubject
% -------------------------------------------
% scatter(all_Lats.Supersubject.Cond, Ycond, 0.5, '.', 'k')


% Spindle latencies Sham
% ----------------------
subplot(2, 1, 2)

y_subjects = -11.5:1:11.5;
y_subjects = y_subjects * 3 + 12;

hold on
for i_subj = 1:size(all_Lats.Condition, 2)
    
    spindle_num = numel(all_Lats.Condition(:, i_subj));
    
    jitterSpan = 25;
    jitterVals = randi([-jitterSpan, jitterSpan], ...
        1, spindle_num);
    jitterVals = jitterVals ./ 200;
    
    y_sub = y_subjects(i_subj);
    
%     scatter(all_Lats.Condition(:, i_subj), ...
%         repmat(i_subj, 1, spindle_num)+jitterVals, 0.5, '.', 'k')
    scatter(all_Lats.Condition(:, i_subj), ...
        repmat(y_sub, 1, spindle_num)+jitterVals, 0.5, '.', 'k')
    
end
% boxplot(all_Lats.Condition, 'Orientation', 'horizontal', ...
%     'OutlierSize', 2, 'Symbol', '.', 'PlotStyle', 'compact', ...
%     'BoxStyle', 'outline', 'Colors', [0, 0, 1], 'MedianStyle', 'line', ...
%     'Jitter', 0)
ylim([0, size(all_Lats.Condition, 2)+1])
% yticklabels({filesSOseries.name})
% hold on

% % Wave form Conditions
% % --------------------
% yyaxis right
% ylabel('Amplitude (uv)')
% plot(v_times_WF, average_WF.Cond, ...
%     'Color',        [0, 0, 0], ...
%     'LineWidth',    2)
% ylim(limitsWF)

xlim([-350, 250])
% xticks([])
ylabel('Subjects')
% xlabel('Spindle latencies (samples)')

