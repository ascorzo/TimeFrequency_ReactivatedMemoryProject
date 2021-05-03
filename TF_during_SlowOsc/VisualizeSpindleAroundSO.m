% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

%       - File paths
SOseriespath            = ['D:\germanStudyData\datasetsSETS\', ...
                            'Ori_CueNight\preProcessing\', ...
                            'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
                            '05-Mar-2021_Cue\SO_timeSeries\'];
TFseriespath            = ['D:\germanStudyData\datasetsSETS\', ...
                            'Ori_CueNight\preProcessing\', ...
                            'EEGLABFilt_Mastoids_Off_On_200Hz_Oct_NEW\', ...
                            '05-Mar-2021_Cue\SO_timeSeries\TF_matrices\'];
%       - General
PM.Conditions           = {'ShamOn', 'OdorOn'};
% Off periods ignored since baseline done during TF
PM.ClustOI              = 'all';
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
        
        SOseries = load([SOseriespath, filesSOseries(i_subj).name]);
        % Has to be loaded into new variable, otherwise PM will get
        % overwritten.
        TFseries = load([TFseriespath, filesTFseries(i_subj).name]);
        file_loaded = 1; %#ok<NASGU>
    end
    
    
    
    %% Define channels to go through
    %  --------------------------------------------------------------------
    
    if ~strcmp(PM.ClustOI, 'all')
        Cluster = PM.Clust.(PM.ClustOI);
    else
        Cluster = fieldnames(SOseries.SO_timeSeries)';
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

v_times_WF    = TFseries.PM.cfg_seldat.latency(1) * ...
                  TFseries.PM.Info.TrialParameters.s_fs  + 1 : 1 : ...
                  TFseries.PM.cfg_seldat.latency(2) * ...
                  TFseries.PM.Info.TrialParameters.s_fs;

minTF         = min([average_TF.Sham(:); average_TF.Odor(:)]);
maxTF         = max([average_TF.Sham(:); average_TF.Odor(:)]);
limitsTF      = [-max([abs(minTF), abs(maxTF)]), ...
                    max([abs(minTF), abs(maxTF)])];

minWF         = min([average_WF.Sham(:); average_WF.Odor(:)]);
maxWF         = max([average_WF.Sham(:); average_WF.Odor(:)]);
limitsWF      = [-max([abs(minWF), abs(maxWF)]), ...
                    max([abs(minWF), abs(maxWF)])];
                
% Scatter plot: x are latencies, y are "jittered" values around 0
Yodor = rand(numel(all_Lats.Supersubject.Odor), 1) * ...
    2*limitsWF(2) - limitsWF(2);
Ysham = rand(numel(all_Lats.Supersubject.Sham), 1) * ...
    2*limitsWF(2) - limitsWF(2);

figure


% Time-frequency Sham
% -------------------
subplot(2, 2, 1)

yyaxis left
pcolor(v_times_TF, v_frequencies, average_TF.Sham);
shading interp
colorbar;
set(gca, 'clim', limitsTF)
ylabel('Distance from spindle peak (Hz)')
xlabel('Time (samples)')

hold on

% Wave form Sham
% --------------
yyaxis right
plot(v_times_WF, average_WF.Sham, ...
    'Color',        [1, 0, 0], ...
    'LineWidth',    2)
ylim(limitsWF)

% Overlay spindle latencies from supersubject
% -------------------------------------------
scatter(all_Lats.Supersubject.Sham, Ysham, 0.5, '.', 'k')

title('Vehicle')

% Spindle latencies Sham
% ----------------------
subplot(2, 2, 3)

yyaxis left
yticks([])
yyaxis right
boxplot(all_Lats.Sham, 'Orientation', 'horizontal', ...
    'OutlierSize', 2, 'Symbol', '.')
xlim([-400, 400])
xlabel('Spindle latencies (samples)')
yticklabels({filesSOseries.name})


% Time-frequency Odor
% -------------------
subplot(2, 2, 2)

yyaxis left
pcolor(v_times_TF, v_frequencies, average_TF.Odor);
shading interp
colorbar;
set(gca, 'clim', limitsTF)
ylabel('Distance from spindle peak (Hz)')
xlabel('Time (samples)')

hold on

% Wave form Odor
% --------------
yyaxis right
plot(v_times_WF, average_WF.Odor, ...
    'Color',        [1, 0, 0], ...
    'LineWidth',    2)
ylim(limitsWF)

% Overlay spindle latencies from supersubject
% -------------------------------------------
scatter(all_Lats.Supersubject.Odor, Yodor, 0.5, '.', 'k')

title('Cue')


% Spindle latencies Odor
% ----------------------
subplot(2, 2, 4)

yyaxis left
yticks([])
yyaxis right
boxplot(all_Lats.Odor, 'Orientation', 'horizontal', ...
    'OutlierSize', 2, 'Symbol', '.')
xlim([-400, 400])
xlabel('Spindle latencies (samples)')
yticklabels({filesSOseries.name})

