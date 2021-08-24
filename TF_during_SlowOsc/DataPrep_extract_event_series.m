% Script that extracts time windows around event center points.

% -------------------------------------------------------------------------
% ---------------------------               -------------------------------
% ---------------------------   I N P U T   -------------------------------
% ---------------------------               -------------------------------
% -------------------------------------------------------------------------

% 1. EEGLAB-structured datasets:
%       - Continuous / Whole recordings
%       - Trial-baseed / Epoched extracted time series
%
% 2. An event detection file that contains information about latencies of
%    detected events
%
%       Continuous time series and trial time series
%       where the latter have been extracted from the continuous ones. 
%       Continuous datasets contain non-valid trials and do not contain 
%       latencies in a by-trial manner while Trial based datasets contain 
%       only valid triggers. Event detection files contain latencies for 
%       events that are comprised between 0 and end latencies of each trial
%       where they have been detected, hence the importance of Trial-based
%       datasets.

% 3. Parameters:
%    - File paths
filePath_CONTINUOUS = ['D:\germanStudyData\datasetsSETS\Ori_CueNight\', ...
                        'preProcessing\NREM'];
filePath_TRIALS     = ['D:\germanStudyData\datasetsSETS\Ori_CueNight\', ...
                        'preProcessing\TRIALS'];
eventfilePath       = ['D:\germanStudyData\datasetsSETS\Ori_CueNight\', ...
                        'preProcessing\EEGLABFilt_Mastoids_Off_On_', ...
                        '200Hz_Oct_NEW\12-Jun-2021_Cue\'];
eventFile           = '12-Jun-2021_17-08-46_AllData.mat';
% eventFile           = '13-Jun-2021_15-33-48_AllData.mat';
%    - Channels to include
chanOfInterest      = 'all';
%    - Conditions to include
PM.Conditions       = {'OdorOn', 'ShamOn'};
%    - Size of desired time window to extract
PM.s_timeWindow     = 6;
%    - Time stamps of start and end of slow oscillations
PM.SOlocations      = [2, 4];
%    - Time stamp to center time window around: Can be based on angles from
%      hilbert transform or time stamps of the delta band time series
%      ({'timeseries', 3} or {'angles', -90} for example)
% PM.SOcenter        = {'angles', 0}; % (0 is pos. peak, -90 = inflection)
PM.SOcenter        = {'timeseries', 7};
% 1) = To which time bin corresponds that particular event
% 2) = startTime
% 3) = midTime (SO: down-up zero crossing)
% 4) = endTime
% 5) = duration
% 6) = maxTime
% 7) = minTime
% 8) = minAmp
% 9) = maxAmp
% 10)= p2pAmp
% 11)= p2pTime
% 12)= Power
% 13)= Frequency
%    - Time stamp of spindle to couple to SO
PM.SSCoupling       = 6;
%    - Allowed range during odor stimulation (seconds)
PM.allowed_range    = 7.5;
%    - Various parameters
fileNames           = {...
    'RC_051_sleep','RC_091_sleep','RC_121_sleep','RC_131_sleep',...
    'RC_141_sleep','RC_161_sleep','RC_171_sleep','RC_201_sleep',...
    'RC_241_sleep','RC_251_sleep','RC_261_sleep','RC_281_sleep',...
    'RC_291_sleep','RC_301_sleep',...
    'RC_392_sleep','RC_412_sleep','RC_442_sleep','RC_452_sleep',...
    'RC_462_sleep','RC_472_sleep','RC_482_sleep','RC_492_sleep',...
    'RC_512_sleep'};
% fileNames           = {...
%     'RC_052_sleep','RC_092_sleep','RC_122_sleep','RC_132_sleep',...
%     'RC_142_sleep','RC_162_sleep','RC_172_sleep','RC_202_sleep',...
%     'RC_242_sleep','RC_252_sleep','RC_262_sleep','RC_282_sleep',...
%     'RC_292_sleep','RC_302_sleep',...
%     'RC_391_sleep','RC_411_sleep','RC_441_sleep','RC_451_sleep',...
%     'RC_461_sleep','RC_471_sleep','RC_481_sleep','RC_491_sleep',...
%     'RC_511_sleep'};
% saveFolder            = [eventfilePath, 'SO_timeSeries', filesep];
saveFolder            = [eventfilePath, 'SO_timeSeries_MinTime_7.5s', filesep];
PM.stimulation_seq    = 'OFF_ON';

% -------------------------------------------------------------------------
% --------------------------                 ------------------------------
% --------------------------   O U T P U T   ------------------------------
% --------------------------                 ------------------------------
% -------------------------------------------------------------------------

% 1. Fieldtrip-structured time series (trials) centered around event time
%    stamps.
%    Each dataset (subject) will be saved separately.
% 2. Spindle latencies relative to Slow oscillation center latencies for
%    each channel for each condition subject-wise (in samples)



%% Set up userland
%  ------------------------------------------------------------------------

eeglab nogui
addpath(genpath('D:\MATLAB\fieldtrip-20200831'))

% Load file containing slow oscillation time stamps in trials
if ~exist('events_loaded', 'var') || events_loaded == 0
    load([eventfilePath eventFile])
    events_loaded = 1;
end

% Prepare storage
if exist(saveFolder, 'dir') ~= 7
    mkdir(saveFolder)
end

% Set up channels of interest
if strcmp(chanOfInterest, 'all')
    PM.Channels = Info.ROIs.str_chans;
else
    PM.Channels = chanOfInterest;
end

% This will point to continuous datasets as well as epoched ones.
% This way, we can just search for event numbers in the EEG structure in an
% easy way:
% 1. Get the time stamp of a slow osc. in a given trial.
% 2. Look up the event number (EEG structure) inside the epoched datasets.
%    This will allow us to search for this event number in the continuous
%    dataset.
% 3. Look at latency of event number in the continuous dataset and adapt so
%    time stamp accordingly.
fileNames_CONTINUOUS = strcat(fileNames, '_NREM');
fileNames_TRIALS     = strcat(fileNames, '_TRIALS');



%% Separate events, before checking SO
%  ------------------------------------------------------------------------

for i_subj = 1:numel(fileNames)
    
    clearvars SO_timeSeries
    disp(strcat('Subject:', {' '}, fileNames(i_subj)))
    
    if exist([saveFolder, char(fileNames(i_subj)), '.mat'], 'file') == 2
        disp('... skipped because exsists already')
        continue 
    end
    
    
    % Check for compatibility
    if ~contains(Info.Subjects(i_subj), fileNames(i_subj))
        % Checkpoint
        error('Incompatible datasets')
    end
    
    % Load datasets:
    % 1. Continuous dataset will serve for extraction of time window
    % around slow oscillations.
    % 2. Trial datasets will serve as indicator of event in EEG structure
    % based on trials containing slow oscillation
    if ~exist('data_loaded', 'var') || data_loaded == 0
        dataType    = '.set';
        [EEGWhole] = f_load_data(...
            strcat(char(fileNames_CONTINUOUS(i_subj)), dataType), ...
            filePath_CONTINUOUS, dataType);
        
        dataType    = '.set';
        [EEGTrialCue] = f_load_data(...
            strcat(char(fileNames_TRIALS(i_subj)), ...
            '_', PM.stimulation_seq, '_Odor', ...
            dataType), filePath_TRIALS, dataType);
        
        dataType    = '.set';
        [EEGTrialVehicle] = f_load_data(...
            strcat(char(fileNames_TRIALS(i_subj)), ...
            '_', PM.stimulation_seq, '_Sham',...
            dataType), filePath_TRIALS, dataType);
        
        data_loaded = 1;
    end
    
    
    
    for i_chan = 1:numel(PM.Channels)
        % /!\ Each channel has to be done separately since detected slow
        % oscillations vary between channels
        channel             = char(PM.Channels(i_chan));
        
        
        % Select data according to current channel
        idx_chan2rejWhole   = ...
            find(~strcmp({EEGWhole.chanlocs.labels}, channel));
        
        idx_chan2rejCue     = ...
            find(~strcmp({EEGTrialCue.chanlocs.labels}, channel));
        idx_chan2rejVehicle = ...
            find(~strcmp({EEGTrialVehicle.chanlocs.labels}, channel));
        
        
        if any(~ismember(idx_chan2rejWhole, idx_chan2rejCue)) || ...
                any(~ismember(idx_chan2rejWhole, idx_chan2rejVehicle))
            % Checkpoint
            error('Incompatible datasets')
        end
        
        [EEGWholeNative]    = ...
            pop_select( EEGWhole, 'nochannel', idx_chan2rejWhole);
        [EEGTrialCueIN]     = ...
            pop_select( EEGTrialCue, 'nochannel', idx_chan2rejCue);
        [EEGTrialVehicleIN] = ...
            pop_select( EEGTrialVehicle, 'nochannel', idx_chan2rejVehicle);
        
        
        
        %% Search for valid trials inside continuous time series
        %  ----------------------------------------------------------------
        
        % The field 'mffkey_gidx' is constant and did not change after 
        % trial rejection. This allows for searching of valid trials inside
        % continuous dataset
        % Triggers will be trial midpoints Off [-15s 0s] and On [0s 15s]
        [~, ~, Trig_OI_Cue]     = intersect(...
            {EEGTrialCueIN.event.mffkey_gidx}, ...
            {EEGWholeNative.event.mffkey_gidx});
        
        [~, ~, Trig_OI_Vehicle] = intersect(...
            {EEGTrialVehicleIN.event.mffkey_gidx}, ...
            {EEGWholeNative.event.mffkey_gidx});
        
        Trig_OI_Cue     = sort(Trig_OI_Cue);
        Trig_OI_Vehicle = sort(Trig_OI_Vehicle);
        
        
        % Find latency of period onsets in the original file (samples)
        Latencies.OdorOn   = [EEGWholeNative.event(Trig_OI_Cue).latency];
        Latencies.ShamOn   = [EEGWholeNative.event(Trig_OI_Vehicle).latency];
        
        Latencies.OdorOff  = Latencies.OdorOn - ...
            Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs;
        Latencies.ShamOff  = Latencies.ShamOn - ...
            Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs;
        
        
        
        %% Prepare trial time vectors
        %  ----------------------------------------------------------------
        
        timeOff = 1:Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs;
        timeOn  = Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs + 1: ...
            (Info.TrialParameters.s_TimeBeforeZero + ...
            Info.TrialParameters.s_TimeAfterZero) * ...
            Info.TrialParameters.s_fs;
        
        
        
        %% Add Slow oscillation event latencies to continuous time series
        %  ----------------------------------------------------------------
                
        for condition = PM.Conditions
                        
            SO_events = OverallSlowOsc.(...
                channel).(char(condition))(:, i_subj);
%             SO_events = OverallSpindles.(...
%                 channel).(char(condition))(:, i_subj);
            
            % Important to reset here! Otherwise, by working directly with
            % EEGWholeNATIVE, leftover time stamps from previous condition 
            % loops will be taken into account as well, thus epoching the 
            % EEG set into epochs belonging to current condition but also 
            % to previous ones.
            EEGWholeIN = EEGWholeNative; 
            
            
            % We will also extract relative latencies of spindles according
            % to slow oscillation center points.
            SS_events           = OverallSpindles.(...
                                    channel).(char(condition))(:, i_subj);
            coupled_spindles    = [];
            
            
            OriginalStartTime   = [];
            OriginalCentTime    = [];
            OriginalendTime     = [];
            Lat_diff            = [];
            for trial = 1:length(SO_events)
                
                if isempty(SO_events{trial})
                    if trial - 1 ~= numel(Latencies.(char(condition)))
                        % Checkpoint
                        error('Incompatible trials')
                    end
                    break
                end
                if size(SO_events{trial}, 1) == 1 && ...
                        all(isnan(SO_events{trial}))
                    continue
                end
                
                for SO = 1:size(SO_events{trial}, 1)
                    
                    
                    
                    %% Extract time points according to axis crossings
                    %  -----------------------------------------------
                    
                    % Time of SO occurence in the trial (in samples)
                    startTime = SO_events{trial}(SO, PM.SOlocations(1));
                    endTime   = SO_events{trial}(SO, PM.SOlocations(2));
                    
                    
                    % Just to filter out trial chunks
                    % -------------------------------
                    accepted_range = PM.allowed_range;
                    
                    if endTime > accepted_range * Info.TrialParameters.s_fs
                        continue
                    end
                    
                    
                    if strcmp(PM.SOcenter(1), 'timeseries')
                        
                        centTime  = SO_events{trial}(SO, PM.SOcenter{2});
                        
                    elseif strcmp(PM.SOcenter(1), 'angles')
                                                
                        % Search for the angle of interest of the signal
                        % between start end end of slow oscillation. For
                        % this, we:
                        % 1. Extract the time series of the current
                        %    trial and the current condition
                        % 2. Filter it accordingly to parameters used 
                        %    during SO detection
                        % 3. Determine the angles of the slow oscillation
                        % 4. Seach for the angle of interest
                        
                        centTime  = NaN; % Make sure to reset
                        
                        % 1. Extract trial condition series
                        if contains(condition, 'OdorOn')
                            DataIN = ...
                                EEGTrialCueIN.data(1, timeOn, trial);
                        elseif contains(condition, 'OdorOff')
                            DataIN = ...
                                EEGTrialCueIN.data(1, timeOff, trial);
                        elseif contains(condition, 'ShamOn')
                            DataIN = ...
                                EEGTrialVehicleIN.data(1, timeOn, trial);
                        elseif contains(condition, 'ShamOff')
                            DataIN = ...
                                EEGTrialVehicleIN.data(1, timeOff, trial);
                        end
                        
                        % 2. Filter the signal into delta band using 
                        %    parameters of SO detection (stored in
                        %    Info.SOparameters): For general Slow osc. 
                        %    detection, signal was filtered in whole delta
                        %    range; for phase-coupling (PC), signal was 
                        %    filtered into low delta range.
                        DataFilt    = f_filter_deltaband(DataIN', ...
                                        Info.SOparameters.butter_order, ...
                                        Info.TrialParameters.s_fs, ...
                                        Info.SOparameters.LowPassFr, ...
                                        Info.SOparameters.HighPassFr);
                        DataFiltPC  = f_filter_deltaband(DataIN', ...
                                        Info.SOparameters.butter_order, ...
                                        Info.TrialParameters.s_fs, ...
                                        Info.SOparameters.LowPassPhCpl, ...
                                        Info.SOparameters.HighPassFr);
                        
                        % 3. Determination of SO angles
                        DataHil     = hilbert(DataFiltPC);
                        DataAngles  = angle(DataHil);
                        
                        SO_Dat      = DataFiltPC(startTime:endTime)';
                        SO_Hil      = real(DataHil(startTime:endTime))';
                        SO_Ang      = DataAngles(startTime:endTime)';
                        
                        % 4. Angle point determination
                        % We take first derivative and look where in it the
                        % difference between data points and zero is
                        % smallest.
                        % v_grad        = gradient(SO_Hil) ./ ...
                        %                   gradient(1:numel(SO_Hil));
                       	% v_sub         = zeros(1, numel(v_grad));
                        % [~, i_delt]	= min(diff([v_grad; v_sub]));
                        
                        % Above method is rejected since it is working with
                        % time series data and not angle data. Furthermore,
                        % it is only working if angle of interest is -90.
                        % We should probably work with angles directly 
                        % since they have been used in phase-shift analysis
                        v_grad      = SO_Ang;
                        v_sub       = repmat(deg2rad(PM.SOcenter{2}), ...
                                        1, numel(v_grad));
                        [~, i_delt]	= min(abs(diff([v_grad; v_sub])));
                        
                        % My brain explodes, therefore we visualize!
                        % run p_visualize_time_stamps.m
                        % close all
                        
                        centTime        = startTime + i_delt - 1;
                        
                        % Difference to zero-crossing of delta signal?
                        Lat_diff = [Lat_diff, ...
                            SO_events{trial}(SO, 3) - centTime];
                        
                    end
                    
                    
                    % Some Slow oscillations, when filtered in low delta
                    % band opposed to whole delta band do not have
                    % appropriate wave forms to determine angle inflection
                    % points. See example of flawed SO detection as figure.
                    % We skip these here:
                    if isnan(centTime)
                        continue
                    end
                    
                    
                    % Time of SO occurence in the whole recording
                    OriginalStartTime   = cat(1, OriginalStartTime, ...
                        Latencies.(char(condition))(trial) + ...
                        startTime);
                    OriginalCentTime    = cat(1, OriginalCentTime, ...
                        Latencies.(char(condition))(trial) + ...
                        centTime);
                    OriginalendTime     = cat(1, OriginalendTime, ...
                        Latencies.(char(condition))(trial) + ...
                        endTime);
                    
                    
                    
                    %% Extract Spindles that are coupled to this very SO
                    %  ----------------------------------------------------
                    
                    for SS = 1:size(SS_events{trial}, 1)
                        curr_spindle = SS_events{trial}(SS, :);
                        
                        % Verify whether spindle time stamp of choice
                        % falls between onset (2) and end latency (4) of SO
                        if curr_spindle(PM.SSCoupling) > ...
                                SO_events{trial}(SO, PM.SOlocations(1)) ...
                                && curr_spindle(PM.SSCoupling) < ...
                                SO_events{trial}(SO, PM.SOlocations(2))
                            
                            % If so, we catenate for this condition for
                            % this subject for this channel the spindle
                            % latency relative to center point of SO
                            rel_latency = ...
                                curr_spindle(PM.SSCoupling) - centTime;
                            
                            coupled_spindles = [coupled_spindles, ...
                                rel_latency];
                        end
                        
                    end
                    
                    
                end
            end
            
            
            
            %% Add event latencies of Slow osc. to EEG.event structure
            %  ------------------------------------------------------------
            
            Lat_diff
            
            s_last_event = numel(EEGWholeIN.event);
            for SO = 1:numel(OriginalStartTime)
                % Add triggers to the original EEG
                EEGWholeIN.event(SO+s_last_event).type      = 'SO_start';
                EEGWholeIN.event(SO+s_last_event).latency   = ...
                    OriginalStartTime(SO);
                EEGWholeIN.event(SO+s_last_event).duration  = 0;
                EEGWholeIN.event(SO+s_last_event).label     = 'SO_start';
            end
            
            s_last_event = numel(EEGWholeIN.event);
            for SO = 1:numel(OriginalCentTime)
                % Add triggers to the original EEG
                EEGWholeIN.event(SO+s_last_event).type      = 'SO_Cent';
                EEGWholeIN.event(SO+s_last_event).latency   = ...
                    OriginalCentTime(SO);
                EEGWholeIN.event(SO+s_last_event).duration  = 0;
                EEGWholeIN.event(SO+s_last_event).label     = 'SO_Cent';
            end
            
            s_last_event = numel(EEGWholeIN.event);
            for SO = 1:numel(OriginalendTime)
                % Add triggers to the original EEG
                EEGWholeIN.event(SO+s_last_event).type      = 'SO_End';
                EEGWholeIN.event(SO+s_last_event).latency   = ...
                    OriginalendTime(SO);
                EEGWholeIN.event(SO+s_last_event).duration  = 0;
                EEGWholeIN.event(SO+s_last_event).label     = 'SO_End';
            end
            
            
            
            %% Epoch Data
            %  ------------------------------------------------------------
        
            clearvars EEGWholeOUT outFT
            
            if ~isempty(OriginalCentTime)
                % Epoch the data around center point of slow osc.
                EEGWholeOUT = pop_epoch(EEGWholeIN, {'SO_Cent'}, ...
                    [-PM.s_timeWindow/2 PM.s_timeWindow/2]);
                
                % Reject huge trailing data that will not be used
                if isfield(EEGWholeOUT, 'rejecteddata')
                    EEGWholeOUT = rmfield(EEGWholeOUT, 'rejecteddata');
                end
                if isfield(EEGWholeOUT, 'rejectedchanlocs')
                    EEGWholeOUT = rmfield(EEGWholeOUT, 'rejectedchanlocs');
                end
                
                
                
                %% Transform EEG structure to Fieldtrip
                %  ------------------------------------------------------------
                
                outFT = eeglab2fieldtrip(EEGWholeOUT, 'raw');
                SO_timeSeries.(channel).(char(condition)) = outFT;
            else
                SO_timeSeries.(channel).(char(condition)) = struct();
                EEGWholeOUT.lst_changes = [];
            end
            
            
            
            %% Store SO-relative latencies of coupled spindles
            %  ------------------------------------------------------------
            
            SS_latencies.(channel).(char(condition)) = coupled_spindles;
            
        end
    end
    
    
    %% Save subjects datasets
    %  --------------------------------------------------------------------
    PM.Name         = char(fileNames(i_subj));
    PM.Info         = Info;
    PM.ListChanges  = EEGWholeOUT.lst_changes;
        
    save(strcat(saveFolder, char(fileNames(i_subj)), '.mat'), ...
        'PM', 'SO_timeSeries', 'SS_latencies', '-v7.3')
    
    
    data_loaded = 0; % Load next subject
    
end
