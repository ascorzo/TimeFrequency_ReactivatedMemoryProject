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
                        '200Hz_Oct_NEW\', '05-Mar-2021_Cue_NEW\'];
eventFile           = '05-Mar-2021_22-23-41_AllData.mat';
%    - Channels to include
chanOfInterest      = 'all';
%    - Conditions to include
PM.Conditions       = {'OdorOn', 'ShamOn', 'OdorOff', 'ShamOff'};
%    - Size of desired time window to extract
PM.s_timeWindow     = 6;
%    - Time stamp to center time window around
PM.SOlocations      = [2, 3, 4]; % See list below
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
saveFolder            = [eventfilePath, 'SO_timeSeries', filesep];
PM.stimulation_seq    = 'OFF_ON';

% -------------------------------------------------------------------------
% --------------------------                 ------------------------------
% --------------------------   O U T P U T   ------------------------------
% --------------------------                 ------------------------------
% -------------------------------------------------------------------------

% 1. Fieldtrip-structured time series (trials) centered around event time
%    stamps.
%    Each dataset (subject) will be saved separately.



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
    dataType    = '.set';
    [EEGWhole] = f_load_data(...
        strcat(char(fileNames_CONTINUOUS(i_subj)), dataType), ...
        filePath_CONTINUOUS, dataType);
    
    dataType    = '.set';
    [EEGTrialCue] = f_load_data(...
        strcat(char(fileNames_TRIALS(i_subj)), ...
        '_', PM.stimulation_seq, '_Odor',...
        dataType), filePath_TRIALS, dataType);
    
    dataType    = '.set';
    [EEGTrialVehicle] = f_load_data(...
        strcat(char(fileNames_TRIALS(i_subj)), ...
        '_', PM.stimulation_seq, '_Sham',...
        dataType), filePath_TRIALS, dataType);
    
    
    
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
        
        [EEGWholeIN]        = ...
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
            {EEGWholeIN.event.mffkey_gidx});
        
        [~, ~, Trig_OI_Vehicle] = intersect(...
            {EEGTrialVehicleIN.event.mffkey_gidx}, ...
            {EEGWholeIN.event.mffkey_gidx});
        
        Trig_OI_Cue     = sort(Trig_OI_Cue);
        Trig_OI_Vehicle = sort(Trig_OI_Vehicle);
        
        
        % Find latency of period onsets in the original file (samples)
        Latencies.OdorOn   = [EEGWholeIN.event(Trig_OI_Cue).latency];
        Latencies.ShamOn   = [EEGWholeIN.event(Trig_OI_Vehicle).latency];
        
        Latencies.OdorOff  = Latencies.OdorOn - ...
            Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs;
        Latencies.ShamOff  = Latencies.ShamOn - ...
            Info.TrialParameters.s_TimeBeforeZero * ...
            Info.TrialParameters.s_fs;
        
        
        
        %% Add Slow oscillation event latencies to continuous time series
        %  ----------------------------------------------------------------
                
        for condition = PM.Conditions
                        
            SO_events = OverallSlowOsc.(...
                channel).(char(condition))(:, i_subj);
            
            OriginalStartTime   = [];
            OriginalCentTime    = [];
            OriginalendTime     = [];
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
                    
                    % Time of SO occurence in the trial (in samples)
                    StartTime   = SO_events{trial}(SO, PM.SOlocations(1));
                    centTime    = SO_events{trial}(SO, PM.SOlocations(2));
                    endTime     = SO_events{trial}(SO, PM.SOlocations(3));
                    
                    % Time of SO occurence accorging to trigger in the
                    % whole recording (in miliseconds)
                    TrialStartTime  = StartTime;
                    TrialCentTime   = centTime;
                    TrialendTime    = endTime;
                    
                    % Time of SO occurence in the whole recording
                    OriginalStartTime   = cat(1, OriginalStartTime, ...
                        Latencies.(char(condition))(trial) + ...
                        TrialStartTime);
                    OriginalCentTime    = cat(1, OriginalCentTime, ...
                        Latencies.(char(condition))(trial) + ...
                        TrialCentTime);
                    OriginalendTime     = cat(1, OriginalendTime, ...
                        Latencies.(char(condition))(trial) + ...
                        TrialendTime);
                    
                end
            end
            
            
            
            %% Add event latencies of Slow osc. to EEG.event structure
            %  ------------------------------------------------------------
            
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
            % Epoch the data around center point of slow osc.
            EEGWholeOUT = pop_epoch(EEGWholeIN, {'SO_Cent'}, ...
                [-PM.s_timeWindow/2 PM.s_timeWindow/2]);
            
            % Reject huge trailing data that will not be used
            EEGWholeOUT = rmfield(EEGWholeOUT, 'rejecteddata');
            EEGWholeOUT = rmfield(EEGWholeOUT, 'rejectedchanlocs');
            
            
            
            %% Transform EEG structure to Fieldtrip
            %  ------------------------------------------------------------
            
            outFT = eeglab2fieldtrip(EEGWholeOUT, 'raw');
            SO_condition.(char(condition)) = outFT;
        end
        
        SO_timeSeries.(channel) = SO_condition;
    end
    
    
    %% Save subjects datasets
    %  --------------------------------------------------------------------
    PM.Name         = char(fileNames(i_subj));
    PM.Info         = Info;
    PM.ListChanges  = EEGWholeOUT.lst_changes;
        
    save(strcat(saveFolder, char(fileNames(i_subj)), '.mat'), ...
        'PM', 'SO_timeSeries', '-v7.3')
    
end
