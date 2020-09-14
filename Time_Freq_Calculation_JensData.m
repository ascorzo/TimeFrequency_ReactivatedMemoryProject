
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    CUE NIGHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = 'D:\GermanData\DATA\CleanByJens\Cue\';
files = dir(strcat(filepath,'*.mat'));

%Datasubj.Events = Events;
%addpath(genpath('C:\Users\lanan\Documents\MATLAB\fieldtrip-20190828'))
ft_defaults



for subj = 1%:numel(files)
    disp(strcat('Sujeto:',num2str(subj)))
    filename = files(subj).name;
    Datasubj = load(strcat(filepath,filename));
    trigger_index  = find(strcmp({Datasubj.Events.value}, 'DIN1')); % index of the triggers in which odor valve was switched on
    trigger_info   = [Datasubj.Events(trigger_index).sample];% sample in which the odor was switched on (either Cue or Vehicle)
    
    epoched_data = []; trigInfo = [];
    preStim      = 15*Datasubj.data.fsample; % samples to include before trigger
    postStim      = 15*Datasubj.data.fsample;% samples to include after trigger

    for trial=1:size(Datasubj.data.sampleinfo,1)
        %----------------------------------------------------------------------
        % check if the moment when the trigger was presented coincides with
        % a non rejected recording
        %----------------------------------------------------------------------
        sampleinfo = Datasubj.data.sampleinfo(trial,:); %start and end sample of each trial
        trigInSeg = find(trigger_info>=sampleinfo(1) & trigger_info<=sampleinfo(2)); % triggers in segment
        
        if isempty(trigInSeg)
            continue
        end
        trigId    = trigInSeg; % Trigger ids that are in the segment
        trigInSeg = trigger_info(trigInSeg); % samples in which each trigger in segment is presented
        data = Datasubj.data.trial{trial}; % temporal EEG data with the data in the n segment
        
        for ti = 1:length(trigInSeg)
            trigIdx     = trigInSeg(ti) - sampleinfo(1); %relative sample of the trigger according to the start of the segment
            
            try
                x            = data(:,trigIdx-preStim:trigIdx+postStim-1); %data selected for each epoch around the trigger time
                epoched_data = cat(3, x, epoched_data); %all data divided in epochs
                trigInfo     = [trigInfo trigId(ti)]; %vector saving the trigger ids that are in the non-rejected segments
            catch
                disp('Exceeding index value,skipping epoch');
                continue
            end
        end
    end
    %%
    %--------------------------------------------------------------------------
    % Identify which epochs are Vehicle(sham) odor, and which ones are cue odor
    %--------------------------------------------------------------------------
    get_cidx = [];
    for chan = 1:numel(trigInfo)
        get_cidx{chan} = ...
            Datasubj.Events(trigger_index(trigInfo(chan))).orig.keys(2).key.data.data;
    end
    
    Sham_trigs = find(mod(str2double(get_cidx),2)==0);% index of the triggers that are Vehicle switched on
    Odor_trigs = find(mod(str2double(get_cidx),2)> 0);% index of the triggers that are Cue Odor switched on
    
    Sham_Epochs = epoched_data(:,:,Sham_trigs);
    Odor_Epochs = epoched_data(:,:,Odor_trigs);
    
    %% Importing to Fieldtrip
    time = 1:size(epoched_data,2);
    time = time/Datasubj.data.fsample - 15;
    ft_warning off
    
    %--------------------------------------------------------------------------
    % Cue data
    %--------------------------------------------------------------------------
    cue_data        = [];
    cue_data.label  = Datasubj.data.label;
    cue_data.time   = repmat({time}, [size(Odor_Epochs,3),1])';
    
    for trial = 1:size(Odor_Epochs,3); cue_data.trial{trial} = squeeze(Odor_Epochs(:,:,trial)); end
    
    cue_data.fsample = Datasubj.data.fsample;
    cue_data.trialinfo = repmat([1 length(time)], [size(Odor_Epochs,3) 1]);
    
    %check everything is working
    cfg = [];
    cue_data = ft_preprocessing(cfg, cue_data);
    
    
    %--------------------------------------------------------------------------
    % Sham data
    %--------------------------------------------------------------------------
    sham_data        = [];
    sham_data.label  = Datasubj.data.label;
    sham_data.time   = repmat({time}, [size(Sham_Epochs,3),1])';
    
    for trial = 1:size(Sham_Epochs,3); sham_data.trial{trial} = squeeze(Sham_Epochs(:,:,trial)); end
    
    sham_data.fsample = Datasubj.data.fsample;
    sham_data.trialinfo = repmat([1 length(time)], [size(Sham_Epochs,3) 1]);
    
    %check everything is working
    cfg = [];
    sham_data = ft_preprocessing(cfg, sham_data);
    
    %% -- PSD
    
    %------- Parameters -----------------------------------------------
    
    cfg                 = [];
    cfg.method          = 'mtmconvol';
    cfg.taper           = 'hanning';
    cfg.output          = 'pow';
    cfg.pad             = 'nextpow2';
    cfg.foi             =  0:0.25:20;                        
    cfg.t_ftimwin       = ones(length(cfg.foi),1);          
    cfg.toi             = -15:0.25:15;
    cfg.keeptrials      = 'yes';
    
    %------- Cue data -----------------------------------------------
    Time_Freq_Cue{subj}   = ft_freqanalysis(cfg, cue_data);

    %------- Sham data -----------------------------------------------
    Time_Freq_Sham_CN{subj}  = ft_freqanalysis(cfg, sham_data);
    
    
end

%------------------------------------------------------------------------
% Find common channels between Subjects
%------------------------------------------------------------------------

commonChans = Time_Freq_Cue{1}.label;
for subj = 1:length(Time_Freq_Cue)
    commonChans = intersect(commonChans, Time_Freq_Cue{subj}.label);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    PLACEBO  NIGHT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filepath = 'D:\GermanData\DATA\CleanByJens\Placebo\';
files = dir(strcat(filepath,'*.mat'));



for subj = 1:numel(files)
    disp(strcat('Sujeto:',num2str(subj)))
    filename = files(subj).name;
    Datasubj = load(strcat(filepath,filename));
    trigger_index  = find(strcmp({Datasubj.Events.value}, 'DIN1')); % index of the triggers in which odor valve was switched on
    trigger_info   = [Datasubj.Events(trigger_index).sample];% sample in which the odor was switched on (either Cue or Vehicle)
    
    epoched_data = []; trigInfo = [];
    preStim      = 15*Datasubj.data.fsample; % samples to include before trigger
    postStim      = 15*Datasubj.data.fsample;% samples to include after trigger

    for trial=1:size(Datasubj.data.sampleinfo,1)
        %----------------------------------------------------------------------
        % check if the moment when the trigger was presented coincides with
        % a non rejected recording
        %----------------------------------------------------------------------
        sampleinfo = Datasubj.data.sampleinfo(trial,:); %start and end sample of each trial
        trigInSeg = find(trigger_info>=sampleinfo(1) & trigger_info<=sampleinfo(2)); % triggers in segment
        
        if isempty(trigInSeg)
            continue
        end
        trigId    = trigInSeg; % Trigger ids that are in the segment
        trigInSeg = trigger_info(trigInSeg); % samples in which each trigger in segment is presented
        data = Datasubj.data.trial{trial}; % temporal EEG data with the data in the n segment
        
        for ti = 1:length(trigInSeg)
            trigIdx     = trigInSeg(ti) - sampleinfo(1); %relative sample of the trigger according to the start of the segment
            
            try
                x            = data(:,trigIdx-preStim:trigIdx+postStim-1); %data selected for each epoch around the trigger time
                epoched_data = cat(3, x, epoched_data); %all data divided in epochs
                trigInfo     = [trigInfo trigId(ti)]; %vector saving the trigger ids that are in the non-rejected segments
            catch
                disp('Exceeding index value,skipping epoch');
                continue
            end
        end
    end
    %%
    %--------------------------------------------------------------------------
    % Identify which epochs are Vehicle(sham) odor, and which ones are cue odor
    %--------------------------------------------------------------------------
    get_cidx = [];
    for chan = 1:numel(trigInfo)
        get_cidx{chan} = ...
            Datasubj.Events(trigger_index(trigInfo(chan))).orig.keys(2).key.data.data;
    end
    
    Sham_trigs = find(mod(str2double(get_cidx),2)==0);% index of the triggers that are Vehicle switched on
    Odor_trigs = find(mod(str2double(get_cidx),2)> 0);% index of the triggers that are Cue Odor switched on
    
    Sham_Epochs = epoched_data(:,:,Sham_trigs);
    Odor_Epochs = epoched_data(:,:,Odor_trigs);
    
    %% Importing to Fieldtrip
    time = 1:size(epoched_data,2);
    time = time/Datasubj.data.fsample - 15;
    ft_warning off
    
    %--------------------------------------------------------------------------
    % Cue data
    %--------------------------------------------------------------------------
    Placebo_data        = [];
    Placebo_data.label  = Datasubj.data.label;
    Placebo_data.time   = repmat({time}, [size(Odor_Epochs,3),1])';
    
    for trial = 1:size(Odor_Epochs,3); Placebo_data.trial{trial} = squeeze(Odor_Epochs(:,:,trial)); end
    
    Placebo_data.fsample = Datasubj.data.fsample;
    Placebo_data.trialinfo = repmat([1 length(time)], [size(Odor_Epochs,3) 1]);
    
    %check everything is working
    cfg = [];
    Placebo_data = ft_preprocessing(cfg, Placebo_data);
    
    
    %--------------------------------------------------------------------------
    % Sham data
    %--------------------------------------------------------------------------
    sham_data        = [];
    sham_data.label  = Datasubj.data.label;
    sham_data.time   = repmat({time}, [size(Sham_Epochs,3),1])';
    
    for trial = 1:size(Sham_Epochs,3); sham_data.trial{trial} = squeeze(Sham_Epochs(:,:,trial)); end
    
    sham_data.fsample = Datasubj.data.fsample;
    sham_data.trialinfo = repmat([1 length(time)], [size(Sham_Epochs,3) 1]);
    
    %check everything is working
    cfg = [];
    sham_data = ft_preprocessing(cfg, sham_data);
    
    %% -- PSD
    
    %------- Parameters -----------------------------------------------
    
    cfg                 = [];
    cfg.method          = 'mtmconvol';
    cfg.taper           = 'hanning';
    cfg.output          = 'pow';
    cfg.pad             = 'nextpow2';
    cfg.foi             =  0:0.25:20;                        
    cfg.t_ftimwin       = ones(length(cfg.foi),1);         
    cfg.toi             = -15:0.25:15;
    cfg.keeptrials      = 'yes';
    
    %------- Cue data -----------------------------------------------
    Time_Freq_Placebo{subj}   = ft_freqanalysis(cfg, Placebo_data);

    %------- Sham data -----------------------------------------------
    Time_Freq_Sham_PN{subj}  = ft_freqanalysis(cfg, sham_data);
    
    
end

%------------------------------------------------------------------------
% Find common channels between Subjects
%------------------------------------------------------------------------
for subj = 1:length(Time_Freq_Placebo)
    commonChans = intersect(commonChans, Time_Freq_Placebo{subj}.label);
end