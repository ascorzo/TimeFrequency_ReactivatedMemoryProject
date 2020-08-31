% p_SigExplore.m
clear all
filepath='D:\GermanData\DATA\CleanByJens\Cue\';
files = dir(strcat(filepath,'*.mat'));

%load('m_GaborCuetotal_JensData.mat')

s_MinFreqHz = 0.3;
s_MaxFreqHz = 20;
s_FreqSeg = 500;
s_StDevCycles = 3;
s_Magnitudes = 1;
s_SquaredMag = 0;
s_TimeStep = 0.5;
srate = 1000;
preStimTime = 15;
postStimTime = 15;

for subj = 1:numel(files)
    
    disp(strcat('Subj:',num2str(subj)))

    m_GaborTrials = [];
    v_TimeTrials = [];
    m_GaborCue = [];
    m_GaborSham = [];
    
    filename = files(subj).name;
    Datasubj = load(strcat(filepath,filename));
    
    disp(strcat('ntrials:',num2str(numel(Datasubj.data.trial))))
    
    trigger_index  = find(strcmp({Datasubj.Events.value}, 'DIN1')); % index of the triggers in which odor valve was switched on
    trigger_info   = [Datasubj.Events(trigger_index).sample];% sample in which the odor was switched on (either Cue or Vehicle)
    
    epoched_data = []; trigInfo = []; v_TimeAxisTotal = [];
    preStim      = preStimTime*Datasubj.data.fsample; % samples to include before trigger
    postStim      = postStimTime*Datasubj.data.fsample;% samples to include after trigger
    
    for trial = 1:numel(Datasubj.data.trial)
        disp(strcat('trial',num2str(trial)))
        DataTemp = Datasubj.data.trial{1,trial};
        m_GaborChans = [];
        for nchan = 1:size(DataTemp,1)
            
            %% Gabor Wavelet analysis
            [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
                f_GaborAWTransformMatlab(...
                DataTemp(nchan,:), ...
                srate, ...
                s_MinFreqHz, ...
                s_MaxFreqHz, ...
                s_FreqSeg, ...
                s_StDevCycles, ...
                s_Magnitudes, ...
                s_SquaredMag, ...
                [], [], ...
                s_TimeStep);
            
            m_GaborChans(:,:,nchan) = m_GaborWT;
        end
        m_GaborTrials = cat(2,m_GaborTrials,m_GaborChans);
        v_TimeAxisTotal = cat(2,v_TimeAxisTotal,...
                v_TimeAxis+(Datasubj.data.sampleinfo(trial,1))/Datasubj.data.fsample);
    end
    
    clear DataTemp m_GaborChans

    %----------------------------------------------------------------------
    % Identify events in non-rejected segments
    %----------------------------------------------------------------------
    for n=1:size(Datasubj.data.sampleinfo,1)
        
        sampleinfo = Datasubj.data.sampleinfo(n,:); %start and end sample of each trial
        trigInSeg = ...
            find(trigger_info>=sampleinfo(1) & trigger_info<=sampleinfo(2)); % triggers in segment
        if isempty(trigInSeg)
            continue
        end
        trigId    = trigInSeg; % Trigger ids that are in the segment
        trigInSeg = trigger_info(trigInSeg); % samples in which each trigger in segment is presented            
        for ti = 1:length(trigInSeg)
            start_sample = trigInSeg(ti) - preStim;
            end_sample = trigInSeg(ti) + postStim;
            if start_sample>=sampleinfo(1) && end_sample<=sampleinfo(2)
                trigInfo = [trigInfo trigId(ti)]; %vector saving the trigger ids that are in the non-rejected segments
            end
        end
    end
   
     
    %--------------------------------------------------------------------------
    % Identify which epochs are Vehicle(sham) odor, and which ones are cue odor
    %--------------------------------------------------------------------------
    for i = 1:numel(trigInfo)
        get_cidx{i} = ...
            Datasubj.Events(trigger_index(trigInfo(i))).orig.keys(2).key.data.data;
    end
    
    Sham_trigs = find(mod(str2double(get_cidx),2)==0);% index of the triggers that are Vehicle switched on (within the trigInfo vector)
    Cue_trigs = find(mod(str2double(get_cidx),2)> 0);% index of the triggers that are Cue Odor switched on (within the trigInfo vector)
    
    %----------------------------------------------------------------------
    % cut events
    %----------------------------------------------------------------------
    
    %--- for Cue epochs
    for ct = 1:numel(Cue_trigs)
        trigSample = ...
            Datasubj.Events(trigger_index(trigInfo(Cue_trigs(ct)))).sample;
        trigTime = trigSample/Datasubj.data.fsample;
        
        [~,trigGaborTime] = min(abs(v_TimeAxisTotal-trigTime));
        
        trigInterval = [trigGaborTime-preStimTime/s_TimeStep,...
            trigGaborTime+postStimTime/s_TimeStep];
        
        epoch = ct;
        m_GaborCue(:,:,:,epoch) = ...
            m_GaborTrials(:,trigInterval(1):trigInterval(2),:);
    end
    
    %--- for Sham epochs
    
    for st = 1:numel(Sham_trigs)
        trigSample = ...
            Datasubj.Events(trigger_index(trigInfo(Sham_trigs(st)))).sample;
        trigTime = trigSample/Datasubj.data.fsample;
        
        [~,trigGaborTime] = min(abs(v_TimeAxisTotal-trigTime));
        
        trigInterval = [trigGaborTime-preStimTime/s_TimeStep,...
            trigGaborTime+postStimTime/s_TimeStep];
        
        epoch = st;
        m_GaborSham(:,:,:,epoch) = ...
            m_GaborTrials(:,trigInterval(1):trigInterval(2),:);
    end
    
    name = filename(1:4);
    m_GaborCue_total.(name) = m_GaborCue;
    m_GaborSham_total.(name) = m_GaborSham;
    
    
    clear m_GaborCue m_GaborSham  Datasubj m_GaborTrials get_cidx
end

save('m_GaborCuetotal_JensData_500seg.mat','m_GaborCue_total','m_GaborSham_total','v_FreqAxis','-v7.3')

% %% Fill the Missing channels with NaNs 
% clear all
% filepath='D:\GermanData\DATA\CleanByJens\Placebo\';
% files = dir(strcat(filepath,'*.mat'));
% 
% load('m_GaborPlacebototal_JensData.mat');
% m_GaborPlacebototal = m_GaborPlacebo_total;
% m_GaborShamtotal = m_GaborSham_total;
% m_GaborPlacebo_total = [];
% m_GaborSham_total = []; includedChannels = [];
% Allchans = {'E1';'E2';'E3';'E4';'E5';'E6';'E7';'E8';'E9';'E10';'E11';...
%     'E12';'E13';'E14';'E15';'E16';'E17';'E18';'E19';'E20';'E21';'E22';...
%     'E23';'E24';'E25';'E26';'E27';'E28';'E29';'E30';'E31';'E32';...
%     'E33';'E34';'E35';'E36';'E37';'E38';'E39';'E40';'E41';'E42';...
%     'E43';'E44';'E45';'E46';'E47';'E48';'E49';'E50';'E51';'E52';...
%     'E53';'E54';'E55';'E56';'E57';'E58';'E59';'E60';'E61';'E62';...
%     'E63';'E64';'E65';'E66';'E67';'E68';'E69';'E70';'E71';'E72';...
%     'E73';'E74';'E75';'E76';'E77';'E78';'E79';'E80';'E81';'E82';...
%     'E83';'E84';'E85';'E86';'E87';'E88';'E89';'E90';'E91';'E92';...
%     'E93';'E94';'E95';'E96';'E97';'E98';'E99';'E100';'E101';'E102';...
%     'E103';'E104';'E105';'E106';'E107';'E108';'E109';'E110';'E111';...
%     'E112';'E113';'E114';'E115';'E116';'E117';'E118';'E119';'E120';...
%     'E121';'E122';'E123';'E124';'E125';'E126';'E127';'E128';'E129'};
%     
% 
% for subj = 1:numel(files)
%     name = files(subj).name(1:6);
%     load(strcat(filepath,files(subj).name));
%     channels = data.label;
%     clear data Events
%     m_GaborPlacebo_total.(name) = ...
%         NaN(size(m_GaborPlacebototal.(name(1:4)),1),...
%         size(m_GaborPlacebototal.(name(1:4)),2),...
%         numel(Allchans),...
%         size(m_GaborPlacebototal.(name(1:4)),4));
%     
%     m_GaborSham_total.(name) = ...
%         NaN(size(m_GaborShamtotal.(name(1:4)),1),...
%         size(m_GaborShamtotal.(name(1:4)),2),...
%         numel(Allchans),...
%         size(m_GaborShamtotal.(name(1:4)),4));
%     
%     for chan = 1:numel(channels)
%         includedChannels(chan) = find(strcmp(channels(chan),Allchans));
%     end
%     
%     m_GaborPlacebo_total.(name)(:,:,includedChannels,:) = ...
%         m_GaborPlacebototal.(name(1:4));
%     
%     m_GaborSham_total.(name)(:,:,includedChannels,:) = ...
%         m_GaborShamtotal.(name(1:4));
%     
%     clear channels includedChannels
% end
% 
% save('m_GaborPlacebototal_JensData.mat','m_GaborPlacebo_total','m_GaborSham_total','v_FreqAxis','-v7.3')
