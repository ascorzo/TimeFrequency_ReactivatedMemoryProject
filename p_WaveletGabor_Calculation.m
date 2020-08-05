% p_SigExplore.m
eeglab nogui

filepath='D:\GermanData\DATA\CustomFiltered\';
%filepath='/home/sleep/Documents/DAVID/Datasets/Ori_PlaceboNight/preProcessing/CustomFiltered/';
files = dir(strcat(filepath,'*.set'));


for subj = 1:numel(files)
    m_GaborCue = [];
    m_GaborSham = [];
    filename = files(subj).name;
    EEG = pop_loadset(strcat(filepath,filename));
    
    %     %% Separar por DIN1 y DIN2
    All_DIN1 = find(strcmp({EEG.event.code},'DIN1'));
    
    All_DIN2 = find(strcmp({EEG.event.code},'DIN2'));
    
    %% Separar por pares e impares
    
    get_cidx= {EEG.event.mffkey_cidx};
    
    Sham_Epochs = find(mod(str2double(get_cidx),2)==0);
    Cue_Epochs = find(mod(str2double(get_cidx),2)~= 0);
    
    [ShamOn] = intersect(All_DIN1,Sham_Epochs);
    [CueOn] = intersect(All_DIN1,Cue_Epochs);
    
    ShamOnLat = [EEG.event(ShamOn).latency];
    CueOnLat = [EEG.event(CueOn).latency];
    
    for nchan = 1:numel(EEG.chanlocs)-1
        
        %% Gabor Wavelet analysis
        
        s_MinFreqHz = 0.3;
        s_MaxFreqHz = 20;
        s_FreqSeg = 40;
        s_StDevCycles = 3;
        s_Magnitudes = 1;
        s_SquaredMag = 0;
        s_TimeStep = 0.5;
        
        [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
            f_GaborAWTransformMatlab(...
            EEG.data(nchan,:), ...
            EEG.srate, ...
            s_MinFreqHz, ...
            s_MaxFreqHz, ...
            s_FreqSeg, ...
            s_StDevCycles, ...
            s_Magnitudes, ...
            s_SquaredMag, ...
            [], [], ...
            s_TimeStep);
        
        %% cut by events
        
        CueOnLatSec = CueOnLat/EEG.srate;
        ShamOnLatSec = ShamOnLat/EEG.srate;
        
        CueOnLatRound = floor(CueOnLatSec) + ...
            floor( (CueOnLatSec-floor(CueOnLatSec))/0.5) * 0.5;
        CueOnStart = CueOnLatRound-15;
        CueOnEnd = CueOnLatRound+15;
        
        
        for epoch = 1:numel(CueOnLatRound)
            [~,posStart]= find(v_TimeAxis==CueOnStart(epoch));
            [~,posEnd]= find(v_TimeAxis==CueOnEnd(epoch));
            
            if isempty(posEnd)
                continue
            end
            
            m_GaborCue(:,:,nchan,epoch) = m_GaborWT(:,posStart:posEnd);
        end
        
        ShamOnLatRound = floor(ShamOnLatSec) + ...
            floor( (ShamOnLatSec-floor(ShamOnLatSec))/0.5) * 0.5;
        ShamOnStart = ShamOnLatRound-15;
        ShamOnEnd = ShamOnLatRound+15;
        
        
        for epoch = 1:numel(ShamOnLatRound)
            [~,posStart]= find(v_TimeAxis==ShamOnStart(epoch));
            [~,posEnd]= find(v_TimeAxis==ShamOnEnd(epoch));
            
            if isempty(posEnd)
                continue
            end
            
            m_GaborSham(:,:,nchan,epoch) = m_GaborWT(:,posStart:posEnd);
        end
    end
    
    name = strcat('RC ',filename(4:6));
    m_GaborPlacebo_total.(name) = m_GaborCue;
    m_GaborSham_total.(name) = m_GaborSham; 
end