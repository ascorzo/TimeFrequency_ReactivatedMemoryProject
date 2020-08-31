% p_SigExplore.m
eeglab nogui

filepath='D:\GermanData\DATA\RawData\preProcessing\DataChans\Cue\';
%filepath='/home/sleep/Documents/DAVID/Datasets/Ori_PlaceboNight/preProcessing/CustomFiltered/';
files = dir(strcat(filepath,'*.set'));

s_MinFreqHz = 0.3;
s_MaxFreqHz = 20;
s_FreqSeg = 40;
s_StDevCycles = 3;
s_Magnitudes = 1;
s_SquaredMag = 0;
s_TimeStep = 0.5;

for subj = 1%2:numel(files)
    
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
    
    OUTEEG_CueOn = pop_epoch( EEG,[],[-15 15],'eventindices',CueOn);
    OUTEEG_ShamOn = pop_epoch( EEG,[],[-15 15],'eventindices',ShamOn);
    
    for nchan = 1:numel(EEG.chanlocs)
        for epoch = 1:size(OUTEEG_CueOn.data,3)
            %% Gabor Wavelet analysis
            
            %For Cue condition
            Data = squeeze(OUTEEG_CueOn.data(nchan,:,epoch));
            
            [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
                f_GaborAWTransformMatlab(...
                Data, ...
                OUTEEG_CueOn.srate, ...
                s_MinFreqHz, ...
                s_MaxFreqHz, ...
                s_FreqSeg, ...
                s_StDevCycles, ...
                s_Magnitudes, ...
                s_SquaredMag, ...
                [], [], ...
                s_TimeStep);
            
            m_GaborCue(:,:,nchan,epoch) = m_GaborWT;
        end
        
        for epoch = 1:size(OUTEEG_ShamOn.data,3)
            
            %For Sham condition
            Data = squeeze(OUTEEG_ShamOn.data(nchan,:,epoch));
            
            [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
                f_GaborAWTransformMatlab(...
                Data, ...
                OUTEEG_ShamOn.srate, ...
                s_MinFreqHz, ...
                s_MaxFreqHz, ...
                s_FreqSeg, ...
                s_StDevCycles, ...
                s_Magnitudes, ...
                s_SquaredMag, ...
                [], [], ...
                s_TimeStep);
            
            m_GaborSham(:,:,nchan,epoch) = m_GaborWT;
        end
    end
    name = strcat('RC ',filename(4:6));
    m_GaborCue_total.(name) = m_GaborCue;
    m_GaborSham_total.(name) = m_GaborSham;
end