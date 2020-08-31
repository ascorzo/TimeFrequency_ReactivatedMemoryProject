% p_SigExplore.m
eeglab nogui

filepath='D:\GermanData\DATA\TrialGroups_Off_On\';
%filepath='/home/sleep/Documents/DAVID/Datasets/Ori_PlaceboNight/preProcessing/CustomFiltered/';
filesOdor = dir(strcat(filepath,'*Odor*.set'));
filesSham = dir(strcat(filepath,'*Sham*.set'));

s_MinFreqHz = 0.3;
s_MaxFreqHz = 20;
s_FreqSeg = 40;
s_StDevCycles = 3;
s_Magnitudes = 1;
s_SquaredMag = 0;
s_TimeStep = 0.5;

for subj = 1:numel(filesOdor)
    
    m_GaborCue = [];
    filenameOdor = filesOdor(subj).name;
    
    EEG_Odor = pop_loadset(strcat(filepath,filenameOdor));
    
    for nchan = 1:numel(EEG_Odor.chanlocs)
        for epoch = 1:size(EEG_Odor.data,3)
        
        %% Gabor Wavelet analysis
        Data = squeeze(EEG_Odor.data(nchan,:,epoch));

        [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
            f_GaborAWTransformMatlab(...
            Data, ...
            EEG_Odor.srate, ...
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
    end
    
    name = strcat('RC ',filenameOdor(4:6));
    m_GaborCue_total.(name) = m_GaborCue;
end

for subj = 2:numel(filesSham)
    
    m_GaborSham = [];
    filenamSham = filesSham(subj).name;
    EEG_Sham = pop_loadset(strcat(filepath,filenamSham));
    
    
    for nchan = 1:numel(EEG_Sham.chanlocs)
        for epoch = 1:size(EEG_Sham.data,3)
        
        %% Gabor Wavelet analysis
        Data = squeeze(EEG_Sham.data(nchan,:,epoch));

        [m_GaborWT, v_TimeAxis, v_FreqAxis, v_StDevArray] = ...
            f_GaborAWTransformMatlab(...
            Data, ...
            EEG_Sham.srate, ...
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
    
    name = strcat('RC ',filenamSham(4:6));
    m_GaborSham_total.(name) = m_GaborSham; 
end