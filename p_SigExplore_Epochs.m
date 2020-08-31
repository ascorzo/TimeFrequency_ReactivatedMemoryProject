% p_SigExplore.m
%addpath('/home/andrea/Documents/MatlabFunctions/eeglab2019_1/');

eeglab nogui

filepath='D:\GermanData\DATA\CustomFiltered\';
%filepath='/home/sleep/Documents/DAVID/Datasets/Ori_PlaceboNight/preProcessing/CustomFiltered/';
files = dir(strcat(filepath,'*.set'));
%chanOI = 36;%104


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
    
    %% Power Spectral Density
    
    %     [spectra,freqs] = ...
    %         spectopo(EEG.data, 0, EEG.srate,'freq',0:1:20,'chanlocs',EEG.chanlocs);
    %
    %     figure
    %     plot(freqs,spectra)
    %     title(filename(1:6))
    
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


%% Calculate and plot the mean of all Subjects for each Channel
subjects = fieldnames(m_GaborCue_total);
baselineTime = 15;
s_TimeStep = 0.5;
v_time = -15:0.5:15;

for nchan = 1:128
    for subj = 1:numel(subjects)
        
        DataSubj = m_GaborCue_total.(string(subjects(subj)));
        DataSubj = DataSubj(:,:,nchan,:);
        DataSubj = reshape(DataSubj,[size(DataSubj,1),...
            size(DataSubj,2),size(DataSubj,4)]);
        
        %baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:),2);
        baseline = mean(DataSubj,2);
        DataSubjBasCorr = DataSubj-baseline;
        
        MeanCueEpochs(subj,:,:) = mean(DataSubjBasCorr,3);
        
        
        DataSubj = m_GaborSham_total.(string(subjects(subj)));
        DataSubj = DataSubj(:,:,nchan,:);
        DataSubj = reshape(DataSubj,[size(DataSubj,1),...
            size(DataSubj,2),size(DataSubj,4)]);
        
        %baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:),2);
        baseline = mean(DataSubj,2);
        DataSubjBasCorr = DataSubj-baseline;
        
        MeanShamEpochs(subj,:,:) = mean(DataSubjBasCorr,3);
        
    end
    
%     AllMeanCueEpochs = squeeze(mean(MeanCueEpochs,1));
%     AllMeanShamEpochs = squeeze(mean(MeanShamEpochs,1));
%     
     
%     vlims = [min(min(AllMeanCueEpochs(:),AllMeanShamEpochs(:))),...
%         max(max(AllMeanCueEpochs(:),AllMeanShamEpochs(:)))];
%     subplot(2,1,1)
%     f_ImageMatrix(AllMeanCueEpochs,v_time,v_FreqAxis, vlims, [], [], 0);
%     title('Cue Epochs')
%     colorbar
%     subplot(2,1,2)
%     f_ImageMatrix(AllMeanShamEpochs,v_time,v_FreqAxis, vlims, [], [], 0);
%     title('Sham Epochs')
%     colorbar
%     suptitle(strcat('Chan E',num2str(nchan)))
%     
%     saveas(gcf,strcat('Chan E',num2str(nchan),'.png'))
    
    s_LastInd = find(v_FreqAxis <=1, 1);
    s_FirstInd = find(v_FreqAxis <=2, 1);
    m_TFFreqBands = ...
        mean(MeanCueEpochs(:,s_FirstInd:s_LastInd, :),2);
    m_TFFreqBands = reshape(m_TFFreqBands,...
        [size(m_TFFreqBands,1),size(m_TFFreqBands,3)]);
    plot(v_time,m_TFFreqBands,'k','linewidth',0.5)
    hold on
    b1 = plot(v_time,mean(m_TFFreqBands,1),'k','linewidth',3,...
        'DisplayName','Cue Mean');

    
    
    m_TFFreqBands = ...
        mean(MeanShamEpochs(:,s_FirstInd:s_LastInd, :),2);
    m_TFFreqBands = reshape(m_TFFreqBands,...
        [size(m_TFFreqBands,1),size(m_TFFreqBands,3)]);
    plot(v_time,m_TFFreqBands,'r','linewidth',0.5)
    hold on
    b2 = plot(v_time,mean(m_TFFreqBands,1),'r','linewidth',3,...
        'DisplayName','Sham Mean');
    legend ([b1,b2])
    title(strcat('Chan E',num2str(nchan)))
    
    hold off
    saveas(gcf,strcat('Chan_E',num2str(nchan),'_Theta.png'))  
end


%% Calculate the mean of all Subjects for all channels

subjects = fieldnames(m_GaborCue_total);
baselineTime = 15;
s_TimeStep = 0.5;

for subj = 1:numel(subjects)
    
    DataSubj = m_GaborCue_total.(string(subjects(subj)));
    
    baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:,:),2);
    DataSubjBasCorr = DataSubj-baseline;
    
    MeanCueChannels = mean(DataSubjBasCorr,3);
    MeanCueChannels = reshape(MeanCueChannels,...
        [size(MeanCueChannels,1),size(MeanCueChannels,2),...
        size(MeanCueChannels,4)]);
    MeanCueEpochs(subj,:,:) = mean(MeanCueChannels,3);
    
    
    DataSubj = m_GaborSham_total.(string(subjects(subj)));
    
    baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:,:),2);
    DataSubjBasCorr = DataSubj-baseline;
   
    MeanShamChannels = mean(DataSubjBasCorr,3);
    MeanShamChannels = reshape(MeanShamChannels,...
        [size(MeanShamChannels,1),size(MeanShamChannels,2),...
        size(MeanShamChannels,4)]);
    MeanShamEpochs(subj,:,:) = mean(MeanShamChannels,3);
    
end

MeanCueSubjects = squeeze(mean(MeanCueEpochs,1));
MeanShamSubjects = squeeze(mean(MeanShamEpochs,1));

v_time = -15:0.5:15;
%% Statistical analysis
for i = 1:size(MeanCueEpochs,2)
    for j = 1:size(MeanCueEpochs,3)
        cond1 = MeanCueEpochs(:,i,j);
        cond2 = MeanShamEpochs(:,i,j);
        x = [cond1,cond2];
        [h,p,sigPairs] = ttest_bonf(x,[1 2],0.05,0);
        m_Statistic(i,j) = p;
    end  
end

m_StatisticSignificant = (m_Statistic<=0.01);
m_StatisticSignificant = m_StatisticSignificant.*m_Statistic;
f_ImageMatrix(m_StatisticSignificant,v_time,v_FreqAxis, [], [], [], 0);
colorbar
title('p-Values')

%% Plot the mean of all Subjects for all channels

vlims = [min(min(MeanCueEpochs(:),MeanShamEpochs(:))),...
    max(max(MeanCueEpochs(:),MeanShamEpochs(:)))]*0.1;
subplot(2,1,1)
f_ImageMatrix(MeanCueSubjects,v_time,v_FreqAxis, vlims, [], [], 0);
title('Cue Epochs')
colorbar
subplot(2,1,2)
f_ImageMatrix(MeanShamSubjects,v_time,v_FreqAxis, vlims, [], [], 0);
title('Sham Epochs')
colorbar
suptitle('AllChans Cue Night')
%saveas(gcf,'AllChans_Cue_Night.png') 


figure
s_LastInd = find(v_FreqAxis <=8, 1);
s_FirstInd = find(v_FreqAxis <=12, 1);
m_TFFreqBands = ...
    mean(MeanCueEpochs(:,s_FirstInd:s_LastInd, :),2);
m_TFFreqBands = reshape(m_TFFreqBands,...
    [size(m_TFFreqBands,1),size(m_TFFreqBands,3)]);
plot(v_time,m_TFFreqBands,'k','linewidth',0.5)
hold on
CuePlot = mean(m_TFFreqBands,1);

m_TFFreqBands = ...
    mean(MeanShamEpochs(:,s_FirstInd:s_LastInd, :),2);
m_TFFreqBands = reshape(m_TFFreqBands,...
    [size(m_TFFreqBands,1),size(m_TFFreqBands,3)]);
plot(v_time,m_TFFreqBands,'r','linewidth',0.5)
hold on
b1 = plot(v_time,CuePlot,'k','linewidth',3,...
    'DisplayName','Cue Mean');
b2 = plot(v_time,mean(m_TFFreqBands,1),'r','linewidth',3,...
    'DisplayName','Sham Mean');
legend ([b1,b2])
title('AllChan Cue Night')

hold off
%saveas(gcf,'AllChan_Cue_Night_Delta.png')

%% Get analysis for each subject (including statistics)

for subj = 1:numel(subjects)
    vlims = [min(min(MeanCueEpochs(:),MeanShamEpochs(:))),...
        max(max(MeanCueEpochs(:),MeanShamEpochs(:)))]*0.3;
    subplot(2,1,1)
    f_ImageMatrix(squeeze(MeanCueEpochs(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Cue Epochs')
    colorbar
    subplot(2,1,2)
    f_ImageMatrix(squeeze(MeanShamEpochs(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Sham Epochs')
    colorbar
    suptitle('AllChans Cue Night')
    %saveas(gcf,'AllChans_Cue_Night.png')
    
end

%% Statistical analysis by subject, all channels
v_time = -15:0.5:15;
subjects = fieldnames(m_GaborCue_total);
baselineTime = 15;
s_TimeStep = 0.5;

for subj = 1:numel(subjects)
    
    DataSubj = m_GaborCue_total.(string(subjects(subj)));
    
    baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:,:),2);
    DataSubjBasCorr = DataSubj-baseline;
    
    MeanCueChannels = mean(DataSubjBasCorr,3);
    MeanCueChannels = reshape(MeanCueChannels,...
        [size(MeanCueChannels,1),size(MeanCueChannels,2),...
        size(MeanCueChannels,4)]);

    DataSubj = m_GaborSham_total.(string(subjects(subj)));
    
    baseline = mean(DataSubj(:,1:(baselineTime/s_TimeStep),:,:),2);
    DataSubjBasCorr = DataSubj-baseline;
   
    MeanShamChannels = mean(DataSubjBasCorr,3);
    MeanShamChannels = reshape(MeanShamChannels,...
        [size(MeanShamChannels,1),size(MeanShamChannels,2),...
        size(MeanShamChannels,4)]);
    
    for i = 1:size(MeanCueChannels,1)
        for j = 1:size(MeanCueChannels,2)
            cond1 = MeanCueChannels(i,j,:);
            cond2 = MeanShamChannels(i,j,1:floor(size(MeanShamChannels,3)/2));
            minlength = min(size(cond1,3),size(cond2,3));
            x = [cond1(:,:,1:minlength),cond2(:,:,1:minlength)];
            [h,p,sigPairs] = ttest_bonf(x,[1 2],0.05,0);
            m_Statistic(i,j) = p;
        end
    end
    
    m_StatisticSignificant = (m_Statistic<=0.05);
    %m_StatisticSignificant = m_StatisticSignificant.*m_Statistic;
    figure
    f_ImageMatrix(m_StatisticSignificant,v_time,v_FreqAxis, [], [], [], 0);
    colorbar
    title('p-Values First Trials')
    
%     for i = 1:size(MeanCueChannels,1)
%         for j = 1:size(MeanCueChannels,2)
%             cond1 = MeanCueChannels(i,j,floor(size(MeanCueChannels,3)/2)+1:end);
%             cond2 = MeanShamChannels(i,j,floor(size(MeanShamChannels,3)/2)+1:end);
%             minlength = min(size(cond1,3),size(cond2,3));
%             x = [cond1(:,:,1:minlength),cond2(:,:,1:minlength)];
%             [h,p,sigPairs] = ttest_bonf(x,[1 2],0.05,0);
%             m_Statistic(i,j) = p;
%         end
%     end
%     
%     m_StatisticSignificant = (m_Statistic<=0.05);
%     %m_StatisticSignificant = m_StatisticSignificant.*m_Statistic;
%     figure
%     f_ImageMatrix(m_StatisticSignificant,v_time,v_FreqAxis, [], [], [], 0);
%     colorbar
%     title('p-Values Last Trials')
end



