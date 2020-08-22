%Load time-frequency transforms
m_Gabor.CueNight = load('m_GaborCuetotal_JensData.mat');
m_Gabor.PlaceboNight = load('m_GaborPlacebototal_JensData.mat');
load('Channels.mat')

subjectsCue = fieldnames(m_Gabor.CueNight.m_GaborCue_total);
subjectsPlac = fieldnames(m_Gabor.PlaceboNight.m_GaborPlacebo_total);
baselineTime = 15;
s_TimeStep = 0.5;
v_time = -15:0.5:15;
v_FreqAxis = m_Gabor.CueNight.v_FreqAxis;

v_GoodLearners = [1 3 9 10 14 15 18 19 22];
v_BadLearners = [2 4 5 6 7 8 11 12 13 16 17 20 21 23];

%-------------------------------------------------------------------------
% create complete arrays for each condition
% size: Subj x Freq x TimePnts x Channels x Trials
%-------------------------------------------------------------------------
s_MinEpochsOdor = NaN;
s_MinEpochsSham = NaN;

%--- Get min number of epochs to create equal arrays for Cue Night -------
for subj = 1:numel(subjectsCue)
    s_MinEpochsOdor = min(s_MinEpochsOdor,...
        size(m_Gabor.CueNight.m_GaborCue_total.(string(subjectsCue(subj))),4));
    s_MinEpochsSham = min(s_MinEpochsSham,...
        size(m_Gabor.CueNight.m_GaborSham_total.(string(subjectsCue(subj))),4));
end

%----- Fill the arrays for Cue Night --------
for subj = 1:numel(subjectsCue) 
    %Cue Night Cue Odor
    m_GaborCue(subj,:,:,:,:) = ...
        m_Gabor.CueNight.m_GaborCue_total.(string(subjectsCue(subj)))(:,:,:,1:s_MinEpochsOdor);
    
    %Cue Night Sham Odor
    m_GaborSham_CN(subj,:,:,:,:) = ...
        m_Gabor.CueNight.m_GaborSham_total.(string(subjectsCue(subj)))(:,:,:,1:s_MinEpochsSham);
end

s_MinEpochsOdor = NaN;
s_MinEpochsSham = NaN;
%--- Get min number of epochs to create equal arrays for Placebo Night ---
for subj = 1:numel(subjectsPlac)
    s_MinEpochsOdor = min(s_MinEpochsOdor,...
        size(m_Gabor.PlaceboNight.m_GaborPlacebo_total.(string(subjectsPlac(subj))),4));
    s_MinEpochsSham = min(s_MinEpochsSham,...
        size(m_Gabor.PlaceboNight.m_GaborSham_total.(string(subjectsPlac(subj))),4));
end

%----- Fill the arrays for Placebo Night --------
for subj = 1:numel(subjectsPlac)
    %Placebo Night Placebo Odor
    m_GaborPlacebo(subj,:,:,:,:) = ...
        m_Gabor.PlaceboNight.m_GaborPlacebo_total.(string(subjectsPlac(subj)))(:,:,:,1:s_MinEpochsOdor);
    
    %Placebo Night Sham Odor
    m_GaborSham_PN(subj,:,:,:,:) = ...
        m_Gabor.PlaceboNight.m_GaborSham_total.(string(subjectsPlac(subj)))(:,:,:,1:s_MinEpochsSham); 
end
clear m_Gabor
%% 
%--------------------------------------------------------------------------
% Baseline Correction
%--------------------------------------------------------------------------

%Cue Night Cue Odor
m_GaborCue_BasCorr = ...
    baselinecorrection(m_GaborCue,baselineTime,s_TimeStep);

%Cue Night Sham Odor
m_GaborSham_CN_BasCorr = ...
    baselinecorrection(m_GaborSham_CN,baselineTime,s_TimeStep);

%Placebo Night Placebo Odor
m_GaborPlacebo_BasCorr = ...
    baselinecorrection(m_GaborPlacebo,baselineTime,s_TimeStep);

%Placebo Night Sham Odor
m_GaborSham_PN_BasCorr = ...
    baselinecorrection(m_GaborSham_PN,baselineTime,s_TimeStep);

%% 
%--------------------------------------------------------------------------
% Mean of channels for all epochs for each subject
%--------------------------------------------------------------------------
s_MinSize = NaN;
%Cue Night Cue Odor
MeanChannels_Cue = squeeze(nanmean(m_GaborCue_BasCorr,4));
s_MinSize = min(s_MinSize,size(MeanChannels_Cue,4));

%Cue Night Sham Odor
MeanChannels_Sham_CN = squeeze(nanmean(m_GaborSham_CN_BasCorr,4));
s_MinSize = min(s_MinSize,size(MeanChannels_Sham_CN,4));

%Equal number of Epochs
MeanChannels_Cue = MeanChannels_Cue(:,:,:,1:s_MinSize);
MeanChannels_Sham_CN = MeanChannels_Sham_CN(:,:,:,1:s_MinSize);

s_MinSize = NaN;
%Placebo Night Placebo Odor
MeanChannels_Placebo = squeeze(nanmean(m_GaborPlacebo_BasCorr,4));
s_MinSize = min(s_MinSize,size(MeanChannels_Placebo,4));

%Placebo Night Sham Odor
MeanChannels_Sham_PN= squeeze(nanmean(m_GaborSham_PN_BasCorr,4));
s_MinSize = min(s_MinSize,size(MeanChannels_Sham_PN,4));

%Equal number of Epochs
MeanChannels_Placebo = MeanChannels_Placebo(:,:,:,1:s_MinSize);
MeanChannels_Sham_PN = MeanChannels_Sham_PN(:,:,:,1:s_MinSize);

%--------------------------------------------------------------------------
% Mean of epochs and channels for each subject
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Ch_Cue = squeeze(nanmean(MeanChannels_Cue,4));

%Cue Night Sham Odor
Mean_Ep_Ch_Sham_CN = squeeze(nanmean(MeanChannels_Sham_CN,4));

%Placebo Night Placebo Odor
Mean_Ep_Ch_Placebo = squeeze(nanmean(MeanChannels_Placebo,4));

%Placebo Night Sham Odor
Mean_Ep_Ch_Sham_PN = squeeze(nanmean(MeanChannels_Sham_PN,4));

%--------------------------------------------------------------------------
% Mean of epochs only
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Cue = squeeze(nanmean(m_GaborCue_BasCorr,5));

%Cue Night Sham Odor
Mean_Ep_Sham_CN = squeeze(nanmean(m_GaborSham_CN_BasCorr,5));

%Placebo Night Placebo Odor
Mean_Ep_Placebo = squeeze(nanmean(m_GaborPlacebo_BasCorr,5));

%Placebo Night Sham Odor
Mean_Ep_Sham_PN = squeeze(nanmean(m_GaborSham_PN_BasCorr,5));
%--------------------------------------------------------------------------
% Mean of epochs and subjects for each channel
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Subj_Cue = squeeze(nanmean(Mean_Ep_Cue,1));

%Cue Night Sham Odor
Mean_Ep_Subj_Sham_CN = squeeze(nanmean(Mean_Ep_Sham_CN,1));

%Placebo Night Placebo Odor
Mean_Ep_Subj_Placebo = squeeze(nanmean(Mean_Ep_Placebo,1));

%Placebo Night Sham Odor
Mean_Ep_Subj_Sham_PN = squeeze(nanmean(Mean_Ep_Sham_PN,1));

%--------------------------------------------------------------------------
% Mean All (All Subjects,epochs and Channels)
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_All_Cue = squeeze(nanmean(Mean_Ep_Ch_Cue,1));

%Cue Night Sham Odor
Mean_All_Sham_CN = squeeze(nanmean(Mean_Ep_Ch_Sham_CN,1));

%Placebo Night Placebo Odor
Mean_All_Placebo = squeeze(nanmean(Mean_Ep_Ch_Placebo,1));

%Placebo Night Sham Odor
Mean_All_Sham_PN = squeeze(nanmean(Mean_Ep_Ch_Sham_PN,1));

%% Statistical analysis for each subject across epochs (mean channels)
% Condition1 and Condition2 should be in a shape of freq x Time x Epochs
for subj = 1:size(MeanChannels_Cue,1) 
    %----- For Cue Night----
    Cond1 = squeeze(MeanChannels_Cue(subj,:,:,:));
    Cond2 = squeeze(MeanChannels_Sham_CN(subj,:,:,:));
    [m_Statis_pval_Ep_Ch_Cue(subj,:,:),...
        m_SigP_Ep_Ch_Cue(subj,:,:),...
        m_SigPvalue_Ep_Ch_Cue(subj,:,:)] = ...
        GetStatistics(Cond1,Cond2);
    
    %----- For Placebo Night----
    Cond1 = squeeze(MeanChannels_Placebo(subj,:,:,:));
    Cond2 = squeeze(MeanChannels_Sham_PN(subj,:,:,:));
    [m_Statis_pval_Ep_Ch_Placebo(subj,:,:),...
        m_SigP_Ep_Ch_Placebo(subj,:,:),...
        m_SigPvalue_Ep_Ch_Placebo(subj,:,:)] = ...
        GetStatistics(Cond1,Cond2);
end

%% Statistical analysis across subjects (mean channels, mean epochs)
% Condition1 and Condition2 should be in a shape of freq x Time x Subj

%----- For Cue Night----
[m_Statis_pval_All_Cue, m_SigP_All_Cue, m_SigPvalue_All_Cue] = ...
    GetStatistics(permute(Mean_Ep_Ch_Cue,[2,3,1]),...
    permute(Mean_Ep_Ch_Sham_CN,[2,3,1]));

%----- For Placebo Night----
[m_Statis_pval_all_Placebo, m_SigP_All_Placebo, m_SigPvalue_All_Placebo] = ...
    GetStatistics(permute(Mean_Ep_Ch_Placebo,[2,3,1]),...
    permute(Mean_Ep_Ch_Sham_PN,[2,3,1]));

%% Statistical analysis for each channel across subjects (mean epochs)
% Condition1 and Condition2 should be in a shape of freq x Time x Subj
for chan = 1:size(Mean_Ep_Cue,4) 
    %----- For Cue Night----
    Cond1 = squeeze(Mean_Ep_Cue(:,:,:,chan));
    Cond2 = squeeze(Mean_Ep_Sham_CN(:,:,:,chan));
    [m_Statis_pval_Ep_Subj_Cue(chan,:,:),...
        m_SigP_Ep_Subj_Cue(chan,:,:),...
        m_SigPvalue_Ep_Subj_Cue(chan,:,:)] = ...
        GetStatistics(permute(Cond1,[2,3,1]),...
        permute(Cond2,[2,3,1]));
    
    %----- For Placebo Night----
    Cond1 = squeeze(Mean_Ep_Placebo(:,:,:,chan));
    Cond2 = squeeze(Mean_Ep_Sham_PN(:,:,:,chan));
    [m_Statis_pval_Ep_Subj_Placebo(chan,:,:),...
        m_SigP_Ep_Subj_Placebo(chan,:,:),...
        m_SigPvalue_Ep_Subj_Placebo(chan,:,:)] = ...
        GetStatistics(permute(Cond1,[2,3,1]),...
        permute(Cond2,[2,3,1]));
end

%% Analysis only on the Delta frequency band
%---------------------------------------------------------------------
% For each subject across epochs (mean channels)
%---------------------------------------------------------------------
m_Delta_MeanChannels_Cue = [];
m_Delta_MeanChannels_Sham_CN = [];
m_Delta_MeanChannels_Placebo = [];
m_Delta_MeanChannels_Sham_PN = [];

for subj = 1:size(MeanChannels_Cue,1)
    %----- For Cue Night----
    DataSubjCue = squeeze(MeanChannels_Cue(subj,:,:,:));
    DataSubjCue = permute(DataSubjCue,[3,1,2]);
    m_Delta_MeanChannels_Cue(subj,:,:) = ...
        DeltaBandCalculation(DataSubjCue,v_FreqAxis);
    DataSubjSham = squeeze(MeanChannels_Sham_CN(subj,:,:,:));
    DataSubjSham = permute(DataSubjSham,[3,1,2]);
    m_Delta_MeanChannels_Sham_CN(subj,:,:) = ...
        DeltaBandCalculation(DataSubjSham,v_FreqAxis);
    
    %----- For Placebo Night----
    DataSubjPlacebo = squeeze(MeanChannels_Placebo(subj,:,:,:));
    DataSubjPlacebo = permute(DataSubjPlacebo,[3,1,2]);
    m_Delta_MeanChannels_Placebo(subj,:,:) = ...
        DeltaBandCalculation(DataSubjPlacebo,v_FreqAxis);
    DataSubjSham = squeeze(MeanChannels_Sham_PN(subj,:,:,:));
    DataSubjSham = permute(DataSubjSham,[3,1,2]);
    m_Delta_MeanChannels_Sham_PN(subj,:,:) = ...
        DeltaBandCalculation(DataSubjSham,v_FreqAxis);
end

%---------------------------------------------------------------------
% Across subjects (mean channels, mean epochs)
%---------------------------------------------------------------------
%% ----- For Cue Night ----
m_Delta_MeanAll_Cue = ...
    DeltaBandCalculation(Mean_Ep_Ch_Cue,v_FreqAxis);
m_Delta_MeanAll_Sham_CN = ...
    DeltaBandCalculation(Mean_Ep_Ch_Sham_CN,v_FreqAxis);

%----- For Placebo Night ----
m_Delta_MeanAll_Placebo = ...
    DeltaBandCalculation(Mean_Ep_Ch_Placebo,v_FreqAxis);
m_Delta_MeanAll_Sham_PN = ...
    DeltaBandCalculation(Mean_Ep_Ch_Sham_PN,v_FreqAxis);
%%
%---------------------------------------------------------------------
% For each channel Across subjects (mean epochs)
%---------------------------------------------------------------------
m_Delta_MeanEpochs_Cue = [];
m_Delta_MeanEpochs_Sham_CN = [];
m_Delta_MeanEpochs_Placebo = [];
m_Delta_MeanEpochs_Sham_PN = [];

for chan = 1:size(Mean_Ep_Cue,4)
    %----- For Cue Night----
    Cond1 = squeeze(Mean_Ep_Cue(:,:,:,chan));
    Cond2 = squeeze(Mean_Ep_Sham_CN(:,:,:,chan));
    m_Delta_MeanEpochs_Cue(chan,:,:) = ...
        DeltaBandCalculation(Cond1,v_FreqAxis);
    m_Delta_MeanEpochs_Sham_CN(chan,:,:) = ...
        DeltaBandCalculation(Cond2,v_FreqAxis);
    
    %----- For Placebo Night----
    Cond1 = squeeze(Mean_Ep_Placebo(:,:,:,chan));
    Cond2 = squeeze(Mean_Ep_Sham_PN(:,:,:,chan));
    m_Delta_MeanEpochs_Placebo(chan,:,:) = ...
        DeltaBandCalculation(Cond1,v_FreqAxis);
    m_Delta_MeanEpochs_Sham_PN(chan,:,:) = ...
        DeltaBandCalculation(Cond2,v_FreqAxis);
end

%% Plot
PlotTimeFrequency()
PlotStatistics()
plotDeltaBandCalculation()
PerformWilcoxon()

%% Wilcoxon test for Delta Band analyses
function PerformWilcoxon()
%%
v_time2 = -15:2.5:15-2.5;
eeglab nogui
%---------------------------------------------------------------------
% For each subject across epochs (mean channels)
%---------------------------------------------------------------------
for subj = 1:size(m_Delta_MeanChannels_Cue,1)
    %----- For Cue Night----
    f_WilcTest(strcat('Delta Cue Night Subj',num2str(subj)),...
        'Time(sec)',' ','Cue','Sham',...
        squeeze(m_Delta_MeanChannels_Cue(subj,:,:)),...
        squeeze(m_Delta_MeanChannels_Sham_CN(subj,:,:)),...
        v_time,'-r')
    saveas(gcf,strcat('Cue_Night_Wilcoxon_Subj',num2str(subj),'.png'))
    
    %----- For Placebo Night----
    f_WilcTest(strcat('Delta Placebo Night Subj',num2str(subj)),...
        'Time(sec)',' ','Placebo','Sham',...
        squeeze(m_Delta_MeanChannels_Placebo(subj,:,:)),...
        squeeze(m_Delta_MeanChannels_Sham_PN(subj,:,:)),...
        v_time,'-b')
    saveas(gcf,strcat('Placebo_Night_Wilcoxon_Subj',num2str(subj),'.png'))
    close all
end
%---------------------------------------------------------------------
% Across subjects (mean channels, mean epochs)
%---------------------------------------------------------------------

%----- For Cue Night----
stats_m_Delta_MeanAll_cue = ...
    f_WilcTest('Delta Band (0.5 - 4Hz) for Cue Night','Time(sec)',' ','Cue','Sham',...
    m_Delta_MeanAll_Cue,m_Delta_MeanAll_Sham_CN,v_time,'-r');
ylim([-30 30])
%saveas(gcf,strcat('Cue_Night_Wilcoxon.png'))

%----- For Placebo Night----
stats_m_Delta_MeanAll_placebo = ...
    f_WilcTest('Delta Band (0.5 - 4Hz) for Placebo Night','Time(sec)',' ','Placebo','Sham',...
    m_Delta_MeanAll_Placebo,m_Delta_MeanAll_Sham_PN,v_time,'-b');
ylim([-30 30])
%saveas(gcf,strcat('Placebo_Night_Wilcoxon.png'))
close all

%% ---------------------------------------------------------------------
% For each channel Across subjects (mean epochs)
%---------------------------------------------------------------------

for chan = 1:129%size(m_Delta_MeanEpochs_Cue,1)%
    
    try 
    %----- For Cue Night----
    stats_Delta_MeanEpochs_Cue = ...
        f_WilcTest(strcat('Delta Band (0.5 - 4Hz) for Cue Night Channel ',channels(chan)),...
        'Time(sec)',' ','Cue','Sham',...
        squeeze(m_Delta_MeanEpochs_Cue(chan,:,:)),...
        squeeze(m_Delta_MeanEpochs_Sham_CN(chan,:,:)),...
        v_time,'-r');
    ylim([-30 30])
    %hold on
    %yyaxis right
    %plot(v_time2,[stats_Delta_MeanEpochs_Cue.pval],'-g','DisplayName','p-value')
    %ylim([0 1])
    saveas(gcf,strcat('Cue_Night_Wilcoxon_Chan',num2str(chan),'.png'))
    
    %----- For Placebo Night----
    stats_Delta_MeanEpochs_Placebo = ...
        f_WilcTest(strcat('Delta Band (0.5 - 4Hz) for Placebo Night Channel ',channels(chan)),...
        'Time(sec)',' ','Placebo','Sham',...
        squeeze(m_Delta_MeanEpochs_Placebo(chan,:,:)),...
        squeeze(m_Delta_MeanEpochs_Sham_PN(chan,:,:)),...
        v_time,'-b');
    ylim([-30 30])
    %hold on
    %yyaxis right
    %plot(v_time2,[stats_Delta_MeanEpochs_Placebo.pval],'-g','DisplayName','p-value')
    %ylim([0 1])
    saveas(gcf,strcat('Placebo_Night_Wilcoxon_Chan',num2str(chan),'.png'))
    close all
    catch
        continue
    end
end
end

%% PLOTS!!!! Time-Frequency
function PlotTimeFrequency()
%--------------------------------------------------------------------------
% Plot Mean of epochs and channels for each subject
%--------------------------------------------------------------------------
for subj = 1:size(Mean_Ep_Ch_Cue,1)
    vlims = [min(min(Mean_Ep_Ch_Cue(:),Mean_Ep_Ch_Sham_CN(:))),...
        max(max(Mean_Ep_Ch_Cue(:),Mean_Ep_Ch_Sham_CN(:)))]*0.3;
    
    %----- For Cue Night----
    figure
    subplot(2,1,1)
    f_ImageMatrix(squeeze(Mean_Ep_Ch_Cue(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Cue Epochs')
    colorbar
    subplot(2,1,2)
    f_ImageMatrix(squeeze(Mean_Ep_Ch_Sham_CN(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Sham Epochs')
    colorbar
    suptitle(strcat('AllChans Cue Night Subj',num2str(subj)))
    saveas(gcf,strcat('AllChans_Cue_Night_Subj',num2str(subj),'.png'))
    
    
    %----- For Placebo Night----
    figure
    subplot(2,1,1)
    f_ImageMatrix(squeeze(Mean_Ep_Ch_Placebo(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Placebo Epochs')
    colorbar
    subplot(2,1,2)
    f_ImageMatrix(squeeze(Mean_Ep_Ch_Sham_PN(subj,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Sham Epochs')
    colorbar
    suptitle(strcat('AllChans Placebo Night Subj',num2str(subj)))
    saveas(gcf,strcat('AllChans_Placebo_Night_Subj',num2str(subj),'.png'))
    
    close all
end

%--------------------------------------------------------------------------
% Plot Mean of epochs and subjects for each channel
%--------------------------------------------------------------------------
for chan = 1:size(Mean_Ep_Subj_Cue,3)
    vlims = [min(min(Mean_Ep_Subj_Cue(:),Mean_Ep_Subj_Sham_CN(:))),...
        max(max(Mean_Ep_Subj_Cue(:),Mean_Ep_Subj_Sham_CN(:)))]*0.3;
    
    %----- For Cue Night----
    figure
    subplot(2,1,1)
    f_ImageMatrix(squeeze(Mean_Ep_Subj_Cue(:,:,chan)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Cue Epochs')
    colorbar
    subplot(2,1,2)
    f_ImageMatrix(squeeze(Mean_Ep_Subj_Sham_CN(:,:,chan)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Sham Epochs')
    colorbar
    suptitle(strcat('AllSubj Cue Night Channel',channels(chan)))
    saveas(gcf,strcat('AllSubj_Cue_Night_Chan',channels(chan),'.png'))
   
    %----- For Placebo Night----
    figure
    subplot(2,1,1)
    f_ImageMatrix(squeeze(Mean_Ep_Subj_Placebo(:,:,chan)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Placebo Epochs')
    colorbar
    subplot(2,1,2)
    f_ImageMatrix(squeeze(Mean_Ep_Subj_Sham_PN(:,:,chan)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Sham Epochs')
    colorbar
    suptitle(strcat('AllSubj Placebo Night Channel',channels(chan)))
    saveas(gcf,strcat('AllSubj_Placebo_Night_Chan',channels(chan),'.png'))
    
    close all
end

%--------------------------------------------------------------------------
% Plot Mean All (All Subjects,epochs and Channels)
%--------------------------------------------------------------------------

vlims = [min(min(Mean_All_Cue(:),Mean_All_Sham_CN(:))),...
    max(max(Mean_All_Cue(:),Mean_All_Sham_PN(:)))]*0.3;

%----- For Cue Night----
figure
subplot(2,1,1)
f_ImageMatrix(Mean_All_Cue,...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Cue Epochs')
colorbar
subplot(2,1,2)
f_ImageMatrix(Mean_All_Sham_CN,...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Sham Epochs')
colorbar
suptitle('Mean All Cue Night')
saveas(gcf,'All_Cue_Night.png')


%----- For Placebo Night----
figure
subplot(2,1,1)
f_ImageMatrix(Mean_All_Placebo,...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Placebo Epochs')
colorbar
subplot(2,1,2)
f_ImageMatrix(Mean_All_Sham_PN,...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Sham Epochs')
colorbar
suptitle('Mean All Placebo Night')
saveas(gcf,'All_Placebo_Night.png')

close all
end
function PlotStatistics()
%% Plot Statistics
vlims = [0.95 1];
%--------------------------------------------------------------------------
% Plot Statistical analysis for each subject across epochs (mean channels)
%--------------------------------------------------------------------------
for subj = 1:size(m_SigPvalue_Ep_Ch_Cue,1)
    
    figure
    %----- For Cue Night----
    subplot(2,1,1)
    f_ImageMatrix(squeeze(m_SigPvalue_Ep_Ch_Cue(subj,:,:)),...
        v_time,v_FreqAxis,vlims , [], [], 0);
    title('Cue Night p-values')
    colorbar
    %----- For Placebo Night----
    subplot(2,1,2)
    f_ImageMatrix(squeeze(m_SigPvalue_Ep_Ch_Placebo(subj,:,:)),...
        v_time,v_FreqAxis,vlims, [], [], 0);
    title('Placebo Night p-values')
    colorbar
    suptitle(strcat('p-values Subj',num2str(subj)))
    saveas(gcf,strcat('p-values_Subj',num2str(subj),'.png'))
end
close all
%--------------------------------------------------------------------------
%Plot Statistical analysis across subjects (mean channels, mean epochs)
%--------------------------------------------------------------------------
figure
%----- For Cue Night----
subplot(2,1,1)
f_ImageMatrix(squeeze(m_SigPvalue_All_Cue),...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Cue Night p-values')
colorbar
%----- For Placebo Night----
subplot(2,1,2)
f_ImageMatrix(squeeze(m_SigPvalue_All_Placebo),...
    v_time,v_FreqAxis, vlims, [], [], 0);
title('Placebo Night p-values')
colorbar
suptitle('p-values All')
saveas(gcf,'p-values.png')

%--------------------------------------------------------------------------
% Plot Statistical analysis for each channel across subjects (mean epochs)
%--------------------------------------------------------------------------
for chan = 1:size(Mean_Ep_Subj_Cue,3)
    
    figure
    %----- For Cue Night----
    subplot(2,1,1)
    f_ImageMatrix(squeeze(m_SigPvalue_Ep_Subj_Cue(chan,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Cue Night p-values')
    colorbar
    %----- For Placebo Night----
    subplot(2,1,2)
    f_ImageMatrix(squeeze(m_SigPvalue_Ep_Subj_Placebo(chan,:,:)),...
        v_time,v_FreqAxis, vlims, [], [], 0);
    title('Placebo Night p-values')
    colorbar
    suptitle(strcat('p-values chan',channels(chan)))
    saveas(gcf,strcat('p-values_chan',channels(chan),'.png'))
end
close all
end
function plotDeltaBandCalculation()
%% Plot Delta band Calculation
% --------------------------------------------------------------------------
% Plot Delta band Calculation for each subject across epochs (mean channels)
%--------------------------------------------------------------------------
for subj = 1:size(m_Delta_MeanChannels_Cue,1)
    %----- For Cue Night ----
    Cond1 = squeeze(m_Delta_MeanChannels_Cue(subj,:,:));
    Cond2 = squeeze(m_Delta_MeanChannels_Sham_CN(subj,:,:));
    
    plot(v_time,Cond1,'r','linewidth',0.5)
    hold on
    CuePlot = nanmean(Cond1,1);
    
    plot(v_time,Cond2,'k','linewidth',0.5)
    hold on
    b1 = plot(v_time,CuePlot,'r','linewidth',3,...
        'DisplayName','Cue Mean');
    b2 = plot(v_time,nanmean(Cond2,1),'k','linewidth',3,...
        'DisplayName','Sham Mean');
    legend ([b1,b2])
    title(strcat('Subj',num2str(subj),'Cue Night'))
    hold off
    saveas(gcf,strcat('Subj',num2str(subj),'_Cue_Night_Delta.png'))
    
    %----- For Placebo Night ----
    Cond1 = squeeze(m_Delta_MeanChannels_Placebo(subj,:,:));
    Cond2 = squeeze(m_Delta_MeanChannels_Sham_PN(subj,:,:));
    
    plot(v_time,Cond1,'b','linewidth',0.5)
    hold on
    PlaceboPlot = nanmean(Cond1,1);
    
    plot(v_time,Cond2,'k','linewidth',0.5)
    hold on
    b1 = plot(v_time,PlaceboPlot,'b','linewidth',3,...
        'DisplayName','Placebo Mean');
    b2 = plot(v_time,nanmean(Cond2,1),'k','linewidth',3,...
        'DisplayName','Sham Mean');
    legend ([b1,b2])
    title(strcat('Subj',num2str(subj),'Placebo Night'))
    hold off
    saveas(gcf,strcat('Subj',num2str(subj),'_Placebo_Night_Delta.png'))
end
close all

%---------------------------------------------------------------------
% Plot Delta band Calculation Across subjects (mean channels, mean epochs)
%---------------------------------------------------------------------
%----- For Cue Night ----
plot(v_time,m_Delta_MeanAll_Cue,'r','linewidth',0.5)
hold on
CuePlot = nanmean(m_Delta_MeanAll_Cue,1);

plot(v_time,m_Delta_MeanAll_Sham_CN,'k','linewidth',0.5)
hold on
b1 = plot(v_time,CuePlot,'r','linewidth',3,...
    'DisplayName','Cue Mean');
b2 = plot(v_time,nanmean(m_Delta_MeanAll_Sham_CN,1),'k','linewidth',3,...
    'DisplayName','Sham Mean');
legend ([b1,b2])
title('Cue Night')
ylim([-50 50])
hold off
saveas(gcf,'Cue_Night_Delta.png')

%----- For Placebo Night ----
plot(v_time,m_Delta_MeanAll_Placebo,'b','linewidth',0.5)
hold on
PlaceboPlot = nanmean(m_Delta_MeanAll_Placebo,1);

plot(v_time,m_Delta_MeanAll_Sham_PN,'k','linewidth',0.5)
hold on
b1 = plot(v_time,PlaceboPlot,'b','linewidth',3,...
    'DisplayName','Placebo Mean');
b2 = plot(v_time,nanmean(m_Delta_MeanAll_Sham_PN,1),'k','linewidth',3,...
    'DisplayName','Sham Mean');
legend ([b1,b2])
title('Placebo Night')
ylim([-50 50])
hold off
saveas(gcf,'Placebo_Night_Delta.png')
close all

% --------------------------------------------------------------------------
% Plot Delta band Calculation for each channel Across subjects (mean epochs)
%--------------------------------------------------------------------------
for chan = 1:size(Mean_Ep_Cue,4)
    %----- For Cue Night ----
    Cond1 = squeeze(m_Delta_MeanEpochs_Cue(chan,:,:));
    Cond2 = squeeze(m_Delta_MeanEpochs_Sham_CN(chan,:,:));
    
    plot(v_time,Cond1,'r','linewidth',0.5)
    hold on
    CuePlot = nanmean(Cond1,1);
    
    plot(v_time,Cond2,'k','linewidth',0.5)
    hold on
    b1 = plot(v_time,CuePlot,'r','linewidth',3,...
        'DisplayName','Cue Mean');
    b2 = plot(v_time,nanmean(Cond2,1),'k','linewidth',3,...
        'DisplayName','Sham Mean');
    legend ([b1,b2])
    title(strcat('Channel',channels(chan),'Cue Night'))
    hold off
    saveas(gcf,strcat('Chan',channels(chan),'_Cue_Night_Delta.png'))
    
    %----- For Placebo Night ----
    Cond1 = squeeze(m_Delta_MeanEpochs_Placebo(chan,:,:));
    Cond2 = squeeze(m_Delta_MeanEpochs_Sham_PN(chan,:,:));
    
    plot(v_time,Cond1,'b','linewidth',0.5)
    hold on
    PlaceboPlot = nanmean(Cond1,1);
    
    plot(v_time,Cond2,'k','linewidth',0.5)
    hold on
    b1 = plot(v_time,PlaceboPlot,'b','linewidth',3,...
        'DisplayName','Placebo Mean');
    b2 = plot(v_time,nanmean(Cond2,1),'k','linewidth',3,...
        'DisplayName','Sham Mean');
    legend ([b1,b2])
    title(strcat('Channel',channels(chan),'Placebo Night'))
    hold off
    saveas(gcf,strcat('Chan',channels(chan),'_Placebo_Night_Delta.png'))
end
close all
end
%% 
function m_TFFreqBands = DeltaBandCalculation(TimeFrequencyMatrix,v_FreqAxis)

%TimeFrequencyMatrix shape shoudl be (var x freq x time)
%Establish the frequencies of interest
MinFreq = 0.5;
MaxFreq = 4;

% Determine the position of the frequencies in the Frequency vector
s_LastInd = find(v_FreqAxis <= MinFreq, 1);
s_FirstInd = find(v_FreqAxis <= MaxFreq, 1);

%Calculate the mean of power in that band of frequency
m_TFFreqBands = ...
    squeeze(nanmean(TimeFrequencyMatrix(:,s_FirstInd:s_LastInd, :),2));
end
%%
function Baseline_Corr = baselinecorrection(Data,baselineTime,s_TimeStep)
baseline = nanmean(Data(:,:,1:(baselineTime/s_TimeStep),:,:),3);
Baseline_Corr = Data-baseline;
end
%%
% Condition1 and Condition2 should be given in a shape of freq x Time x
% Variation
function [m_Statistics, m_sigp, m_sigpvalue] = ...
    GetStatistics(Condition1, Condition2)

m_Statistics = [];
for i = 1:size(Condition1,1)
    for j = 1:size(Condition1,2)
        cond1 = Condition1(i,j,:);
        cond2 = Condition2(i,j,:);
        minlength = min(size(cond1,3),size(cond2,3));
        x = [squeeze(cond1(:,:,1:minlength)),squeeze(cond2(:,:,1:minlength))];
        [~,p,~] = ttest_bonf(x,[1 2],0.05,0);
        m_Statistics(i,j) = p;
    end
end
%-- Get significant p-values (<=0.05)
m_sigp = (m_Statistics <=0.05);

%-- Multiply the boolean of significant p-values by the p-value
m_sigpvalue = m_sigp.*(1-m_Statistics);
end


%%
% v_time2 = -15+1:2.5:15;
% for freq = 1:size(Mean_Ep_Ch_Cue,2)
%     name = strcat('Band',num2str(freq));
%     For cue Night
%     c1 = squeeze(Mean_Ep_Ch_Cue(:,freq,:));
%     c2 = squeeze(Mean_Ep_Ch_Sham_CN(:,freq,:));
%     statsCue = test_wilcoxon_cvar(c1, c2,'Cue','sham','signrank',1,1);
%     pvalCue(freq,:) = statsCue.pval;
%     
%     For Placebo Night
%     For cue Night
%     c1 = squeeze(Mean_Ep_Ch_Placebo(:,freq,:));
%     c2 = squeeze(Mean_Ep_Ch_Sham_PN(:,freq,:));
%     statsPlacebo = test_wilcoxon_cvar(c1, c2,'Placebo','sham','signrank',1,1);
%     pvalPlacebo(freq,:) = statsPlacebo.pval;
% end
% figure
% subplot(2,1,1)
% f_ImageMatrix((pvalCue<=0.05).*(0.1-pvalCue),...
%     v_time2,v_FreqAxis,[], [], [], 0);
% xlabel('Time(sec)')
% title('Cue Night')
% colorbar
% subplot(2,1,2)
% f_ImageMatrix((pvalPlacebo<=0.05).*(0.1-pvalPlacebo),...
%     v_time2,v_FreqAxis,[], [], [], 0);
% xlabel('Time(sec)')
% title('Placebo Night')
% colorbar
% suptitle(strcat('p-values Wilcoxon test All'))