%Load time-frequency transforms
%m_Gabor.CueNight = load('m_Gabor_CueNight_AllChans.mat');
%m_Gabor.PlaceboNight = load('m_Gabor_PlaceboNight_AllChans.mat');

subjectsCue = fieldnames(m_Gabor.CueNight.m_GaborCue_total);
subjectsPlac = fieldnames(m_Gabor.PlaceboNight.m_GaborPlacebo_total);
baselineTime = 15;
s_TimeStep = 0.5;
v_time = -15:0.5:15;
v_FreqAxis = m_Gabor.CueNight.v_FreqAxis;

%----------------------------------------------------------------
% create complete arrays for each condition
% size: Subj x Freq x TimePnts x Channels x Trials
%----------------------------------------------------------------
s_MinEpochsOdor = NaN;
s_MinEpochsSham = NaN;

%--- Get min number of epochs to create equal arrays for Cue Night -------
for subj = 1:numel(subjectsCue)
    s_MinEpochsOdor = min(s_MinEpochsOdor,...
        size(m_Gabor.CueNight.m_GaborCue_total.(string(subjectsCue(subj))),4));
    s_MinEpochsSham = min(s_MinEpochsSham,...
        size(m_Gabor.CueNight.m_GaborSham_total.(string(subjectsCue(subj))),4));
end

%----- Fill the arrays for Cue Night--------
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
%--- Get min number of epochs to create equal arrays for Placebo Night ----
for subj = 1:numel(subjectsCue)
    s_MinEpochsOdor = min(s_MinEpochsOdor,...
        size(m_Gabor.PlaceboNight.m_GaborPlacebo_total.(string(subjectsPlac(subj))),4));
    s_MinEpochsSham = min(s_MinEpochsSham,...
        size(m_Gabor.PlaceboNight.m_GaborSham_total.(string(subjectsPlac(subj))),4));
end

%----- Fill the arrays for Placebo Night--------
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
%---------------------------------------------------------
% Baseline Correction
%---------------------------------------------------------

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
% Mean by channels for each subject
%--------------------------------------------------------------------------
%Cue Night Cue Odor
MeanChannels_Cue = squeeze(mean(m_GaborCue_BasCorr,4));

%Cue Night Sham Odor
MeanChannels_Sham_CN = squeeze(mean(m_GaborSham_CN_BasCorr,4));

%Placebo Night Placebo Odor
MeanChannels_Placebo = squeeze(mean(m_GaborPlacebo_BasCorr,4));

%Placebo Night Sham Odor
MeanChannels_Sham_PN= squeeze(mean(m_GaborSham_PN_BasCorr,4));
%--------------------------------------------------------------------------
% Mean of epochs and channels for each subject
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Ch_Cue = squeeze(mean(MeanChannels_Cue,4));

%Cue Night Sham Odor
Mean_Ep_Ch_Sham_CN = squeeze(mean(MeanChannels_Sham_CN,4));

%Placebo Night Placebo Odor
Mean_Ep_Ch_Placebo = squeeze(mean(MeanChannels_Placebo,4));

%Placebo Night Sham Odor
Mean_Ep_Ch_Sham_PN = squeeze(mean(MeanChannels_Sham_PN,4));

%--------------------------------------------------------------------------
% Mean of epochs only
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Cue = squeeze(mean(m_GaborCue_BasCorr,5));

%Cue Night Sham Odor
Mean_Ep_Sham_CN = squeeze(mean(m_GaborSham_CN_BasCorr,5));

%Placebo Night Placebo Odor
Mean_Ep_Placebo = squeeze(mean(m_GaborPlacebo_BasCorr,5));

%Placebo Night Sham Odor
Mean_Ep_Sham_PN = squeeze(mean(m_GaborSham_PN_BasCorr,5));
%--------------------------------------------------------------------------
% Mean of epochs and subjects for each channel
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_Ep_Subj_Cue = squeeze(mean(Mean_Ep_Cue,1));

%Cue Night Sham Odor
Mean_Ep_Subj_Sham_CN = squeeze(mean(Mean_Ep_Sham_CN,1));

%Placebo Night Placebo Odor
Mean_Ep_Subj_Placebo = squeeze(mean(Mean_Ep_Placebo,1));

%Placebo Night Sham Odor
Mean_Ep_Subj_Sham_PN = squeeze(mean(Mean_Ep_Sham_PN,1));

%--------------------------------------------------------------------------
% Mean All (All Subjects,epochs and Channels)
%--------------------------------------------------------------------------
%Cue Night Cue Odor
Mean_All_Cue = squeeze(mean(Mean_Ep_Ch_Cue,1));

%Cue Night Sham Odor
Mean_All_Sham_CN = squeeze(mean(Mean_Ep_Ch_Sham_CN,1));

%Placebo Night Placebo Odor
Mean_All_Placebo = squeeze(mean(Mean_Ep_Ch_Placebo,1));

%Placebo Night Sham Odor
Mean_All_Sham_PN = squeeze(mean(Mean_Ep_Ch_Sham_PN,1));

%% Statistical analysis for each subject across epochs (mean channels)
for subj = 1:numel(subjects) 
    %----- For Cue Night----
    [m_Statis_pval_Ep_Ch_Cue, m_SigP_Ep_Ch_Cue, m_SigPvalue_Ep_Ch_Cue] = ...
        GetStatistics(Mean_Ep_Ch_Cue, Mean_Ep_Ch_Sham_CN);
    
    %----- For Placebo Night----
    [m_Statis_pval_Ep_Ch_Placebo, m_SigP_Ep_Ch_Placebo, m_SigPvalue_Ep_Ch_Placebo] = ...
        GetStatistics(Mean_Ep_Ch_Placebo, Mean_Ep_Ch_Sham_PN);
end

%% Statistical analysis across subjects (mean channels, mean epochs)
%----- For Cue Night----
[m_Statis_pval_All_Cue, m_SigP_All_Cue, m_SigPvalue_All_Cue] = ...
    GetStatistics(Mean_All_Cue, Mean_All_Sham_CN);

%----- For Placebo Night----
[m_Statis_pval_all_Placebo, m_SigP_All_Placebo, m_SigPvalue_All_Placebo] = ...
    GetStatistics(Mean_All_Placebo, Mean_All_Sham_PN);

%% Statistical analysis for each channel across subjects (mean epochs)

%%
function Baseline_Corr = baselinecorrection(Data,baselineTime,s_TimeStep)
baseline = mean(Data(:,:,1:(baselineTime/s_TimeStep),:,:),3);
Baseline_Corr = Data-baseline;
end

%%
function [m_Statistics, m_sigp, m_sigpvalue] = ...
    GetStatistics(Condition1, Condition2)

for i = 1:size(Condition1,1)
    for j = 1:size(Condition1,2)
        cond1 = Condition1(i,j,:);
        cond2 = Condition2(i,j,:);
        minlength = min(size(cond1,3),size(cond2,3));
        x = [cond1(:,:,1:minlength),cond2(:,:,1:minlength)];
        [~,p,~] = ttest_bonf(x,[1 2],0.05,0);
        m_Statistics(i,j) = p;
    end
end
%-- Get significant p-values (<=0.05)
m_sigp = (m_Statistic<=0.05);

%-- Multiply the boolean of significant p-values by the p-value
m_sigpvalue = m_sigp.*m_Statistics;
end
