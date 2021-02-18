addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/')

ft_defaults
addpath('/gpfs01/born/group/Andrea/fieldtrip-20200828/qsub')

ft_warning off
%--------------------------------------------------------------------------
% Time frequency parameters
%--------------------------------------------------------------------------
s_tstep = 0.05; % try with 0.05
s_fstep = 0.2; % 0.05

cfg_Tf                      = [];
cfg_Tf.method               = 'wavelet'; %'mtmconvol';
cfg_Tf.output               = 'pow';
cfg_Tf.foi                  =  0.5:s_fstep:20; % Reduce to 20
cfg.width                   = 12;
cfg_Tf.toi                  = -25:s_tstep:55;
cfg_Tf.keeptrials           = 'yes';


%--------------------------------------------------------------------------
% Parameters for Baseline Correction 
%--------------------------------------------------------------------------
cfg_Bas                     = [];
cfg_Bas.baseline            = [-15 45];
cfg_Bas.baselinetype        = 'zscore';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For Declarative Associated Odor Night
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%filepath = 'D:\GermanData\DATA\RawData\preProcessing\Epoched_90SecTrial_MastoidRef_Interp\OdorD_Night\';

filepath = '/gpfs01/born/group/Andrea/ReactivatedConnectivity/Time-Frequency_FT/TF_OdorD_Night/';

%files = dir(strcat(filepath,'*.set'));
files = dir(strcat(filepath,'*.mat'));


for subj = 1:numel(files)
    
    %______________________________________________________________________
    %
    %Cue Odor
    %______________________________________________________________________
    
    disp(strcat('Sujeto: ',num2str(subj)))
    
%     %--------- Load Data --------------------------------------------------
%     addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
%     EEGOdor = pop_loadset(strcat(filepath,files(subj).name));
%     dataOdor = eeglab2fieldtrip(EEGOdor,'raw');
%     rmpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1\'))
    
    
    %-----------Time-Frequency Calculation---------------------------------
%     Time_Freq_DA_Temp  = ft_freqanalysis(cfg_Tf, dataOdor);
%     
%     %--------------Select time of interest---------------------------------
%     cfg = [];
%     cfg.latency = [-15 45];
%     Time_Freq_DA{subj}  = ft_selectdata(cfg, Time_Freq_DA_Temp);


    load(strcat(filepath,files(subj).name));
    %-----------baseline Correction ---------------------------------------
    
    Time_Freq_DA_Baseline = ft_freqbaseline(cfg_Bas,Time_Freq );
    
    %-----------Mean of all trials, each channel --------------------------
    cfg = [];
    cfg.avgoverrpt = 'yes';
    Time_Freq_Cue_baseline2 = rmfield(Time_Freq_DA_Baseline,  'trialinfo');
    Time_Freq_DA_Mean = ft_selectdata(cfg, Time_Freq_Cue_baseline2);
    
    %-----------Separate conditions -----------------------  
    cfg = [];
    cfg.latency = [-15 15];
    Time_Freq_Cue = ft_selectdata(cfg, Time_Freq_DA_Mean);
    
    cfg = [];
    cfg.latency = [15 45];
    Time_Freq_Vehicle = ft_selectdata(cfg, Time_Freq_DA_Mean);
    Time_Freq_Vehicle.time = Time_Freq_DA_Mean.time;
    %----------------------------------------------------------------------
    %           Combine Mean of channels by clusters
    %----------------------------------------------------------------------
    p_clustersOfInterest
    
    clusters = fieldnames(Clust);

    for cluster = 1:numel(clusters)
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Cue.label);
        Time_Freq_Cue_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Cue.powspctrm(ind2,:,:),1));
        
        [~,~,ind2] = intersect(Clust.(clusters{cluster}),Time_Freq_Vehicle.label);
        Time_Freq_Vehicle_clust.(clusters{cluster})(subj,:,:) = ...
            squeeze(mean(Time_Freq_Vehicle.powspctrm(ind2,:,:),1));
        
    end

end



%% Calculation of power change in specific bands

SpindleBand     = [10 18];
DeltaBand       = [0.5 4];
ThetaBand       = [4 8];
BetaBand        = [18 30];
time = Time_Freq.time ;
subjects = 1:numel(files);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
BetaIdx = find(Time_Freq.freq>=BetaBand(1) &...
    Time_Freq.freq<=BetaBand(2));

for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Cue.(clusters{cluster}) = squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_Beta_Vehicle.(clusters{cluster}) = squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster})(subjects,BetaIdx,:),2));

end
%% plot CUE Night

addpath(genpath('/gpfs01/born/group/Andrea/eeglab2019_1/'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All

for cluster = 1:numel(clusters)
    
    y_lims = [];%[-0.1 0.1];
    x_lims_parcial = [-15 15];
    time_parcial = Time_Freq.time;
    %z_lims = [min(Time_Freq_Cue_Frontal(:)),max(Time_Freq_Cue_Frontal(:))];
    
    figure
    total_subplots = 6;
    count = 1;
    frequencies = Time_Freq.freq;
    
    subplot(total_subplots,1,count)
    f_ImageMatrix(squeeze(mean(Time_Freq_Cue_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(pink)
    %colorbar
    title('TF Odor')
 
    count = count+1;
    subplot(total_subplots,1,count)
    f_ImageMatrix(squeeze(mean(Time_Freq_Vehicle_clust.(clusters{cluster}),1)),time_parcial,frequencies,y_lims)
    xlim(x_lims_parcial)
    colormap(pink)
    %colorbar
    title('TF Vehicle')
    
    v_time = Time_Freq_Cue.time ;
    ValidTimeCue = v_time(~isnan(Time_Freq_Cue_clust.(clusters{cluster})(1,1,:)));
    ValidTimeVehicle = v_time(~isnan(Time_Freq_Vehicle_clust.(clusters{cluster})(1,1,:)));
    
    ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
    start_sample = ceil(ValidTime(1)-v_time(1)/s_tstep);
    end_sample = floor(ValidTime(end)-v_time(1)/s_tstep);
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_WilcTest('Beta',...
        ' ',' ','Odor D','Vehicle',...
        v_Beta_Cue.(clusters{cluster}),...
        v_Beta_Vehicle.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample)
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_WilcTest('Spindle',...
        ' ',' ','Odor D','Vehicle',...
        v_Spindle_Cue.(clusters{cluster}),...
        v_Spindle_Vehicle.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample)
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_WilcTest('Theta',...
        ' ',' ','Odor D','Vehicle',...
        v_Theta_Cue.(clusters{cluster}),...
        v_Theta_Vehicle.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample)
    
    count = count+1;
    subplot(total_subplots,1,count)
    f_WilcTest('Delta',...
        'Time(sec)',' ','Odor D','Vehicle',...
        v_Delta_Cue.(clusters{cluster}),...
        v_Delta_Vehicle.(clusters{cluster}),...
        v_time,'-r',start_sample,end_sample)
    
    
    set(gcf,'position',[81,35,1252,926])
    
    l = suptitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    %saveas(gcf,strcat(clusters{cluster},'.png'))

end

 