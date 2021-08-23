addpath('C:\Users\lanan\Documents\Github\TimeFrequency_ReactivatedMemoryProject\')

p_clustersOfInterest
clusters = fieldnames(Clust);
clusters{length(clusters)+1}  = 'all';

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Clust_ToPlot_DNight.mat')
load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Clust_ToPlot_MNight.mat')

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_MeanAllchans_ToPlot_DNight.mat')
load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_MeanAllchans_ToPlot_MNight.mat')

load('G:\Mi unidad\2021\AnalysisTemp\Using\Time_Freq_Parms')

Time_Freq_OdorD_clust.all = Time_Freq_OdorD_All;
Time_Freq_VehicleD_clust.all = Time_Freq_VehicleD_All;

Time_Freq_OdorM_clust.all = Time_Freq_OdorM_All;
Time_Freq_VehicleM_clust.all = Time_Freq_VehicleM_All;

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       Calculation of power change in specific bands
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time_Freq = Time_Freq_Parms;
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

SpindleBand     = [12 16];
DeltaBand       = [1 4];
Around_2_5      = [2 8];
ThetaBand       = [4 8];
SWBand          = [0.5 2];
%BetaBand       = [18 30];
time            = Time_Freq.time;
subjects        = 1:23; %numel(filesOdor);


SpindleIdx = find(Time_Freq.freq>=SpindleBand(1) &...
    Time_Freq.freq<=SpindleBand(2));
DeltaIdx = find(Time_Freq.freq>=DeltaBand(1) &...
    Time_Freq.freq<=DeltaBand(2));
ThetaIdx = find(Time_Freq.freq>=ThetaBand(1) &...
    Time_Freq.freq<=ThetaBand(2));
Around_2_5Idx = find(Time_Freq.freq>=Around_2_5(1) &...
    Time_Freq.freq<=Around_2_5(2));
SWIdx = find(Time_Freq.freq>=SWBand(1) &...
    Time_Freq.freq<=SWBand(2));

%clusters = {'All'};
%--------------------------------------------------------------------------
% For Declarative Associated Odor Night
%--------------------------------------------------------------------------
for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorD.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_VehicleD.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

end

%--------------------------------------------------------------------------
% For Motor Associated Odor Night
%--------------------------------------------------------------------------

for cluster = 1:numel(clusters)
    % -----------------For Cue----------------------------------
    v_Spindle_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_OdorM.(clusters{cluster}) = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

    % -----------------For Vehicle----------------------------------
    v_Spindle_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SpindleIdx,:),2));
    v_Delta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,DeltaIdx,:),2));
    v_Theta_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,ThetaIdx,:),2));
    v_SW_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,SWIdx,:),2));
    v_Around_2_5_VehicleM.(clusters{cluster}) = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster})(subjects,Around_2_5Idx,:),2));

end


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1'))
addpath('./Scripts_Wilc')

%p_plot_TimeFreq_All
tf = [];
% Time_Freq.time = -1:0.05:5;
% Time_Freq.freq = 0.5:0.05:20;

for cluster = 1:numel(clusters)
    
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = Time_Freq_Parms.time;
    total_subplots_row      = 5;
    total_subplots_column   = 2;
    count                   = 1;
    frequencies             = Time_Freq_Parms.freq;
    v_xlim                  = [-5 25];
    
    TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
    TF_VehicleD = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster}),1));
    
    TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
    TF_VehicleM = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD-TF_VehicleD),...
        min(TF_OdorM- TF_VehicleM)])*0.6;
    top_TF     = max([max(TF_OdorD(:)- TF_VehicleD(:)),...
        max(TF_OdorM-TF_VehicleM)])*0.6;
    
    % -- ylims Spindles
    bottom_Spindles  = min([min(v_Spindle_OdorD.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        min(v_Spindle_OdorM.(clusters{cluster})(:)),...
        min(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Spindles     = max([max(v_Spindle_OdorD.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleD.(clusters{cluster})(:)),...
        max(v_Spindle_OdorM.(clusters{cluster})(:)),...
        max(v_Spindle_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Theta
    bottom_Theta  = min([min(v_Theta_OdorD.(clusters{cluster})(:)),...
        min(v_Theta_VehicleD.(clusters{cluster})(:)),...
        min(v_Theta_OdorM.(clusters{cluster})(:)),...
        min(v_Theta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Theta     = max([max(v_Theta_OdorD.(clusters{cluster})(:)),...
        max(v_Theta_VehicleD.(clusters{cluster})(:)),...
        max(v_Theta_OdorM.(clusters{cluster})(:)),...
        max(v_Theta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims Delta
    bottom_Delta  = min([min(v_Delta_OdorD.(clusters{cluster})(:)),...
        min(v_Delta_VehicleD.(clusters{cluster})(:)),...
        min(v_Delta_OdorM.(clusters{cluster})(:)),...
        min(v_Delta_VehicleM.(clusters{cluster})(:))])*0.2;
    top_Delta     = max([max(v_Delta_OdorD.(clusters{cluster})(:)),...
        max(v_Delta_VehicleD.(clusters{cluster})(:)),...
        max(v_Delta_OdorM.(clusters{cluster})(:)),...
        max(v_Delta_VehicleM.(clusters{cluster})(:))])*0.12;
    
    % -- ylims SW
    bottom_SW  = min([min(v_SW_OdorD.(clusters{cluster})(:)),...
        min(v_SW_VehicleD.(clusters{cluster})(:)),...
        min(v_SW_OdorM.(clusters{cluster})(:)),...
        min(v_SW_VehicleM.(clusters{cluster})(:))])*0.2;
    top_SW     = max([max(v_SW_OdorD.(clusters{cluster})(:)),...
        max(v_SW_VehicleD.(clusters{cluster})(:)),...
        max(v_SW_OdorM.(clusters{cluster})(:)),...
        max(v_SW_VehicleM.(clusters{cluster})(:))])*0.12;
    
    
    start_sample = 1;
    end_sample = length(time_parcial);
    
    
    %-------------
    zmin = bottom_TF;
    zmax = top_TF;
    s_step                      = (zmax-zmin)/1000;
    mapStep                     = 0:s_step:1-s_step;
    
    red                         = [];
    blue                        = [];
    for i = 1:numel(mapStep)
        s_color                 = s_step*i-s_step;
        red(end+1:end+i, :)     = repmat([1, s_color, s_color], i, 1);
        blue(end+1:end+i, :)    = repmat([s_color, s_color, 1], i, 1);
    end
    s_color                     = s_step*i;
    red(end+1:end+i, :)         = repmat([1, s_color, s_color], i, 1);
    blue(end+1:end+i, :)        = repmat([s_color, s_color, 1], i, 1);
    
    s_filling                   = round(1 * size(red, 1));
    filling                     = repmat([1, 1, 1], s_filling, 1);
    
    custommap                   = [blue; filling; red(end:-1:1, :)];
    
   

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
    
    % -- D Night ----
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorD-TF_VehicleD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)')
    colormap(custommap)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor D - Vehicle')
    xlabel('')
    
    % -- M Night ----
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorM-TF_VehicleM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(custommap)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor M - Vehicle')
    xlabel('')
    a = colorbar('Position',[0.928478964586473,0.778617710514844,0.012297734442653,0.112311010815339]);
    a.Label.String = 'Zscore';
    
    
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),'-',num2str(SpindleBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Vehicle',...
        v_Spindle_OdorD.(clusters{cluster}),...
        v_Spindle_VehicleD.(clusters{cluster}),...
        time_parcial,'-r','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Spindle (',...
        num2str(SpindleBand(1)),'-',num2str(SpindleBand(2)),' Hz)'),...
        'Time(s)',' ','Odor M','Vehicle',...
        v_Spindle_OdorM.(clusters{cluster}),...
        v_Spindle_VehicleM.(clusters{cluster}),...
        time_parcial,'-b','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Spindles top_Spindles])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    
    %----------------------------------------------------------------------
    % Plot Theta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),'-',num2str(ThetaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Vehicle',...
        v_Theta_OdorD.(clusters{cluster}),...
        v_Theta_VehicleD.(clusters{cluster}),...
        time_parcial,'-r','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    ylabel('Zscore','position',[-8.333333381017056,-0.118067079381196,-1])
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Theta (',...
        num2str(ThetaBand(1)),'-',num2str(ThetaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor M','Vehicle',...
        v_Theta_OdorM.(clusters{cluster}),...
        v_Theta_VehicleM.(clusters{cluster}),...
        time_parcial,'-b','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Theta top_Theta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot Delta Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),'-',num2str(DeltaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Vehicle',...
        v_Delta_OdorD.(clusters{cluster}),...
        v_Delta_VehicleD.(clusters{cluster}),...
        time_parcial,'-r','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('Delta (',...
        num2str(DeltaBand(1)),'-',num2str(DeltaBand(2)),' Hz)'),...
        'Time(s)',' ','Odor M','Vehicle',...
        v_Delta_OdorM.(clusters{cluster}),...
        v_Delta_VehicleM.(clusters{cluster}),...
        time_parcial,'-b','-k',start_sample,end_sample);
    legend(tf(count),'off')
    xlim(v_xlim)
    ylim([bottom_Delta top_Delta])
    xlabel('')
    xticks([-5 0 5 10 15 20 25])
    
    %----------------------------------------------------------------------
    % Plot SW Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),'-',num2str(SWBand(2)),' Hz)'),...
        'Time(s)',' ','Odor D','Vehicle',...
        v_SW_OdorD.(clusters{cluster}),...
        v_SW_VehicleD.(clusters{cluster}),...
        time_parcial,'-r','-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.363996762828148,0.035676324060544,0.100728156265703,0.040631750283952]);
    
    xticks([-5 0 5 10 15 20 25])
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_nonParametricTest(strcat('SW (',...
        num2str(SWBand(1)),'-',num2str(SWBand(2)),' Hz)'),...
        'Time(s)',' ','Odor M','Vehicle',...
        v_SW_OdorM.(clusters{cluster}),...
        v_SW_VehicleM.(clusters{cluster}),...
        time_parcial,'-b','-k',start_sample,end_sample);
    xlim(v_xlim)
    ylim([bottom_SW top_SW])
    legend('Position',...
        [0.804773461857274,0.035676324224359,0.100728156265703,0.040631750283952]);
    xticks([-5 0 5 10 15 20 25])

    %----------------------------------------------------------------------
    % Adjust Plot
    %----------------------------------------------------------------------
    linkaxes(tf,'x')
    set(gcf,'position',[453.8,41.8,824,740.8])
    
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
    
    savepath = 'G:\Mi unidad\2021\AnalysisTemp\Figures\TF_Publication\1st_Reduction\';
    saveas(gcf,strcat(savepath,clusters{cluster},'.png'))
    close all
end

