p_clustersOfInterest
clusters = fieldnames(Clust);

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
    time_parcial            = Time_Freq_Odor.time;
    total_subplots_row      = 6;
    total_subplots_column   = 2;
    count                   = 1;
    frequencies             = Time_Freq_Odor.freq;
    v_xlim                  = [-5 25];
    
    TF_OdorD = squeeze(mean(Time_Freq_OdorD_clust.(clusters{cluster}),1));
    TF_VehicleD = squeeze(mean(Time_Freq_VehicleD_clust.(clusters{cluster}),1));
    
    TF_OdorM = squeeze(mean(Time_Freq_OdorM_clust.(clusters{cluster}),1));
    TF_VehicleM = squeeze(mean(Time_Freq_VehicleM_clust.(clusters{cluster}),1));
    
    bottom_TF  = min([min(TF_OdorD(:)),min(TF_VehicleD(:)),...
        min(TF_OdorM(:)),min(TF_VehicleM(:))])*0.6;
    top_TF     = max([max(TF_OdorD(:)),max(TF_VehicleD(:)),...
        max(TF_OdorM(:)),max(TF_VehicleM(:))])*0.6;
    
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
    

    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
    
    % -- D Night ----
    tf(count) = subplot(total_subplots_row,total_subplots_column,1);
    f_ImageMatrix(TF_OdorD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)','position',[-7.724637712257495,-7.696978993229502,1])
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor D')
    xlabel('')
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,3);
    f_ImageMatrix(TF_VehicleD,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Vehicle')
    xlabel('')
    
    % -- M Night ----  
    
    tf(count) = subplot(total_subplots_row,total_subplots_column,2);
    f_ImageMatrix(TF_OdorM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor M')
    xlabel('')
    
    
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,4);
    f_ImageMatrix(TF_VehicleM,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Vehicle')
    xlabel('')
    a = colorbar('Position',...
    [0.923139158749658,0.659827213822894,0.011812297561022,0.23110151187905]);
    a.Label.String = 'Zscore';
      
    %----------------------------------------------------------------------
    % Plot Spindle Power
    %----------------------------------------------------------------------
    % -- D Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,5);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,6);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,7);
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
    ylabel('Zscore','position',[-7.724637728843142,-0.18949351753476,-1])
    
    % -- M Night ----
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,8);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,9);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,10);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,11);
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
    tf(count) = subplot(total_subplots_row,total_subplots_column,12);
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
    
    savepath = 'G:\Mi unidad\2021\AnalysisTemp\Figures\TF_Publication\';
    saveas(gcf,strcat(savepath,clusters{cluster},'.png'))
    close all
end

