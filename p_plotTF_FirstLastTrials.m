%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('C:\Users\lanan\Documents\MATLAB\eeglab2019_1'))
addpath('./Scripts_Wilc')

tf = [];

clusters = fieldnames(TF_OdorFirstTrial_All);

clusters = {'All'};

for cluster = 1:numel(clusters)
    
    figure
    y_lims                  = [];
    x_lims_parcial          = [-5 25];
    time_parcial            = Time_Freq_Odor.time;
    total_subplots_row      = 2;
    total_subplots_column   = 1;
    count                   = 1;
    frequencies             = Time_Freq_Odor.freq;
    v_xlim                  = [-5 25];
    
    TF_OdorD_FirstTrial = ...
        squeeze(mean(TF_OdorFirstTrial_All.(clusters{cluster}),1));
    
    TF_OdorD_LastTrial = ...
        squeeze(mean(TF_OdorLastTrial_All.(clusters{cluster}),1));
    
    TF_VehicleD_FirstTrial = ...
        squeeze(mean(TF_VehicleFirstTrial_All.(clusters{cluster}),1));
    
    TF_VehicleD_LastTrial = ...
        squeeze(mean(TF_VehicleLastTrial_All.(clusters{cluster}),1));
    
    
    bottom_TF  = min([min(TF_OdorD_FirstTrial(:)),...
        min(TF_OdorD_LastTrial(:)),min(TF_VehicleD_FirstTrial(:)),...
        min(TF_VehicleD_LastTrial(:))])*0.8;
    
    top_TF     = max([max(TF_OdorD_FirstTrial(:)),...
        max(TF_OdorD_LastTrial(:)),max(TF_VehicleD_FirstTrial(:)),...
        max(TF_VehicleD_LastTrial(:))])*0.8;
    
    start_sample = 1;
    end_sample = length(time_parcial);
    
    
    %----------------------------------------------------------------------
    % Plot TF Odor
    %----------------------------------------------------------------------
    
    % -- first trial----
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorD_LastTrial-TF_VehicleD_LastTrial,...
        time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)')
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor - Vehicle Last stimulation')
    xlabel('')
    colorbar
    
    
    % -- Odor D last trial---- 
    count = count+1;
    tf(count) = subplot(total_subplots_row,total_subplots_column,count);
    f_ImageMatrix(TF_OdorD_FirstTrial-TF_VehicleD_FirstTrial,time_parcial,frequencies,y_lims)
    xlim(v_xlim)
    ylabel('Frequency(Hz)')
    colormap(parula)
    caxis manual
    caxis([bottom_TF top_TF]);
    title('Odor - Vehicle First stimulation')
    xlabel('')
    colorbar
    
   
    l = sgtitle(clusters{cluster});
    set(l, 'Interpreter', 'none')
end