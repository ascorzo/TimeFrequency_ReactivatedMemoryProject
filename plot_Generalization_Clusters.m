addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light


p_clustersOfInterest
clusters = fieldnames(Clust);

for cluster = 1:numel(clusters)
    
%     for subj = 1:numel(cf_Generalization.(clusters{cluster}))
%         figure
%         %F = Time_Freq_Odor.freq;
%         mv_plot_2D(cf_Generalization.(clusters{cluster}){subj}, 'x', F, 'y', F)
%         xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
%         l = title(strcat('Subj_',num2str(subj),{' '},...
%             'Frequency generalization',{' '},...
%             clusters{cluster},{' '},'Cluster'));
%         set(l, 'Interpreter', 'none')
%         colorbar('off')
%         colorbar
%         caxis([min(cf_Generalization.(clusters{cluster}){subj}(:)),...
%             max(cf_Generalization.(clusters{cluster}){subj}(:))])
%         
% %         saveas(gcf,strcat(...
% %             'C:\Users\lanan\Desktop\Temp\Figures\OdorMvsVehicleGeneralization\',...
% %             'Subj_',num2str(subj),'_',clusters{cluster},'.png'))
% %         
% %         close all
%     end
    %
    result_average_Generalization.(clusters{cluster}) = ...
        mv_combine_results(resultGeneralization.(clusters{cluster}), 'average');
    
    figure
    %F = Time_Freq_Odor.freq;
    mv_plot_2D(result_average_Generalization.(clusters{cluster}).perf{1}, 'x', F, 'y', F)
    xlabel('Test Frequency (Hz)'), ylabel('Train Frequency (Hz)')
    l = title(strcat('Average_',...
        'Frequency generalization',{' '},...
        clusters{cluster},{' '},'Cluster'));    set(l, 'Interpreter', 'none')
    colorbar('off')
    a = colorbar;
    a.Label.String = 'AUC';
    caxis([min(result_average_Generalization.(clusters{cluster}).perf{1}(:)),...
        max(result_average_Generalization.(clusters{cluster}).perf{1}(:))])

    
    
    saveas(gcf,strcat(...
        'G:\Mi unidad\2021\AnalysisTemp\Figures\Odor_DvsM_Generalization\',...
        clusters{cluster},'Average.png'))
    close all
    
    
    %% Plot Above threshold
    
    s_threshold = 0.6;
    clim = [0 1];
    
    for subj = 1:numel(cf_Generalization.(clusters{cluster}))
        cf_Generalization_threshold.(clusters{cluster}){subj}=...
            (cf_Generalization.(clusters{cluster}){subj}>=s_threshold).*...
            cf_Generalization.(clusters{cluster}){subj};
        
    end
%     
%     for subj = 1:numel(cf_Generalization_threshold.(clusters{cluster}))
%         figure
%         %    F = Time_Freq_Odor.freq;
%         mv_plot_2D(cf_Generalization_threshold.(clusters{cluster}){subj}, 'x', F, 'y', F)
%         xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
%         l = title(strcat('Subj_',num2str(subj),{' '},...
%             'Frequency generalization',{' '},...
%             clusters{cluster},{' '},'Cluster'));        set(l, 'Interpreter', 'none')
%         colorbar('off')
%         colorbar
%         caxis([min(cf_Generalization_threshold.(clusters{cluster}){subj}(:)),...
%             max(cf_Generalization_threshold.(clusters{cluster}){subj}(:))])
%         
% %         saveas(gcf,strcat(...
% %             'C:\Users\lanan\Desktop\Temp\Figures\OdorMvsVehicleGeneralization\AboveThreshold\',...
% %             'Subj',num2str(subj),'_',clusters{cluster},'.png'))
% %         
% %         close all
%     end
    
    result_average_Generalization.(clusters{cluster}) = ...
        mv_combine_results(resultGeneralization.(clusters{cluster}), 'average');
    
    
    averageThreshold.(clusters{cluster}) = ...
        result_average_Generalization.(clusters{cluster}).perf{1}.*...
        (result_average_Generalization.(clusters{cluster}).perf{1}>=0.55);
    
%     figure
%     % F = Time_Freq_Odor.freq;
%     mv_plot_2D(averageThreshold.(clusters{cluster}), 'x', F, 'y', F)
%     xlabel('Test Frequency'), ylabel('Train Frequency')
%     l = title(strcat('Average_',...
%         'Frequency generalization',{' '},...
%         clusters{cluster},{' '},'Cluster'));    set(l, 'Interpreter', 'none')
%     colorbar('off')
%     a = colorbar;
%     a.Label.String = 'AUC';
%     caxis([0 0.55])
    
%     saveas(gcf,strcat(...
%         'C:\Users\lanan\Desktop\Temp\Figures\OdorMvsVehicleGeneralization\AboveThreshold\',...
%         clusters{cluster},'Average.png'))
%     close all
    
end