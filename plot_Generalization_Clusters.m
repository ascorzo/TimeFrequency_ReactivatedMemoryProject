addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light

for subj = 1:numel(cf_Generalization)
    
    figure
    %F = Time_Freq_Odor.freq;
    mv_plot_2D(cf_Generalization{subj}, 'x', F, 'y', F)
    xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
    title('Frequency generalization using channels-x-times as features')
    colorbar('off')
    colorbar
    caxis([min(cf_Generalization{subj}(:)),...
        max(cf_Generalization{subj}(:))])
    
    saveas(gcf,strcat(...
        'C:\Users\lanan\Desktop\Temp\Figures\OdorDvsVehicleGeneralization\',...
        'Subj',num2str(subj),'.png'))
    
    close all
end
% 
result_average_Generalization = mv_combine_results(resultGeneralization, 'average');

figure
%F = Time_Freq_Odor.freq;
mv_plot_2D(result_average_Generalization.perf{1}, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Frequency generalization using channels x time as features')
colorbar('off')
colorbar
caxis([min(result_average_Generalization.perf{1}(:)),...
        max(result_average_Generalization.perf{1}(:))])

saveas(gcf,strcat(...
        'C:\Users\lanan\Desktop\Temp\Figures\OdorDvsVehicleGeneralization\',...
        'Avegare.png'))
close all


%% Plot Above threshold

s_threshold = 0.6;
clim = [0 1];

for subj = 1:numel(cf_Generalization)   
    cf_Generalization_threshold{subj}=...
        (cf_Generalization{subj}>=s_threshold).*cf_Generalization{subj};
     
end

for subj = 1:numel(cf_Generalization_threshold)
    figure
%    F = Time_Freq_Odor.freq;
    mv_plot_2D(cf_Generalization_threshold{subj}, 'x', F, 'y', F)
    xlabel('Test Frequency [Hz]'), ylabel('Train Frequency [Hz]')
    title('Frequency generalization using channels-x-times as features')
    colorbar('off')
    colorbar
    caxis([min(cf_Generalization_threshold{subj}(:)),...
        max(cf_Generalization_threshold{subj}(:))])
    
    saveas(gcf,strcat(...
        'C:\Users\lanan\Desktop\Temp\Figures\OdorDvsVehicleGeneralization\AboveThreshold\',...
        'Subj',num2str(subj),'.png'))
    
    close all
end

result_average_Generalization = mv_combine_results(resultGeneralization, 'average');


averageThreshold = result_average_Generalization.perf{1}.*...
    (result_average_Generalization.perf{1}>=0.55);

figure
% F = Time_Freq_Odor.freq;
mv_plot_2D(averageThreshold, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Frequency generalization using channels x time as features')
colorbar('off')
colorbar
caxis([0 0.55])

saveas(gcf,strcat(...
        'C:\Users\lanan\Desktop\Temp\Figures\OdorDvsVehicleGeneralization\AboveThreshold\',...
        'Avegare.png'))
close all


%% statistics

cfg_stat = [];
cfg_stat.metric          = 'auc';
cfg_stat.test            = 'permutation';
cfg_stat.correctm        = 'cluster';  % correction method is cluster
cfg_stat.n_permutations  = 100;

cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.alpha           = 0.05; % use standard significance threshold of 5%

cfg_stat.design          = 'within';

cfg_stat.statistic       = 'wilcoxon';
cfg_stat.null            = 0.5;


cfg_stat.clustercritval  = 1.65;
% z-val = 1.65 corresponds to uncorrected p-value = 0.1
% z-val = 1.96 corresponds to uncorrected p-value = 0.05
% z-val = 2.58 corresponds to uncorrected p-value = 0.01

stat_level2 = mv_statistics(cfg_stat, resultGeneralization);

figure
mv_plot_2D(stat_level2.mask_with_cluster_numbers, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Cluster permutation')
colorbar('off')
colorbar

