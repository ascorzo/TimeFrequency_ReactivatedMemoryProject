%% Group statistics: within-subject cluster permutation test for AUC
% For each subject and every time point, we have calculated AUC values. We
% will now perform a cluster permutation test. See Maris & Oostenveld's
% paper or the FieldTrip tutorials (https://www.fieldtriptoolbox.org/tutorial/cluster_permutation_timelock/) 
% for explanations on cluster permutation tests
cfg_stat = [];
cfg_stat.metric          = 'auc';
cfg_stat.test            = 'permutation';
cfg_stat.correctm        = 'cluster';  % correction method is cluster
cfg_stat.n_permutations  = 1000;

cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.alpha           = 0.05; % use standard significance threshold of 5%

cfg_stat.design          = 'within';

cfg_stat.statistic       = 'wilcoxon';
cfg_stat.null            = 0.5;


cfg_stat.clustercritval  = 1.96;
% z-val = 1.65 corresponds to uncorrected p-value = 0.1
% z-val = 1.96 corresponds to uncorrected p-value = 0.05
% z-val = 2.58 corresponds to uncorrected p-value = 0.01

stat_level2 = mv_statistics(cfg_stat, result_freq);


%result_merge = mv_combine_results(result_freq, 'merge');
%result_merge = mv_select_result(result_merge, 'auc'); % select AUC only
%mv_plot_result(result_merge)

% Instead of plotting all subjects together, we can also calculate the
% grand average across subjects and plot this. To this end, we only need to
% replace 'merge' by 'average' when calling mv_combine_results. Note that
% the shaded area is now the standard deviation across subjects.
result_average = mv_combine_results(result_freq, 'average');
%result_average = mv_select_result(result_average, 'auc');
%mv_plot_result(result_average)

% plot the grand average result again and indicate the cluster in bold
%mv_plot_result(result_average,Time_Freq_Odor{subj}.time, 'mask', stat_level2.mask)

