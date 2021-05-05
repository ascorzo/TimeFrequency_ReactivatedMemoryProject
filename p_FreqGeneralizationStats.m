%% compare generalization D night vs M night

D_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorDvsVehicle_FreqGeneralization.mat');
M_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorMvsVehicle_FreqGeneralization.mat');

F = D_Generalization.F;

subjects = numel(D_Generalization.cf_Generalization);

for subj = 1:subjects
    
    for freq = 1:numel(F)
        AllSubj_GenDiagonal_DNight(subj,freq) = ...
            D_Generalization.cf_Generalization{subj}(freq,freq);   
        
        AllSubj_GenDiagonal_MNight(subj,freq) = ...
            M_Generalization.cf_Generalization{subj}(freq,freq); 
    end
end
% 
% result_average_DGeneralization = ...
%     mv_combine_results(D_Generalization.resultGeneralization, 'average');
% result_average_MGeneralization = ...
%     mv_combine_results(M_Generalization.resultGeneralization, 'average');
% 
% for freq = 1:numel(F)
%     MeanSubj_Generalization(freq) = ...
%         result_average_Generalization.perf{1}(freq,freq);
% end

%plot(F,AllSubj_Generalization)
%%
f_nonParametric_FreqGen('Frequency Generalizatioon',...
        'Time(sec)',' ','D Night','M Night',...
        AllSubj_GenDiagonal_DNight,...
        AllSubj_GenDiagonal_MNight,...
        F,'-r','-b',1,numel(F));
xticks(F(1:1/0.1:end))
%ylim([0.4 0.65])

%% stats on complete TF


for subj = 1:subjects
    
    AllSubj_Generalization_DNight(:,:,subj) = ...
        D_Generalization.cf_Generalization{subj};
    
    AllSubj_Generalization_MNight(:,:,subj) = ...
        M_Generalization.cf_Generalization{subj};
    
end

data = {AllSubj_Generalization_DNight; AllSubj_Generalization_MNight};
paired = 'on';

method = 'bootstrap';
perms = 100;

[stats, ~, pvals, ~] = statcond( data,'paired', paired,...
    'method',method,'naccu',perms);

%[~, pmask] = fdr(pvals,0.05,'nonParametric');

pmask = pvals<=0.05;

figure
mv_plot_2D((1-pvals).*pmask, 'x', F, 'y', F)
xlabel('Test Frequency'), ylabel('Train Frequency')
title('Cluster permutation')
colorbar('off')
colorbar
caxis([0.98 1])