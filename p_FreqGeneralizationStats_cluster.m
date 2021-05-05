%% compare generalization D night vs M night
% 

addpath('C:\Users\lanan\Documents\MATLAB\MVPA-Light\startup\')
startup_MVPA_Light

%For Frequency generalization plot
% D_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorDvsVehicle_FreqGeneralization_clusters.mat');
% M_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorMvsVehicle_FreqGeneralization_clusters.mat');


%For Time generalization plot
D_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorDvsVehicle_TimeGeneralization_to20_clusters.mat');
M_Generalization = load('C:\Users\lanan\Desktop\Temp\OdorMvsVehicle_TimeGeneralization_to20_clusters.mat');

%For Frequency generalization plot
% F = D_Generalization.F;

%For Time generalization plot
t = D_Generalization.t;


subjects = numel(D_Generalization.cf_Generalization.central);

p_clustersOfInterest
clusters = fieldnames(Clust);

for cluster = 1:numel(clusters)
    for subj = 1:subjects
        
        %for Dim = 1:numel(F)
        for Dim = 1:numel(t)
            AllSubj_GenDiagonal_DNight.(clusters{cluster})(subj,Dim) = ...
                D_Generalization.cf_Generalization.(clusters{cluster}){subj}(Dim,Dim);
            
            AllSubj_GenDiagonal_MNight.(clusters{cluster})(subj,Dim) = ...
                M_Generalization.cf_Generalization.(clusters{cluster}){subj}(Dim,Dim);
        end
    end

    %%
    eeglab nogui
    
%     %For Frequency generalization plot
%     figure
%     f_nonParametric_FreqGen(strcat('Frequency Generalization',{' '},(clusters{cluster})),...
%         'Frequency(Hz)',' ','D Night','M Night',...
%         AllSubj_GenDiagonal_DNight.(clusters{cluster}),...
%         AllSubj_GenDiagonal_MNight.(clusters{cluster}),...
%         F,'-r','-b',1,numel(F));
%     xticks(F(1:1/0.1:end))
%     set(gcf,'position',[3.4,342,1528,420])
%     %ylim([0.4 0.65])
%     l = title(strcat('Frequency Generalization',{' '},(clusters{cluster})));
%     set(l, 'Interpreter', 'none')
%     
%     
%     saveas(gcf,strcat(...
%         'C:\Users\lanan\Desktop\Temp\Figures\Freq_Generalization_DvsM_Night\',...
%         clusters{cluster},'.png'))
    
    
    %For Time generalization plot
    figure
    f_nonParametric_FreqGen(strcat('Time Generalization',{' '},(clusters{cluster})),...
        'Time(s)',' ','D Night','M Night',...
        AllSubj_GenDiagonal_DNight.(clusters{cluster}),...
        AllSubj_GenDiagonal_MNight.(clusters{cluster}),...
        t,'-r','-b',1,numel(t));
    set(gcf,'position',[3.4,342,1528,420])
    %ylim([0.4 0.65])
    l = title(strcat('Time Generalization',{' '},(clusters{cluster})));
    set(l, 'Interpreter', 'none')
    
    
%     saveas(gcf,strcat(...
%         'C:\Users\lanan\Desktop\Temp\Figures\Time_Generalization_DvsM_Night\',...
%         clusters{cluster},'.png'))
    
    
    %% stats on complete TF
    
    
    for subj = 1:subjects   
        AllSubj_Generalization_DNight.(clusters{cluster})(:,:,subj) = ...
            D_Generalization.cf_Generalization.(clusters{cluster}){subj};  
        AllSubj_Generalization_MNight.(clusters{cluster})(:,:,subj) = ...
            M_Generalization.cf_Generalization.(clusters{cluster}){subj};
    end
    
    data = {AllSubj_Generalization_DNight.(clusters{cluster}); ...
        AllSubj_Generalization_MNight.(clusters{cluster})};
    paired = 'on';
    
    method = 'bootstrap';
    perms = 100;
    
    [stats, ~, pvals, ~] = statcond( data,'paired', paired,...
        'method',method,'naccu',perms);
    
    [~, pmask] = fdr(pvals,0.05,'nonParametric');
    
    %pmask = pvals<=0.05;
    
    
%     %For Frequency generalization plot
%     figure
%     mv_plot_2D((1-pvals).*pmask, 'x', F, 'y', F)
%     xlabel('Test Frequency'), ylabel('Train Frequency')
%     l = title(strcat('Cluster permutation',{' '},(clusters{cluster})));
%     set(l, 'Interpreter', 'none')
%     colorbar('off')
%     colorbar
%     caxis([0.98 1])
%     
%     saveas(gcf,strcat(...
%         'C:\Users\lanan\Desktop\Temp\Figures\Freq_Generalization_DvsM_Night\',...
%         clusters{cluster},'_Matrix','.png'))
    
    %For Frequency generalization plot
    figure
    mv_plot_2D((1-pvals).*pmask, 'x', t, 'y', t)
    xlabel('Test Time(s)'), ylabel('Train Time(s)')
    l = title(strcat('Cluster permutation',{' '},(clusters{cluster})));
    set(l, 'Interpreter', 'none')
    colorbar('off')
    colorbar
    caxis([0.98 1])
    
%     saveas(gcf,strcat(...
%         'C:\Users\lanan\Desktop\Temp\Figures\Time_Generalization_DvsM_Night\',...
%         clusters{cluster},'_Matrix','.png'))
%     
%     close all
    
end