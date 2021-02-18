y_lims = [-0.1 0.1];
x_lims_total = [-20 50];
x_lims_parcial = [-15 15];
time_parcial = Time_Freq_Cue.time;
z_lims = [min(Time_Freq_Cue_Frontal(:)),max(Time_Freq_Cue_Frontal(:))];
%----plot frontal channels----

% TF plot
figure
total_subplots = 3;
count = 1;
frequencies = Time_Freq.freq;
 
subplot(total_subplots,1,count)
f_ImageMatrix(squeeze(mean(Time_Freq_Cue_Central,1)),time_parcial,frequencies,y_lims)
xlim(x_lims_parcial)
colormap(pink)
%colorbar
title('TF Odor Central channels')



count = count+1;
subplot(total_subplots,1,count)
f_ImageMatrix(squeeze(mean(Time_Freq_Vehicle_Central,1)),time_parcial,frequencies,y_lims)
xlim(x_lims_parcial)
colormap(pink)
%colorbar
title('TF Vehicle Central channels')

v_time = Time_Freq_Cue.time ;
ValidTimeCue = v_time(~isnan(v_Delta_Cue_Frontal(1,:)));
ValidTimeVehicle = v_time(~isnan(v_Delta_Vehicle_Frontal(1,:)));

ValidTime = intersect(ValidTimeCue,ValidTimeVehicle);
start_sample = ceil(ValidTime(1)-v_time(1)/s_tstep);
end_sample = floor(ValidTime(end)-v_time(1)/s_tstep);

count = count+1;
subplot(total_subplots,1,count)
f_WilcTest('Delta Central',...
     'Time(sec)',' ','Odor D','Vehicle',...
     v_Delta_Cue_Central,...
     v_Delta_Vehicle_Central,...
     v_time,'-r',start_sample,end_sample)
 

set(gcf,'position',[86,37,613,948])
