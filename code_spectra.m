% TFR
tfrlimits						= [0.5 20];
cycles							= 10;

cfg                             = [];
cfg.method                      = 'wavelet';
cfg.output                      = 'pow';
cfg.width                       = cycles;
cfg.toi                         = 'all'; 
cfg.foi                         = tfrlimits(1):0.005:tfrlimits(2);
data_freq						= ft_freqanalysis(cfg, SelectedData);

cfg								= [];
cfg.baseline					= [-15 0]; % the tone should start at about -2.8
cfg.baselinetype				= 'relchange';
cfg.zlim						= [-2 2];
cfg.showlabels					= 'no';
cfg.colormap					= flipud(brewermap(64, 'RdBu'));
ft_singleplotTFR(cfg, data_freq)

% Spectrum
cfg								= [];
cfg.toilim						= [2 15];
data_pos						= ft_redefinetrial(cfg, SelectedData);
cfg.toilim						= [-15 2];
data_pre						= ft_redefinetrial(cfg, SelectedData);

cfg								= [];
cfg.length						= 4; 
cfg.overlap						= .9;
data_pos						= ft_redefinetrial(cfg, data_pos);
data_pre						= ft_redefinetrial(cfg, data_pre);

% Estimate fractal component
freqs							= 2:.05:45;

cfg								= [];
cfg.foi							= freqs;
cfg.method						= 'irasa';
cfg.pad							= 'nextpow2'; 
frac_pos						= ft_freqanalysis(cfg, data_pos);
frac_pre						= ft_freqanalysis(cfg, data_pre);

% Estimate mixed (normal) power spectrum 
cfg.method 						= 'mtmfft';
cfg.taper 						= 'hanning';
mix_pos							= ft_freqanalysis(cfg, data_pos);
mix_pre							= ft_freqanalysis(cfg, data_pre);

% Calculate the oscillatory component by subtracting the fractal from the
% mixed component
cfg								= [];
cfg.parameter					= 'powspctrm';
cfg.operation					= 'subtract';
osci_pos						= ft_math(cfg, mix_pos, frac_pos); 
osci_pre						= ft_math(cfg, mix_pre, frac_pre); 

% Use percent change for even more obvious peaks
cfg.operation					= 'divide';
perc_pos						= ft_math(cfg, osci_pos, frac_pos); % calculate percent change
perc_pre						= ft_math(cfg, osci_pre, frac_pre); % calculate percent change

figure
sgtitle('Pre-cue (black) vs. Post-cue (red) spectra')

subplot(4,1,1)
plot(mix_pre.freq, mix_pre.powspctrm, 'k--'), hold on
plot(mix_pos.freq, mix_pos.powspctrm, 'r'), xlim([freqs(1) freqs(end)]), ylim([0 5])
title('Raw spectrum')

subplot(4,1,2)
plot(frac_pre.freq, frac_pre.powspctrm, 'k--'), hold on
plot(frac_pos.freq, frac_pos.powspctrm, 'r'), xlim([freqs(1) freqs(end)]), ylim([0 5]) 
title('Fractal component only')

subplot(4,1,3)
plot(osci_pre.freq, osci_pre.powspctrm, 'k--'), hold on
plot(osci_pos.freq, osci_pos.powspctrm, 'r'), xlim([freqs(1) freqs(end)])
title('Oscillatory component only')

subplot(4,1,4)
plot(perc_pre.freq, perc_pre.powspctrm, 'k--'), hold on
plot(perc_pos.freq, perc_pos.powspctrm, 'r'), xlim([freqs(1) freqs(end)])
title('Oscillatory relative to fractal component ')

