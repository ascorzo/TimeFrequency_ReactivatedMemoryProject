cfg = [];
%cfg.trials = 17;
cfg.channel = 'E36';

SelectedData = ft_selectdata(cfg, cue_data);

cfg                 = [];
cfg.method          = 'mtmconvol';
cfg.taper           = 'hanning';
cfg.output          = 'pow';
cfg.pad             = 'nextpow2';
cfg.foi             =  0:0.005:20;
cfg.t_ftimwin       = ones(length(cfg.foi),1);
cfg.toi             = -15:0.005:15;
cfg.keeptrials      = 'yes';

%------- Cue data -----------------------------------------------
Time_Freq_Cue   = ft_freqanalysis(cfg, SelectedData);

cfg                     = [];
cfg.channel             = 'E36'; % C3
cfg.baseline            = [-15 0];
cfg.baselinetype        = 'relative';
cfg.zlim = [min(Time_Freq_Cue.powspctrm(:)),...
    max(Time_Freq_Cue.powspctrm(:))]*0.002;
ft_singleplotTFR(cfg, Time_Freq_Cue);
colormap(jet)

%% From Forcato

% Time-Frequency Transformation
	cfg                             = [];
	cfg.method                      = 'wavelet';                        
	cfg.output                      = 'pow';
	cfg.width                       = 10;                             
	cfg.toi                         = -15:0.005:15;
	cfg.foi                         = 0.5:0.005:20;    
	data_freq						= ft_freqanalysis(cfg, SelectedData); % FOR EVERY TRAIL Warning: output frequencies are different from input frequencies, multiples of the same bin were requested but not given
 
	
	% Plotting	
	cfg						= [];
	cfg.baseline			= [-15 0]; % the tone should start at about -2.8
	cfg.baselinetype		= 'relative';
	
	cfg.zlim				= [0 2.5];
	%cfg.ylim				= [3 20];
	cfg.xlim				= [-15 15];
    
    ft_singleplotTFR(cfg, data_freq);
    colormap(jet)
