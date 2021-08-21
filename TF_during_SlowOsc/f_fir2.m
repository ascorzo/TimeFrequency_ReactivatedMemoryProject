function [FilteredEEGEnvelope, FilteredEEG, FrBand] = ...
    f_fir2(signal, srate, min_freq, max_freq)
    
global isoctave

% min_freq = round(min_freq,2);
% max_freq = round(max_freq, 2);
revAdded    = 2.5 * srate; % 2.5 for Customfiltered data was enough
% Taking the first and last time points (numel = order) from original
% signal and add it to left/right edge respectively in reverse order
% --> Filtered time points towards the edges will take into account a
% "mirrored" signal and generate this way time points from the original
% signal.
signal      = [signal(revAdded:-1:1) signal signal(end:-1:end-revAdded+1)]';  
ro          = min(0.5,mean([ min_freq 0])); % rollover

if isempty(min_freq)
    res_f   = [0 max_freq max_freq+ro srate/2] / (srate/2);
    res     = [1 1 0 0];
elseif isempty(max_freq)
    res_f   = [0 min_freq-ro min_freq srate/2] / (srate/2);
    res     = [0 0 1 1];
else
    res_f   = [0 min_freq-ro min_freq max_freq max_freq+ro srate/2];
    res_f   =  res_f / (srate/2);
    res     = [0 0 1 1 0 0];
end


%% Here are introduced differences between Octave and Matlab:
% - xF gives same shape of curve but with slightly different amplitudes and
%   slopes of curve. Reading through online forums, this is neglectable and
%   comes from differences in how the programs Octave and Matlab are
%   dealing with the hardware. Whatever that means...

N           = length(signal)- (1-rem(length(signal), 2)); % the order 60000 --> 59999
F           = fir2(N-1,res_f,res); % filter (time)
xF          = abs(fft(F)); % filter (freq)

sign        = fft(double(signal'), N, 2);
sign        = sign .* xF;

% Up to here, no real differences between Matlab and Octave. Slight
% differences in values but VERY small. Matlab will produce a double
% array, while Octave will keep the double array with complex numbers.
sign        = ifft(sign);
% Solution is to force Octave to convert the array to real().
if exist('OCTAVE_VERSION', 'builtin') ~= 0
  sign = real(sign);
end


FilteredEEG = sign';
FilteredEEG = FilteredEEG(revAdded+1:end-revAdded);

FilteredEEGEnvelope = abs(hilbert(FilteredEEG));
FrBand      = [min_freq, max_freq];

end

