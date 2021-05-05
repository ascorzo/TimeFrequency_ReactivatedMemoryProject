function FiltereData = f_filter_deltaband(v_data, s_order, fsample, ...
    low_pass, high_pass)

[d, e]          = butter(s_order, 2 * low_pass / fsample, 'low');
FiltereData     = filtfilt(d, e, v_data); %Filter Signal
[d, e]          = butter(s_order, 2 * high_pass / fsample, 'high');
FiltereData     = filtfilt(d, e, FiltereData);

end