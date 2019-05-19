function s_vol = vy_vol_thresh(s_vol,projthresh)


% set threshold
% projthresh = 0.5;
val = s_vol.anatomy;
val(abs(val) < projthresh*max(abs(val(:)))) = NaN;
val = val./max(abs(val(:)));
s_vol.anatomy = val;