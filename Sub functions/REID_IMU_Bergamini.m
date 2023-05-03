%% REID_IMU_Bergamini

% Bergamini et al. 2012
% Found that the IC and TC corresponded to maxima and minima in
% the second derivative of resultant angular velocity

function [timings, stances, segmented] = REID_IMU_Bergamini(data, Fs, max_step_freq, min_stance_t)

% Re-sample to 200 Hz
Fs_Bergamini = 200;
w_resample = scale_time(data(:,5:7),1,size(data,1),round(size(data,1)/Fs*Fs_Bergamini));
w_res = vecnorm(w_resample')';
% They used Luo 2006's wavelet-mediated differentiation method to obtain angular jerk
% This differentiation function is available at:
% Jianwen Luo (2022). Numerical differentiation based on wavelet transforms (https://www.mathworks.com/matlabcentral/fileexchange/13948-numerical-differentiation-based-on-wavelet-transforms)
% In the Luo 2006 paper a quadratic wavelet function was used
% No other info provided by Bergamini et al. so we assume the same was used there
% To validate this assumption we digitized their data and were able to
% reproduce their Fig1C (angular jerk) from their Fig1B (angular velocity)
% using the following combination of parameters...
% Add the quadratic spline function
wavemngr('add','Quadratic Spline','spl',4,'','spline_scale_function',[-1 1]);
% Use the Luo method to double differentiate resultant angular velocity
j_res = derivative_dwt(derivative_dwt(w_res','spl',4,1/Fs_Bergamini),'spl',4,1/Fs_Bergamini)';
% Find max in j_res
pkdist = round(Fs_Bergamini/max_step_freq);
[mag_1 jerk_max_ind] = findpeaks(j_res,"MinPeakDistance",pkdist);
% Burn the first and last peak as they can correspond to erroneous features
IC = jerk_max_ind(2:end-1);
% Look through IC-IC and find the minima between them
for step_count = 1:size(IC,1)-1
    [mag_2 jerk_min_ind] = min(j_res(IC(step_count,1):IC(step_count+1,1)));
    TC(step_count,1) = jerk_min_ind + IC(step_count,1) - 1;
end % for step_count
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
[IC,TC] = REID_IMU_crash_catch(min_stance_t*Fs_Bergamini/1000,IC,TC);
% Create tables and structures
timings = table;
stances = table;
segmented = struct;
% Convert back into original timestamps
% Recalling that data were resampled
IC = IC*Fs/Fs_Bergamini;
TC = TC*Fs/Fs_Bergamini;
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact = TC + data(1,1) - 1;
% Fill stance table with timestamps and zeros (no stance)
timings.left_stance(1:size(IC,1)) = 0.5;
stances.time = data(:,1);
stances.left_stance = zeros(size(data,1),1);
stances.right_stance = zeros(size(data,1),1);
% Loop through ICs and extract data segmented by gait events
for step_count = 1:size(IC,1)-1
    stances.left_stance(IC(step_count,1):TC(step_count,1)) = 0.5; % No method is reported for stance side determination so call both sides 0.5
    stances.right_stance(IC(step_count,1):TC(step_count,1)) = 0.5;
    segmented.(strcat("stance_",num2str(step_count),"_unknown")) = data(IC(step_count):TC(step_count),:);
    segmented.(strcat("swing_",num2str(step_count),"_unknown")) = data(TC(step_count)+1:IC(step_count + 1)-1,:);
end % for step_count

end % function