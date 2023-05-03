%% REID_IMU_Purcell

% Purcell et al. 2006
% Report identifying IC using a minima in posterior acceleration (i.e., a max in anterior) (+x) that corresponds in time with a maxima in resultant
% However, the figure they present--and what I was able to successfully implement here--suggest that a minima in x was actually used
% Likely that their figure describing coordinate system was mislabeled

function [timings, stances, segmented] = REID_IMU_Purcell(data, location, Fs, max_step_freq, min_stance_t)

pkdist = round(Fs/max_step_freq)*2;
% Find minima in direction of progression (+x)
[mag_1 x_min_ind_ic] = findpeaks(-data(:,2));
% Calculate resultant acceleration
a_res = vecnorm(data(:,2:4)')';
% Find peaks in resultant
[mag_2 res_max_ind] = findpeaks(a_res,'MinPeakDistance',pkdist);
% Find the x_ic_peak nearest to the res_peak and call it IC
for res_peak_count = 1:size(res_max_ind,1)
    [mag_3 x_ic_ind] = min(abs(x_min_ind_ic - res_max_ind(res_peak_count)));
    IC(res_peak_count,1) = x_min_ind_ic(x_ic_ind);
end % for res_peaks
% Report identifying TC as mean index between local minima in posterior and maxima in medial
% Given error in coordinate definition this converts to maxima in posterior (-x) and maxima in medial (-z)
% Since we are using the ISB +z equals right convention
% We need to multiply z values by -1 if using a left-shank mounted IMU
if strcmp(location,'Left shank') == 1
    data(:,4) = -data(:,4); % convert from z is +right to z is +lateral
end % if left side
% Look from IC-IC to find peaks
for IC_count = 1:size(IC,1)-1
    % The -z peak is not very pronounced so introduce a new constraint
    % Assume TC happens somewhere between 20-60% of IC-IC time
    plus20 = round((IC(IC_count+1,1)-IC(IC_count,1))/5);
    plus60 = 3*plus20;
    % Find maxima in medial (-z) (minima in lateral)
    [mag_4 z_min_ind] = findpeaks(-data(IC(IC_count,1)+plus20:IC(IC_count,1)+plus60,4));
    % Find maxima in posterior (-x) (minima in anterior)
    [mag_5 x_min_ind_tc] = findpeaks(-data(IC(IC_count,1)+plus20:IC(IC_count,1)+plus60,2),"SortStr","descend");
    if isempty(x_min_ind_tc)
        [mag_6 x_min_ind_tc] = max(-data(IC(IC_count,1)+plus20:IC(IC_count,1)+plus60,2));
    end
    % Find differences between z_min_inds and x_min_ind_tcs
    [mag_7 z_closest_to_x] = min(abs(z_min_ind - x_min_ind_tc(1)));
    % Call the point between the x_tc_min and its closest z_tc_min TC
    TC(IC_count,1) = round(mean([x_min_ind_tc(1) z_min_ind(z_closest_to_x)])) + IC(IC_count,1) + plus20;
end % for IC_count
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
[IC,TC] = REID_IMU_crash_catch(min_stance_t*Fs/1000,IC,TC);
% Create tables and structures
timings = table;
stances = table;
segmented = struct;
% Correct timings back to original timestamps
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact = TC + data(1,1) - 1;
% Set left side boolean to true or false based on location (1 = left; 0 = right)
if strcmp(location,'Left shank') == 1
    timings.left_stance(1:size(IC,1)) = 1;
    side = 'left';
elseif strcmp(location,'Right shank') == 1
    timings.left_stance(1:size(IC,1)) = 0;
    side = 'right';
else
    error('Unrecognized wearable placement location')
end
% Fill stance table with timestamps and zeros (no stance)
stances.time = data(:,1);
stances.left_stance = zeros(size(data,1),1);
stances.right_stance = zeros(size(data,1),1);
% Loop through ICs and extract data segmented by gait events
for step_count = 1:size(IC,1) - 1
    segmented.(strcat("stance_",num2str(step_count),'_',side)) = data(IC(step_count):TC(step_count),:);
    segmented.(strcat("swing_",num2str(step_count),'_',side)) = data(TC(step_count)+1:IC(step_count + 1)-1,:);
    if strcmp(side,'left')
        stances.left_stance(IC(step_count,1):TC(step_count,1)) = 1;
    else
        stances.right_stance(IC(step_count,1):TC(step_count,1)) = 1;
    end % if left
end % for step_count

end % function