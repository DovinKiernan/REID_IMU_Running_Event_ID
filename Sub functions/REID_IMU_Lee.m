%% REID_IMU_Lee

% Lee et al., 2010
% Identified IC as large maxima in AP acceleration
% TC, was identified as second, lower magnitude, maxima in AP acceleration
% Left and right stance were identified by large ML accelerations
% coinciding in time with IC (just after IC)

function [timings, stances, segmented] = REID_IMU_Lee(data, Fs, max_step_freq, min_stance_t)

% Identify peaks in AP acceleration
% Lee et al do not fully describe their coordinate system
% Their Fig 1 is labelled "forward acceleration," however, the data they
% present are consistent with ~posterior positive, we have executed this
% code based on that assumption
pkdist = round(Fs/max_step_freq);
[mag_1, ap_min_ind] = findpeaks(-data(:,2),"MinPeakDistance",pkdist);
% Burn the first and last peak
ap_min_ind(1) = [];
ap_min_ind(end) = [];
% Call ind1 IC
IC = ap_min_ind;
% Look between 2 ICs for a smaller AP acceleration maxima
% There can be a negative AP peak immediately after IC that is larger
% than the target peak, to avoid that peak burn min_stance_time - 10 ms
for step_count = 1:size(IC,1)-1
    [mag_2 ap_secondary_min_ind] = findpeaks(-data(IC(step_count,1)+round(Fs/1000*min_stance_t)-round(Fs/100):IC(step_count+1,1),2),"SortStr","descend","NPeaks",1);
    if ~isempty(ap_secondary_min_ind)
        TC(step_count,1) = ap_secondary_min_ind + IC(step_count,1) + round(Fs/1000*min_stance_t) - round(Fs/100) - 1;
    else
        TC(step_count,1) = NaN;
    end
end % for step_count
% Burn the last IC as there is no corresponding TC
IC(end) = [];
% If any TC is NaN-flagged, delete the corresponding IC
IC(isnan(TC)) = [];
TC(isnan(TC)) = [];
% Look just after the IC in time for a peak in the ML (we'll call it 1/3rd min stance time) 
% Direction of this peak is not fully described but assuming it is supposed to be in direction of contact side
% Their sample rate is 100 (1/10th) what was used to pilot this code
% Lower sample frequency appears to be missing a transient that creates a rapid large magnitude peak in opposite of described direction on many steps
% % % Approach 1 -- IDed stance correctly 63% of time in testing
% % % So we're going to grab both positive and negative peaks and then look at timing
% % % The transient that Lee does not describe should come first
% % for step_count = 1:size(IC,1)
% %     [ML_pos_mag ML_pos_ind] = findpeaks(data(IC(step_count,1):IC(step_count,1)+round(Fs/1000*min_stance_t/3),4),"SortStr","descend","NPeaks",1);
% %     [ML_neg_mag ML_neg_ind] = findpeaks(-data(IC(step_count,1):IC(step_count,1)+round(Fs/1000*min_stance_t/3),4),"SortStr","descend","NPeaks",1);
% %     % Find which of these two large peaks occurred first (i.e., is the transient not described by Lee)
% %     % Recalling that our coordinate system is +right
% %     ML_two_largest = sortrows([ML_pos_ind, -ML_pos_mag; ML_neg_ind, ML_neg_mag]);
% %     if ~isempty(ML_two_largest) && ML_two_largest(1,2) > 0
% %         left_stance(step_count,1) = 0;
% %     elseif ~isempty(ML_two_largest) && ML_two_largest(1,2) < 0
% %         left_stance(step_count,1) = 1;
% %     else
% %         % In case it equals 0 or doesn't exist, call it 0.5
% %         left_stance(step_count,1) = 0.5;
% %     end % if pos or neg first
% % end % step_count
% Approach 2 -- IDed stance correctly 83% of time in testing
% We can use filtered ML data to get a waveform more in line with what Lee described
% Low-pass filter ML data
order = 4; % order
Fc = 10; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
ML_filt = filtfilt(b1,b2,data(:,4));
for step_count = 1:size(IC,1)
    [ML_pos_mag ML_pos_ind] = findpeaks(ML_filt(IC(step_count,1):TC(step_count,1)),"SortStr","descend","NPeaks",1);
    [ML_neg_mag ML_neg_ind] = findpeaks(-ML_filt(IC(step_count,1):TC(step_count,1)),"SortStr","descend","NPeaks",1);
    % If you don't find a value call it 0
    if isempty(ML_pos_mag)
        ML_pos_mag = 0;
    end
    if isempty(ML_neg_mag)
        ML_neg_mag = 0;
    end
    % Recall coordinate system is +right so...
    if abs(ML_pos_mag) > abs(ML_neg_mag)
        left_stance(step_count,1) = 0;
    elseif abs(ML_pos_mag) < abs(ML_neg_mag)
        left_stance(step_count,1) = 1;
    else
        % In case it equals 0 or doesn't exist, call it 0.5
        left_stance(step_count,1) = 0.5;
    end % if pos or neg first
end % step_count
% Use mean of every other row to check for misidentified stances
if mean(left_stance(1:2:end)) >= mean(left_stance(2:2:end))
    left_stance(1:2:end) = 1;
    left_stance(2:2:end) = 0;
elseif mean(left_stance(1:2:end)) < mean(left_stance(2:2:end))
    left_stance(1:2:end) = 0;
    left_stance(2:2:end) = 1;
else
    % If this happens it will crash when it tries to assign values to the output variables
    disp('ERROR: cannot reliably determine stance side')
end
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
[IC,TC,left_stance] = REID_IMU_crash_catch(min_stance_t*Fs/1000,IC,TC,left_stance);
% Create tables and structures
timings = table;
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact = TC + data(1,1) - 1;
timings.left_stance = left_stance;
stances = table;
stances.time = data(:,1);
stances.left_stance = zeros(size(data,1),1);
stances.right_stance = zeros(size(data,1),1);
segmented = struct;
for step_count = 1:size(IC,1)-1
    if left_stance(step_count) == 1
        stances.left_stance(IC(step_count,1):TC(step_count,1)) = 1;
        segmented.(strcat("stance_",num2str(step_count),"_left")) = data(IC(step_count):TC(step_count),:);
        segmented.(strcat("swing_",num2str(step_count),"_left")) = data(TC(step_count):IC(step_count + 1)-1,:);
    else
        stances.right_stance(IC(step_count,1):TC(step_count,1)) = 1;
        segmented.(strcat("stance_",num2str(step_count),"_right")) = data(IC(step_count):TC(step_count),:);
        segmented.(strcat("swing_",num2str(step_count),"_right")) = data(TC(step_count):IC(step_count + 1)-1,:);
    end
end % for step_count

end % function