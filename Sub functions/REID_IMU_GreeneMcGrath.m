%% REID_IMU_GreeneMcGrath

% Greene et al. 2010 developed an adapative algorithm to identify gait events
% during walking from angular rotation about the medial-lateral axis
% This algorithm was later applied to running by McGrath et al. 2012

function [timings, stances, segmented] = REID_IMU_GreeneMcGrath(data, location, Fs, min_stance_t)

% 5th order Butterworth with a 5-Hz cut-off frequency
order = 5; % order
Fc = 5; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
ML_filt = filtfilt(b1,b2,data(:,7));

% Midswing identified as large positive peak
% Peaks must be at least t1 == 0.5 s apart
t1 = round(Fs/2);
[potential_midswing_mag, potential_midswing_ind] = findpeaks(ML_filt,"MinPeakDistance",t1);

% This potential midswing peak needs to be th2 == 0.8*mean of all data points greater than the mean angular velocity
th2 = 0.8*mean(ML_filt(ML_filt>mean(ML_filt)));
% Delete any potential peaks below this theshold
potential_midswing_ind(potential_midswing_mag<th2) = [];
potential_midswing_mag(potential_midswing_mag<th2) = [];

% And needs to have a preceding minima at least th1 == 0.6*max deg/s less than the maximum
th1 = 0.6*potential_midswing_mag;
[minima_mag, minima_ind] = findpeaks(-ML_filt);
minima_mag = -minima_mag;

% Count backwards through potential midswing peaks (to not mess with count as we delete)
for potential_midswing_count = size(potential_midswing_ind):-1:1
    % Find the index of the minima preceding the potential midswing peak
    closest_preceding_min_ind = potential_midswing_ind(potential_midswing_count) - min(potential_midswing_ind(potential_midswing_count)-minima_ind(potential_midswing_ind(potential_midswing_count)-minima_ind>0));
    % Keep the preceding min if it exists and is at least th1 less than the potential midswing peak
    if ~isempty(closest_preceding_min_ind) &&...
            minima_mag(minima_ind==closest_preceding_min_ind) <= potential_midswing_mag(potential_midswing_count) - th1(potential_midswing_count)
        % If keeping, store it's value
        closest_preceding_min_val = minima_mag(minima_ind==closest_preceding_min_ind); % Don't need, reverse the logic and just delete
    else
        % Otherwise
        potential_midswing_ind(potential_midswing_count) = [];
        potential_midswing_mag(potential_midswing_count) = [];
    end
end
% TC and IC are minima pre- and proceeding this maxima
[maxima_mag, maxima_ind] = findpeaks(ML_filt);
% The IC should have a preceding maximum at least th3 greater than the potential IC
% where th3 == 0.8*|mean of all points less than the mean of the angular velocity|
th3 = 0.8*abs(mean(ML_filt(ML_filt<mean(ML_filt))));
% And should be less than th5 == the mean value of the angular velocity
th5 = mean(ML_filt);
% The TC should be less than th4 == 0.8*mean of all points less than the mean of the angular velocity
th4 = 0.8*mean(ML_filt(ML_filt<mean(ML_filt)));
% And should have a preceding maximum at least th6 greater than the potential TC
% where th6 = 2*th3
% NOTE: Using our own data to pilot this code no maxima preceding potentials TCs ever met this condition
% Therefore, we digitized Greene's data using https://apps.automeris.io/wpd/ and executed this approach on their own data
% Again, this approach was unable to identify any TCs due to no preceding maxes exceeding th6 greater than the potential TC
% Thus, we believe this was a typo in their paper and that this condition
% is supposed to be applied to the maxima PRO-ceeding a potential TC
% We use that approach here
th6 = 2*th3;
% Can sometimes fail to find any potential midswing inds
% Skip if empty
if ~isempty(potential_midswing_ind)
    % Preallocate ICs
    IC(1:size(potential_midswing_ind,1),1) = NaN;
    % Count backward through steps (confirmed midswing peaks)
    for step_count = size(potential_midswing_ind,1):-1:1
        % Find indices of minima proceding the midswing peak (potential ICs) within window t2
        % TC and IC minima must be within t2 == 1.5 s of the maxima
        t2 = round(1.5*Fs);
        % However, if this exceeds the size of the data set just look to the end
        % if ~isempty(closest_preceding_min_ind) &&...
        %         potential_midswing_ind(step_count)+t2 > size(ML_filt,1)
        %     t2_temp = size(ML_filt,1) - potential_midswing_ind(step_count);
        % else
        %     t2_temp = t2;
        % end
        if potential_midswing_ind(step_count)+t2 > size(data,1)
            [potential_IC_mag potential_IC_ind] = findpeaks(-ML_filt(potential_midswing_ind(step_count):end));
        else
            [potential_IC_mag potential_IC_ind] = findpeaks(-ML_filt(potential_midswing_ind(step_count):potential_midswing_ind(step_count)+t2));
        end
        if isempty(potential_IC_ind)
            if potential_midswing_ind(step_count)+t2 > size(data,1)
                [potential_IC_mag potential_IC_ind] = max(-ML_filt(potential_midswing_ind(step_count):end));
            else
                [potential_IC_mag potential_IC_ind] = max(-ML_filt(potential_midswing_ind(step_count):potential_midswing_ind(step_count)+t2));
            end
        end
        potential_IC_mag = -potential_IC_mag;
        potential_IC_ind = potential_IC_ind + potential_midswing_ind(step_count) - 1;
        % Delete magnitude >= th5
        potential_IC_ind(potential_IC_mag>=th5) = [];
        potential_IC_mag(potential_IC_mag>=th5) = [];
        % With a preceding maxima at least th3 greater than the potential IC
        if ~isempty(potential_IC_ind)
            for potential_IC_count = size(potential_IC_ind,1):-1:1
                % Find the index of the maxima preceding the potential IC
                closest_max_preceding_IC_ind = potential_IC_ind(potential_IC_count) - min(potential_IC_ind(potential_IC_count)-maxima_ind(potential_IC_ind(potential_IC_count)-maxima_ind>0));
                % Call the potential IC the true IC if it has a preceding max that is at least th3 greater than it's own value
                if ~isempty(closest_max_preceding_IC_ind) &&...
                        maxima_mag(maxima_ind==closest_max_preceding_IC_ind) >= potential_IC_mag(potential_IC_count) + th3
                    % This will overwrite if you find a valid IC closer to the peak
                    IC(step_count,1) = potential_IC_ind(potential_IC_count);
                else
                    % Not a valid IC (either no preceding max, or max is too small)
                    potential_IC_ind(potential_IC_count) = [];
                    potential_IC_mag(potential_IC_count) = [];
                end
            end % for potential IC count
        end % isempty
        % Find indices of the minima preceding the peak (potential TCs) within window t2
        % If this exceeds the size of the data set just look from the start
        if potential_midswing_ind(step_count)-t2 < 1
            t2_temp = potential_midswing_ind(step_count) - 1;
        else
            t2_temp = t2;
        end
        [potential_TC_mag potential_TC_ind] = findpeaks(-ML_filt(potential_midswing_ind(step_count)-t2_temp:potential_midswing_ind(step_count)));
        if isempty(potential_TC_ind)
            [potential_TC_mag potential_TC_ind] = max(-ML_filt(potential_midswing_ind(step_count)-t2_temp:potential_midswing_ind(step_count)));
        end
        potential_TC_mag = -potential_TC_mag;
        potential_TC_ind = potential_TC_ind + potential_midswing_ind(step_count) - t2_temp - 1;
        % Magnitude < th4
        potential_TC_ind(potential_TC_mag>=th4) = [];
        potential_TC_mag(potential_TC_mag>=th4) = [];
        % Greene state that the potential TC must be PRE-ceded by a maximum at least th6 bigger than itself
        % However, as stated above, this approach does not work (even on Greene's own data) and we believe this is a typo
        % Thus, we apply the threshold to the maxima PRO-ceeding the potential TC
        for potential_TC_count = 1:size(potential_TC_ind,1)
            if ~isempty(potential_TC_ind)
                % Find the index of the maxima preceding the potential TC
                closest_proceding_max_ind = potential_TC_ind(potential_TC_count) + min(maxima_ind(maxima_ind-potential_TC_ind(potential_TC_count)>0)-potential_TC_ind(potential_TC_count));
                % Call the potential TC the true TC if it has a proceding max at least th6 greater than it's own value
                if ~isempty(closest_proceding_max_ind) &&...
                        maxima_mag(maxima_ind==closest_proceding_max_ind) >= potential_TC_mag(potential_TC_count) + th6
                    % This will overwrite if you find a valid TC closer to the peak
                    TC(step_count,1) = potential_TC_ind(potential_TC_count);
                else
                    % Not a valid TC so NaN-flag (either no preceding max, or max is too small)
                    potential_TC_ind(potential_TC_count) = NaN;
                    potential_TC_mag(potential_TC_count) = NaN;
                end
            else
                % Not a valid TC so NaN-flag (either no preceding max, or max is too small)
                potential_TC_ind(potential_TC_count) = NaN;
                potential_TC_mag(potential_TC_count) = NaN;
            end
        end
        % Delete NaN-flags
        potential_TC_ind(isnan(potential_TC_ind)) = [];
        potential_TC_mag(isnan(potential_TC_mag)) = [];
    end % step_count
end % empty potential midswing ind
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
% Correct IC timings back to original timestamps
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