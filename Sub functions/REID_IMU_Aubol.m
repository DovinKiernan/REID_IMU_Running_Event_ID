%% REID_IMU_Aubol

% Aubol (2020)
% Uses resultant acceleration and jerk from an accelerometer on the
% anteromedial distal tibia to identify IC
% No TC identification

function [timings, stances, segmented] = REID_IMU_Aubol(data, location, Fs, max_step_freq)

% Filter at 70 Hz
order = 4; % order
Fc = 70; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
data_filt = filtfilt(b1,b2,data);
% Calculate resultant acceleration
a_res = vecnorm(data_filt(:,2:4)')';
% Calculate resultant jerk
jerk_res = diff(a_res);
% Peaks in resultant acceleration
% Added a pkdist constraint based on expected timings during running
pkdist = round(Fs/max_step_freq)*2;
[a_res_max_mag, a_res_max_ind] = findpeaks(a_res,'MinPeakDistance',pkdist,'SortStr','descend');
% Peaks in resultant jerk within 150 ms of a peak in resultant acceleration
[mag_1, jerk_res_max_ind] = findpeaks(jerk_res);
% Count backwards throught the jerk maxes and delete ones that aren't within 150 ms
% Then, find minima in resultant acceleration within 75 ms before resultant jerk peak
% with a minimum prominence of 0.2*the third highest resultant acceleration peak
% If there are 3 or more a_res peaks, take the third biggest
if size(a_res_max_mag,1) >= 3
    min_prom = 0.2*a_res_max_mag(3);
else
    % Otherwise, just take the smallest
    min_prom = 0.2*min(a_res_max_mag);
end
% Pre-allocate for mins
a_res_min_ind = [];
for jerk_count = size(jerk_res_max_ind,1):-1:1
    % Take the absolute difference between the current jerk max and ALL
    % a_res maxes. If the jerk max isn't within 150 ms of ANY of the a_res
    % maxes then delete it
    % Or, if it's less than 3 frames from the data start (can't find peaks)
    if (~any(abs(jerk_res_max_ind(jerk_count)-a_res_max_ind) < round(Fs*150/1000))) ||...
            (jerk_res_max_ind(jerk_count) < 3)
        jerk_res_max_ind(jerk_count) = [];
    else
        % For surviving jerk_maxes, look 75 ms back and find a_res minima
        % meeting the prominence threshold
        % Check the 75 ms window is within our data
        if jerk_res_max_ind(jerk_count) - round(Fs*75/1000) > 0
            [mag_2, a_res_min_ind_temp] = findpeaks(-a_res(jerk_res_max_ind(jerk_count) - round(Fs*75/1000):jerk_res_max_ind(jerk_count)),'MinPeakProminence',min_prom);
            % Correct time stamp
            a_res_min_ind_temp = a_res_min_ind_temp + jerk_res_max_ind(jerk_count) - round(Fs*75/1000) - 1;
            % Add to min matrix
            a_res_min_ind = [a_res_min_ind; a_res_min_ind_temp];
        else
            % Otherwise just look from the start of the data set
            [mag_2, a_res_min_ind_temp] = findpeaks(-a_res(1:jerk_res_max_ind(jerk_count)),'MinPeakProminence',min_prom);
            % Add to min matrix
            a_res_min_ind = [a_res_min_ind; a_res_min_ind_temp];
        end % data within search window
    end % within 150 ms
end % jerk_count
% Eliminate redundant mins
a_res_min_ind = unique(a_res_min_ind);
% Eliminate minima that have a 40 Hz or faster a_res peak preceding them
% First approach here was to find peaks peaks with a max width of 1000/40 = 25 ms
% That preceded mins and then delete the min
% Overly aggressive and deleted too many mins
% Thus, moved to a second solution: Simply looking for ANY peak occurring
% in a 25 ms window preceding a min
% Much more liberal, seems to work
for min_count = size(a_res_min_ind,1):-1:1
    % Check the window is within our data
    if a_res_min_ind(min_count) - round(Fs*25/1000) > 0
        [mag_3 fast_peak_ind] = findpeaks(a_res(a_res_min_ind(min_count)-round(Fs*25/1000):a_res_min_ind(min_count)));
    else
        [mag_3 fast_peak_ind] = findpeaks(a_res(1:a_res_min_ind(min_count)));
    end % search window is within data
    % If we find a peak meeting these criteria then delete the min
    if exist('fast_peak_ind','var') && ~isempty(fast_peak_ind)
        a_res_min_ind(min_count) = [];
        clear fast_peak_ind
    end % found a fast peak
end % mincount
% Accept the earliest occuring minima that meets these criteria as the IC
% Reorder by index
a_res_max_ind = sort(a_res_max_ind,'ascend');
% Look from a_res max to next a_res max
% Find all a_res_mins between them
% Take the earliest occurring a_res_min as IC
for step_count = 1:size(a_res_max_ind,1)-1
    potential_mins = double(a_res_min_ind >= a_res_max_ind(step_count)).*double(a_res_min_ind <= a_res_max_ind(step_count+1));
    if exist('potential_mins','var') && sum(potential_mins) > 0
        IC(step_count) = min(a_res_min_ind(logical(potential_mins)));
    else
        IC(step_count) = NaN;
    end
end % step_count
% If no IC is found, NaN-flag
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
IC = REID_IMU_crash_catch(NaN,IC);
% Create tables and structures
timings = table;
segmented = struct;
% No output for stances with this method, so create dummy
stances = NaN;
% Correct timings back to original timestamps
timings.initial_contact = IC + data(1,1) - 1;
% No TC for this method, so fill with NaNs
timings.terminal_contact(1:size(IC,1)) = NaN;
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
% Loop through ICs and extract data segmented by gait events
for step_count = 1:size(IC,1) - 1
    segmented.(strcat("step_",num2str(step_count),'_',side)) = data(IC(step_count):IC(step_count+1)-1,:);
end % for step_count

end % function