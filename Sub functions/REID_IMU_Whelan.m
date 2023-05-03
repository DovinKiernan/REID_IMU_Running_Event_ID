%% REID_IMU_Whelan

% Whelan et al. 2015
% Report identifying IC as the peak acceleration in the A-P axis
% Their figure, however, suggests that IC was actually identified as a
% local maximum immediately preceding a larger max
% We have followed this assumption here
% No TC identification

function [timings, stances, segmented] = REID_IMU_Whelan(data, location, Fs, max_step_freq)

% Filter at 10 Hz
order = 4; % order not specified, guessed 4th
Fc = 10; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
AP_filt = filtfilt(b1,b2,data(:,2));
% Identify peaks in unfiltered
% Added a minimum separation constraint
pkdist = round(Fs/max_step_freq)*2;
[mag_1, AP_max_ind] = findpeaks(AP_filt,'MinPeakDistance',pkdist);
% Whelan et al report using the "time the peak acceleration occurred" to
% define IC, however, their figure suggests they used a local maximum
% immediately preceding the peak acceleration
% Following this assumption, look back in the data and find the local
% maxima immediately preceding our peaks
% They do not provide a time window to look for this preceding peak and
% their figure's x-axis is not labelled; thus, we developed a novel approach: 
% First, starting from the second peak, look between peaks and find the minimum
% Second, look from the minimum to the proceeding peak and take the local max
if size(AP_max_ind,1) > 1
    for step_count = 1:size(AP_max_ind,1)-1
        % Find largest minimum between two AP peaks
        [mag_2, AP_min_ind] = findpeaks(-AP_filt(AP_max_ind(step_count):AP_max_ind(step_count+1)),'SortStr','descend','NPeaks',1);
        % Correct the timestamp
        AP_min_ind = AP_min_ind + AP_max_ind(step_count) + 1;
        % Find the largest max between the minimum and the AP peak
        [mag_3, IC_temp] = findpeaks(AP_filt(AP_min_ind:AP_max_ind(step_count+1)),'SortStr','descend','NPeaks',1);
        % Correct the timestamp
        if ~isempty(IC_temp)
            IC(step_count,1) = IC_temp + AP_min_ind + 1;
        else
            IC(step_count,1) = NaN;
        end
    end % step_count
else
    % If there aren't two AP peaks then make a dummy IC
    IC(step_count,1) = NaN;
end % If there are more than 2 AP peaks
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
% Correct IC timings back to original timestamps
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