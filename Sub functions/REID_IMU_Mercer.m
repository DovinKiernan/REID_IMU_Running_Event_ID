%% REID_IMU_Mercer

% Mercer et al. 2003
% Mercer, I think, first describes this method, but it is based off
% previous work by Shorten 1992; Hamill 1995; Derrick 1998; Derrick 2002
% Initial contact is defined as minimum shank acceleration before peak acceleration
% Terminal contact is defined as minimum shank acceleration after a second local maximum

function [timings, stances, segmented] = REID_IMU_Mercer(data, location, Fs, max_step_freq, min_stance_t)

% Filter at 15 Hz
order = 4; % order
Fc = 15; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
data_filt = filtfilt(b1,b2,data);
% Identify peaks in unfiltered
pkdist = round(Fs/max_step_freq)*2;
[mag_1, unfilt_peak_ind] = findpeaks(data(:,3),"MinPeakDistance",pkdist);
% Identify peak in filtered
pkdist = round(Fs/max_step_freq)*2;
[mag_2, filt_peak_ind] = findpeaks(data_filt(:,3));
% Error if no peaks are found
if isempty(unfilt_peak_ind) || isempty(filt_peak_ind)
    warning('Data does not contain any peaks recognizable by the Mercer function')
end % isempty
% First and last IC may not correspond to consistent features so delete
if size(unfilt_peak_ind,1) > 2
    unfilt_peak_ind(end) = []; unfilt_peak_ind(1) = []; filt_peak_ind(end) = []; filt_peak_ind(1) = [];
end
% For the remaining peaks find the filtered peak closest to the unfiltered
% Then walk back from peak to point in filtered data when values start increasing
for step_count = 1:size(unfilt_peak_ind,1)
    [mag_3, closest_filt_ind(step_count)] = min(abs(unfilt_peak_ind(step_count)-filt_peak_ind));
    ticker = 0;
    while data_filt(filt_peak_ind(closest_filt_ind(step_count)) - ticker,3) - data_filt(filt_peak_ind(closest_filt_ind(step_count))-ticker-1,3) > 0
        ticker = ticker + 1;
    end
    IC(step_count,1) = filt_peak_ind(closest_filt_ind(step_count)) - ticker;
end % step_count
for step_count = 1:size(IC,1) - 1
    % Look between peak tib acceleration and IC and find the point where
    % shank acceleration begins increasing 
    % between the 2nd and 3rd maxima in the filtered data
    [mag_4, max_acc_loc] = findpeaks(data_filt(filt_peak_ind(closest_filt_ind(step_count)):IC(step_count+1),3),"NPeaks",4);
    max_acc_loc = max_acc_loc + filt_peak_ind(closest_filt_ind(step_count)) - 1;
    % Then walk forward from the 3rd maxima until values start increasing
    ticker = 0;
    if size(max_acc_loc,1) > 1
        while data_filt(max_acc_loc(2)+ticker+1,3) - data_filt(max_acc_loc(2)+ticker,3) < 0
            ticker = ticker + 1;
        end
        TC(step_count,1) = max_acc_loc(2) + ticker;
    elseif size(max_acc_loc,1) == 1
        while data_filt(max_acc_loc(1)+ticker+1,3) - data_filt(max_acc_loc(1)+ticker,3) < 0
            ticker = ticker + 1;
        end
        TC(step_count,1) = max_acc_loc(1) + ticker;
    else
        TC(step_count,1) = NaN;
    end
end % for step_count
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