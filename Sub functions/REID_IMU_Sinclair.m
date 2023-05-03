%% REID_IMU_Sinclair

% Sinclair (2013)
% Uses ~longitudinal acceleration from a shank-mounted accelerometer to define IC and TC
% First filters data then finds peaks (we added a constraint that they must be separated by a minimum time)
% Then walks back and finds a zero-crossing with 20 ms of positive data
% following it, defines this as IC
% Then finds the largest local maxima between longitudinal peaks and finds
% a plateau where the data magnitude is NOT decaying by 2% of the current
% signal magnitude per second for at least 5% of minimum stance t
% Calls the first frame of that plateau TC

function [timings, stances, segmented] = REID_IMU_Sinclair(data, location, Fs, max_step_freq, min_stance_t)

% Filter at 60 Hz
order = 4; % order
Fc = 60; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
data_filt = filtfilt(b1,b2,data);

% Find maxima in y-axis acceleration ("peak tibial shock")
% Added a pkdist constraint based on expected timings during running
pkdist = round(Fs/max_step_freq)*2;
[mag_1, y_max_ind] = findpeaks(data_filt(:,3),'MinPeakDistance',pkdist);

% From the y-axis maxima, walk back to the zero-crossing
for y_max_count = 1:size(y_max_ind,1)
    tick = 1;
    zero_found = 0;
    % While there's still data to look through and we haven't found the zero-cross...
    while y_max_ind(y_max_count) - tick > 0 &&...
            zero_found == 0
        % If the signal crosses below zero ensure the 20 ms of data following the zero-crossing are all positive
        % Also ensure the window exists within our data
        if data_filt(y_max_ind(y_max_count)-tick,3) < 0 &&...
            (y_max_ind(y_max_count)-tick+round(20*Fs/1000) < size(data_filt,1)) && ...
            all(data_filt(y_max_ind(y_max_count)-tick+1:y_max_ind(y_max_count)-tick+round(20*Fs/1000),3) > 0)
            
            % Flag the last frame that was above zero as the zero-crossing
            zero_cross(y_max_count) = y_max_ind(y_max_count) - tick + 1;
            zero_found = 1;
        else
            % Otherise, increment up to keep searching for the zero-crossing
            tick = tick + 1;
        end % if crossed below 0
    end
    % If all data has been exhausted and no zero-crossing has been found, NaN-flag
    if y_max_ind(y_max_count) - tick == 0 &&...
            zero_found == 0
        zero_cross(y_max_count) = NaN;
    end
end % for y_max_count
% Language is somewhat unclear for TC identification
% "[TC] was determined using target pattern recognition with a 2% tolerance 
% as the first plateau in the descent phase of the second peak"
% I believe they're saying they called a plateau after a local max between peak axial accelerations the TC
% and that they are defining a plateau as a period in the data where the value is decreasing by less than 2%
% We first interpreted this as a 2% change in magnitude frame-frame (i.e., frame n is 2% of frame n's magnitude off frame n-1; at 1000 Hz)
% However, that constraint was very low and too many frames were accepted as plateaus
% Then tried 2% of local max value change frame-frame, still too low
% Then tried 2% per second change which seems to work better
% Therefore, going with...
% Find a local maximum between y-axis maxima
for y_max_count = 1:size(y_max_ind,1)-1
    % If valid peaks exist...
    if ~isnan(y_max_ind(y_max_count)) && ~isnan(y_max_ind(y_max_count+1))
        [mag_2, local_max_ind] = findpeaks(data_filt(y_max_ind(y_max_count):y_max_ind(y_max_count+1),3),'SortStr','descend','NPeaks',1);
        % Correct frame
        local_max_ind = local_max_ind + y_max_ind(y_max_count) + 1;
        % Take the difference between frames between the local max and the next peak axial acceleration
        % Then find frames where the magnitude is changing by <= -2% of the signal's current magnitude per second
        delta_2 = diff(data_filt(local_max_ind:y_max_ind(y_max_count+1),3))./data_filt(local_max_ind+1:y_max_ind(y_max_count+1),3)*Fs;
        delta_2(delta_2>=-0.02) = 1;
        delta_2(delta_2<-0.02) = 0;
        % The minimum duration of plateau was not specified in original paper
        % However, we digitized their data and were able to determine that the plateau presented in their figure was ~6% of stance t
        % Based on this, we'll look for the first delta_2 with at least 3% stance t duration occurring after the local max in time
        potential_plat = find(delta_2==1);
        plat_found = 0;
        plat_tick = 1;
        while plat_found == 0 &&...
                plat_tick < (size(potential_plat,1)-round(0.05*min_stance_t*Fs/1000))
            if all(diff(potential_plat(plat_tick:plat_tick+round(0.05*min_stance_t*Fs/1000)))==1)
                plat_start(y_max_count) = potential_plat(plat_tick) + local_max_ind + 1;
                plat_found = 1;
            else
                plat_tick = plat_tick + 1;
            end
        end % plat not found
        % If all data has been exhausted and no plateau has been found
        if plat_found == 0 &&...
                plat_tick >= (size(potential_plat,1)-round(0.05*min_stance_t*Fs/1000))
            plat_start(y_max_count) = NaN;
        end
    else % no zero-crossing was found so we can't find the plateau
        plat_start(y_max_count) = NaN;
    end % no NaNs
end % for y_max_count
IC = zero_cross';
TC = plat_start';
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