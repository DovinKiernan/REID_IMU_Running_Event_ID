%% REID_IMU_Benson

% Benson et al. 2019
% Benson provides a number of well-documented functions to identify IC and
% TC points as well as stance side from the acceleration profiles of a
% sacrum-mounted accelerometer
% Here we have minimally adapted their code to work with the same inputs and
% provide the outputs as other methods in this package

function [timings, stances, segmented] = REID_IMU_Benson(data, Fs, min_stance_t)

% Preallocate for results
IC = [];
TC = [];
all_steps_norm = [];
% 10 Hz low-pass Butterworth filter all accelerometer data
% Filter parameters
order = 4; % order
Fc = 10; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
% Run filter
a_filt = filtfilt(b1,b2,data(:,2:4));
% Find peak times and locations in vertical (y) data
[pos_v_mags, pos_v_locs] = findpeaks(a_filt(:,2),Fs,'MinPeakDistance',0.25);
% Convert locations back to frames and round to whole number
pos_v_peaks = ceil((pos_v_locs(:,1)+1/Fs)*Fs);
% For each positive peak starting with the second one
for pos_v_peak_count = 2:size(pos_v_peaks,1)
    % Find negative peaks in AP signals within window of two positive peaks in the vertical signal
    [neg_ap_mags, neg_ap_locs] = findpeaks(-a_filt(pos_v_peaks(pos_v_peak_count-1,1):pos_v_peaks(pos_v_peak_count,1),1));
    % Adjust -ve AP peak locations to match with the full data set length
    neg_ap_locs = neg_ap_locs + pos_v_peaks(pos_v_peak_count-1,1) - 1;
    % Put in descending order so that peak closest to end of window is first
    [neg_ap_locs, new_order] = sort(neg_ap_locs, 'descend');
    % Resort the magnitudes by the order the peaks occur in
    neg_ap_mags = neg_ap_mags(new_order);
    % Depending on how many negative AP peaks there are...
    if isempty(neg_ap_locs)
        % If there is no -ve AP peak
        % Instead find minimum of AP data
        [~, min_loc] = min(a_filt(pos_v_peaks(pos_v_peak_count-1,1):pos_v_peaks(pos_v_peak_count,1),1));
        % Call this minima the IC
        IC(pos_v_peak_count,1) = min_loc + pos_v_peaks(pos_v_peak_count-1,1) - 1;
        % Find max of AP data from preceding peak to IC (AP minima)
        [~, max_idx] = max(a_filt(pos_v_peaks(pos_v_peak_count-1,1):IC(pos_v_peak_count,1),1));
        max_ap_idx = max_idx(1,1) + pos_v_peaks(pos_v_peak_count-1,1) - 1;
        % Look from max AP preceding IC to IC and find biggest slope in the ap acceleration (i.e., a-p jerk)
        jerk = diff(a_filt(max_ap_idx:IC(pos_v_peak_count,1),1))/(1/Fs);
        % Depending on size of jerk array...
        if size(jerk,1) < 3
            % If the jerk array is too short just find the highest value in jerk
            TC(pos_v_peak_count-1,1) = max_ap_idx;
        else
            % If the jerk array is long enough find the peak jerk
            [~, jerk_peak] = findpeaks(jerk, 'SortStr', 'descend', 'NPeaks', 1);
            % If there is no peak in jerk
            if isempty(jerk_peak)
                % Find number of frames equal to 10% of the size of the window between max of AP and IC; at least one frame
                ignore_end = ceil(size(jerk,1)/10);
                % If the jerk is 20 frames or shorter (ignore end is 2 frames)
                if size(jerk,1) <= 2*ignore_end
                    % Find max jerk
                    [~, b] = max(jerk(:,1));
                    jerk_peak = b(1,1);
                else
                    % Find max jerk (instead of findpeaks) ignoring end points
                    [~, b] = max(jerk(ignore_end+1:size(jerk,1)-ignore_end, 1));
                    jerk_peak = b(1,1)+ignore_end;
                end
            end
            % Correct time
            TC(pos_v_peak_count-1,1) = jerk_peak + max_ap_idx - 1;
        end % Slope size
    elseif size(neg_ap_locs,1) == 1
        % There is a single negative AP peak
        IC(pos_v_peak_count,1) = neg_ap_locs(1,1);
        % Find max of AP data between first positive peak in vertical data and IC
        [~, max_idx] = max(a_filt(pos_v_peaks(pos_v_peak_count-1,1):IC(pos_v_peak_count,1),1));
        % Get index in full reference frame
        max_ap_idx = max_idx(1,1) + pos_v_peaks(pos_v_peak_count-1,1) - 1;
        % Find jerk between max of AP and IC
        jerk = diff(a_filt(max_ap_idx:IC(pos_v_peak_count,1),1))/(1/Fs);
        % If not at least three points
        if size(jerk,1) < 3
            % Set TC for previous IC at location of max of AP data between first positive peak in vertical acceleration and IC
            TC(pos_v_peak_count-1,1) = max_ap_idx;
        else
            % Find highest peak in jerk
            [~, jerk_peak] = findpeaks(jerk, 'SortStr', 'descend', 'NPeaks', 1);
            % If no peak in slope
            if isempty(jerk_peak)
                % Find number of frames equal to 10% of the size of the window between max of AP and IC; at least one frame
                ignore_end = ceil(size(jerk,1)/10);
                if size(jerk,1) <= 2*ignore_end
                    % Find max slope ignoring end points
                    [~, b] = max(jerk(:,1));
                    jerk_peak = b(1,1);
                else
                    % Find max slope ignoring end points
                    [~, b] = max(jerk(ignore_end+1:size(jerk,1)-ignore_end, 1));
                    jerk_peak = b(1,1)+ignore_end;
                end
            end
            TC(pos_v_peak_count-1,1) = jerk_peak + max_ap_idx - 1;
        end
    else
        % If there are multiple negative AP peaks
        % Get order of peaks by height
        [~, pk_order] = sort(neg_ap_mags, 'descend');
        % Get order of peaks by timing
        [~, idx_order] = sort(neg_ap_locs, 'descend');
        % Get a reference ranking based on the number of peaks
        ref_rank = 1:size(neg_ap_mags,1);
        % for each peak
        for peak_count = 1:size(neg_ap_mags,1)
            % Find rank based on height
            pk_rank(peak_count,1) = ref_rank(pk_order == peak_count);
            % Find rank based on position
            idx_rank(peak_count,1) = ref_rank(idx_order == peak_count);
        end % peak_count
        % Find mean rank based on height and position
        mean_rank = mean([pk_rank, idx_rank],2);
        % Find lowest rank
        [~,low_rank] = min(mean_rank);
        % IC is lowest ranked peak (tallest & latest)
        IC(pos_v_peak_count,1) = neg_ap_locs(low_rank,1);
        % Remove IC peak and any peaks after from consideration (recall neg ap locs was sorted with latest as first entry)
        neg_ap_locs(1:low_rank,:) = [];
        mean_rank(1:low_rank,:) = [];
        if ~isempty(neg_ap_locs)
            prev_pk = 1;
            while prev_pk <= size(neg_ap_locs,1)
                % If a previous peak is less than 0.1s *BEFORE* the IC
                if abs(neg_ap_locs(prev_pk,1) - IC(pos_v_peak_count-1,1))/Fs < 0.1
                    % Remove from TC contention
                    neg_ap_locs(prev_pk,:) = [];
                else
                    prev_pk = prev_pk + 1;
                end
            end
        end
        if isempty(neg_ap_locs)
            % Find max of AP data
            [~, max_idx] = max(a_filt(pos_v_peaks(pos_v_peak_count-1,1):IC(pos_v_peak_count,1),1));
            max_ap_idx = max_idx(1,1) + pos_v_peaks(pos_v_peak_count-1,1) - 1;
            % Find slope between max of AP and IC
            jerk = diff(a_filt(max_ap_idx:IC(pos_v_peak_count,1),2))/(1/Fs);
            if size(jerk,1) < 3
                TC(pos_v_peak_count-1,1) = max_ap_idx;
            else
                % Find highest peak in slope (closest to zero)
                [~, jerk_peak] = findpeaks(jerk, 'SortStr', 'descend', 'NPeaks', 1);
                % If no peak in slope
                if isempty(jerk_peak)
                    ignore_end = ceil(size(jerk,1)/10);
                    if size(jerk,1) <= 2*ignore_end
                        % Find max slope ignoring end points
                        [~, b] = max(jerk(:,1));
                        jerk_peak = b(1,1);
                    else
                        % Find max slope ignoring end points
                        [~, b] = max(jerk(ignore_end+1:size(jerk,1)-ignore_end, 1));
                        jerk_peak = b(1,1)+ignore_end;
                    end
                end
                TC(pos_v_peak_count-1,1) = jerk_peak + max_ap_idx - 1;
            end
        else
            % TC is next negative peak closest to end
            TC(pos_v_peak_count-1,1) = max(neg_ap_locs);
        end
    end
end % pos_peak_count
if ~isempty(IC)
    IC(IC==0,:) = [];
end
if ~isempty(TC)
    TC(TC==0,:) = [];
end
% Ensure that an IC occurs first
while TC(1) < IC(1)
    TC(1) = [];
end
% Ensure that a TC occurs last
while IC(end) > TC(end)
    IC(end) = [];
end
% Loop through each gait event (IC and TC pair)
% Determine whether it is a left or right step
for gait_event_idx = 1:size(IC,1)
    ML_pos_pk = [];  % initialize as empty
    ML_neg_pk = [];  % initialize as empty
    side_identified = false;  % initialize as false
    % Find magnitude and location of largest positive ML peak within stance
    [ML_pos_pk, ML_pos_loc] = findpeaks(a_filt(IC(gait_event_idx):TC(gait_event_idx),3), 'SortStr', 'descend', 'NPeaks', 1);
    % Find magnitude and location of largest negative ML peak within stance
    [ML_neg_pk, ML_neg_loc] = findpeaks(-a_filt(IC(gait_event_idx):TC(gait_event_idx)-1,3), 'SortStr', 'descend', 'NPeaks', 1);
    % While side has not been identified
    while ~side_identified
        % If no positive AP peak
        if isempty(ML_pos_pk)
            % Set as left step
            left_step(gait_event_idx,1) = 1;
            % Set as true
            side_identified = true;
            % If no negative AP peak
        elseif isempty(ML_neg_pk)
            % Set as right step
            left_step(gait_event_idx,1) = 0;
            % Set as true
            side_identified = true;
            % If both positive and negative peaks exist check if negative peak is later than positive peak
        elseif ML_neg_loc > ML_pos_loc
            % If negative peak is closer to TC than to positive peak (correct)
            if abs(ML_neg_loc - (TC(gait_event_idx)-IC(gait_event_idx)+1)) < abs(ML_neg_loc-ML_pos_loc)
                % Set as right step
                left_step(gait_event_idx,1) = 0;
                % Set as true
                side_identified = true;
                % Positive peak is likely the initial peak right before the major negative peak (incorrect)
            else
                % Retain location of incorrect positive peak
                shift_loc = ML_pos_loc;
                if shift_loc < 0.015*Fs
                    % Find magnitude and location of largest positive peak after incorrect positive peak
                    [ML_pos_pk, ML_pos_loc] = findpeaks(a_filt(IC(gait_event_idx)+shift_loc:TC(gait_event_idx),3), 'SortStr', 'descend', 'NPeaks', 1);
                    % If a new positive peak was found
                    if ~isempty(ML_pos_pk)
                        % Get location of new positive peak in same frame reference as negative peak; while loop will be repeated
                        ML_pos_loc = ML_pos_loc + shift_loc;
                    end
                else
                    % Set as right step
                    left_step(gait_event_idx,1) = 0;
                    % Set as true
                    side_identified = true;
                end
            end
            % Positive peak is later than negative peak
        else
            % If positive peak is closer to TC than to negative peak (correct)
            if abs(ML_pos_loc - (TC(gait_event_idx)-IC(gait_event_idx)+1)) < abs(ML_pos_loc-ML_neg_loc)
                % Set as left step
                left_step(gait_event_idx,1) = 1;
                % Set as true
                side_identified = true;
                % Negative peak is likely the initial peak right before the major positive peak (incorrect)
            else
                % Retain location of incorrect negative peak
                shift_loc = ML_neg_loc;
                if shift_loc < 0.015*Fs
                    % Find magnitude and location of largest negative peak after incorrect negative peak
                    [ML_neg_pk, ML_neg_loc] = findpeaks(-a_filt(IC(gait_event_idx)+shift_loc:TC(gait_event_idx),3), 'SortStr', 'descend', 'NPeaks', 1);
                    % If a new negative peak was found
                    if ~isempty(ML_neg_pk)
                        % Get location of new negative peak in same frame reference as positive peak; while loop will be repeated
                        ML_neg_loc = ML_neg_loc + shift_loc;
                    end
                else
                    % Set as left step
                    left_step(gait_event_idx,1) = 1;
                    % Set as true
                    side_identified = true;
                end
            end
        end
    end
end
% This method has never errored out for me so I don't think it actually
% requires the crash_catch function but I put it here for consistency
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
[IC,TC,left_step] = REID_IMU_crash_catch(min_stance_t*Fs/1000,IC,TC,left_step);
% Create tables and structures
timings = table;
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact = TC + data(1,1) - 1;
timings.left_stance = left_step;
stances = table;
stances.time = data(:,1);
stances.left_stance = zeros(size(data,1),1);
stances.right_stance = zeros(size(data,1),1);
segmented = struct;
for step_count = 1:size(IC,1)-1
    if left_step(step_count) == 1
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