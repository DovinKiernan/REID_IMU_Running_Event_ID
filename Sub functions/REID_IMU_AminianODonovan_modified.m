%% REID_IMU_AminianODonovan_modified

% Aminian et al. 2002; O'Donovan et al. 2009
% Modified as described in our Supplementary Materials
% To identify IC and TC during running, O'Donovan 2009 used a method
% that was originally proposed by Aminian 2002 to identify IC and TC during walking
% This method enters ML (z) angular velocity into a wavelet multiresolution analysis (MRA)
% Neither Aminian nor O'Donovan define the +ve direction in their coordinate system
% However, in Aminian's discussion they describe details consistent with +right
% Figures are also consistent with this interpretation
% Aminian report a 10 level MRA
% Each MRA level represents a frequency band from (1/(2^(k+1)*1/Fs)):(1/(2^k*1/Fs))
% Thus, the sampling frequency (Fs) determines the frequency band of each level
% Aminian report a 200 Hz sampling frequency so detail level 1 should represent 50:100 Hz
% In contrast to this, Aminian state that they "only consider" up to 36 Hz
% Unclear where this discrepancy originates from
% However, to be consistent with the math and sampling frequency reported in their paper,
% and obtain equivalent frequency bands for each level, we resample to 200 Hz before performing the MRA
% Note: ODonovan also report angular velocities close to 200 rad/s (11,000 deg/s) for walking
% They don't report gyro range but values unreasonable for walking (or running)
% Perhaps mislabelled as rad/s when it was actually deg/s?

function [timings, stances, segmented] = REID_IMU_AminianODonovan_modified(data, location, Fs, max_step_freq, min_swing_t, max_swing_t, min_stance_t)

% Aminian report MRA frequency bands consistent with a 72 Hz collection frequency but explicitly state a collection frequency of 200 Hz
% Resample data to their sampling rate in order to achieve identical frequency bands from the MRA
Fs_Aminian = 200;
data_resample = scale_time(data(:,1:7),1,size(data,1),round(size(data,1)/Fs*Fs_Aminian));
% Check that data is long enough to run the MRA
% Needs to be at least 2^10 samples long to run the MRA because floor(log2(length(data))) is max level
if size(data_resample,1) < 2^10
    % If length is too short, zero-pad the data
    pad_length = 2^10-size(data_resample,1);
    data_resample = [data_resample; zeros(pad_length,7)];
    warning('Data set is too short to run the Aminian/ODonovan multiresolution analysis and has been zero-padded')
end
% 5th order Coiflet wavelet with 10 levels
% Rows 1:10 represent details of levels 1:10 (D_2^1:D_2^10)
% Row 11 represents approximation 10 (A_2^10)
mra_data_resample = modwtmra(modwt(data_resample(:,7),'coif5',10),'coif5');
% s_a is the sum of d_2^1:d_2^9
s_a = sum(mra_data_resample(1:9,:),1);
% New 10-level MRA then performed on s_a
mra_s_a = modwtmra(modwt(s_a,'coif5',10),'coif5');
% Two new approximations obtained
% First is designed to enhance the IC component by subtracting approximation 9 from 1
% Where approximation k is equal to A_2^10 + D_2^10 + ... + D_2^k+1
% i.e., row 11:k+1
a_1 = sum(mra_s_a(2:11,:),1);
a_9 = sum(mra_s_a(10:11,:),1); % Extremely low magnitude signal, subtracting it doesn't have a large effect
a_1_9 = a_1 - a_9;
% Second is designed to enhance the TC component by subtracting approximation 9 from 3
a_3 = sum(mra_s_a(4:11,:),1);
a_3_9 = a_3 - a_9;
% Trim off zero-padding if necessary
if exist('pad_length') && pad_length > 0
    a_1_9((size(data_resample,1)-pad_length+1):end) = [];
    a_3_9((size(data_resample,1)-pad_length+1):end) = [];
end
% Aminian report finding global maxima from each approximation
% It is not clear how this was done
% Therefore, we will do it by using find peaks and timings based off previously published running data
pkdist = round(Fs_Aminian/max_step_freq)*2;
[ms_a_1_9_mag ms_a_1_9_ind] = findpeaks(a_1_9,'MinPeakDistance',pkdist);
[ms_a_3_9_mag ms_a_3_9_ind] = findpeaks(a_3_9,'MinPeakDistance',pkdist);
% It seems unlikely given the similarity in the a_1_9 and a_3_9 signals
% but if ms_a_1_9 and ms_a_3_9 are different sizes it could crash the code
% so force them to have the same number of maxima
if size(ms_a_1_9_ind,2) ~= size(ms_a_3_9_ind,2)
    % Cross-correlate the indices of the identified maxima
    [r shift] = xcorr(ms_a_1_9_ind,ms_a_3_9_ind);
    [r_max_mag r_max_ind] = max(r);
    shift = shift(r_max_ind);
    % Negative values indicate ms_a_3_9 needs to shift left (has erroneous values at start)
    % Positive values indicate ms_a_1_9 needs to shift left (has erroneous values at start)
    if shift < 0
        ms_a_3_9_ind(1:abs(shift)) = [];
        ms_a_3_9_mag(1:abs(shift)) = [];
    elseif shift > 0
        ms_a_1_9_ind(1:abs(shift)) = [];
        ms_a_1_9_mag(1:abs(shift)) = [];
    end % if shift <> 0
end % if uneven sizes
% % % Burn the first and last peaks to ensure we're grabbing consistent signal features
% % ms_a_1_9_val(1) = []; ms_a_1_9_val(end) = [];
% % ms_a_1_9_ind(1) = []; ms_a_1_9_ind(end) = [];
% % ms_a_3_9_val(1) = []; ms_a_3_9_val(end) = [];
% % ms_a_3_9_ind(1) = []; ms_a_3_9_ind(end) = [];
% Aminian state these peaks ~correspond to mid swing
% Look forward and back from these peaks for IC and TC, respectively
for step_count = 1:size(ms_a_1_9_ind,2)
    % Aminian find local minima +0.25 to +2.00s after ms_a_1_9_ind
    % Although Aminian developed these timings based off walking,
    % O'Donovan do not report any adjustments to the timings for their
    % application to running
    % Our results suggest that these timing windows do not work for running
    % The +0.25 s delay from midswing consistently misses the IC
    % The -2 s window from midswing grabs TCs from 1-3 steps prior to the current step and often grabs ICs
    % Therefore, we have modified these timings based on running
    % Aminian state that the identified maxes correspond to midswing
    % Therefore we will look min_swing_t/1000/2 to max_swing_t/1000/2
    % When stepping forward in time ensure we don't go outside data window
    if size(a_1_9,2) - (ms_a_1_9_ind(step_count)+round(min_swing_t/1000/2*Fs_Aminian)) < 3
        % The window doesn't exist in our data so NaN-flag
        local_min_a_1_9_mag{step_count} = NaN;
        local_min_a_1_9_ind{step_count} = NaN;
    elseif size(a_1_9,2) - (ms_a_1_9_ind(step_count)+round(max_swing_t/1000/2*Fs_Aminian)) >= 0
        % The entire window exists in our data
        [local_min_a_1_9_mag{step_count} local_min_a_1_9_ind{step_count}] = findpeaks(-a_1_9(ms_a_1_9_ind(step_count)+round(min_swing_t/1000/2*Fs_Aminian):ms_a_1_9_ind(step_count)+round(max_swing_t/1000/2*Fs_Aminian)),'MinPeakHeight',0);
        % If there's no peaks
        if isempty(local_min_a_1_9_ind{step_count})
            % Then find the lowest value (max to be consistent with -signal used in find peaks)
            [local_min_a_1_9_mag{step_count} local_min_a_1_9_ind{step_count}] = max(-a_1_9(ms_a_1_9_ind(step_count)+round(min_swing_t/1000/2*Fs_Aminian):ms_a_1_9_ind(step_count)+round(max_swing_t/1000/2*Fs_Aminian)));
        end % if no peaks
    else
        % The start of the window exists in our data but the end point does not
        [local_min_a_1_9_mag{step_count} local_min_a_1_9_ind{step_count}] = findpeaks(-a_1_9(ms_a_1_9_ind(step_count)+round(min_swing_t/1000/2*Fs_Aminian):end),'MinPeakHeight',0);
        % If there's no peaks
        if isempty(local_min_a_1_9_ind{step_count})
            % Then find the lowest value (max to be consistent with -signal used in find peaks)
            [local_min_a_1_9_mag{step_count} local_min_a_1_9_ind{step_count}] = max(-a_1_9(ms_a_1_9_ind(step_count)+round(min_swing_t/1000/2*Fs_Aminian):end));
        end % if no peaks
    end % If window is outside data range
    local_min_a_1_9_ind{step_count} = local_min_a_1_9_ind{step_count} + ms_a_1_9_ind(step_count) + round(min_swing_t/1000/2*Fs_Aminian) - 1;
    local_min_a_1_9_mag{step_count} = local_min_a_1_9_mag{step_count}*-1;
    % Aminian find local minima -2.00 to -0.05 s before ms_a_3_9_ind
    % Again we modify this to min_swing_t/1000/2 to max_swing_t/1000/2
    % When stepping back in time, ensure we don't go outside data window
    if ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian) < 1 ||...
            ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian) < 3
        % The window doesn't exist in our data so NaN-flag
        local_min_a_3_9_mag{step_count}= NaN; local_min_a_3_9_ind{step_count} = NaN;
    elseif ms_a_3_9_ind(step_count)-round(max_swing_t/1000/2*Fs_Aminian) < 1
        % The start exists but the end doesn't
        [local_min_a_3_9_mag{step_count} local_min_a_3_9_ind{step_count}] = findpeaks(-a_3_9(1:ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian)),'MinPeakHeight',0);
        % If there's no peaks
        if isempty(local_min_a_3_9_ind{step_count})
            % Then find the lowest value (max to be consistent with -signal used in find peaks)
            [local_min_a_3_9_mag{step_count} local_min_a_3_9_ind{step_count}] = max(-a_3_9(1:ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian)));
        end % if no peaks
    else
        % The entire window exists
        [local_min_a_3_9_mag{step_count} local_min_a_3_9_ind{step_count}] = findpeaks(-a_3_9(ms_a_3_9_ind(step_count)-round(max_swing_t/1000/2*Fs_Aminian):ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian)),'MinPeakHeight',0);
        % If there's no peaks
        if isempty(local_min_a_3_9_ind{step_count})
            % Then find the lowest value (max to be consistent with -signal used in find peaks)
            [local_min_a_3_9_mag{step_count} local_min_a_3_9_ind{step_count}] = max(-a_3_9(ms_a_3_9_ind(step_count)-round(max_swing_t/1000/2*Fs_Aminian):ms_a_3_9_ind(step_count)-round(min_swing_t/1000/2*Fs_Aminian)));
        end % if no peaks
        local_min_a_3_9_ind{step_count} = local_min_a_3_9_ind{step_count} + ms_a_3_9_ind(step_count) - round(max_swing_t/1000/2*Fs_Aminian) - 1;
    end % If window is outside data range
    local_min_a_3_9_mag{step_count} = local_min_a_3_9_mag{step_count}*-1;
    % Look forward +0 to +0.075s from each local_min_a_3_9 to find the minima in original signal
    for TC_possibility = 1:size(local_min_a_3_9_ind{step_count},2)
        % Check for NaN
        if isnan(local_min_a_3_9_ind{step_count}(TC_possibility))
            local_min_data_resample_mag{step_count}(TC_possibility) = NaN; local_min_data_resample_ind{step_count}(TC_possibility) = NaN;
        elseif size(data_resample,1) - (local_min_a_3_9_ind{step_count}(TC_possibility)+round(0.075*Fs_Aminian)) > 0
            [local_min_data_resample_mag{step_count}(TC_possibility) local_min_data_resample_ind{step_count}(TC_possibility)] = min(data_resample(local_min_a_3_9_ind{step_count}(TC_possibility):local_min_a_3_9_ind{step_count}(TC_possibility)+round(0.075*Fs_Aminian),7));
            local_min_data_resample_ind{step_count}(TC_possibility) = local_min_data_resample_ind{step_count}(TC_possibility) + local_min_a_3_9_ind{step_count}(TC_possibility) - 1;
        else
            [local_min_data_resample_mag{step_count}(TC_possibility) local_min_data_resample_ind{step_count}(TC_possibility)] = min(data_resample(local_min_a_3_9_ind{step_count}(TC_possibility):end,7));
            local_min_data_resample_ind{step_count}(TC_possibility) = local_min_data_resample_ind{step_count}(TC_possibility) + local_min_a_3_9_ind{step_count}(TC_possibility) - 1;
        end % if isnan
    end % TC_possibility
end % for step_count
% Look back from the first potential IC to the preceding TCs
% Accept the first TC that meets the condition 0.1 s < (IC - TC) < 2.5 s
% i.e., swing phase is between 0.1 and 2.5s long
% Again modified for running to accept 200 to 600 ms
% If no TC meets this condition, then iterate to the next potential IC if possible
for step_count = 1:size(local_min_a_1_9_ind,2)
    % Prepopulate with NaNs in case no pairings satisfy conditions
    IC(step_count) = NaN;
    TC(step_count) = NaN;
    % Ensure there's more than just NaNs in the potential IC and TC arrays
    if any(~isnan(local_min_data_resample_ind{step_count})) && any(~isnan(local_min_a_1_9_ind{step_count}))
        potential_IC_tick = 1;
        IC_identified = 0;
        while IC_identified == 0 && potential_IC_tick <= size(local_min_a_1_9_ind{step_count},2)
            potential_TC_tick = size(local_min_data_resample_ind{step_count},2);
            TC_identified = 0;
            while TC_identified == 0 && potential_TC_tick > 0
                if min_swing_t/1000*Fs_Aminian < local_min_a_1_9_ind{step_count}(potential_IC_tick) - local_min_data_resample_ind{step_count}(potential_TC_tick) &&...
                        local_min_a_1_9_ind{step_count}(potential_IC_tick) - local_min_data_resample_ind{step_count}(potential_TC_tick) < max_swing_t/1000*Fs_Aminian
                    % Call these IC and TC
                    IC(step_count) = local_min_a_1_9_ind{step_count}(potential_IC_tick);
                    TC(step_count) = local_min_data_resample_ind{step_count}(potential_TC_tick);
                    IC_identified = 1;
                    TC_identified = 1;
                else
                   	% Otherwise increment down to look at the next potential TC
                    potential_TC_tick = potential_TC_tick - 1;
                end % min swing t < (IC - TC) < max swing t
            end % TC_identified = 0
            potential_IC_tick = potential_IC_tick + 1;
        end % IC_identified = 0
    end % data exists in potential IC and TC arrays
end % step_count
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Convert back to original frequency and transpose
IC = IC'*Fs/Fs_Aminian;
TC = TC'*Fs/Fs_Aminian;
% Due to rounding, could end up with frames outside of data range so check for and correct it
IC(IC > size(data,1)) = size(data,1);
TC(TC > size(data,1)) = size(data,1);
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