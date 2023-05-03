%% REID_IMU_Auvinet

% Auvinet et al. 2002
% Auvinet identified IC as occuring prior to a maxima in vertical
% acceleration and minima in anterior-posterior acceleration
% TC was identified as the point after the vertical maxima when
% acceleration stopped decreasing
% Left and right were identified by a large maxima in medial-lateral
% acceleration on the stance side

function [timings, stances, segmented] = REID_IMU_Auvinet(data, Fs, max_step_freq, min_stance_t)

% Identify peaks in vertical acceleration
pkdist = round(Fs/max_step_freq);
[mag_1, vert_max_ind] = findpeaks(data(:,3),'MinPeakDistance',pkdist);
% Burn the first and last peak (they often capture erroneous points)
vert_max_ind(1) = [];
vert_max_ind(end) = [];
% Loop through and look between the remaining peaks...
for step_count = 1:size(vert_max_ind,1)-1
    [mag_2 vert_min_before_max(step_count,1)] = min(data(vert_max_ind(step_count):vert_max_ind(step_count+1),3));
    % Correct index
    vert_min_before_max(step_count,1) = vert_min_before_max(step_count,1) + vert_max_ind(step_count) - 1;
    % Should be a minima in a-p that occurs immediately after vert_min_before_max
    [mag_3, ap_min_ind(step_count,1)] = min(data(vert_min_before_max(step_count,1):vert_max_ind(step_count+1,1),2));
    % Correct index
    ap_min_ind(step_count,1) = ap_min_ind(step_count,1) + vert_min_before_max(step_count,1) - 1;
    % Now walk back to the point when values start decreasing
    ticker = 0;
    while data(ap_min_ind(step_count,1) - ticker,2) < data(ap_min_ind(step_count,1) - ticker - 1,2)
        ticker = ticker + 1;
    end
    ap_max_before_min(step_count,1) = ap_min_ind(step_count,1) - ticker;
    % Now take the mean of vert_min_before_max and ap_max_before_min, call it IC
    IC(step_count,1) = round(mean([ap_max_before_min(step_count,1) vert_min_before_max(step_count,1)],2));
    % Auvient 2002 states that TC occurs at end of loading phase
    % Per, Figure 2 in their paper this is ~ where a-p acceleration crosses
    % zero and begins to increase again
    % However, with processing approaches that subtract the gravity
    % component from acceleration this would correspond to ~ - 1
    % Therefore, we look for where VERT acceleration first goes below -1
    % i.e., segment is in ~freefall or no longer "loaded"
    first_vert_offload = find(data(vert_max_ind(step_count+1,1):end,3) <= -1,1,'first');
    first_vert_offload = first_vert_offload + vert_max_ind(step_count+1,1) - 1;
    % Call this TC
    if ~isempty(first_vert_offload)
        TC(step_count,1) = first_vert_offload;
    else
        TC(step_count,1) = NaN;
    end
    % Determine stance side by looking at ML acceleration at IC
    % "sharp ipsilateral acceleration at the IC"
    % Chose a range of ~10 ms on each side of IC
    if mean(data(IC(step_count,1)-round(Fs/100):IC(step_count,1)+round(Fs/100),4)) < 0
        % Sharp spike in ML acceleration is left sense
        % Call stance left
        left_stance(step_count,1) = 1;
    elseif mean(data(IC(step_count,1)-round(Fs/100):IC(step_count,1)+round(Fs/100),4)) > 0
        % Sharp spike in ML acceleration is right sense
        % Call stance right
        left_stance(step_count,1) = 0;
    else
        % In case of parity, call it 0.5
        left_stance(step_count,1) = 0.5;
    end
end % for step_count
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