%% REID_IMU_Wixted

% Wixted et al. 2010
% Identified IC as aproximately co-occuring with an AP acceleration braking spike
% Identified TC as aproximately co-occuring with vertical acceleration negative zero crossing
% Did not propose a method to identify stance side

function [timings, stances, segmented] = REID_IMU_Wixted(data, Fs, max_step_freq, min_stance_t)

pkdist = round(Fs/max_step_freq);
[mag_1 ap_min_ind] = findpeaks(-data(:,2),"MinPeakDistance",pkdist);
% Burn the first and last peaks, call ind1 IC
IC = ap_min_ind(2:end-1);
% Look for where vert acceleration first goes below 0
% We added a constraint here to start looking at min_stance_t - 10 ms (~85 ms)
for step_count = 1:size(IC,1)
    if size(data,1) > IC(step_count,1) + round(min_stance_t*Fs/1000) + 1
        first_vert_offload = find(data(IC(step_count,1)+round(min_stance_t*Fs/1000):end,3) <= 0,1,'first');
        first_vert_offload = first_vert_offload + IC(step_count,1) + round(min_stance_t*Fs/1000) - 1;
    end
    % Call this TC
    if exist('first_vert_offload','var') && ~isempty(first_vert_offload)
        TC(step_count,1) = first_vert_offload;
    else
        TC(step_count,1) = NaN;
    end
    clear first_vert_offload
end % step_count
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
% Convert back into original timestamps
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact = TC + data(1,1) - 1;
% Fill stance table with timestamps and zeros (no stance)
timings.left_stance(1:size(IC,1)) = 0.5;
stances.time = data(:,1);
stances.left_stance = zeros(size(data,1),1);
stances.right_stance = zeros(size(data,1),1);
% Loop through ICs and extract data segmented by gait events
for step_count = 1:size(IC,1)-1
    stances.left_stance(IC(step_count,1):TC(step_count,1)) = 0.5; % No method is reported for stance side determination so call both sides 0.5
    stances.right_stance(IC(step_count,1):TC(step_count,1)) = 0.5;
    segmented.(strcat("stance_",num2str(step_count),"_unknown")) = data(IC(step_count):TC(step_count),:);
    segmented.(strcat("swing_",num2str(step_count),"_unknown")) = data(TC(step_count)+1:IC(step_count + 1)-1,:);
end % for step_count

end % function