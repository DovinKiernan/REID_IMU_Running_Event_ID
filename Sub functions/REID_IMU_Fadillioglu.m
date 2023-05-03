%% REID_IMU_Fadillioglu

% Fadillioglu et al 2020
% Uses angular velocity about the ML axis to identify IC and TC
% Adapted with premission from code provided by Cagla Fadillioglu on 22/10/17
% Minimal changes made to the original code; simply adapted to work with
% the same inputs and give the same output as other methods included in
% this package

function [timings, stances, segmented] = REID_IMU_Fadillioglu(data, location, Fs, min_stance_t)

% Filter signal
order = 2;
Fc = 15;
[b1, b2] = butter(order,Fc/(Fs/2));
GYROUS_15 = filtfilt(b1, b2, data(:,7));
% Complementary Signal for TC detection
Fc = 10;
[b1, b2] = butter(order, Fc/(Fs/2));
ABS1 = data(:,7) - GYROUS_15;
comp_sig = filtfilt(b1, b2, ABS1);
% Swing Phase prediction
[~, LOC] = findpeaks(GYROUS_15, 'MinPeakHeight', deg2rad(250), 'MinPeakDistance', round(Fs/3));
% IC detection
IC_index = 1;
GYROUS_15_detrend = detrend(GYROUS_15);
for LOC_count = 1:numel(LOC) - 1
    if LOC(LOC_count) + Fs < size(data,1) % added to prevent error from looking outside available data
        IC_rel = find(GYROUS_15_detrend(LOC(LOC_count):LOC(LOC_count) + Fs) < 0, 1, 'first');
        if ~isempty(IC_rel)
            IC(IC_index, 1) = IC_rel + LOC(LOC_count) - 1;
        else
            IC(IC_index, 1) = NaN;
        end
        IC_index = IC_index + 1;
    else % otherwise, just look to the end of the data set
        IC_rel = find(GYROUS_15_detrend(LOC(LOC_count):end) < 0, 1, 'first');
        if ~isempty(IC_rel)
            IC(IC_index, 1) = IC_rel + LOC(LOC_count) - 1;
        else
            IC(IC_index, 1) = NaN;
        end
        IC_index = IC_index + 1;
    end % if size < Data
end % for LOC_count
% TC detection
TC_index = 1;
for LOC_count = 2:numel(LOC)
    if LOC(LOC_count) - LOC(LOC_count - 1) < Fs*2 % to check if the peak is just in the beginning
        % Time between midswing maxima
        Y = LOC(LOC_count) - LOC(LOC_count - 1);
        % Finds a minimum for the second search in the following step
        [~, loc_min1] = min(GYROUS_15(LOC(LOC_count) - round(Y/2):LOC(LOC_count) - round(Y/10)));
        index_temp = loc_min1 + LOC(LOC_count) - round(Y/2) - 1;
        if Y < Fs % to distunguish between fast and slow events
            [~, loc_min2] = max(comp_sig(LOC(LOC_count) - round(Y/2):index_temp));
        else
            [~, loc_min2] = min(comp_sig(LOC(LOC_count) - round(Y/2):index_temp + round(Y/10)));
        end % if Y < Fs
        TC(TC_index, 1) = loc_min2 + LOC(LOC_count) - round(Y/2) - 1;
        TC_index = TC_index + 1;
    end % if location < 2 s
end % for LOC_count
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