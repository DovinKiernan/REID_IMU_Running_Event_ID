%% REID_IMU_Reenalda

% Reenalda et al. 2021
% Used the peak downward velocity of the pelvis to determine IC
% Used signal output from proprietary algorithms provided by XSens
% Given proprietary and undisclosed data processing used in the paper 
% we instead used the approach described in our paper to calculate the 
% coordinate systems and filter data
% No additional filtering or data manipulation occurs within this function

function [timings, stances, segmented] = REID_IMU_Reenalda(data, Fs, max_step_freq)

% Calculate velocity by integrating y-axis acceleration
y_vel = cumtrapz(data(:,3));
% Find peak times and locations in vertical (y) data
% Added a constraint that peaks must be separated by a minimum distance
pkdist = round(Fs/max_step_freq);
[mag_1, IC] = findpeaks(-y_vel,'MinPeakDistance',pkdist);
% No method is reported for TC or stance side determination
% If no IC is found, NaN-flag
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
IC = REID_IMU_crash_catch(NaN,IC);
% Create tables and structures
timings = table;
stances = NaN;
segmented = struct;
% Convert back into original timestamps
timings.initial_contact = IC + data(1,1) - 1;
timings.terminal_contact(1:size(IC,1)) = NaN;
timings.left_stance(1:size(IC,1)) = 0.5;
for step_count = 1:size(IC,1)-1
    segmented.(strcat("step_",num2str(step_count),"_unknown")) = data(IC(step_count):IC(step_count+1),:);
end % for step_count

end % function