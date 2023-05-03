%% REID_IMU_Mizrahi

% Mizrahi et al. 2000 method
% Defines initial contact as peak shank axial acceleration

function [timings, stances, segmented] = REID_IMU_Mizrahi(data, location, Fs, max_step_freq)

% Find peaks
pkdist = round(Fs/max_step_freq)*2;
[mag_1 peak_axial_acc_ind] = findpeaks(data(:,3),'MinPeakDistance',pkdist);
% Error if no peaks are found
if isempty(peak_axial_acc_ind)
    warning('Data does not contain any peaks recognizable by the Mizrahi function')
end % isempty
% First and last IC may not correspond to consistent features so delete
if size(peak_axial_acc_ind,1) > 2
    peak_axial_acc_ind(end) = []; peak_axial_acc_ind(1) = [];
end
IC = peak_axial_acc_ind;
% This method had no crashes with the data we used to test it so probably
% doesn't need the crash_catch function but here for consistency
% If no IC is found, NaN-flag
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
IC = REID_IMU_crash_catch(NaN,IC);
% Create tables and structures
timings = table;
segmented = struct;
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
for stride_count = 1:size(IC,1) - 1
    segmented.(strcat('stride_',num2str(stride_count),'_',side)) = data(IC(stride_count):IC(stride_count+1)-1,:);
end % for stride_count
% No output for stances with this method (doesn't have TC), so create dummy
stances = NaN;

end % function