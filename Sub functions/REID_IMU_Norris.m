%% REID_IMU_Norris

% Norris et al. 2015
% ML accelerometer data low-pass filtered at 2 Hz
% IC identified as positive zero-crossing
% No TC identification proposed

function [timings, stances, segmented] = REID_IMU_Norris(data, location, Fs)

% Low-pass filter ML data
order = 2; % order
Fc = 2; % cut-off
[b1, b2] = butter(order,Fc/(Fs/2),'low');
ML_filt = filtfilt(b1,b2,data(:,4));
% Look for positive zero crossings and call them IC
IC_ticker = 1;
for i = 2:size(ML_filt,1)
    zerocross(i-1) = ML_filt(i)*ML_filt(i-1) < 0;
    if zerocross(i-1) == 1 && ML_filt(i) > 0
        IC(IC_ticker,1) = i;
        IC_ticker = IC_ticker + 1;
    end
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
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
stances = NaN;
segmented = struct;
% Correct values
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

end % function