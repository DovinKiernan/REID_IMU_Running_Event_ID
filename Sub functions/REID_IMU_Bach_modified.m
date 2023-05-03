%% REID_IMU_Bach_modified

% Bach et al. 2022
% Adapted from code posted online by Bach et al. at https://github.com/marlow17/PredictingGroundReactionForces
% to work with the same inputs/give the same outputs as other methods in
% this package
% Uses an Echo State Network (sparsely randomly connected reservoir) to estimate vertical ground reaction forces
% Then uses estimated vertical ground reaction forces to estimate IC and TC timings
% The Echo State Network has been previously trained and validated using 75 and 25%, respectively, of Bach 2022's originally provided running data
% Per Bach, 100 networks were trained. We then retained the single network best able to estimate IC and TC events (lowest mean absolute error in validation) for use here
% This trained network and all necessary functions are included with this code as...
% % Bach_ESN.mat
% % test_ESN.m
% % estimate_events.m
% % scale_time.m
% % compute_state_matrix.m
 % Due to the poor performance of the method when reproduced faithfully, we introduce two modifications to originally described code in this version

function [timings, stances, segmented] = REID_IMU_Bach_modified(data, location, Fs, min_stance_t)

% Check if Echo State Network is loaded
if ~exist('Bach_ESN','var')
    % If the Echo State Network isn't loaded yet, then load it
    load('Bach_ESN.mat');
    % Set it to global so it exists outside the function and doesn't need
    % to load every time the function is called
    global Bach_ESN
end
% Process data for use by the Echo State Network
% Resample to match Bach paper Fs = 2000/14 ~ 142.9 Hz
% Higher sampling rate of our test data was introducing higher frequency signal features (e.g., impact transients) that were degrading performance
% If an ESN was trained on data sampled at a higher frequency this problem may be ameliorated
Fs_Bach = 2000/14;
a_resample = scale_time(data(:,2:4),1,size(data,1),round(size(data,1)/Fs*Fs_Bach));
% Calculate Principal Components of acceleration
[x, a_PC] = pca(a_resample,'algorithm','eig');
% Retain only the 1st PC
a_PC = a_PC(:,1);
% If negative, multiply by -1
if x(2,1) < 0
    a_PC = -a_PC;
end
% Normalize by standard deviation
a_PC = a_PC/std(a_PC);
% Estimate "velocity" and "position" (integrals of 1st PC) after removing the DC-component with high pass filtering
order = 2;
Fc = 1;
[b1, b2] = butter(order, Fc/(Fs/2), 'high');
a_PC_int1 = filtfilt(b1, b2, cumtrapz(1/Fs, a_PC));
a_PC_int2 = filtfilt(b1, b2, cumtrapz(1/Fs, a_PC_int1));
% Arrange data in proper format for Echo State Network
input{1} = [a_PC, a_PC_int1, a_PC_int2];
input = cellfun(@(s)[ones(size(s,1),1),s], input, 'UniformOutput', false);
% Feed into the Echo State Network
% Requires the test_esn function
% Originally found in ESNToolbox provided by Manta Lukoševičius at https://www.ai.rug.nl/minds/uploads/ESNToolbox.zip
% Outputs are vertical ground reaction forces (vGRFs)
transient = round(250/1000*Fs_Bach);
esn_out = test_esn(input{1}, Bach_ESN, transient);
% Note, we introduced two modifications here:
% First, we observed a lot of negative values between characteristic vGRF waveforms
% Not sure if these values originate from signals that the ESN wasn't trained on
% In any case, these values do not make sense
% Thus, we zero out all negative values in esn_out
esn_out(esn_out<0) = 0;
% Second, the "estimate_events" function inverts GRF
% This inversion, however, negatively impacts the estimation of gait events
% It is unclear why this was originally included in the code
% However, given the negative impact of the inversion we feed in an inverted esn_out to cancel it out
% These vertical GRFs can be used to estimate IC and TC using the estimateEvents function provided by Bach 2022
[IC, TC, vGRF_predicted] = estimate_events(-esn_out, Fs_Bach);
% Convert back into original frames
IC = round(IC*Fs/Fs_Bach);
TC = round(TC*Fs/Fs_Bach);
% If no IC or TC are found, NaN-flag
if ~exist('TC','var') || isempty(TC)
    TC = NaN;
end
if ~exist('IC','var') || isempty(IC)
    IC = NaN;
end
% Correction for common crashes
[IC,TC] = REID_IMU_crash_catch(min_stance_t*Fs/1000,IC,TC);
% Create structures and tables
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