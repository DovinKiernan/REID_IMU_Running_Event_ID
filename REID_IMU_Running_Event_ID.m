%% REID_IMU -- Running Event ID

% Dovin Kiernan
% University of California Davis Human Performance Lab
% 23/03/31

%% LICENSING

% Original content distributed under GNU GENERAL PUBLIC LICENSE v3
% For licensing of dependencies, consult original sources
% For licensing of redistributed, non-original, code (Fadillioglu, Schmidt, Bach, and Benson methods) consult original sources

%% PURPOSE

% This function takes accelerometer and/or gyroscope data from either the
% shank or low-back/sacrum and enters it into the user's choice of one of several
% possible gait event ID methods. The function then returns timings of
% initial and terminal contact events for the left and/or right foot.

%% INPUTS

% Data 
% % Matrix with each row corresponding to a single frame and...
% % Column 1 -- Time stamps (in ms)
% % Columns 2:4 -- Linear acceleration x, y, and z in g (~9.81 m/s^2 depending on location data were collected)
% % Columns 5:7 -- Angular velocity x, y, and z in rad/s

% Note on data conventions
% % All coordinate systems per ISB and as described in Wu & Cavanagh 1995 and our paper
% % +x is anterior (direction of progression)
% % +y is proximal (or vertical)
% % +z is right

% Method
% % String specifying the method to be used (last name of first author)
% % Shank methods
% % % 'Mizrahi' (Mizrahi et al. 2000)
% % % 'Mercer' (Mercer et al. 2003)
% % % 'Purcell' (Purcell et al. 2006)
% % % 'AminianODonovan' (Aminian et al. 2002; ODonovan et al. 2009)
% % % 'AminianODonovan_modified' (Aminian et al. 2002; ODonovan et al. 2009; modified per our Supplemental Material)
% % % 'GreeneMcGrath' (Greene et al. 2010; McGrath et al. 2012)
% % % 'GreeneMcGrath_modified' (Greene et al. 2010; McGrath et al. 2012; modified per our Supplemental Material)
% % % 'Sinclair' (Sincalir et al. 2013)
% % % 'Whelan' (Whelan et al. 2015)
% % % 'Norris' (Norris et al. 2015)
% % % 'Schmidt' (Schmidt et al. 2016)
% % % 'Aubol' (Aubol et al. 2020)
% % % 'Fadillioglu' (Fadillioglu et al. 2020)
% % % 'Bach' (Bach et al. 2022)
% % % 'Bach_modified' (Bach et al. 2022; modified per our Supplemental Material)
% % Low-back/sacrum methods
% % % 'Auvinet' (Auvinet et al. 2002)
% % % 'Lee' (Lee et al. 2010)
% % % 'Wixted' (Wixted et al. 2010)
% % % 'Bergamini' (Bergamini et al. 2012)
% % % 'Benson' (Benson et al. 2019)
% % % 'Reenalda' (Reenalda et al. 2021)

% Location
% String specifying wearable location
% % 'Left shank'
% % 'Right shank'
% % For sacrum/low-back methods no location input is necessary

% Supervision
% % Only for Schmidt method
% % 'Supervised' (Requires UI to define thresholds for each trial)
% % 'Unsupervised' (Uses default values, no UI)

%% OUTPUTS

% To the extent possible for any given method...
% timings -- Table with columns containing:
% % timings of initial contact 
% % timing of terminal contact
% % boolean for left stance (1 = left stance true, 0 = left stance false/right stance true)
% stances -- Table with columns containing
% % timestamps
% % left stance true (1) or false (0)
% % right stance true (1) or false (0)
% segmented -- Structure with fields for each stance and swing (or, if not possible, for each stride or step)
% % Labelled by side and occurence (from 1 ... n)
% % Data are in original matrix form with timestamps, x, y, z acceleration and x, y, z angular velocity

%% DEPENDENCIES

% This main function calls on sub-functions for each method
% Some of these method sub-functions are themselves dependent on supporting functions and data
% These sub- and supporting-functions are included with the download

% List of supporting functions and data:
% % MATLAB's Wavelet Toolbox -- published at https://www.mathworks.com/products/wavelet.html
% % MATLAB's Signal Processing Toolbox -- published at https://www.mathworks.com/products/signal.html
% % MATLAB's Statistics and Machine Learning Toolbox -- published at https://www.mathworks.com/products/statistics.html
% % MATLAB's Deep Learning Toolbox -- published at https://www.mathworks.com/products/deep-learning.html
% % scaletime -- published by Jan at https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime
% % Bach_ESN.mat -- Echo State Network trained and validated on running data published by Bach et al. 2022 at https://github.com/marlow17/PredictingGroundReactionForces
% % estimate_events.m -- published by Bach et al. 2022 at https://github.com/marlow17/PredictingGroundReactionForces
% % compute_state_matrix.m -- found in ESNToolbox published by Manta Lukoševičius at https://www.ai.rug.nl/minds/uploads/ESNToolbox.zip
% % test_ESN.m -- found in ESNToolbox published by Manta Lukoševičius at https://www.ai.rug.nl/minds/uploads/ESNToolbox.zip
% % wavelet-mediated differentiation method published by Jianwen Luo at https://www.mathworks.com/matlabcentral/fileexchange/13948-numerical-differentiation-based-on-wavelet-transforms, specifically including:
% % % derivative_dwt.m
% % % spline_scale_function.m
% % % spline_function.m
% % REID_IMU_crash_catch.m -- original function to ensure ICs and TCs from each method are standardized and can be used to generate the desired output and that common crashes are avoided

%% CODE

% This wrapper function...
% (1) finds the necessary sub- and supporting-functions
% (2) defines some constants then 
% (3) calls the sub-function for the specified method

function [timings, stances, segmented] = REID_IMU_Running_Event_ID(data, method, location, supervision)

% UI select the path where the sub- and supporting-functions are stored
subfunction_path = uigetdir('Path to gait event ID subfunctions folder');
addpath(subfunction_path);
% % addpath(genpath('SUBFUNCTIONPATHSTRING')) % get rid of UI

% Automation constraints
% To provide novel automation for some methods, constraints were introduced based on
% human running values previously reported across Cavanagh & LaFortune 1980;
% Munro et al., 1987; Cavanagh & Kram 1989; Williams et al, 1991;  DeWit et al., 2000;
% Weyand et al., 2000; Leskinen et al., 2009; Weyand et al., 2010; Meardon et al., 2011;
max_step_freq = 4.75; % maximum 4.75 steps per second
min_stance_t = 95; % minimum 95 ms stance time
max_stance_t = 270; % maximum 270 ms stance time
min_swing_t = 200; % minimum 200 ms swing time
max_swing_t = 600; % maximum 600 ms swing time

% Calculate sample frequency from time stamps
Fs = mean(diff(data(:,1)))*1000;
% % % % Alternatively, UI specify your sample frequency
% % % Fs = cell2mat(inputdlg('Enter sample frequency'));

% Choose which method to call
switch method
    case 'Mizrahi'
        [timings, stances, segmented] = REID_IMU_Mizrahi(data, location, Fs, max_step_freq);
    case 'Mercer'
        [timings, stances, segmented] = REID_IMU_Mercer(data, location, Fs, max_step_freq, min_stance_t);
    case 'Purcell'
        [timings, stances, segmented] = REID_IMU_Purcell(data, location, Fs, max_step_freq, min_stance_t);
    case 'AminianODonovan'
        [timings, stances, segmented] = REID_IMU_AminianODonovan(data, location, Fs, max_step_freq, min_stance_t);
    case 'AminianODonovan_modified'
        [timings, stances, segmented] = REID_IMU_AminianODonovan_modified(data, location, Fs, max_step_freq, min_swing_t, max_swing_t, min_stance_t);
    case 'GreeneMcGrath'
        [timings, stances, segmented] = REID_IMU_GreeneMcGrath(data, location, Fs, min_stance_t);
    case 'GreeneMcGrath_modified'
        [timings, stances, segmented] = REID_IMU_GreeneMcGrath_modified(data, location, Fs, max_step_freq, max_swing_t, min_stance_t);
    case 'Sinclair'
        [timings, stances, segmented] = REID_IMU_Sinclair(data, location, Fs, max_step_freq, min_stance_t);
    case 'Whelan'
        [timings, stances, segmented] = REID_IMU_Whelan(data, location, Fs, max_step_freq);
    case 'Norris'
        [timings, stances, segmented] = REID_IMU_Norris(data, location, Fs);
    case 'Schmidt'
        [timings, stances, segmented] = REID_IMU_Schmidt(data, location, Fs, supervision, min_stance_t);
    case 'Aubol'
        [timings, stances, segmented] = REID_IMU_Aubol(data, location, Fs, max_step_freq);
    case 'Fadillioglu'
        [timings, stances, segmented] = REID_IMU_Fadillioglu(data, location, Fs, min_stance_t);
    case 'Bach'
        [timings, stances, segmented] = REID_IMU_Bach(data, location, Fs, min_stance_t);
    case 'Bach_modified'
        [timings, stances, segmented] = REID_IMU_Bach_modified(data, location, Fs, min_stance_t);
    case 'Auvinet'
        [timings, stances, segmented] = REID_IMU_Auvinet(data, Fs, max_step_freq, min_stance_t);
    case 'Lee'
        [timings, stances, segmented] = REID_IMU_Lee(data, Fs, max_step_freq, min_stance_t);
    case 'Wixted'
        [timings, stances, segmented] = REID_IMU_Wixted(data, Fs, max_step_freq, min_stance_t);
    case 'Bergamini'
        [timings, stances, segmented] = REID_IMU_Bergamini(data, Fs, max_step_freq, min_stance_t);
    case 'Benson'
        [timings, stances, segmented] = REID_IMU_Benson(data, Fs, min_stance_t);
    case 'Reenalda'
        [timings, stances, segmented] = REID_IMU_Reenalda(data, Fs, max_step_freq);
end % switch method

end % function