REID_IMU -- Running Event ID

Written by Dovin Kiernan
University of California Davis Human Performance Lab
23/03/31

Submitted to Sensors May 2023 as "Unsupervised gait event detection using a single wearable accelerometer and/or gyroscope: 
A comparison of methods across running speeds, surfaces, and foot strike patterns"

LICENSING

Original content distributed under GNU GENERAL PUBLIC LICENSE
For licensing of dependencies, consult original sources
For licensing of redistributed, non-original, code (Fadillioglu, Schmidt, Bach, and Benson methods) consult original sources

PURPOSE

This function takes accelerometer and/or gyroscope data from either the
shank or low-back/sacrum and enters it into the user's choice of one of several
possible gait event ID methods. The function then returns timings of
initial and terminal contact events for the left and/or right foot.

INPUTS

Data 
% Matrix with each row corresponding to a single frame and...
% Column 1 -- Time stamps (in ms)
% Columns 2:4 -- Linear acceleration x, y, and z in g (~9.81 m/s^2 depending on location data were collected)
% Columns 5:7 -- Angular velocity x, y, and z in rad/s

Note on data conventions
% All coordinate systems per ISB and as described in Wu & Cavanagh 1995 and our paper
% +x is anterior (direction of progression)
% +y is proximal (or vertical)
% +z is right

Method
% String specifying the method to be used (last name of first author)
% Shank methods
% % 'Mizrahi' (Mizrahi et al. 2000)
% % 'Mercer' (Mercer et al. 2003)
% % 'Purcell' (Purcell et al. 2006)
% % 'AminianODonovan' (Aminian et al. 2002; ODonovan et al. 2009)
% % 'AminianODonovan_modified' (Aminian et al. 2002; ODonovan et al. 2009; modified per our Supplemental Material)
% % 'GreeneMcGrath' (Greene et al. 2010; McGrath et al. 2012)
% % 'GreeneMcGrath_modified' (Greene et al. 2010; McGrath et al. 2012; modified per our Supplemental Material)
% % 'Sinclair' (Sincalir et al. 2013)
% % 'Whelan' (Whelan et al. 2015)
% % 'Norris' (Norris et al. 2015)
% % 'Schmidt' (Schmidt et al. 2016)
% % 'Aubol' (Aubol et al. 2020)
% % 'Fadillioglu' (Fadillioglu et al. 2020)
% % 'Bach' (Bach et al. 2022)
% % 'Bach_modified' (Bach et al. 2022; modified per our Supplemental Material)
% Low-back/sacrum methods
% % 'Auvinet' (Auvinet et al. 2002)
% % 'Lee' (Lee et al. 2010)
% % 'Wixted' (Wixted et al. 2010)
% % 'Bergamini' (Bergamini et al. 2012)
% % 'Benson' (Benson et al. 2019)
% % 'Reenalda' (Reenalda et al. 2021)

Location
% String specifying wearable location
% % 'Left shank'
% % 'Right shank'
% % For sacrum/low-back methods no location input is necessary

Supervision
% Only for Schmidt method
% % 'Supervised' (Requires UI to define thresholds for each trial)
% % 'Unsupervised' (Uses default values, no UI)

OUTPUTS

To the extent possible for any given method...

timings -- Table with columns containing:
% timings of initial contact 
% timing of terminal contact
% boolean for left stance (1 = left stance true, 0 = left stance false/right stance true)

stances -- Table with columns containing
% timestamps
% left stance true (1) or false (0)
% right stance true (1) or false (0)

segmented -- Structure with fields for each stance and swing (or, if not possible, for each stride or step)
% Labelled by side and occurence (from 1 ... n)
% Data are in original matrix form with timestamps, x, y, z acceleration and x, y, z angular velocity

DEPENDENCIES

This main function calls on sub-functions for each method
Some of these method sub-functions are themselves dependent on supporting functions and data
These sub- and supporting-functions are included with the download

List of supporting functions and data:
% MATLAB's Wavelet Toolbox -- published at https://www.mathworks.com/products/wavelet.html
% MATLAB's Signal Processing Toolbox -- published at https://www.mathworks.com/products/signal.html
% MATLAB's Statistics and Machine Learning Toolbox -- published at https://www.mathworks.com/products/statistics.html
% MATLAB's Deep Learning Toolbox -- published at https://www.mathworks.com/products/deep-learning.html
% scaletime -- published by Jan at https://www.mathworks.com/matlabcentral/fileexchange/25463-scaletime
% Bach_ESN.mat -- Echo State Network trained and validated on running data published by Bach et al. 2022 at https://github.com/marlow17/PredictingGroundReactionForces
% estimate_events.m -- published by Bach et al. 2022 at https://github.com/marlow17/PredictingGroundReactionForces
% compute_state_matrix.m -- found in ESNToolbox published by Manta Lukoševičius at https://www.ai.rug.nl/minds/uploads/ESNToolbox.zip
% test_ESN.m -- found in ESNToolbox published by Manta Lukoševičius at https://www.ai.rug.nl/minds/uploads/ESNToolbox.zip
% wavelet-mediated differentiation method published by Jianwen Luo at https://www.mathworks.com/matlabcentral/fileexchange/13948-numerical-differentiation-based-on-wavelet-transforms, specifically including:
% % derivative_dwt.m
% % spline_scale_function.m
% % spline_function.m
% REID_IMU_crash_catch.m -- original function to ensure ICs and TCs from each method are standardized and can be used to generate the desired output and that common crashes are avoided
