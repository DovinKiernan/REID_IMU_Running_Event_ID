function Yi = scale_time(Y, Ti, Tf, Tn)
% Linear interpolation of a matrix over first dimension
% The matrix Y is assumed to be a time-series of data with the first index is
% the time in equidistant steps. These data are interpolated at the new time
% steps Ti. Ti must not be smaller than 1 or greater than number of rows of Y.
% This is a fat-free version of:
%   interp1(1:size(Y, 1), Y, Ti, '*linear') with Y is a real double matrix, no
%   extrapolation, no handling of NaN values.
%
% Method 1:  Yi = ScaleTime(Y, Ti)
%   Y:  Real double matrix, which is interpolated over the 1st dimension.
%       Y must have at least 2 elements.
%       As usual in Matlab a row vector is accepted also.
%   Ti: Time steps for interpolation as vector of doubles. If an element of Ti
%       has no integer value, the corresponding Yi is:
%         Yi = Y(floor(Ti)) * (1 - Ti + floor(Ti)) + ...
%              Y(ceil(Ti))  * (Ti - floor(Ti))
%       Ti need not be monotonically increasing, but just the first and last
%       frame are checked to be inside the limits.
%       Tip: N equidistant time steps between iFrame and nFrame:
%         iFrame:((fFrame - iFrame) / (N - 1)):fFrame
%   Yi: Interpolated values as [length(Ti) x size(Y, 2)] double array.
%
% Method 2:  Yi = ScaleTime(Y, initialT, finalT, numberT)
%   Y:  See above.
%   initialT, finalT, numberT: 3 scalar numbers for equidistant equally spaced
%       interpolation steps. This can be much faster than creating the vector Ti
%       explicitely.
%   Yi: Interpolated values as [numberT x size(Y, 2)] double array.
%
% EXAMPLE:
% Interpolate SIN and COS:
%   T  = linspace(0, 2*pi, 20);  Ti = linspace(1, length(T), length(T) * 3);
%   Y  = [sin(T)', cos(T)'];
%   Yi = ScaleTime(Y, Ti);
%   plot(T, Y, 'ob'); hold('on'); plot(ScaleTime(T, Ti), Yi, 'xr');
%
% SPEED:
% - Usually ScaleTime.m is 2 to 3 times faster than Matlab's INTERP1('*linear'),
%   but for large Y (e.g. [100000 x 10]) and short Ti (e.g. 1:10) the speed gain
%   can reach the factor 30.
% - The MEX is 10 to 125 times faster than INTERP1. The index method (4 inputs)
%   is faster than the vector method (2 input), if Y has less than 100 columns.
%   E.g. Y := [1000 x 1], Ti := [10000 x 1], numbers are iterations/sec ov a
%   1.5GHz Pentium-M, WinXP:
%                            Matlab6/BCC5.5   Matlab7.8/LCC3.8
%     ScaleTime.mex(index):       4900              3715
%     ScaleTime.mex(vector):       550              1150
%     ScaleTime.m:                 230               383
%     INTERP1:                     100               174
%     lininterp1f:                  29                30
%     qinterp1:                    104               174
% - Look at the speed table created by TestScaleTime.
%
% COMPILE:
%   mex -O ScaleTime.c
%
% Tested: Matlab 6.5, 7.7, 7.8
% Author: Jan Simon, Heidelberg, (C) 2009 J@n-Simon.De
%
% Inspired by INTERP1 of Matlab5.3 written by:
%   Clay M. Thompson 7-4-91
%
% See also INTERP1.

% $JRev: R0j V:035 Sum:30BD430B Date:30-Sep-2009 01:00:47 $
% $File: User\JSim\Published\ScaleTime\ScaleTime.m $

% ------------------------------------------------------------------------------
% Get input as matrix or column vector:
[nRow, nCol] = size(Y);
doTranspose  = (nRow == 1 && nCol > 1);
if doTranspose
   nRow = nCol;
   nCol = 1;
   Y    = Y(:);
end

% Get interpolation steps:
switch nargin
   case 4
      % Check for out of range values of Ti:
      if and(Ti < 1, Tf > nRow)
         error(['*** ', mfilename, ': Interpolation frames out of range.']);
      end
      
      if Tn < 1.0
         Yi = [];
         return;
      elseif Tn > 1
         Ti = Ti:((Tf - Ti) / (Tn - 1)):Tf;
      elseif Tf < Ti
         Ti = [];
      end
      
   case 2
      % Check for out of range values of Ti:
      if any(Ti < 1)  % Faster to do 2 tests then using OR - surprising
         error(['*** ', mfilename, ': All [Ti] must be >= 1.']);
      end
      if any(Ti > nRow)
         error(['*** ', mfilename, ...
               ': All [Ti] must be <= number of rows of Y.']);
      end
      
   otherwise
      error(['*** ', mfilename, ': 2 or 4 inputs required.']);
end

% Return the empty matrix, if no interpolation steps are wanted:
if isempty(Ti)
   Yi = [];
   return;
end

% The Matlab method needs at least 2 original time points. Althouh the MEX
% implementation could handle a scalar input also, it is not worth to consider
% this pathological case:
if nRow < 2
   error(['*** ', mfilename, ': Y must have at least 2 rows.']);
end

% Shape Ti:
Ti = Ti(:);

% Interpolation parameters:
Si = Ti - floor(Ti);
Ti = floor(Ti);

% Shift frames on boundary:
d     = (Ti == nRow);
Ti(d) = Ti(d) - 1;
Si(d) = 1;           % Was: Si(d) + 1;

% Now interpolate:
if nCol > 1
   Si = Si(:, ones(1, nCol));  % Expand Si
end
Yi = Y(Ti, :) .* (1 - Si) + Y(Ti + 1, :) .* Si;

%Slower:
%Yi = Y(Ti, :) + (Y(Ti + 1, :) - Y(Ti, :)) .* Si(:, ones(1, nCol));

if doTranspose
   Yi = reshape(Yi, 1, []);
end

return;
