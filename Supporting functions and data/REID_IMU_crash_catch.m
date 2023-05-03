%% Function to correct common errors from the REID_IMU functions and ensure common output

function [IC,TC,left_side] = REID_IMU_crash_catch(min_stance_frames,IC,TC,left_side)

% If the IC only has a NaN-flag just skip the below
if all(isnan(IC))
    IC = NaN;
    % NaN-flag TCs too
    if nargin > 2
        TC = NaN;
    end
    % NaN-flag side too
    if nargin > 3
        left_side = NaN;
    end
else
    % Each method outputs at least IC
    % We will delete any NaN-flags
    % This will possibly delete associated TC or side values bute take this 
    % conservative approach to ensure that we don't get stances that are 
    % actually composed of multiple stance and swing phases
    if nargin > 3
        left_side(isnan(IC)) = [];
    end
    IC(isnan(IC)) = [];
    % We also don't want any repeat values
    % This will also sort into temporal order
    IC = unique(IC);
    % If a TC is also output, then ensure that each TC is matched to the appropriate IC
    if nargin > 2
        % If the TC only has NaN-flags just skip the below
        if all(isnan(TC))
            TC = NaN;
            % Very conservative -- if a method is supposed to deliver paired IC-TC but doesn't deliver TC, then don't accept its IC or side
            IC = NaN;
            if nargin > 3
                left_side = NaN;
            end
        else
            % The match should be the closest TC in time AFTER the IC
            % If a method IDed multiple consecutive TCs without any ICs following them
            % this may cause incorrect stance times and TCs
            % (i.e., if a method prematurely IDed TC during midstance then
            % IDed another TC closer to the true value)
            % Or, if a method IDed the same frame as the IC and TC
            % Thus, we only accept TC-IC > min_stance_t
            TC_rep = repmat(TC,[1 length(IC)]); % TC values in each row repeated n of ICs times in each column
            TC_IC_diff = TC_rep-IC'; % Subtract IC1...n from each TC (i.e., 1,1 = TC1-IC1; 1,2 = TC1-IC2; 2,1 = TC2-IC1)
            TC_IC_diff(TC_IC_diff < 0.5*min_stance_frames) = NaN; % NaN-flag values below the minimum stance time -- lowered to 50% min_stance_t to allow for methods that could systematically shorten stance times
            [~, closest_ind] = min(TC_IC_diff,[],1); % Find the TC that minimizes stance t for each IC (i.e., columns correspond to each IC and values correspond to the TC that best matches it)
            TC_closest = TC_rep(closest_ind)';
            % If there's an unmatched IC AFTER the last TC
            % Then we will get a column of all NaN values
            % Use that to NaN-flag the TC as non-existent
            % and remove it
            TC_missing = sum(isnan(TC_IC_diff),1)';
            TC_closest(TC_missing == size(TC_IC_diff,1)) = NaN;
            % It is possible for multiple ICs to share the same TC
            % This occurs if the method IDs consecutive ICs with no TC following them
            % In this case, we will accept the *last* IC that shares the same TC (minimizing stance t)
            % Then delete all other ICs and sides (possibly deleting a true IC where the method failed to find a TC)
            [TC TC_unique_ind] = unique(TC_closest,'last');
            IC = IC(TC_unique_ind);
            if nargin > 3
                left_side = left_side(TC_unique_ind);
            end
            % We also want to ensure that only the last TC is a NaN
            if any(isnan(TC(1:end-1)))
                IC(isnan(TC(1:end-1))) = [];
                if nargin > 3
                    left_side(isnan(TC(1:end-1))) = [];
                end
                TC(isnan(TC(1:end-1))) = [];
            end
        end % only NaN-flags in the TC
    end % if TC is an input
end % only NaN-flags in the IC
% Last check to ensure not empty
if isempty(IC)
    IC = NaN;
end
if nargin > 2 && isempty(TC)
    TC = NaN;
end
if nargin > 3 && isempty(left_side)
    left_side = NaN;
end

end % function