%% REID_IMU_Schmidt

% Schmidt et al. 2016
% Uses angular velocity about the M-L axis and acceleration in the P-D axis to
% identify IC and TC events during sprinting
% Was developed and validated specifically for sprinting and is not intended
% for a broader range of running speeds
% Code is adapted with permission from R code provided by Marcus Schmidt 22/09/23
% Other than porting to MATLAB and ensuring that it works with same
% inputs/gives same outputs as other methods in this package changes to
% the code were minimal
% The original code has user-defined thresholds and windows
% To facilitate this a UI interface was added 
% However, given the goal of this project was to evaluate unsupervised
% methods of gait event ID, an option has also been added to run the code
% without supervision

function [timings, stances, segmented] = REID_IMU_Schmidt(data, location, Fs, supervision, min_stance_t)

% Loop to re-run function until user is satisfied with results (if supervising)
% If not supervising this loop will be exited on it's first iteration
accepted = 0;
while accepted == 0
    % Code relies on 3 individually scaleable values...
    switch supervision
        % In the original code these values were defined by the user
        % Prompt user to define values...
        case 'Supervised'
            prompt = {'Critical acceleration threshold (g):','Post-IC deadtime (ms):','Min acceleration window (ms):'};
            dlgtitle = 'Individually scalable parameters';
            dims = [1 75];
            definput = {'5','90','150'};
            UI_inputs = inputdlg(prompt,dlgtitle,dims,definput);
            crit_a_thresh = str2num(UI_inputs{1});
            post_IC_dead = str2num(UI_inputs{2});
            min_a_window = str2num(UI_inputs{3});
            % For our purposes, however, we wanted an unsupervised method to identify gait events
            % So, we also have the option to simply use default values...
        case 'Unsupervised'
            crit_a_thresh = 5;
            post_IC_dead = 90;
            min_a_window = 150;
    end % switch supervision
    % Convert ms to frames
    post_IC_dead = round(Fs/(1000/post_IC_dead));
    min_a_window = round(Fs/(1000/min_a_window));
    % Find minimum in the angular velocity about the ML axis
    Gz_min = min(data(:,7));
    % Look through frames
    help_pos_ticker = 1;
    for frame_count = 1:size(data,1)-100
        % Find index and value for local ML angular velocity minima 100 ms following current frame
        [local_Gz_min local_Gz_min_pos] = min(data(frame_count:frame_count+round(Fs/10),7));
        % Find differential of angular velocity for 50 ms following current frame
        descent = diff(data(frame_count:frame_count+round(Fs/20),7));
        % If the local angular velocity is at a global minima OR all differentials in angular velocity are negative and the acceleration exceeds the critical threshold
        if local_Gz_min == Gz_min ||...
                all(descent < 0) && data(frame_count,3) >= crit_a_thresh
            % Flag the current frame to use to idenitfy IC
            help_pos_start(help_pos_ticker) = local_Gz_min_pos + frame_count - 1;
        else
            % Otherwise put dummy data in
            help_pos_start(help_pos_ticker) = NaN;
        end % if
        % Increment up the ticker
        help_pos_ticker = help_pos_ticker + 1;
    end % for frame_count
    % Remove NaNs
    help_pos_start(isnan(help_pos_start)) = [];
    % Retain only unique values
    % This should now contain an array of indices where the IC may occur
    points_of_interest = unique(help_pos_start);
    % Remove points of interest too close to trial start or end to avoid crashes
    points_of_interest(points_of_interest > size(data,1)-round(Fs/10)) = [];
    points_of_interest(points_of_interest < round(Fs/50)) = [];
    % Should now have points_of_interest populated with continuous indices
    % around which IC and TC events are occuring, i.e., each continuous series
    % of indices is associated with one step. Discrete jumps between continuous
    % series of indices therefore represent a jump in time between steps
    % Thus, we can know how many steps are in our dataset based on how many
    % continuous series of indices we have
    is_not_continuous = [true,diff(points_of_interest) ~= 1];
    continuous_indices = points_of_interest(~is_not_continuous);
    continuous_sections = transpose([0; transpose(cumsum(diff(continuous_indices) > 1))]);
    continuous_sections = continuous_sections + 1;
    num_steps = max(continuous_sections);
    % Look through all points of interest
    contact_start_ticker = 1;
    for poi_count = points_of_interest
        % Find vertical acceleration max near point of interest
        [mag_1 local_max_point_Ay] = max(data(poi_count-round(Fs/50):poi_count+round(Fs/10),3));
        % Find vertical acceleration min near point of interest but before max
        [mag_2 local_min_point_Ay] = min(data(poi_count-round(Fs/50):poi_count-round(Fs/50)+local_max_point_Ay-1,3));
        % Adjust timings back into trial frame number
        local_min_point_Ay = local_min_point_Ay + poi_count - round(Fs/50) - 1;
        % If the local min we found occurs after the current point of interest
        if local_min_point_Ay > poi_count
            % Call the timing of that local min the contact start
            contact_start(contact_start_ticker) = local_min_point_Ay;
        else
            % Otherwise, call the current point of interest the contact start
            contact_start(contact_start_ticker) = poi_count;
        end % if
        % Increment up the ticker
        contact_start_ticker = contact_start_ticker + 1;
    end % for poi_count
    % Look through all contact starts
    if exist('contact_start','var') && ~isempty(contact_start)
        contact_end_ticker = 1;
        for contact_start_count = 1:size(contact_start,2)
            % Find the min in angular velocity about the ML axis near the contact start, flag it to help identify TC
            % Make sure the window is within the existing data
            if size(data,1) - (contact_start(contact_start_count)+post_IC_dead+min_a_window) > 0
                [mag_3 help_pos_end] = min(data(contact_start(contact_start_count)+post_IC_dead:contact_start(contact_start_count)+post_IC_dead+min_a_window,7));
                % Otherwise, if the start of the window exists in the data but the end doesn't,
                % Just look until the end of the data
            elseif size(data,1) - (contact_start(contact_start_count)+post_IC_dead) > 0
                [mag_3 help_pos_end] = min(data(contact_start(contact_start_count)+post_IC_dead:end,7));
                % Otherwise, NaN flag the error
            else
                help_pos_end = NaN;
            end
            % Adjust timings back into trial frame number
            help_pos_end = help_pos_end + contact_start(contact_start_count) + post_IC_dead - 1;
            % Find the min in vertical acceleration near the flagged frame
            % If the window is within the existing data
            if ~isnan(help_pos_end) && size(data,1) - (help_pos_end+round(Fs/20)) > 0
                [mag_4 contact_end(contact_end_ticker)] = min(data(help_pos_end:help_pos_end+round(Fs/20),3));
                % Otherwise, just look to the end of the data
            elseif ~isnan(help_pos_end)
                [mag_4 contact_end_temp] = min(data(help_pos_end:end,3));
                % If a temporary contact end was found
                if exist('contact_end_temp','var') && ~isempty(contact_end_temp)
                    contact_end(contact_end_ticker) = contact_end_temp;
                    clear contact_end_temp
                else
                    contact_end(contact_end_ticker) = NaN;
                end
            else
                contact_end(contact_end_ticker) = NaN;
            end
            % Adjust timings back into trial frame number
            if ~isempty(help_pos_end)
                contact_end(contact_end_ticker) = contact_end(contact_end_ticker) + help_pos_end - 1;
            else
                contact_end(contact_end_ticker) = NaN;
            end
            % Increment up the ticker
            contact_end_ticker = contact_end_ticker + 1;
        end % for contact_start_count
    end
    % Should now have two arrays: contact_start and contact_end
    % Each of these is populated with indices where IC and TC are occurring
    % This gives a "stair" appearance where each stair jumps to the estimated
    % contact_start or contact_end timing for a given step
    % HOWEVER, sometimes the code will find multiple estimated timings for a
    % given step. To solve this I have added code to quantify the number of
    % steps based on continuous points of interest (corresponding to a step)
    % Based on this, we can loop through the estimated number of steps and take
    % the mean estimated timing for the continuous section associated with that
    % step
    for step_count = 1:num_steps
        % Get the indices from points_of_interest that correspond to the frames of the current step
        % Then take the estimated timings in contact_start and _end that correspond to those frames
        if ~isempty(continuous_indices) && any(continuous_sections == step_count)
            indices_of_interest = continuous_indices(continuous_sections == step_count);
            IC(step_count,1) = round(mean(contact_start(ismember(points_of_interest, indices_of_interest))));
            TC(step_count,1) = round(mean(contact_end(ismember(points_of_interest, indices_of_interest))));
        else
            IC(step_count,1) = NaN;
            TC(step_count,1) = NaN;
        end
    end
    % Option to display output and adjust UI-defined parameters if supervising
    switch supervision
        case 'Supervised'
            figure('units','normalized','outerposition',[0 0 1 1],'Color', [1 1 1]);
            plot(data(:,3))
            ylabel('acceleration')
            hold on
            fig_prop = gca;
            line(fig_prop.XLim,[crit_a_thresh crit_a_thresh],'LineStyle','--','Color','k')
            yyaxis right
            plot(data(:,7))
            set(gca,'YColor','k')
            ylabel('angular velocity')
            yyaxis left
            scatter(IC, crit_a_thresh-0.1*crit_a_thresh,'v','MarkerEdgeColor','k','MarkerFaceColor','k')
            scatter(TC, crit_a_thresh-0.1*crit_a_thresh,'^','MarkerEdgeColor','k','MarkerFaceColor','k')
            scatter(IC+post_IC_dead, crit_a_thresh+0.1*crit_a_thresh,'o','MarkerEdgeColor','k')
            scatter(IC+post_IC_dead+min_a_window, crit_a_thresh+0.1*crit_a_thresh,'diamond','MarkerEdgeColor','k')
            title('blue--acceleration; orange--angular velocity; downward triangle--IC; upward triangle--TC; dashed line--critical acceleration threshold; circle--post IC deadtime; diamond--minimum acceleration window')
            accept_reject = questdlg('Accept data with current thresholds?', ...
                '', ...
                'Yes','No','Yes');
            if strcmp(accept_reject,'Yes') == 1
                accepted = 1;
                close
            else
                accepted = 0;
                clear IC TC
                close
            end
        case 'Unsupervised'
            accepted = 1;
    end
end % accepted while loop
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
% Once results are accepted (or if running unsupervised) write the data
% Create tables and structures
timings = table;
stances = table;
segmented = struct;
% Correct timings back to original timestamps
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