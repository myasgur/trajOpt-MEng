function [T2, y_STM_T2, STM_T2] = findT2(y0, t_0, t_end, h, Io, mu, n) 
% INPUTS:
%         y0 = state vector in Cr3BP coords [6x1]
%        t_0 = start time (normalized) [1x1]
%      t_end = end time (normalized) [1x1]
%          h = time step
%         Io = initial eye() matrix for STM integration
%
% OUTPUTS:
%         T2 = time at half period (t=T/2)
%   y_STM_T2 = state at t = T/2
%     STM_T2 = STM at t = T/2
%
% FUNCTION DESCRIPTION:
%   This function then iterates through using Runge_Kutta method and is
%   "broken" when the y component of the state vector crosses 0 on the
%   "other side" of the HALO orbit.
%   It iterates through this process until a tolerance of delta y < 1e-11 
%   is achieved after the y-crossing.
%% 
counter = 1; maxCount = 5e4; y_it = zeros(1,42);
% Iterates through time 0 -> t_end
while t_0 < t_end && counter < maxCount
    % Sends t_0, and y0; Io to the
    % passed along to the @stateTransitionMatrix function to be solved.
    y_it(counter,:) = Runge_Kutta_Merson(@stateTransitionMatrix,...
                                            t_0, [y0; Io], h, mu, n)';
    % Pulls out and updates Y and sets it to y0 for the next loop.
    y0 = y_it(counter,1:6)'; t_0 = t_0 + h; % Steps forward in time by h. 
    % Pulls out and updates Io (becomes STM)
    Io = y_it(counter, 7:42)';

    % Checks for y-crossing at y=0.
    if (counter > 1)
        if (y0(2) * y_it(counter-1, 2) < 0)
            break
        end % y0(2) <0
    end

    % Iterates through the loop.
    counter = counter + 1;
end
% After y crossing has been detected, T2 is set to be the time t_0 at 
% which the y-crossing occurs.
T2 = t_0 - h;
KeepGoing = 0;
% While loop to iterate T2 down such that the delta y value of the 
% y-crossing is < 1e-11.
    while KeepGoing < 1 && counter < maxCount 
        % Reverts t_0 back 20 time steps (so that you don't have to iterate 
        % through all of the time span every iteration).
        t_0 = T2 - 20*h;
        % Reduces h by factor of 10 to achieve better resolution on T2.
        h = h/10;
        y0 = y_it(counter-20, 1:6)'; % Sets back the value of y0 wrt t_0 set 
        Io = y_it(counter-20, 7:42)';% Sets back the value of Io wrt t_0 set 
        counter = counter-20; y_it = zeros(1,42);
        while t_0 < T2 + .2 
            % Similar to above, sends time and state to R_G and executes
            % the function @sTM.
            y_it(counter,:) = Runge_Kutta_Merson(@stateTransitionMatrix,...
                                                    t_0, [y0; Io], h, mu, n)';
            % Pulls out state vector.
            y0 = y_it(counter,1:6)'; t_0 = t_0 + h; % Step forward in time, h 
            % Pulls out Io for STM calculation
            Io = y_it(counter,7:42)';
            % Gets the full [state vector; STM] vector stored for later
            % calculations.
            ySTM = y_it(counter,:);
            
            % Detects y-crossing at y = 0.
              %if (y0(2) > 0)
            if (y0(2) * y_it(counter-1, 2) < 0)
            
                % Updates the value for T2 with the time of y-crossing
                % detection. (Will be more accurate due to h decrease.
                T2 = t_0;
                KeepGoing = KeepGoing + .1;
                % Calls function @xSTM(), sending the [state vector; STM]
                % at time T2, which sends back the 6x6 STM (at time T2). 
                [y_STM_T2, STM_T2] = xSTM(ySTM(end,:));
                % Checks if the tolerance of delta y at y-crossing has been 
                % met. Tolerance is 1e-11.
                if (abs(y0(2)) < 1e-15)
                    KeepGoing = 1;
                end
                break
            end % y0(2) <0
            counter = counter + 1;
        end % while t_0 < T2 + .2
    end % while keepGoing < 1
    if counter >= maxCount
        T2 = t_0;
        [y_STM_T2, STM_T2] = xSTM(ySTM(end,:));
        fprintf('After %d iterations, the error of delta-y at y-crossing = %d\n',counter,y0(2))
    end
end
