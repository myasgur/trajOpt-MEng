function [Yi, STM] = xSTM(Y)
% INPUT:
%           Y = vector of state and STM (42x1)
%
% OUTPUTS: 
%          Yi = state vector (6x1)
%         STM = state transition matrix (6x6)
%
    Yi = Y(1:6);
    STM(:,1) = Y(7:12);
    STM(:,2) = Y(13:18);
    STM(:,3) = Y(19:24);
    STM(:,4) = Y(25:30);
    STM(:,5) = Y(31:36);
    STM(:,6) = Y(37:42);
end