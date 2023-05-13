function [deltaV, Y_prop] = StationKeep(y0, Io, propOut, t_0, h, j, mu, n, STMsave, Y_xyz, b ) 
    ti = t_0;
    yman = y0.';
    for i = 1:1:propOut
        Y_prop(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y0; Io], h, mu, n)'; 
        y0 = Y_prop(i,1:6)'; Io = Y_prop(i, 7:42)';
        t_0 = t_0 + propOut*h;
    end

    M2 = j;
    T2 = j+ propOut;
    deltaV = zeros(1,6);
    deltaV(4:6) = zeros(1,3);
    if j+propOut < size(Y_xyz,1)
        STM_mod = STMsave(:,:,T2)/(STMsave(:,:,M2)); % M = manuevor, T = target
        A = STM_mod(1:3,1:3); B = STM_mod(1:3,4:6);
        deltaV(4:6) = -inv(B)*A*(yman(1:3)'-Y_xyz(b,1:3)')-(yman(4:6)'-Y_xyz(b,4:6)');
    end 
end