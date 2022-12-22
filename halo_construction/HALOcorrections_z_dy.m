function dz_correction = HALOcorrections_z_dy(Y_T2, STM_T2, mu) 
% INPUT:
%       Y_T2 = state vector at t = T/2 (y=0 crossing) (6x1)
%     STM_T2 = state transition matrix at t = T/2 (6x6)
%         mu = mass parameter
% OUTPUT:
%       dz_correction = correction to z and dy to trend towards
%                       periodic orbit (2x1)
%

% Pull out the state vector components.
 x = Y_T2(1);  y = Y_T2(2);  z = Y_T2(3);
dx = Y_T2(4); dy = Y_T2(5); dz = Y_T2(6);
% Calculate distances for use in later calculations.
r1=sqrt((mu+x)^2+(y)^2+(z)^2); r2=sqrt((x-1+mu)^2+(y)^2+(z)^2);
% Calculate the x and z accelerations (from CRTBP equations)
xdd = 2*dy + x - (1-mu)*(x+mu)/r1^3 - mu*(x+mu-1)/r2^3; 
zdd = -z*((1-mu)/r1^3 + mu/r2^3);
% Setting up the expression to calculate correction values (this comes
% from a paper on 3D HALO Orbit families by Kathleen Howell).
A = ([STM_T2(4,3), STM_T2(4,5); STM_T2(6,3), STM_T2(6,5)]- ...
        (1/dy)*[xdd; zdd]*[STM_T2(2,3),STM_T2(2,5)]);
B = [-dx; -dz];
% Calculates the corrections and sets them to return. 
if abs(det(A)) <= eps(1)
    dz_correction = [0;0];
else
    dz_correction = A\B; % Return del_z and del_dy
end
end
