function GMatrix=G_CRTBP(x, mu)
% INPUT:
%       x  = position vector in CRTBP coordiantes (nx1)
%       mu = mass parameter
% OUTPUT:
%       GMatrix = Hessian of scalar pseudopotential U(x,y,z), (nxn)
%
% Credit to Ethan Geipel (2019)
%

% Distances
r1 = sqrt((x(1)+mu)^2 + x(2)^2 + x(3)^2);
r2 = sqrt((x(1)-(1-mu))^2 + x(2)^2 + x(3)^2);

% The gradient has three components which we will call u1, u2, u3. 
% The differential of the  gradient is the matrix of partials of these 
% functions. These will be denoted by u1_x, u1_y, and so forth.
u1_x = 1 - (1-mu)*(1/(r1^3) - 3*((x(1)+mu)^2)/(r1^5))- mu*(1/(r2^3) - 3*((x(1)-(1- mu))^2)/(r2^5));
u2_y = 1 - (1-mu)*(1/(r1)^3 - 3*x(2)^2/r1^5) - mu*(1/r2^3 - 3*x(2)^2/r2^5);
u3_z=(-1)*(1-mu)*(1/(r1)^3-3*x(3)^2/r1^5)-mu*(1/r2^3-3*x(3)^2/r2^5);
u1_y=3*(1-mu)*x(2)*(x(1)+mu)/r1^5+3*mu*x(2)*(x(1)-(1-mu))/r2^5;
u1_z=3*(1-mu)*x(3)*(x(1)+mu)/r1^5+3*mu*x(3)*(x(1)-(1-mu))/r2^5;
u2_z=3*(1-mu)*x(2)*x(3)/r1^5+3*mu*x(2)*x(3)/r2^5;

% Equality of mixed partials gives (as all the terms are already partials 
% of the potential function);
u3_y=u2_z;
u2_x=u1_y;
u3_x=u1_z;
% Then G is the matrix of partials
GMatrix=[u1_x, u1_y, u1_z;
         u2_x, u2_y, u2_z;
         u3_x, u3_y, u3_z];
end