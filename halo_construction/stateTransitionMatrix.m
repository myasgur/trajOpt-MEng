function Xreturn = stateTransitionMatrix(~, X, mu) 
% INPUT:
%           X = state and STM at time t-1 (42x1)
%          mu = mass parameter for 3-body system
% OUTPUT:
%          Xreturn = state and STM at time t (42x1)
%

Xreturn = zeros(42,1);
% Obtain r and v state vectors from initial expression.
r = X(1:3); v = X(4:6);
% Sets up Phi by reconstructing the 36x1 matrix version
PHIo(:,1) = X(7:12);
PHIo(:,2) = X(13:18);
PHIo(:,3) = X(19:24);
PHIo(:,4) = X(25:30);
PHIo(:,5) = X(31:36);
PHIo(:,6) = X(37:42);
% Used in constructing the G matrix.
Omega = [0,1,0;-1,0,0;0,0,0];
% Calls a function (obtained from James Mireles at FAU)
% the 3x3 matrix with the derivatives necessary to construct the G 
% matrix.
dadr = G_CRTBP([r; v], mu);
% Constructs the G matrix from above terms and generic code.
G = [zeros(3) eye(3); dadr 2*Omega];
% Multiplies G and Phi to iterate through the Phi process.
PHI = G*PHIo;
% Reconstructs the return matrix to 36x1 rather than 6x6. 
Xreturn(7:12) = PHI(:,1);
Xreturn(13:18) = PHI(:,2);
Xreturn(19:24) = PHI(:,3);
Xreturn(25:30) = PHI(:,4);
Xreturn(31:36) = PHI(:,5);
Xreturn(37:42) = PHI(:,6);
% Variables from matrix
 x = X(1);  y = X(2);  z = X(3);
dx = X(4); dy = X(5); dz = X(6);
% Distances, used in iterating through using CRTBP method.
r1 = sqrt((mu+x)^2 + (y)^2 + (z)^2);
r2 = sqrt((x-1+mu)^2 + (y)^2 + (z)^2);
% Implementation of the CRTBP state vector integration process. 
Xreturn(1:6) = [dx;
                dy;
                dz; 
                2*dy + x-(1-mu)*(x+mu)/r1^3 - mu*(x-(1-mu))/r2^3; 
                -2*dx + y - (1-mu)*y/r1^3 - mu*y/r2^3; 
                -(1-mu)*z/r1^3 - mu*z/r2^3];
end
