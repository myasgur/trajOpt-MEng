function y = Runge_Kutta_Merson(func,t,y0,dt, mu, n)
% INPUTS:
%           func = function handle to be integrated
%              t = timespan [kx1];
%             y0 = initial conditions [1xm]
%             dt = time step (TU)
%             mu = mass parameter of 3-body system
%              n = mean motion of 3-body system (1/TU)
% OUTPUT:
%              y = integrated state [kxm]

% Implementation of the Runge-Kutta integration approach
eta0 = y0;
k0   = dt*func(t,eta0, mu)*n;
eta1 = eta0 + k0/3;
k1   = dt*func(t+dt/3,eta1, mu)*n;
eta2 = eta0 + (k0+k1)/6;
k2   = dt*func(t+dt/3,eta2, mu)*n;
eta3 = eta0 + (k0+3*k2)/8;
k3   = dt*func(t+dt/2,eta3, mu)*n;
eta4 = eta0 + (k0-3*k2+4*k3)/2;
k4   = dt*func(t+dt,eta4, mu)*n;

% Final return value
y = eta0 + (k0+4*k3+k4)/6;
end