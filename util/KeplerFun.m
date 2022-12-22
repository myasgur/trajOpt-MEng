function [Afun, Bfun, Mfun, rfun, vfun, rfunv, vfunv] =  KeplerFun()
 %Keplerian Orbit Propagation Helper Functions
 %Returns a set of anonymous functions for use in Keplerian Orbit
 %propagation
 %
 %INPUTS
 %   None
 %
 %OUTPUTS
 %   Afun, Bfun  Thiel-Innes elements for inputs a,e,I,omega,Omega
 %   MFun        Mean anomaly for inputs of n, t, t_p
 %   rfun        Position vector (3x1) for inputs a,e,I,omega,Omega
 %   vfun        Velocity vector (3x1) for inputs A,B,E,e,n
 %   rfunv       Position vectors (3xn) for inputs a,e,I,omega,Omega,E.  All
 %               inputs are assumed scalar except for E, which is assumed to
 %               be nx1. 
 % 
 %NOTES
 %   All functions operate on a standard set of Keplerian elements:
 %   a       Semi-major axis (DU - distance unit)
 %   e       Eccentricity
 %   I       Inclination (rad)
 %   omega   Argument of periapsis (rad)
 %   Omega   Longitude of the ascending node (rad)
 %   n       Mean motion (rad/TU - time unit)
 %   t       Time (TU)
 %   t_p     Time of periapsis passage (TU)
 %   E       Eccentric anomaly
 
 % Thiel-Innes elements
 Afun = @(a,e,I,omega,Omega) ...
     [a.*(cos(Omega).*cos(omega) - sin(Omega).*cos(I).*sin(omega));...
      a.*(sin(Omega).*cos(omega) + cos(Omega).*cos(I).*sin(omega));...
      a.*sin(I).*sin(omega)];
 
 Bfun = @(a,e,I,omega,Omega) ...
     [-a.*sqrt(1-e.^2).*(cos(Omega).*sin(omega) + ...
                sin(Omega).*cos(I).*cos(omega));...
      a.*sqrt(1-e.^2).*(-sin(Omega).*sin(omega) + ...
                cos(Omega).*cos(I).*cos(omega));...
      a.*sqrt(1-e.^2).*sin(I).*cos(omega)];
 
 % Mean anomaly 
 Mfun = @(n,t,t_p) mod(n*(t - t_p),2*pi);
 
 %position and velocity vectors (3x1)
 rfun = @(A,B,E,e) A*(cos(E) - e) + B*sin(E);
 vfun = @(A,B,E,e,n) (n/(1 - e*cos(E)))*(-A*sin(E) + B*cos(E));
 
 %positions for multiple time steps
 rfunv = @(a,e,I,omega,Omega,E) Afun(a,e,I,omega,Omega)*(cos(E) - e).' + ...
     Bfun(a,e,I,omega,Omega)*(sin(E)).';

 %velocities for multiple time steps
 vfunv = @(a,e,I,omega,Omega,E,n) ((-n.*sin(E)./(1 - e.*cos(E))).*Afun(a,e,I,omega,Omega).' ...
     + (n.*cos(E)./(1 - e.*cos(E))).*Bfun(a,e,I,omega,Omega).').';
 end