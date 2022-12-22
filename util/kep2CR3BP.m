function [r0,v0] = kep2CR3BP(a,e,I,w,Om,tp,t0,tf)
% INPUT: 
%       a  = semimajor axis (DU)
%       e  = eccentricity
%       i  = inclination (deg)
%       w  = argument of periapsis (deg)
%       Om = longitude of the asc./des node (deg)
%       tp = time of periapse passage (TU)
%
% OUTPUT:
%       r0 = [x,y,z] = position in CR3BP coords (DU) 
%       v0 = [vx,vy,vz] = velocity in CR3BP coords (VU)      

%% Useful constants and helper functions
mu_s = 2.9591220828559093E-04;   %sun GM in AU^3/day^2
mu_e_s = 2.9591309705483544E-04; %earth+sun GM, AU^3/day^2
mu_e = mu_e_s-mu_s;              %earth GM
R_e = 6378.137;                  %earth radius, km
kmAU = 149597870.700;            %1 AU in km
a_e =  1.000373836656026E+00 ;   %AU
e_e = 0;%1.712127710968187E-02;
I_e = 2.777040607882003E-03*pi/180; %rad
w_e = 3.043573249748720E+02*pi/180;
O_e = 1.596967974767415E+02*pi/180;
t_p_e = 0;%2458853.731945450883; %JD

rotMats = DCMs();
[~, ~, Mfun, ~, ~, rfunv, vfunv] =  KeplerFun();

%% Setup CR3BP
mu = mu_e/mu_e_s; %=3.0035e6
mu1 = 1-mu;
mu2 = mu;

% calculate mean motion
n_e = sqrt(mu_e_s/a_e^3); % 1/day

DU = a_e; %AU, canonical distance unit
TU = 1/n_e; %days, canonical time unit

% convert angular orbital elements to radians
I = I*pi/180; 
w = w*pi/180; 
Om = Om*pi/180; 

% calculate mean motion
n = sqrt(mu_e/a^3)*TU; % rad/TU

ts = linspace(t0,tf,5e4); %integration time

%% Notation definitions
% P = spacecraft
% 1 = sun
% 2 = earth
% O = sun/earth barycenter
% _G = geocentric inertial, ecliptic 
% _I = heliocentric inertial, ecliptic
% _B = barycentric rotating frame, ecliptic
% _P = earth-sun perifocal frame

% earth rotation matrices
PCI = rotMats{3}(w_e)*rotMats{1}(I_e)*rotMats{3}(O_e); % inertial -> perifocal
ICP = PCI.'; % perifocal -> inertial

% pos,vel of sun (m1) relative to sun-earth barycenter, perifocal
r1O_P = [-mu2*cos(ts);-mu2*sin(ts);zeros(size(ts))];
v1O_P = [ mu2*sin(ts);-mu2*cos(ts);zeros(size(ts))];

% pos,vel of earth (m2) relative to sun-earth barycenter, perifocal
r2O_P = [ mu1*cos(ts);mu1*sin(ts);zeros(size(ts))];
v2O_P = [-mu1*sin(ts);mu1*cos(ts);zeros(size(ts))];

% pos,vel of sun (m1) relative to sun-earth barycenter, ecliptic
r1O_I = ICP*r1O_P;
v1O_I = ICP*v1O_P;

% pos,vel of earth (m2) relative to sun-earth barycenter, ecliptic
r2O_I = ICP*r2O_P;
v2O_I = ICP*v2O_P;

% pos,vel of P relative to earth (m2)
rP2_G = rfunv(a,e,I,w,Om,invKepler(Mfun(n,ts,tp),e)); % DU
vP2_G = vfunv(a,e,I,w,Om,invKepler(Mfun(n,ts,tp),e),n); % VU

% pos,vel of P relative to sun-earth barycenter
rPO_I = rP2_G + r2O_I;
vPO_I = vP2_G + v2O_I;

% pos,vel of P relative to sun (m1)
rP1_I = rPO_I - r1O_I; %DU
vP1_I = vPO_I - v1O_I; %VU

%% Determine rotating CR3BP equivalent state
%Find true anomaly at initial time and use as angle offset
nu0 = invKepler(mod(n_e*(ts(1)*TU - t_p_e),2*pi),e_e);
th0 =  w_e + nu0;
%compute rotating frame components
rPO_B = zeros(size(rP1_I));
vPO_B = zeros(size(vP1_I));
for j = 1:length(ts)
    BCI = rotMats{3}((ts(j)-ts(1)) + th0)*rotMats{1}(I_e)*rotMats{3}(O_e);
    rPO_B(:,j) = BCI*rPO_I(:,j);
    vPO_B(:,j) = BCI*vPO_I(:,j) - cross([0;0;1], rPO_B(:,j));
end

%trajectory
x = rPO_B(1,:);  y = rPO_B(2,:);  z = rPO_B(3,:);
xd = vPO_B(1,:); yd = vPO_B(2,:); zd = vPO_B(3,:);

t = 1; % pull the initial condition (t=1)
r0 = [x(t) y(t) z(t)].';
v0 = [xd(t) yd(t) zd(t)].';

end