global MU DU TU VU MscU FU mu
global T_max c epsilon rho R_e
global t0 tf r0 v0 m0 rf vf

%% Useful constants
% mu_s = 2.959122082322128e-04;    % sun GM in AU^3/day^2
% mu_e = 8.887692446706600e-10;    % earth GM
% mu_e_s = mu_s+mu_e;              % earth+sun GM, AU^3/day^2

mu_e = 398600.436;               % earth GM in km^3/s^2
mu_m = 4902.800066;              % moon GM in km^3/s^2
mu_e_m = mu_e+mu_m;              % earth+moon GM, km^3/s^2

% kmAU = 149597870.700;            % 1 AU in km
% GMconv = kmAU^3/86400^2;         % convert GM in AU^3/day^2 to km^3/s^2
R_e = 6378.137;                  % km
% a_e =  1.00000261 ;              % AU, Earth semimajor axis
a_m = 384400;                    % km, moon semimajor axis
e_m = 0.0554;	                 % moon eccentricity
w_m = 318.15*pi/180;	         % rad, moon argument of periapsis
M_m = 135.27*pi/180;	         % rad, moon mean anomaly
i_m = 5.16*pi/180;	             % rad, moon inclination
Om_m = 125.08*pi/180;            % rad, moon longitude of asc. node

rotMats = DCMs();
[Afun, Bfun, Mfun, rfun, vfun, ~, ~] =  KeplerFun();

%% Problem Setup
% normalized M,L,T quantities
% m1 = 1.989 * 10^30; m2 = 5.972 * 10^24; %masses of primaries (sun, earth)

m1 = 5.97219 * 10^24; m2 = 7.34767309 * 10^22; %masses of primaries (earth, moon)
G = 6.6743e-20; %km^3/s^2
% mu = 3.04036e-6; % mass parameter of Earth-Sun in CR3BP

% calculate orbital period of primaries
T_m = 2*pi*sqrt(a_m^3/mu_e_m); % sec

% calculate mean motion
n_m = 2*pi/T_m; % rad/sec

% tp_m = t - M_m/n_m;

% Nondimensionalized Canonical Units
MU = m1+m2;             % mass unit (kg)
DU = a_m;               % distance unit (km)
TU = T_m/(2*pi);        % time unit (sec)
VU = DU/TU;             % velocity unit (km/s)
MscU = 400;             % S/C mass unit (kg)
FU = 1000*MscU*DU/TU^2; % force unit (N)
mu = m2/MU;             % mass parameter

% Nondimensionalized Metric Units

T_max = 6/FU; %N/FU, max thrust of engine
Isp = 3000; % s, specific impulse of engine
g0 = 9.8066e-3; % km/s^2, Earth gravitational acceleration at sea level
c = (Isp*g0/VU); % effective velocity
m_dry = 20; % kg, structural spacecraft mass
m_tot = 1500; % kg, total spacecraft wet mass

%

% two-body orbital elements (equatorial GTO)
rp = (400+R_e)/DU; ra = (35864+R_e)/DU; % DU
a = (ra+rp)/2; % DU
e = (ra-rp)/(ra+rp);
I = 0*pi/180; % rad
w = 245*pi/180; % rad
Om = 0*pi/180; % rad
Tp = 2*pi*sqrt(a^3/(1-mu)); % TU
n = 2*pi/Tp; % rad/TU
tp = 0.0*Tp;% TU

maxdt = 6*86400/TU; % TU (days)
t0 = 0;
tf = t0 + maxdt;


% Calculate A and B matrices for spacecraft (P)
A = Afun(a,e,I,w,Om);
B = Bfun(a,e,I,w,Om);

% Define function to calculate eccentric anomaly for spacecraft
E = @(t) invKepler(Mfun(n,t,tp),e);

% known initial conditions
r0 = rfun(A,B,E(t0),e) + [-mu;0;0];
v0 = vfun(A,B,E(t0),e,n) + [0;-mu;0];

% ri  [-0.019488511458668; -0.016033479812051; 0] LU
% vi  [8.918881923678198; -4.081793688818725; 0] VU


m0 = m_tot/MscU; % 1 spacecraft mass unit = 400 kg

% terminal boundary conditions
halo=load('halo.mat').halo;
itarg = 15004; % index of target state out of saved halo [r,v]

% choose to target either halo (0) or stable manifold (1)
target = 1;

switch target
    case 0
        rf = halo(itarg,1:3).';
        vf = halo(itarg,4:6).';
    case 1
        [~,~,targ_pt]=genManifold(halo,itarg,true);
        rf = targ_pt(1:3);
        vf = targ_pt(4:6);
end


% homotopy variables
epsilon = 1; %energy optimal problem
rho = 1; % smoothing parameter
