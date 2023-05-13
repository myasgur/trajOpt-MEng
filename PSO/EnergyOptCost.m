%% Cost Function for PSO
% based on Eqs (31-33) of Jiang, Baoyin, et al., 2012 - Practical techniques for 
% low-thrust trajectory optimization with homotopic approach

function z = EnergyOptCost(x)
% INPUT: 
%           x = search-variable with N components
% OUTPUT: 
%           z = F(x) = energy optimal cost function 
%
%% Parameters
% global MU DU TU VU MscU FU mu 
global T_max c epsilon 
global t0 tf r0 v0 m0 l0 rf vf

W = diag([10 10 10 1 1 1 1]); %weighting matrix

%% Initial Costate Guess (Search Variables)
lams = psoSearchVars(x);

% Note: l0 changes in psoSearchVars() although not return since global

%% Integrate Trajectory

y0 = [r0;v0;m0;lams]; % initial conditions
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@cr3bp_event); % integration options
[t,y,~,~,~]=ode78(@cr3bp_EOM,[t0 tf],y0,opts);

% turn off warning that integration failed
[~, MSGID] = lastwarn();
warning('off', MSGID)

% Plot
% plot_traj(y)
% pause()

%% Unpack State/Costate
r = y(:,1:3); v = y(:,4:6); m = y(:,7); 
lr = y(:,8:10); lv = y(:,11:13); lm = y(:,14);

%% Determine Throttle for Full Trajectory
u = zeros(length(t),1);

for i = 1:length(t)
    S = 1-c.*norm(lv(i,:))./(l0*m(i)) - lm(i)/l0; % switching function
    
    % u = 0.5*(1-tanh(S./rho)); %hyperbolic tangent smoothing
        
    if S > epsilon
        u(i) = 0;
    elseif S < -epsilon
        u(i) = 1;
    else
        u(i) = 0.5*(1 - S/epsilon);
    end
end

%% Shooting Function 
rdiff = r(end,:).' - rf;
vdiff = v(end,:).' - vf;

phi = [rdiff;vdiff;lm(end)];


%% Cost Function

if length(t) == 1
    z = inf;
else
    z = (phi)'*W*phi;
end

end