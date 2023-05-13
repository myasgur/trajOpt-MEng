%% Cost Function for PSO
% based on Eqs (31-33) of Jiang, Baoyin, et al., 2012 - Practical techniques for 
% low-thrust trajectory optimization with homotopic approach

function z = TimeOptCost(x)
% INPUT: 
%           x = search-variable with N components
% OUTPUT: 
%           z = F(x) = time optimal cost function 
%
%% Parameters
% global MU DU TU VU MscU FU mu    
global T_max c t0 r0 v0 m0 rf vf

pf = 100; % penalty factor (weight)

%% Initial Costate and Final Time Guess (Search Variables)
% lams = psoSearchVars(x(1:7));
lams = x(1:7).';
tf = abs(x(8));

%% Integrate Trajectory

y0 = [r0;v0;m0;lams]; % initial conditions
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@cr3bp_event); % integration options
[tspan,y,~,~,~]=ode113(@cr3bp_EOM_time,[t0 tf],y0,opts);
% tspan = t0:0.0001:tf;
% y = ode4(@cr3bp_EOM_time,tspan,y0);

% Plot
% plot_traj(y)
% pause()

%% Unpack State/Costate
r = y(:,1:3); v = y(:,4:6); m = y(:,7); 
lr = y(:,8:10); lv = y(:,11:13); lm = y(:,14);

%% Determine Throttle and Hamiltonian for Full Trajectory
u = zeros(length(tspan),1);
Ht = zeros(length(tspan),1);

for i = 1:length(tspan)
    S = -c.*norm(lv(i,:))./(m(i)) - lm(i); % switching function
   
    if S > 0
        u(i) = 0;
    elseif S < 0
        u(i) = 1;
    else
        u(i) = 0.5;
    end

    [g,h]=ghFunctions(r(i,:),v(i,:));

    Ht(i) = dot(lr(i,:),v(i,:)) + ...
            dot(lv(i,:),(g + h + u(i)*T_max/m(i))) + ...
            - lm(i)*u(i)*T_max/c + 1;

end

%% Shooting Function 
rdiff = r(end,:).' - rf;
vdiff = v(end,:).' - vf;

phi = [rdiff;vdiff;lm(end);Ht(end)];

%% Cost Function
z = pf*(phi)'*phi; %tf>=2

end