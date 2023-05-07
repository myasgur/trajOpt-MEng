%% Clear
clear all; close all; clc;
%% Problem Setup
problem_setup

%%
%global options for integration
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@cr3bp_event);

% minimize fuel use such that the final state stops on the nominal halo
% and the mass costate is zero at t=tf.
% optimization state = [lr0,lv0,lm0]
fopts = optimoptions('fmincon','Display','iter','EnableFeasibilityMode',true,'Algorithm','sqp');

objfun = @costFun;
confun = @nonlcon;

% upper and lower bounds on search variable
lb = [-100 -100 -100 -10 -10 -10  0];
ub = [ 100  100  100  10  10  10 10];

% number of iterations of optimization
N = 100;
i = 1;

while i<N
    disp(['------------------ ' ...
        'Iteration ',num2str(i),' of ',num2str(N), ...
         ' ------------------'])
    X0 = [unifrnd(-40,40,3,1);unifrnd(-2,2,3,1);unifrnd(0,2)]; % random initial guess on interval [0 1]
    [Xopt,~,exitflag,~] = fmincon(objfun,X0(:),[],[],[],[],lb,ub,confun,fopts);
    i=i+1;
    if exitflag == 1
        break
    end
end
% lam_opt = psoSearchVars(Xopt);


Zopt0 = [r0;v0;m0;Xopt(:)];
% Reconstruct trajectory with optimized ICs
[tOpt,Zopt]=ode78(@cr3bp_EOM,[t0 tf],Zopt0(:),opts);
plot_traj(Zopt)
%% Objective Function
function obj = costFun(x)
global t0 tf r0 v0 m0 rf vf  c  epsilon

% z0 = [r0,v0,m0,lr0,lv0,lm0];
% lams = psoSearchVars(x);
lams=x;
z0 = [r0;v0;m0;lams];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@cr3bp_event);
[t,z,~,~,~] = ode78(@cr3bp_EOM,[t0 tf],z0(:),opts);

u = zeros(length(t),1);

% Unpack State/Costate
r = z(:,1:3); v = z(:,4:6); m = z(:,7); 
lr = z(:,8:10); lv = z(:,11:13); lm = z(:,14);


for i = 1:length(t)
    S = 1-c.*norm(lv(i,:))./(m(i)) - lm(i); % switching function

    % u = 0.5*(1-tanh(S./rho)); %hyperbolic tangent smoothing

    if S > epsilon
        u(i) = 0;
    elseif S < -epsilon
        u(i) = 1;
    else
        u(i) = 0.5*(1 - S/epsilon);
    end
end

% Shooting Function 
rdiff = r(end,:).' - rf;
vdiff = v(end,:).' - vf;

phi = [rdiff;vdiff;lm(end)];
% pf = 1; % terminal constraint weight
W = diag([10 10 10 1 1 1 1]); %weight matrix

if length(t) == 1
    obj = inf;
else
    % obj = T_max/c*trapz(t,u.^2) + pf*(phi)'*phi;
    obj = phi'*W*phi;
end
end

%% Nonlinear Constraints for Optimization

function [c,ceq] = nonlcon(x)
global t0 tf r0 v0 m0 rf vf
%z = [r0,v0,m0,lr0,lv0,lm0];
% lams = psoSearchVars(x);
lams = x;
z0 = [r0;v0;m0;lams];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@cr3bp_event); % integration options
[~,z,~,~,~]=ode78(@cr3bp_EOM,[t0 tf],z0,opts);


rdiff = z(end,1:3).' - rf;
vdiff = z(end,4:6).' - vf;

% r(tf)-r_halo = 0,  v(tf)-v_halo = 0, lambda_m(tf) = 0
ceq = [rdiff;vdiff;z(end,14)];
c = [];

end
