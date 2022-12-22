function varargout = main()
%% Useful constants
mu_s = 2.9591220828559093E-04;   %sun GM in AU^3/day^2
mu_e_s = 2.9591309705483544E-04; %earth+sun GM, AU^3/day^2
mu_e = mu_e_s-mu_s;              %earth GM
kmAU = 149597870.700;            %1 AU in km
GMconv = kmAU^3/86400^2;         % convert GM in AU^3/day^2 to km^3/s^2
R_e = 6378.137;                  % km

%% Problem Setup
% normalized M,L,T quantities
m1 = 1.989 * 10^30; m2 = 5.972 * 10^24; %masses of primaries
G = 6.6743e-20; %km^3/s^2
% mu = 3.04036e-6; % mass parameter of Earth-Sun in CR3BP
MUnit = m1+m2; % nondim mass unit (kg)
TUnit = 3.156e7; % nondim time unit (seconds) 
DUnit = 1.49598e8; % nondim distance unit (km)
VUnit = DUnit/TUnit; % nondim speed unit (km/s)
MscUnit = 400; % nondim S/C mass unit (kg)
units = [MUnit,MscUnit,TUnit,DUnit,VUnit];
mu = m2/MUnit;

T_max = .25; %N, max thrust of engine
Isp = 3000; % s, specific impulse of engine
g0 = 9.8066e-3; % km/s^2, Earth gravitational acceleration at sea level
c = Isp*g0; % km/s, effective velocity
m_dry = 20/MscUnit; % kg, structural spacecraft mass
m_tot = 1; % kg, total spacecraft wet mass 


% two-body orbital elements (equatorial GTO)
rp = (35864+R_e)/DUnit; ra = (35864+R_e)/DUnit;
a = (ra+rp)/2;
e = (ra-rp)/(ra+rp);
I = 0*pi/180; 
w = 0*pi/180; 
Om = 0*pi/180; 
tp = 0;

maxdt = 110*86400./TUnit; % 110 days in seconds
t0 = 0; 
tf = t0 + maxdt; 
tspan = linspace(t0,tf,5e4);

[r0,v0]=GTO2CR3BP(a,e,I,w,Om,tp,t0,tf);
m0 = 1; % 1 spacecraft mass unit = 400 kg


% initial guess of costates (to be optimized)
l = [0 0 0.5 0 -0.0003 -0.0015 0.0001].';
   
%% Load Nominal Halo Orbit and Stable Manifold
% Ls = getLpoints(mu);
% Ax = 250000; Az = 420000; m = 1; Lpoint = 2; plt=0;
% halo = nominalHalo(mu,Lpoint,Ax,Az,m,plt);
halo=load('halo.mat').halo;
itarg = 700; % index of target state out of saved halo [r,v]
% rtarg = halo(itarg,1:3).';
% vtarg = halo(itarg,4:6).';

% generate the stable manifold and choose a point near earth to target
[~,~,targ_pt]=genManifold(halo,itarg,false);
rtarg = targ_pt(1:3);
vtarg = targ_pt(4:6); 

%% Fuel Optimal Integration with Indirect Shooting Method
epsilon = 1; % start with energy optimal
rho = 1; % bang-bang smoothing parameter, vary from 1 to 1e-5
% create struct of parameters
params = struct('c',c,'Tmax',T_max,'mu',mu,'rho',rho,'epsilon',epsilon, ...
    'm1',m1,'m2',m2,'MUnit',MUnit,'MscUnit',MscUnit,'TUnit',TUnit,'DUnit', ...
    DUnit,'VUnit',VUnit);

% opts = odeset('AbsTol',1e-6);
tic

%global options for integration
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

% minimize fuel use such that the final state stops on the nominal halo
% and the mass costate is zero at t=tf.
% optimization state = [lr0,lv0,lm0]
fopts = optimoptions('fmincon','Display','iter','StepTolerance',eps, ...
    'EnableFeasibilityMode',true);

Lopt0 = l; % [costates]
objfun = @maxMass;
confun = @nonlcon;

% upper and lower bounds on costates
lb = [-1e0*ones(3,1);-1e-1*ones(3,1);-1e-15];
ub = [1e0*ones(3,1);1e-1*ones(3,1);1];

[Lopt] = fmincon(objfun,Lopt0(:),[],[],[],[],lb,ub,confun,fopts);

Zopt0 = [[r0;v0;m0;Lopt(:)], eye(14)];
% Reconstruct trajectory with optimized ICs
[tOpt,Zopt]=ode45(@(t,z) varEqsTrajOpt(t,z,params),tspan,Zopt0(:),opts);
toc

% Output
varargout = {[tOpt,Zopt]};

%% Plotting

figure(1)
clf
hold on
plot3(Zopt(:,1),Zopt(:,2),Zopt(:,3),'k') % trajectory
plot3(1-mu,0,0,'b*') %earth
plot3(halo(:,1),halo(:,2),halo(:,3),'-b')
plot3(rtarg(1),rtarg(2),rtarg(3),'.b','MarkerSize',14)
plot3(Zopt0(1,1),Zopt0(2,1),Zopt0(3,1),'.r','MarkerSize',14) % initial point
set(gca,'FontName','Times','FontSize',16)
hold off
legend({'Spacecraft','Earth','Halo'},'Location','best')
xlabel('$\mathbf{\hat{e}}_r \, \rightarrow$','Interpreter','Latex')
ylabel('$\mathbf{\hat{e}}_{\theta} \, \rightarrow$','Interpreter','Latex')
grid on
xlim([0.99 1.02])
axis equal
%% Objective Function
    function mdiff = maxMass(l)
        % z = [r0,v0,m0,lr0,lv0,lm0];
        z = reshape([[r0;v0;m0;l] eye(14)],14*15,1);
        zint=ode4(@(t,z) varEqsTrajOpt(t,z,params),tspan,z(:));
        mdiff = abs(m_tot-zint(end,7));
    end

%% Nonlinear Constraints for Optimization

    function [c,ceq] = nonlcon(l)
        
        %z = [r0,v0,m0,lr0,lv0,lm0,STM];
        
        z = reshape([[r0;v0;m0;l] eye(14)],14*15,1);
        zint=ode4(@(t,z) varEqsTrajOpt(t,z,params),tspan,z(:));

        rdiff = zint(end,1:3).' - rtarg;
        vdiff = zint(end,4:6).' - vtarg;

        % r(tf)-r_halo = 0,  v(tf)-v_halo = 0, lambda_m(tf) = 0
        ceq = [rdiff;vdiff;zint(end,14)];
        c = [];
        
    end

end
