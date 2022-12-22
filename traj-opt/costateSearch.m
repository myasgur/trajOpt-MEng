% Costate search 
function L0 = costateSearch()
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
MUnit = m1+m2;          % nondim mass unit (kg)
TUnit = 3.156e7;        % nondim time unit (seconds) 
DUnit = 1.49598e8;      % nondim distance unit (km)
VUnit = DUnit/TUnit;    % nondim speed unit (km/s)
MscUnit = 400;          % nondim S/C mass unit (kg)
units = [MUnit,MscUnit,TUnit,DUnit,VUnit];
mu = m2/MUnit;

T_max = .25;             % N, max thrust of engine
Isp = 3000;              % s, specific impulse of engine
g0 = 9.8066e-3;          % km/s^2, Earth gravitational acceleration at sea level
c = Isp*g0;              % km/s, effective velocity
m_dry = 20/MscUnit;      % structural spacecraft mass
m_tot = 1;               % total spacecraft wet mass 

% two-body orbital elements (equatorial GEO)
rp = (35864+R_e)/DUnit; ra = (35864+R_e)/DUnit;
a = (rp+ra)/2;
e = (ra-rp)/(ra+rp);
I = 0*pi/180; 
w = 0*pi/180; 
Om = 0*pi/180; 
Tp = 2*pi*sqrt((a*DUnit)^3/(mu_e*GMconv))/TUnit;
tp = 0;

maxdt = 140*86400./TUnit; % 6 months in seconds
t0 = tp; 
tf = t0 + maxdt; 
tspan = linspace(t0,tf,1e4);

[r0,v0]=kep2CR3BP(a,e,I,w,Om,tp,t0,tf);
m0 = m_tot; % 1 spacecraft mass unit = 400 kg

% initial guess of costates (optimized)
lr0 = [0 0 0].';
lv0 = [5e-8 -8e-6 2e-8].';
lm0 = 0.02;
l = [lr0;lv0;lm0];

% l = [-0.00001 0 0 -0.00000001 0 0 0.0001].'; %best manifold shooting case

%% Load Nominal Halo Orbit
% Ls = getLpoints(mu);
% Ax = 250000; Az = 420000; m = 1; Lpoint = 2; plt=0;
% halo = nominalHalo(mu,Lpoint,Ax,Az,m,plt);
halo=load('halo.mat').halo;
itarg = 700; % index of target state out of saved halo [r,v]
% [~,~,targ_pt]=genManifold(halo,itarg,true);
% rtarg = targ_pt(1:3);
% vtarg = targ_pt(4:6); 
rtarg = halo(itarg,1:3).';
vtarg = halo(itarg,4:6).';

%% Fuel Optimal Integration with Indirect Shooting Method
epsilon = 1; % start with energy optimal
rho = 1; % bang-bang smoothing parameter, vary from 1 to 1e-5
% create struct of parameters
params = struct('c',c,'Tmax',T_max,'mu',mu,'rho',rho,'epsilon',epsilon, ...
    'm1',m1,'m2',m2,'MUnit',MUnit,'MscUnit',MscUnit,'TUnit',TUnit,'DUnit', ...
    DUnit,'VUnit',VUnit);

%% Validate number of timesteps
Nt = [1e3 5e3 1e4 5e4 1e5 5e5]; % number of timesteps
linestyles = {'o';'.';'-.';':';'--';'-'};
z0 = reshape([[r0;v0;m0;lr0;lv0;lm0] eye(14)],14*15,1); % initial conditions
figure(1)
clf
for j = 1:length(Nt)
    tspan = linspace(t0,tf,Nt(j));
    zint=ode4(@(t,z) varEqsTrajOpt(t,z,params),tspan,z0(:));
    hold on
    plot3(zint(:,1),zint(:,2),zint(:,3),linestyles{j}) % trajectory
end
plot3(1-mu,0,0,'b*') %earth
plot3(halo(:,1),halo(:,2),halo(:,3),'-b')
plot3(rtarg(1),rtarg(2),rtarg(3),'.b','MarkerSize',14)
hold off
legend({'1k','5k','10k','50k','100k','500k','Earth','Halo'},'Location','best')
xlabel('$\mathbf{\hat{e}}_r \, \rightarrow$','Interpreter','Latex')
ylabel('$\mathbf{\hat{e}}_{\theta} \, \rightarrow$','Interpreter','Latex')
grid on
axis equal
xlim([0.99 1.02])

%% Random Costate Searching
N = 10000; % number of random iterations
ls = zeros(7,N); % initialize array to store all of the random costates
feas= zeros(1,N); % feasibility (infinity norm of equality constraints)

for i = 1:N
    ls(:,i) = [1e-1*randn(3,1);
               1e-3*randn(3,1);
               randn(1)]; % search on the interval [-1,1]

    [~,ceq] = nonlcon(ls(:,i));
    feas(i)= max(ceq);
    disp(feas(i))

    if abs(feas(i))<1e-5
        Z0 = reshape([[r0;v0;m0;ls(:,i)] eye(14)],14*15,1);
        [Z]=ode4(@(t,z) varEqsTrajOpt(t,z,params),tspan,Z0(:));

        CR3BPplot(1,Z)
        xlim([0.998,1.012])
        legend({'','','Earth','','','Starting Point','Target',''})

    end
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
%% Trajectory Plotter in Rotating Frame
    function CR3BPplot(fignum,z)
        figure(fignum)
        hold on
        plot3(z(:,1),z(:,2),z(:,3),'k') % trajectory
        plot3(1-mu,0,0,'.b') %earth
        plot3(halo(:,1),halo(:,2),halo(:,3),'-b')
        plot3(z(1,1),z(1,2),z(1,3),'.r','MarkerSize',14) % initial point
        plot3(rtarg(1),rtarg(2),rtarg(3),'.b','MarkerSize',14)
        hold off
        xlabel('$\mathbf{\hat{e}}_x \, \rightarrow$','Interpreter','Latex')
        ylabel('$\mathbf{\hat{e}}_y \, \rightarrow$','Interpreter','Latex')
        set(gca,'FontName','Times','FontSize',16)
        grid on
        axis equal
        axis([0.99 1.02 -7.5e-3 7.5e-3 -inf inf])

    end


end