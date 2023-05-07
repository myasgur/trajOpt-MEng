function plot_traj(y)

% Global Variables
global rf mu R_e t0 DU

[~, ~, Mfun, ~, ~, rfunv, ~] =  KeplerFun();

% Load Halo
halo=load('halo.mat').halo;

% Generate Initial Osculating Orbit
% two-body orbital elements (equatorial GTO)
rp = (400+R_e)/DU; ra = (35864+R_e)/DU; % DU
a = (ra+rp)/2; % DU
e = (ra-rp)/(ra+rp);
I = 0*pi/180; % rad
w = 245*pi/180; % rad
Om = 0*pi/180; % rad
Tp = 2*pi*sqrt(a^3/(1-mu));
n = 2*pi/Tp;
tp = 0.0*Tp;% TU

ts = linspace(t0,t0+Tp,1e4); % only 1 period of initial orbit

% known initial conditions
E  = invKepler(Mfun(n,ts,tp),e);
rs = rfunv(a,e,I,w,Om,E); % centered on Earth
% vs = vfunv(a,e,I,w,Om,E,n);

% center trajectory on Earth
z = y(:,1:3) - [-mu;0;0].';

% Plotting
figure(1) % Rotating-body frame
clf
hold on
plot3(z(:,1),z(:,2),z(:,3),'k') % trajectory
plot3(rs(1,:),rs(2,:),rs(3,:),'--m') % osculating trajectory
plot3(1-mu,0,0,'b*') %earth
plot3(halo(:,1),halo(:,2),halo(:,3),'-b')
plot3(rf(1),rf(2),rf(3),'.b','MarkerSize',14)
plot3(z(1,1),z(1,2),z(1,3),'.r','MarkerSize',14) % initial point
set(gca,'FontName','Times','FontSize',16)
hold off
legend({'Spacecraft','Starting Orbit','Earth','Halo','','Start'},'Location','best')
xlabel('$\mathbf{\hat{e}}_r \, \rightarrow$','Interpreter','Latex')
ylabel('$\mathbf{\hat{e}}_{\theta} \, \rightarrow$','Interpreter','Latex')
grid on
% xlim([0.99 1.02])
% xlim([0.9996 1.0004])
axis equal
set(gcf,'WindowStyle','Docked')
end