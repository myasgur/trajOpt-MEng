%% Generate Invariant Manifold from Sun-Earth Halo Orbit
function [Ws,Wu,targ_pt] = genManifold(halo,i_h,plt)
% INPUT: 
%       halo = [nh x 6] array of states for halo orbit or struct containing
%              parameters to generate nominal periodic orbit
%              {Lpoint,Ax,Az,m,plt}
%           EX: halo = struct('Lpoint',1,'Ax',3.2e3,'Az',4.9e3,'m', 3,'plt', false); 
%               [Ws,Wu,targ_pt] = genManifold(halo,15004,true);
%
%        i_h = [nh x 1] indices of target point(s) on the halo from which 
%               to generate the manifolds 
%
% OUTPUTS: 
%     Ws, Wu = [nh x N x 6] stable/unstable manifolds to/from given 
%                point(s) on halo
%     targ_pt= [6 x nh] target point(s) for trajectory shooting
%                       ** recommended to select only 1 pt **
%% Useful constants
% mu_s = 2.9591220828559093E-04;   %sun GM in AU^3/day^2
% mu_e_s = 2.9591309705483544E-04; %earth+sun GM, AU^3/day^2
% mu_e = mu_e_s-mu_s;              %earth GM
% kmAU = 149597870.700;            %1 AU in km
% a_e =  1.000373836656026E+00 ;   %AU
mu_e = 398600.436;               % earth GM in km^3/s^2
mu_m = 4902.800066;              % moon GM in km^3/s^2
mu_e_m = mu_e+mu_m;              % earth+moon GM, km^3/s^2
a_m = 384400;                    %km
%% Problem Setup
mu = 0.012153619140872; %=3.0035e6

% calculate mean motion
n_m = sqrt(mu_e_m/a_m^3); % rad/sec

DU = a_m; %km
TU = 1/n_m; %sec

if isstruct(halo)
% Generate the JWST Halo Orbit
    % Lpoint = 2; % L1 or L2
    % Ax = 250000; Az = 420000; % amplitudes
    % m = 1; % northern = 1, southern = 3
    % plt=0; % plotting? Y/N
    halo = nominalHalo(mu,halo.Lpoint,halo.Ax,halo.Az,halo.m,halo.plt);
end

% or load it from a saved file
% halo = load('halo.mat').halo;

nh = length(i_h); % number of starting points
N = 1e5; % integrated points
Ws = nan(nh,N,6); % stable manifold
Wu = Ws; % unstable manifold
tspan = linspace(0,7.5,N); % coast time (140 days nominal)
opts = odeset('AbsTol',1e-16,'RelTol',1e-13,'Events',@manifoldStop);
for i = 1:nh
    Z0 = halo(i_h(i),:).'-[0;0;0;0;1e-8;0]; % small perturbation off halo
    % integrate forwards in time away from the halo...
    [tu,Zu,~,~,~] = ode113(@cr3bp_eom,tspan,Z0,opts);
    Wu(i,1:size(Zu,1),:) = Zu;
    % ...and then backwards in time onto the halo 
    [ts,Zs,~,~,~] = ode113(@cr3bp_eom,fliplr(tspan),Z0,opts);
    Ws(i,size(Zs,1):-1:1,:) = Zs;
end

% For optimization, select the first point on the stable manifold as the
% target for the controlled segment of the trajectory:
targ_pt = reshape(squeeze(Ws(:,1,:)),6,nh);

%% Plotting
if plt
    figure(1)
    clf
    for i = 1:nh
        hold on
        plot3(Ws(i,1:10:length(ts),1),Ws(i,1:10:length(ts),2), ...
            Ws(i,1:10:length(ts),3),'-b')
        plot3(Ws(i,1,1),Ws(i,1,2),Ws(i,1,3),'.b')

        plot3(Wu(i,1:10:length(tu),1),Wu(i,1:10:length(tu),2), ...
            Wu(i,1:10:length(tu),3),'-r')
        plot3(Wu(i,1,1),Wu(i,1,2),Wu(i,1,3),'.r')
    end
    plot3(halo(:,1),halo(:,2),halo(:,3),'-k')
    plot3(1-mu,0,0,'b.','MarkerSize',10)
    plot3( -mu,0,0,'b.','MarkerSize',10)
    hold off
%     axis([-inf 1.3 -0.075 0.075 -0.1 0.1])
%     zlim([-0.1 0.1])
    axis equal
    
    box on
    grid on
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    zlabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')
    legend({'Stable','','Unstable'})
    set(gca,'FontName','Times','FontSize',14)
end
%% CR3BP Equations of Motion
    function DY = cr3bp_eom(~,Y)

        x = Y(1); y = Y(2); z = Y(3);
        v_x = Y(4); v_y = Y(5); v_z = Y(6);

        r1sq = sqrt((mu+x).^2+y.^2+z.^2);
        r2sq = sqrt((mu+x-1).^2+y.^2+z.^2);

        r1_32 = 1./(r1sq.^3);
        r2_32 = 1./(r2sq.^3);

        Ux=x-mu.*(mu+x-1).*r2_32-(1-mu).*(mu+x).*r1_32;

        Uy=y-mu.*y.*r2_32-(1-mu).*y.*r1_32;

        Uz=-mu.*z.*r2_32-(1-mu).*z.*r1_32;

        DY = [v_x;
            v_y;
            v_z;
            2*v_y + Ux;
            -2*v_x + Uy;
            Uz];
    end

%% Manifold Event Function
    function [value,isterminal,dir] = manifoldStop(~,z)
        mu = 0.012153619140872;
        r1sq = sqrt((mu+z(1)).^2+z(2).^2+z(3).^2);
        r2sq = sqrt((mu+z(1)-1).^2+z(2).^2+z(3).^2);

        value = [r1sq-6378.1/DU; r2sq-1738.1/DU]; % stop at body surfaces
        isterminal = [1;1];
        dir = [0;0];
    end
end