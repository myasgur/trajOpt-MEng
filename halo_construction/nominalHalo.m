function [halo] = nominalHalo(mu,Lpoint,Ax,Az,m,plt)
% INPUT:
%         mu = mass parameter of 3-body system
%     Lpoint = Lagrange point number (1,2,or 3) where the halo is centered
%         Ax = halo orbit amplitude in x (km)
%         Az = halo orbit amplitude in z (km)
%          m = parity (either northern = 1, southern = 3)
%        plt = plotting Y/N? (boolean)
%
% OUTPUT:
%      halo = state vectrix for nominal halo over 1 closed orbit [nx6]
%
% EXAMPLE:
%     halo = nominalHalo(3.0035e-6,2,2.5e5,4.2e5,1,1);
%     halo = nominalHalo(0.012153619140872,1,1e4,9e3,1,1);
%%
clc; close all
I6 = eye(6); Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];
fprintf('Cleared, closed, and ready to go.\n')
%% Initial Variable Set Ups
% rotMats = DCMs(); % call function to load rotation matrix functions
% m1 = 1.989 * 10^30; m2 = 5.972 * 10^24; % Actual masses in kg (earth-sun)
m1 = 5.97219 * 10^24; m2 = 7.34767309 * 10^22; %masses of primaries (earth, moon)
G = 6.6743e-20; %km^3/s^2
% ae = 1.000373836656026; % Earth semi-major axis (AU)
am = 384400; % Moon semi-major axis (km)
% kmAU = 149597870.700; % AU to km
DU = am; % distance unit
nm = sqrt((m1+m2)*G/DU^3); n = 1; % SI and nondim mean motion
% TU = 1/ne/86400; % time unit
Ls = getLpoints(mu); % Calculate x-position of collinear Lagrange pts
L = Ls(Lpoint); % Select given Lagrange pt

InitialConditions = HALO_3OA(mu, Ax, Az, nm, m, Lpoint); % Richardson approx.
x0 = InitialConditions(1);
z0 = InitialConditions(3);
dy0 = InitialConditions(5);

% Integrate initial guess to then be corrected
Y_3OA = [L*am+x0; 0; z0; 0; dy0/nm; 0]/am; 
i = 1; Y_final = zeros(1,42); t_0 = 0; h = 0.0005;
[T2, ~, ~] = findT2(Y_3OA, t_0, 10, h, Io, mu, 1);
y_int = Y_3OA;
while t_0 < T2
    Y_final(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, ...
        [y_int; Io], h, mu, 1)';
    y_int = Y_final(i,1:6)'; Io = Y_final(i, 7:42)';
    t_0 = t_0 + h; i = i+1;
end
Y_3OA_traj = Y_final;

%% Differential Correction Process
t_end = 2*T2; y_guess = Y_3OA; maxIter = 50; iter = 1;
Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];
[T2, y_STM_T2, ~] = findT2(y_guess, 0, 3, 0.0005, Io, mu, n);
fprintf('Begining Differential Correction...\n')
while (abs(y_STM_T2(4)) > 1e-10) && (abs(y_STM_T2(6)) > 1e-10) && iter < maxIter
    Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];
    [T2, y_STM_T2, STM_T2] = findT2(y_guess, 0, 3, 0.0005, Io, mu, n);
    if (mod(iter,2) == 1) % z and dy
        fprintf('Fixing z and dy\n')
        dx_correction = HALOcorrections_z_dy(y_STM_T2, STM_T2, mu);
        y_guess = [y_guess(1); y_guess(2); y_guess(3)+ dx_correction(1); ...
            y_guess(4); y_guess(5) + dx_correction(2); y_guess(6)];
    end
    if (mod(iter,2) == 0) % x and dy
        fprintf('Fixing x and dy\n')
        dx_correction = HALOcorrections_x_dy(y_STM_T2, STM_T2, mu);
        y_guess = [y_guess(1)+ dx_correction(1) ; y_guess(2); y_guess(3); ...
            y_guess(4); y_guess(5) + dx_correction(2); y_guess(6)];
    end
    t_end = 2*T2;
    iter = iter+1;
    if iter==maxIter
        disp('Not enough iterations to converge.')
    end
end
fprintf('Differential Corrections Complete!\n')
y_int = y_guess;
t_0 = 0; i = 2; h = 0.001;
Y_final = zeros(1,42);
Y_final(1,:) = [y_int; Io];
while t_0 < t_end + 1*h
    Y_final(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, ...
        [y_int; Io], h, mu, 1).';
    y_int = Y_final(i,1:6).'; t_0 = t_0 + h; Io = Y_final(i, 7:42)';
    i = i+1;
end

%% First Round of Plotting
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    zlabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')
if plt==1
    figure(1) % Planar Views of the Stable (Differentially Corrected) Orbit
    subplot(1,3,1)
    plot(Y_final(:,1)-L, Y_final(:,2), 'k', 'LineWidth', 1.5); hold on
    plot(0,0,'kx')
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    axis square; box off; title('X-Y Plot')
    xlim([min(Y_final(:,1))*0.9959, max(Y_final(:,1))*1.0025]-L);
    ylim([min(Y_final(:,2))*1.1, max(Y_final(:,2))*1.1])

    subplot(1,3,2)
    plot(Y_final(:,1)-L, Y_final(:,3), 'k', 'LineWidth', 1.5); hold on
    plot(0,0,'kx')
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')
    axis square; box off; title('X-Z Plot')
    xlim([min(Y_final(:,1))*0.9959, max(Y_final(:,1))*1.0025]-L);
    ylim([min(Y_final(:,3))*1.1, max(Y_final(:,3))*1.1])

    subplot(1,3,3)
    plot(Y_final(:,2), Y_final(:,3), 'k', 'LineWidth', 1.5); hold on
    plot(0,0,'kx')
    xlabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')
    axis square; box off; title('Y-Z Plot')
    xlim([min(Y_final(:,2))*1.1, max(Y_final(:,2))*1.1])
    ylim([min(Y_final(:,3))*1.1, max(Y_final(:,3))*1.1])


    figure(2) % 3OA_T2 vs Differentially Corrected Stable Orbit
    plot3(Y_3OA_traj(:,1), Y_3OA_traj(:,2), Y_3OA_traj(:,3), ...
        'Color', [0.9290 0.6940 0.1250], 'LineWidth', 0.8);
    hold on
    plot3(Y_final(:,1), Y_final(:,2), Y_final(:,3), ...
        'Color', [0.4660 0.6740 0.1880], 'LineWidth', 0.8);
    plot3(L, 0, 0, 'kx')
    legend('3OA Initial HALO Orbit (Thru T2)', ...
        'Differentially Corrected HALO Orbit', ['L' num2str(Lpoint)])
    set(gca,'YTickLabel',[]); xlabel('X')
    set(gca,'XTickLabel',[]); ylabel('Y')
    set(gca,'ZTickLabel',[]); zlabel('Z')
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    zlabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')

end
%% Nominal Orbit Acquisition
fprintf('Propogating for 1 period (~%3.0f days)..\n', T2*2/nm/24/60/60)
h = T2/15000; % "Time"-step size
t_0 = 0; t_end = 2*T2; Y_nominal = zeros(1,42);
y0 = y_guess;
Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];
Y_nominal(1,1:6) = y_guess; %Y_base = zeros(1,6);
% Iterates state and STM through 1 orbit to aquire the nominal trajectory
for i = 1:1:t_end/h + 10
    Y_nominal(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, ...
        t_0, [y0; Io], h, mu, n)';
    y0 = Y_nominal(i,1:6)'; Io = Y_nominal(i, 7:42)';
    t_0 = t_0 + h;
end % Steps

halo = zeros(2,6); STMsave = zeros(6,6,2);
% Breaks up integration into State Vector and STM
for i = 1:1:length(Y_nominal)
    [Yi, STM] = xSTM(Y_nominal(i,:));
    halo(i,:) = Yi;
    STMsave(:,:,i) = STM;
end

fprintf('Acquired nominal.\n')

% %% Extract JWST Horizons Data
% if plt==true
% 
%     JWST = readHorizons('horizons_results-jwst.txt',1); %wrt to Earth
%     JWt = JWST{1}; % time (JD)
%     JWr = [JWST{3},JWST{4},JWST{5}]; % position (km)
% 
%     % Rotate into Sun-Earth perifocal frame
%     JWrrot = zeros(size(JWr)); 
%     nu0 = 94.6*pi/180; % rad
%     for i = 1:length(JWt)
%         DCM = rotMats{3}((JWt(i)-JWt(1))./TU+nu0);
% 
%         % Inertial positions and velocities in perifocal
%         JWrrot(i,:) = DCM*(JWr(i,:)).'./DU; % sun to particle
%     end
% end

%% Second Round of Plotting
if plt==true
    % Axis Stuff
    xlmax = max(halo(:,1)); xlmin = min(halo(:,1));
    ylmax = max(halo(:,2)); ylmin = min(halo(:,2));
    zlmax = max(halo(:,3)); zlmin = min(halo(:,3));
    xlw = 1.25*(xlmax - xlmin)/2; xlm = (xlmax + xlmin)/2;
    ylw = 1.25*(ylmax - ylmin)/2; ylm = (ylmax + ylmin)/2;
    zlw = 1.25*(zlmax - zlmin)/2; zlm = (zlmax + zlmin)/2;

    figure(3) % Nominal Trajectory Plot
    plot3(halo(:,1), halo(:,2), halo(:,3), 'Color', [0.9 0.1 0.1], ...
        'LineWidth', 1.8)
    hold on;
%     plot3(1+JWrrot(:,1),JWrrot(:,2),JWrrot(:,3),'Color',[0, 0.4470, ...
%         0.7410],'LineWidth',1.4)
    plot3(L, 0, 0,'kx');
    hold off 
    axis square
    grid on
    box on
    xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]);
    xlabel('$\mathbf{\hat{e}}_x$','Interpreter','Latex')
    ylabel('$\mathbf{\hat{e}}_y$','Interpreter','Latex')
    zlabel('$\mathbf{\hat{e}}_z$','Interpreter','Latex')
%     legend({'Nominal-DC','JWST-Horizons'})
end
set(findall(groot,'type','axes'),'FontName','Times','FontSize',14)
save halo.mat halo
end