clc; clear; close all
I6 = eye(6); Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)]; 
fprintf('Cleared, closed, and ready to go.\n')
%% Initial Variable Set Ups
m1 = 1.989 * 10^30; m2 = 5.972 * 10^24; % Actual Masses (earth-sun)
G = (6.67408e-11); mu = 3.04036e-6;
au = 149597870.700;
ne = -sqrt((m1+m2)*G/((au)*10^3)^3); n = 1;
numOrbs = 14; % Number of Periods to Integrate over
Ls = getLpoints(mu);
L1 = Ls(1); L2=Ls(2); L3=Ls(3);
% Howell mu = 0.04 Demo --
% mu = 0.04; au = 1;
% Y_3OA = [0.723268; 0; 0.040000; 0; 0.198019; 0];

% An orbit around the Earth!
% yinit = [au-mu*au+10000; 0; 5000; 0; (6.3)/ne; 0]/au;

% Richardson Demo --
% x0 = 2.452044207*10^5; y0 = 0; z0 = 1.003274912*10^5; 
% dx0 = 0; dy0 = -2.91755907*10^2/1000; dz0 = 0; fprintf('Starting 3rd Order Approximation...\n')
Ax = 205000; Az = 125000;
InitialConditions = HALO_3OA(mu, Ax, Az, ne);
x0 = InitialConditions(1);
z0 = InitialConditions(3);
dy0 = -InitialConditions(5);

Y_3OA = [L1*au+x0; 0; z0; 0; L1*au-(L1*au*ne+dy0)/ne; 0]/au; 
i = 1; Y_final = zeros(1,42); t_0 = 0; h = 0.0005;
[T2, ~, ~] = findT2(Y_3OA, t_0, 10, h, Io, mu, 1);
y_int = Y_3OA;
while t_0 < T2
    Y_final(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y_int; Io], h, mu, 1)';
    y_int = Y_final(i,1:6)'; Io = Y_final(i, 7:42)'; 
    t_0 = t_0 + h; i = i+1;
end
Y_3OA_traj = Y_final;
  
%% Differential Correction Process
t_end = 2*T2; y_guess = Y_3OA; zeroPoint = 30; 
fprintf('Begining Differential Correction...\n') 
for narrowDown = 1:1:20
    Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];
    [T2, y_STM_T2, STM_T2] = findT2(y_guess, 0, 3, 0.0005, Io, mu, n); 
    if (mod(narrowDown,3) == 2) % x and z
        fprintf('Fixing x and z\n')
        dx_correction = HALOcorrections_x_z(y_STM_T2, STM_T2, mu);
        y_guess = [y_guess(1) + dx_correction(1); y_guess(2); y_guess(3)+ dx_correction(2); ...
                   y_guess(4); y_guess(5); y_guess(6)];
    end
    if (mod(narrowDown,3) == 1) % z and dy
        fprintf('Fixing z and dy\n')
        dx_correction = HALOcorrections_z_dy(y_STM_T2, STM_T2, mu);
        y_guess = [y_guess(1); y_guess(2); y_guess(3)+ dx_correction(1); ...
        y_guess(4); y_guess(5) + dx_correction(2); y_guess(6)];
    end
    if (mod(narrowDown,3) == 0) % x and dy
        fprintf('Fixing x and dy\n')
        dx_correction = HALOcorrections_x_dy(y_STM_T2, STM_T2, mu);
        y_guess = [y_guess(1)+ dx_correction(1) ; y_guess(2); y_guess(3); ...
                   y_guess(4); y_guess(5) + dx_correction(2); y_guess(6)];
    end
    t_end = 2*T2;
    if (abs(y_STM_T2(4)) < 1e-10 && abs(y_STM_T2(6)) < 1e-10)
        fprintf('Differential Corrections Complete!\n')
        break 
    end
end
y_int = y_guess;
t_0 = 0; i = 2; h = 0.001;
Y_final = zeros(1,42);
Y_final(1,:) = [y_int; Io];
while t_0 < t_end + 1*h
    Y_final(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y_int; Io], h, mu, 1).';
    y_int = Y_final(i,1:6).'; t_0 = t_0 + h; Io = Y_final(i, 7:42)'; 
    i = i+1;
end

%% First Round of Plotting
figure(1) % Planar Views of the Stable (Differentially Corrected) Orbit
subplot(1,3,1)
plot(Y_final(:,1), Y_final(:,2), 'k', 'LineWidth', 1.5); hold on
plot(L2,0,'kx')
xlabel('X'); ylabel('Y'); axis square; box off; title('X-Y Plot') 
xlim([min(Y_final(:,1))*.9975 max(Y_final(:,1))*1.0025]); 
ylim([-max(Y_final(:,2))*1.1 max(Y_final(:,2))*1.1])

subplot(1,3,2)
plot(L2(1)-Y_final(:,2), Y_final(:,3), 'k', 'LineWidth', 1.5); hold on 
plot(L2,0,'kx')
xlabel('Y'); ylabel('Z'); axis square; box off; title('Y-Z Plot') 
xlim([min(L2(1)+Y_final(:,2))*.9995 max(L2(1)+Y_final(:,2))*1.0005]); 
ylim([min(Y_final(:,3))*1.05 max(Y_final(:,3))*1.05])

subplot(1,3,3)
plot(Y_final(:,1), Y_final(:,3), 'k', 'LineWidth', 1.5); hold on
plot(L2,0,'kx')
xlabel('X'); ylabel('Z'); axis square; box off; title('X-Z Plot') 
xlim([min(Y_final(:,1))*.9975 max(Y_final(:,1))*1.0025]); 
ylim([min(Y_final(:,3))*1.05 max(Y_final(:,3))*1.05])

figure(2) % 3OA_T2 vs Differentially Corrected Stable Orbit
plot3(Y_3OA_traj(:,1), Y_3OA_traj(:,2), Y_3OA_traj(:,3), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 0.8);
hold on
plot3(Y_final(:,1), Y_final(:,2), Y_final(:,3), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 0.8);
plot3(L1, 0, 0, 'kx')
legend('3OA Initial HALO Orbit (Thru T2)', 'Differentially Corrected HALO Orbit', 'L2') 
xlim([0.985 0.995]);
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y')
set(gca,'ZTickLabel',[]); zlabel('Z')


%% Nominal Orbit Acquisition for numOrbs Number of Orbits
fprintf('Propogating for %2.0f periods (~%3.0f days)..\n',numOrbs, -numOrbs*T2*2/ne/24/60/60) 
h = T2/1500; % "Time"-step size
t_0 = 0; t_end = 2*T2; Y_nominal = zeros(1,42); 
y0 = y_guess;
Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)]; 
Y_nominal(1,1:6) = y_guess; Y_base = zeros(1,6);
% Iterates state and STM through 1 orbit to aquire the nominal trajectory 
for i = 1:1:t_end/h + 10
    Y_nominal(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y0; Io], h, mu, n)'; 
    y0 = Y_nominal(i,1:6)'; Io = Y_nominal(i, 7:42)';
    t_0 = t_0 + h;
end % Steps

Y_xyz = zeros(2,6); STM = zeros(6,6,2); STMsave = zeros(6,6,2);
% Breaks up integration into State Vector and STM
for i = 1:1:length(Y_nominal)
    [Yi, STM] = xSTM(Y_nominal(i,:)); 
    Y_xyz(i,:) = Yi;
    STMsave(:,:,i) = STM;
end
% Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];

% Repeats the nominal trajectory numOrbs times so station keeping has the 
% nominal orbit to target. Even stable initial conditions deviates in the 
% co-linear L point HALO orbits, so you can't just integrate the initial 
% condition for 2+ periods.
Y_xyz = repmat(Y_xyz(1:end,:),numOrbs,1); 
STMsave = repmat(STMsave(:,:,1:end),1,1,numOrbs); 
Y_nominal = repmat(Y_nominal(1:end,:),numOrbs,1); 
fprintf('Acquired nominal.\n')

%% Station Keeping Correction Method
% Initial Deviation, either random perturbation or 0 (since L(1-3) orbits
% are unstable anyway).
initialDeviation = [rand(1)*1; rand(1)*10; rand(1)*1; ...
rand(1)*.0001*ne; rand(1)*.001*ne; rand(1)*.0001*ne]/au; % Random Deviation 
% initialDeviation = [0; 0; 0; 0; 0; 0]; % Zero Initial Deviation
y0 = Y_xyz(1,1:6)' + initialDeviation;
Io = [I6(:,1); I6(:,2); I6(:,3); I6(:,4); I6(:,5); I6(:,6)];

% Setting up variables/matrices used in station keeping strategy.
y0_nc = y0; Io_nc = Io; Y_delta_nc = zeros(1,42); Y_delta_nc(1,1:6) = y0; 
t_0 = 0; burnCounter = 0; burnSpot = zeros(1,3);
DeltaVStore = zeros(1,6); DeltaVStoreP = zeros(1,6);
delay = 500; quickFix = 0; % Sets delay time between burns (1000 = 60 days) 
Y_delta = zeros(size(Y_xyz,1),42);
Y_delta(1,1:6) = y0; Y_delta(1,7:end) = Io;
fprintf('Beginning Integration with Course Corrections...\n')
propOut = 675; % Time to propogate out to calculate deltaV's (41 days). 
for i = 2:1:length(Y_xyz)
    % Integrate
    Y_delta(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y0; Io], h, mu, n)'; 
    Y_delta_nc(i,:) = Runge_Kutta_Merson(@stateTransitionMatrix, t_0, [y0_nc; Io_nc], h, mu, n)'; 
    y0 = Y_delta(i,1:6)';
    Io = Y_delta(i, 7:42)';
    y0_nc = Y_delta_nc(i,1:6)';
    Io_nc = Y_delta_nc(i,7:42)';

    if( delay < 1)
        t_0 = t_0-h;
        % Gets closest nominal location
        displaces = Y_delta(i,1:3)'-Y_xyz(:,1:3)';
        displaces = sqrt(displaces(1,:).^2+displaces(2,:).^2+displaces(3,:).^2);
        [a, b] = min(displaces);
        [deltaV, Y_prop] = StationKeep(Y_delta(i,1:6)', Y_delta(i,7:end)', ...
            propOut, t_0, h, i, mu, 1, STMsave, Y_xyz, b);
        delay = 50;

        % Checks if conditions for making a manuever have been met (minimum
        % delta V threshold, 2 m/s max, etc.)
        if ((abs(norm(deltaV)*au*ne*1000) > 1e-1 && ...
                abs(norm(deltaV)*au*ne*1000) < 2) || (quickFix == 1) )
%             burnCounter(end+1) = i;
            delay = 1000;
            fprintf('deltaV = %3.2f m/s ', norm(deltaV)*au*(-ne)*1000) 
            fprintf('occuring at t = %4.0f days\n', -i*h/ne/24/60/60) 
            y0 = y0 + deltaV'; y0_nc = Y_delta_nc(i,1:6)';
            Io = Y_delta(i, 7:42)';
            burnSpot(i,:) = y0(1:3);
            DeltaVStore(i,:) = deltaV;
            DeltaVStoreP(i,:) = deltaV;
            quickFix = 0;
        end
        % Strategy imposed for burns over 2m/s (none for t_t = 41 days).
        if (abs(norm(deltaV)*au*ne*1000) > 2)
            delay = 300;
            quickFix = 1;
        end

    end % delay < 1
    t_0 = t_0 + h; % Step forward in time. 
    delay = delay - 1; % Handles the delay between burns.
end % 2:1:length(Y_xyz)

fprintf('Done with course correction!\n')

%% Process Deviations
% Actually - gets deviation from Orbit (the "right" way)
fprintf('Processing deviation data...\n') 
for i = 1:1:size(Y_delta,1)
    if (mod(i,size(Y_delta,1)/4) == 0)
    fprintf('%1.0f percent..\n',100*i/size(Y_delta,1))
    end
    % Gets closest nominal location
    displaces = Y_xyz(1:size(Y_xyz,1)/numOrbs+10,1:3)'-Y_delta(i,1:3)'; 
    displaces = sqrt(displaces(1,:).^2+displaces(2,:).^2+displaces(3,:).^2); 
    [a, b] = min(displaces);

    delta_xyz(i,:) = Y_xyz(b,1:3) - Y_delta(i,1:3); 
    delta_dxyz(i,:) = Y_xyz(b,4:6) - Y_delta(i,4:6);

end
x_orbits = 0:numOrbs/(length(delta_xyz)):numOrbs; 
fprintf('Finished processing deviation data!\n')


%% Second Round of Plotting
figure(3) % x,y,z Deviations from Nominal (position) 
plot(x_orbits(:,1:end-1), au*delta_xyz(1:end,:), 'LineWidth', 1) 
legend(['x'; 'y'; 'z']); xlim([0 numOrbs]); box off
xlabel('Orbit #'); ylabel('Distance [km]'); title('Deviation from Nominal') 
%
figure(4) % x,y,z Deviations from Nominal (velocity) 
plot(x_orbits(:,1:end-1), au*delta_dxyz(1:end,:)*(-ne)*1000, 'LineWidth', 1) 
legend(['dx'; 'dy'; 'dz']); xlim([0 numOrbs]); box off
xlabel('Orbit #'); ylabel('Velocity [km/s]'); title('Deviation from Nominal')

% Axis Stuff
xlmax = max(Y_xyz(:,1)); xlmin = min(Y_xyz(:,1)); 
ylmax = max(Y_xyz(:,2)); ylmin = min(Y_xyz(:,2)); 
zlmax = max(Y_xyz(:,3)); zlmin = min(Y_xyz(:,3));
xlw = 1.25*(xlmax - xlmin)/2; xlm = (xlmax + xlmin)/2; 
ylw = 1.25*(ylmax - ylmin)/2; ylm = (ylmax + ylmin)/2; 
zlw = 1.25*(zlmax - zlmin)/2; zlm = (zlmax + zlmin)/2;

figure(5) % Trajectory Plots
subplot(1,3,1) % Nominal
plot3(Y_xyz(:,1), Y_xyz(:,2), Y_xyz(:,3), 'Color', [0.4660 0.6740 0.1880], ...
        'LineWidth', 0.9) 
hold on; plot3(au*L1, 0, 0); title('Nominal'); axis square
xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]); 
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y')
set(gca,'ZTickLabel',[]); zlabel('Z')

subplot(1,3,2) % Uncorrected Trajectory
plot3(Y_delta_nc(:,1), Y_delta_nc(:,2), Y_delta_nc(:,3), 'Color', ...
        [0.9290 0.6940 0.1250], 'LineWidth', 0.9)
hold on; plot3(au*L1, 0, 0); title('Not Corrected'); axis square
xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]); 
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y')
set(gca,'ZTickLabel',[]); zlabel('Z')

subplot(1,3,3) % Corrected Trajectory
plot3(Y_delta(:,1), Y_delta(:,2), Y_delta(:,3), 'k', 'LineWidth', 0.9) 
hold on; plot3(au*L1, 0, 0)
xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]); 
title('Corrected'); hold on; axis square 
plot3(burnSpot(:,1),burnSpot(:,2),burnSpot(:,3), 'ro') 
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y')
set(gca,'ZTickLabel',[]); zlabel('Z')
% -- end subplots Figure 5

figure(6) % Nominal and Corrected Side-by-Side 
subplot(1,2,1) % Nominal Trajectory
plot3(Y_xyz(:,1), Y_xyz(:,2), Y_xyz(:,3), 'Color', [0.4660 0.6740 0.1880], ...
    'LineWidth', 0.9) 
hold on; plot3(au*L1, 0, 0); title('Nominal')
xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]); 
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y') 
set(gca,'ZTickLabel',[]); zlabel('Z')

subplot(1,2,2) % Corrected Trajectory w/ Course Correction Burns 
plot3(Y_delta(:,1), Y_delta(:,2), Y_delta(:,3), 'k', 'LineWidth', 0.9) 
hold on; plot3(au*L1, 0, 0); title('Corrected'); hold on;
xlim([xlm-xlw xlm+xlw]); ylim([ylm-ylw ylm+ylw]); zlim([zlm-zlw zlm+zlw]); 
plot3(burnSpot(:,1),burnSpot(:,2),burnSpot(:,3), 'ro', 'LineWidth', 2) 
set(gca,'YTickLabel',[]); xlabel('X')
set(gca,'XTickLabel',[]); ylabel('Y')
set(gca,'ZTickLabel',[]); zlabel('Z')
% -- end subplots Figure 6
fprintf('Plotting Complete!\n')

%% Burn Profile
clear DVStoreP DeltaVval
x_days = x_orbits*(t_end/-ne/24/60/60);

for i = size(DeltaVStore,1):1:size(x_orbits,2) 
    DeltaVStore(i,1:6) = zeros(1,6); 
    DeltaVStoreP(i,1:6) = zeros(1,6);
end
DVStoreP = sqrt(DeltaVStoreP(:,4).^2+DeltaVStoreP(:,5).^2+DeltaVStoreP(:,6).^2); 
counter = 1;
% Iterates through and determines location of actual burns vs empty slots
% in the array.
for i = 1:1:size(x_orbits,2)
    if (DVStoreP(i) == 0)
        DVStoreP(i) = -au;
    end
    if (DVStoreP(i) > 0)
        DeltaVval(counter,1) = DVStoreP(i)*au*abs(ne)*1000;
        DeltaVval(counter,2) = (i/size(x_orbits,2)*numOrbs)*(t_end/-ne/24/60/60);
        counter = counter + 1;
    end 
end
for i = 1:1:size(DeltaVval,1)-1
    dtburn(i) = DeltaVval(i+1,2)-DeltaVval(i,2);
end

% Plots the burn profile.
figure(7)
plot(x_days,DVStoreP*au*abs(ne)*1000,'ro', 'LineWidth', 2); 
hold on
ylim([0 max(DVStoreP)*au*abs(ne)*1000*1.25]) 
plot(DeltaVval(1:end-0,2),DeltaVval(1:end-0,1),'r') 
xlabel('Time [days]')
ylabel('\DeltaV [m/s]')
title('\DeltaV Station Keeping Maneuvors Over Time')
box off

% Calculate total delta V and fuel required.
delVtotal = sum(DeltaVval(:,1)); m0 = 1850;
Isp = 175;
mf = m0/exp(delVtotal/Isp/(9.81)); dm = m0-mf;
fprintf('Total delta V: %1.2f m/s\n', delVtotal) 
fprintf('Fuel used: %1.2f kg\n', dm)
fprintf('Mean coast time: %1.2f days\n', mean(dtburn))
%% Duration Between Burns
tdays = [23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53];
dVdays = [193.4, 53.7, 35.1, 7.92, 7.1, 4.4, 3.9, 5.8, 6.2, 17.65, 57.7];
meandays = [60.5, 60.6, 61.4, 62.2, 66.2, 73.1, 78.7, 75.6, 77.2, 71.4, 68.1];

% Plots duration between burns and total delta V for t_t range. (Hardcoded 
% with results from using above code).
figure(8)
plot(tdays, dVdays/max(dVdays), 'r', 'LineWidth', 2); hold on
plot(tdays, meandays/max(meandays), 'bo', 'LineWidth', 3)
box off; legend('Total Delta V','Mean Time Between Burns') 
title('Normalized Total Delta-V & Mean Time Between Burns vs t_t') 
xlabel('t_t [days]'); ylabel('Normalized')

fprintf('End.\n')
