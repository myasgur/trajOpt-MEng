%% Particle Swarm Optimization
clc;
clear;
close all;

%% Problem Definition
problem_setup;                  % Load Problem Variables

CostFunction = @(x) EnergyOptCost(x); % Cost Function

nVar = 7;                       % Number of Unknown (Decision) Variables

VarSize = [1 nVar];             % Matrix Size of Decision Variables

VarMin = zeros(VarSize); % Lower Bounds of Decision Variables
VarMax = ones(VarSize);  % Upper Bounds of Decision Variables

%% Parameters of PSO

MaxIt = 1000;   % Maximum Number of Iterations
MaxRuns = 100;  % Maximum Number of Runs

nPop = 20;      % Population (Swarm) Size

wmin = 0.4;     % Minimum Inertia Coefficient
wmax = 0.9;     % Maximum Inertia Coefficient

c1min = 0.5;    % Minimum Personal Acceleration Coefficient
c1max = 2.5;    % Maximum Personal Acceleration Coefficient

c2min = 0.5;    % Minimum Social/Global Acceleration Coefficient
c2max = 2.5;    % Maximum Social/Global Acceleration Coefficient

Vmax = 0.8;     % Maximum Velocity

%% Initialization

% Particle Template
empty_particle.Position = []; 
empty_particle.Velocity = []; 
empty_particle.Cost = []; 
empty_particle.Best.Position = [];
empty_particle.Best.Cost = [];


    % Create Population Array
    particle = repmat(empty_particle, nPop, 1);

    % Initialize Global Best
    GlobalBest.Cost = Inf;

    % Initialize Population Members
    for i = 1:nPop

        % Generate Random Solution
        for j = 1:nVar
            particle(i).Position(j) = unifrnd(VarMin(j),VarMax(j),1);
        end

        % Initialize Velocity
        particle(i).Velocity = zeros(VarSize);

        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);

        % Update Personal Best
        particle(i).Best.Position = particle(i).Position;
        particle(i).Best.Cost = particle(i).Cost;

        % Update Global Best
        if particle(i).Best.Cost < GlobalBest.Cost

            GlobalBest = particle(i).Best;

        end
    end

    % Array to Hold Best Cost Value on Each Iteration
    BestCosts = zeros(MaxIt,1);

    % Initialize PSO Parameters
    w = wmax;
    c1 = c1max;
    c2 = c2min;

    %% Main Loop of PSO
for k = 1:MaxRuns
    disp(['------------------ ' ...
        'Run ',num2str(k),' of ',num2str(MaxRuns), ...
        ' ------------------'])
    for it = 1:MaxIt
        for i = 1:nPop

            % Update Velocity
            particle(i).Velocity = w*particle(i).Velocity...
                + c1*rand(VarSize).*(particle(i).Best.Position - particle(i).Position)...
                + c2*rand(VarSize).*(GlobalBest.Position - particle(i).Position);

            % Ensure Velocity Is In Range [-Vmax Vmax]
            particle(i).Velocity = min(max(particle(i).Velocity,-Vmax),Vmax);

            % Update Position
            particle(i).Position = particle(i).Position + particle(i).Velocity;

            % Constrain Position
            for j = 1:nVar
                particle(i).Position(j) = min(max(particle(i).Position(j),VarMin(j)),VarMax(j));
            end

            % Evaluation
            particle(i).Cost = CostFunction(particle(i).Position);

            % Update Personal Best
            if particle(i).Cost < particle(i).Best.Cost

                particle(i).Best.Position = particle(i).Position;
                particle(i).Best.Cost = particle(i).Cost;

                % Update Global Best
                if particle(i).Best.Cost < GlobalBest.Cost

                    GlobalBest = particle(i).Best;

                end
            end


        end

        % Store the Best Cost on Every Iteration
        BestCosts(it) = GlobalBest.Cost;

        % Display Iteration Information
        disp(['Iteration ',num2str(it),': Best Cost = ',num2str(BestCosts(it),8)])

        if it>1 && BestCosts(it)<BestCosts(it-1)
            % Plot the Global Best Trajectory for the Current Iteration
            lams = psoSearchVars(GlobalBest.Position);
            y0 = [r0;v0;m0;lams]; % initial conditions
            opts = odeset('RelTol',1e-10,'AbsTol',1e-10); % integration options
            [t,y]=ode78(@cr3bp_EOM,[t0 tf],y0,opts);
            %         tspan = t0:0.0001:tf;
            %         y = ode4(@cr3bp_EOM_time,tspan,y0);
            plot_traj(y)
            pause(1)
        end

        % Update PSO Parameters
        w  =  wmax + (it/MaxIt)*( wmin -  wmax);
        c1 = c1max + (it/MaxIt)*(c1min - c1max);
        c2 = c2min + (it/MaxIt)*(c2max - c2min);

%         % Quitting Procedure
%         if it>100 && BestCosts(it)==BestCosts(it-40)
%             disp('Did not converge.')
%             break
%         end
    end

    % Reinitialize PSO Parameters
    w = wmax;
    c1 = c1max;
    c2 = c2min;

end


%% Results