%% Particle Swarm Optimization using the particleswarm() Function 
% included with the Global Optimization Toolbox
clc;
clear;
close all;

%% Problem Definition
problem_setup;                  % Load Problem Variables

CostFunction = @(x) EnergyOptCost(x); % Cost Function

nVar = 7;
VarMin = [-100 -100 -100 -10 -10 -10  0];  % Lower Bounds of Decision Variables
VarMax = [ 100  100  100  10  10  10 10];  % Upper Bounds of Decision Variables

options = optimoptions('particleswarm','SwarmSize',250,'Display','iter');

[x,fval,exitflag,output] = particleswarm(CostFunction,nVar,VarMin,VarMax,options);