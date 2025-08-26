clear; close all; clc;
% Problem preparation 
 problem.nVar = 3;             % Number of search variables
 problem.ub   = [0.5, 30, 30]; % Upper bounds of search variables
 problem.lb   = [0.2, 10, 10]; % Lower bounds of search variables
 problem.fobj = @ObjFunc;      % Objective function script
% PSO parameters 
 PopNo        = 50; % Population number (i.e. no. of particles in each swarm)
 IterNo       = 60; % Iteration number
 RunNo        = 5; % Number of independed runs (repetitions of experiment)
 visFlag      = 1; % Set this to 0 if you do not want visualization
 BestSolutions_PSO = zeros(1 , RunNo); % Initialization of solution vector
 [GBEST , cgcurve] = PSO( PopNo, IterNo, problem , RunNo,visFlag );