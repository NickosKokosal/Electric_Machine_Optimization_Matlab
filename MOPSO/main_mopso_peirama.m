clear; close all; clc;
% Problem preparation 
 problem.nVar = 5;             % Number of search variables
 problem.ub   = [0.5, 30, 30, 0.8, 4]; % Upper bounds of search variables
 problem.lb   = [0.2, 10, 10, 0.4, 1]; % Lower bounds of search variables
 problem.fobj = @MObjFunc;      % Objective function script
% PSO parameters 
 PopNo        = 25; % Population number (i.e. no. of particles in each swarm)
 IterNo       = 25; % Iteration number
 RunNo        = 1; % Number of independed runs (repetitions of experiment)
 visFlag      = 1; % Set this to 0 if you do not want visualization
 BestSolutions_MOPSO = zeros(2 , RunNo); % Initialization of solution vector
 [GBEST , cgcurve] = MOPSO( PopNo, IterNo, problem , RunNo,visFlag );