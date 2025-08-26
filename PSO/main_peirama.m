clear; close all; clc;
% Problem preparation 
 problem.nVar = 5;             % Number of search variables
 problem.ub   = [0.55, 24, 28, 0.80, 5.0]; % Upper bounds of search variables
 problem.lb   = [0.30, 8, 10, 0.55, 1.5]; % Lower bounds of search variables
 problem.fobj = @ObjFunc;      % Objective function script
% PSO parameters 
 PopNo        = 25; % Population number (i.e. no. of particles in each swarm)
 IterNo       = 25; % Iteration number
 RunNo        = 1 ; % Number of independed runs (repetitions of experiment)
 visFlag      = 1; % Set this to 0 if you do not want visualization
 BestSolutions_PSO = zeros(1 , RunNo); % Initialization of solution vector
 [GBEST , cgcurve] = PSO( PopNo, IterNo, problem , RunNo,visFlag );
 