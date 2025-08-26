function [GBEST , cgCurve ] = MOPSO ( PopulationNo, IterationsNo,  problem, repetitions, Vis )
format shortG
% Define the details of the objective function
nVar = problem.nVar;
ub   = problem.ub;
lb   = problem.lb;
fobj = problem.fobj;
% Extra variables for data visualization
%%% average_objective = pagetranspose(zeros(repetitions, IterationsNo,2));
% average_objective = zeros(repetitions, IterationsNo);
average_objective(1,1,:)=[0 0];
%cgCurve = zeros(repetitions, IterationsNo);
%CC=[];
 for zr=1:repetitions 
  for zt=1:IterationsNo
  cgCurve(zr,zt,:) =[0 0];
  end
 end
FirstP_D1 = zeros(repetitions , IterationsNo);
Solution_found=zeros(1, repetitions);
position_history = zeros(PopulationNo , IterationsNo , nVar );
tt=0;
% Define the PSO's paramters
for repNo= 1:repetitions
wMax = 0.9;
wMin = 0.2;
c1 = 2;
c2 = 2;
vMax = (ub - lb) .* 0.2;
vMin  = -vMax;
% The PSO algorithm
% Initialize the particles
for swarm = 1 : PopulationNo
    Swarm.Particles(swarm).X = (ub-lb) .* rand(1,nVar) + lb;
    Swarm.Particles(swarm).V = zeros(1, nVar);
    Swarm.Particles(swarm).PBEST.X = zeros(1,nVar);
    Swarm.Particles(swarm).PBEST.O = inf;
    Swarm.GBEST.X = zeros(1,nVar);
    Swarm.GBEST.O = inf;
end
% Main loop
for iter = 1 : IterationsNo
tic    
    % Calcualte the objective value
    for swarm = 1 : PopulationNo
    fprintf('Run: %d/%d, Iteration: %d/%d, Swarm: %d/%d \n',...
            repNo,repetitions,iter,IterationsNo,swarm,PopulationNo);
        currentX = Swarm.Particles(swarm).X;
        position_history(swarm , iter , : ) = currentX;
        
        Swarm.Particles(swarm).O = fobj(currentX);
%%%        average_objective(repNo,iter,:) =  average_objective(repNo,iter,:)  + Swarm.Particles(swarm).O;
        % Update the PBEST
        if Swarm.Particles(swarm).O < Swarm.Particles(swarm).PBEST.O
            Swarm.Particles(swarm).PBEST.X = currentX;
            Swarm.Particles(swarm).PBEST.O = Swarm.Particles(swarm).O;
        end
        % Update the GBEST
        if Swarm.Particles(swarm).O < Swarm.GBEST.O
            Swarm.GBEST.X = currentX;
            Swarm.GBEST.O = Swarm.Particles(swarm).O;
        end
    end
    % Update the X and V vectors
    w = wMax - iter .* ((wMax - wMin) / IterationsNo);
    FirstP_D1(repNo,iter) = Swarm.Particles(1).X(1);
    
    for swarm = 1 : PopulationNo
        Swarm.Particles(swarm).V = w .* Swarm.Particles(swarm).V + c1 .* rand(1,nVar) .* (Swarm.Particles(swarm).PBEST.X - Swarm.Particles(swarm).X) ...
            + c2 .* rand(1,nVar) .* (Swarm.GBEST.X - Swarm.Particles(swarm).X);
        % Check velocities
        index1 = find(Swarm.Particles(swarm).V > vMax);
        index2 = find(Swarm.Particles(swarm).V < vMin);
        
        Swarm.Particles(swarm).V(index1) = vMax(index1);
        Swarm.Particles(swarm).V(index2) = vMin(index2);
        
        Swarm.Particles(swarm).X = Swarm.Particles(swarm).X + Swarm.Particles(swarm).V;
        % Check positions
        index1 = find(Swarm.Particles(swarm).X > ub);
        index2 = find(Swarm.Particles(swarm).X < lb);
        
        Swarm.Particles(swarm).X(index1) = ub(index1);
        Swarm.Particles(swarm).X(index2) = lb(index2);
    end
    
    if Vis == 1
     fprintf('--------------------------------------------\n');
     fprintf('In Iteration : %d \n',iter);
       fprintf('Best ObjF value so far           : %f   \n',Swarm.GBEST.O);
     for q=1:nVar
       fprintf('Best value for variable #%d so far: %f \n',q,Swarm.GBEST.X(q));
     end
     fprintf('--------------------------------------------\n');
    end
    
%   CC(repNo,iter,:)=[CC;Swarm.GBEST.O];
    cgCurve(repNo,iter,:)=Swarm.GBEST.O;
%   cgCurve(repNo,iter) = Swarm.GBEST.O;
%%%    average_objective(repNo,iter) = average_objective(repNo,iter) / PopulationNo;
    
    fileName = ['MO Results after iteration # ',num2str(iter)]; save(fileName);

  et=toc; tt=tt+et;
  fprintf('Iteration Time   : %9.2f sec \n',et);
  fprintf('Total time so far: %9.2f min\n',tt/60);
  fprintf('--------------------------------------------\n');  
end

GBEST(repNo) = Swarm.GBEST;

%if Vis == 1
%    iterations = 1: IterationsNo;
%    semilogy(iterations,cgCurve(repNo,:));
%    grid on; hold on;
%    title('Convergence curve(s)');
%    xlabel('Iterations');
%    ylabel('Objective Function Value');
%end

%-----------------------visualization
if Vis==1
 for zr=1:repetitions
  f1=cgCurve(zr,:,1);
  f2=cgCurve(zr,:,2);
  f1=(1./f1)*100;
  f2=(1./f2)*1000;
  plot(f1,f2);hold on;
  set(gca, 'XDir','reverse')
 end
 grid on;
 title('Pareto front(s)');
 xlabel('f1 - efficiency (%)');
 ylabel('f2 - mass (kg)');
end
%-----------------------visualization

end
saveas(gcf,'Pareto_Fronts.png')
% CALCULATION FOR PARETO 
if Vis==1
    ideal = [100, 0];  %ideal point 100% efficiency and 0 kg mass
    all_best_points = [];

    figure;
    hold on;

    for zr = 1:repetitions
        f1 = (1 ./ cgCurve(zr,:,1)) * 100;
        f2 = (1 ./ cgCurve(zr,:,2)) * 1000;

        distances = sqrt((f1 - ideal(1)).^2 + (f2 - ideal(2)).^2);
        [minDist, idx] = min(distances);

        best_eff = f1(idx);
        best_mass = f2(idx);

        all_best_points = [all_best_points; zr, best_eff, best_mass, minDist];

        plot(f1, f2, 'b'); hold on;
        plot(best_eff, best_mass, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
        set(gca, 'XDir','reverse')
    end

    grid on;
    title('Pareto front(s) with best Euclidean points');
    xlabel('f1 - efficiency (%)');
    ylabel('f2 - mass (kg)');
    legend('Pareto fronts','Best point per run')
    saveas(gcf, 'Pareto_with_best_points.png');

    % SAVE OF THE RESULTS
    T = array2table(all_best_points, ...
        'VariableNames', {'Run','Efficiency','Mass','EuclideanDist'});
    writetable(T, 'Best_Points_Euclidean.csv');
end

end