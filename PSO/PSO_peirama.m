function [GBEST , cgCurve ] = PSO ( PopulationNo, IterationsNo,  problem, repetitions, Vis )
format shortG
% Define the details of the objective function
nVar = problem.nVar;
ub   = problem.ub;
lb   = problem.lb;
fobj = problem.fobj;

cgCurve = zeros(repetitions, IterationsNo);
FirstP_D1 = zeros(repetitions , IterationsNo);
Solution_found=zeros(1, repetitions);
position_history = zeros(PopulationNo , IterationsNo , nVar );
tt=0;

for repNo= 1:repetitions
    wMax = 0.9;
    wMin = 0.2;
    c1 = 2;
    c2 = 2;
    vMax = (ub - lb) .* 0.2;
    vMin  = -vMax;

    for swarm = 1 : PopulationNo
        Swarm.Particles(swarm).X = (ub-lb) .* rand(1,nVar) + lb;
        Swarm.Particles(swarm).V = zeros(1, nVar);
        Swarm.Particles(swarm).PBEST.X = zeros(1,nVar);
        Swarm.Particles(swarm).PBEST.O = inf;
        Swarm.GBEST.X = zeros(1,nVar);
        Swarm.GBEST.O = inf;
    end

    for iter = 1 : IterationsNo
        tic
        for swarm = 1 : PopulationNo
            fprintf('Run: %d/%d, Iteration: %d/%d, Swarm: %d/%d \n',...
                repNo,repetitions,iter,IterationsNo,swarm,PopulationNo);
            currentX = Swarm.Particles(swarm).X;
            position_history(swarm , iter , : ) = currentX;

            [Swarm.Particles(swarm).O, eff_tmp, Mtotal_tmp] = fobj(currentX);
            if Swarm.Particles(swarm).O < Swarm.Particles(swarm).PBEST.O
                Swarm.Particles(swarm).PBEST.X = currentX;
                Swarm.Particles(swarm).PBEST.O = Swarm.Particles(swarm).O;
                
            end
            if Swarm.Particles(swarm).O < Swarm.GBEST.O
                Swarm.GBEST.X = currentX;
                Swarm.GBEST.O = Swarm.Particles(swarm).O;
            end
        end

        w = wMax - iter .* ((wMax - wMin) / IterationsNo);
        FirstP_D1(repNo,iter) = Swarm.Particles(1).X(1);

        for swarm = 1 : PopulationNo
            Swarm.Particles(swarm).V = w .* Swarm.Particles(swarm).V + c1 .* rand(1,nVar) .* (Swarm.Particles(swarm).PBEST.X - Swarm.Particles(swarm).X) ...
                + c2 .* rand(1,nVar) .* (Swarm.GBEST.X - Swarm.Particles(swarm).X);

            index1 = find(Swarm.Particles(swarm).V > vMax);
            index2 = find(Swarm.Particles(swarm).V < vMin);

            Swarm.Particles(swarm).V(index1) = vMax(index1);
            Swarm.Particles(swarm).V(index2) = vMin(index2);

            Swarm.Particles(swarm).X = Swarm.Particles(swarm).X + Swarm.Particles(swarm).V;

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
            fprintf('Motor Efficiency: %.2f %%\n', eff*100);
            fprintf('Motor Weight    : %.4f kg\n', Mtotal);
            fprintf('--------------------------------------------\n');
        end

        cgCurve(repNo,iter,:) = Swarm.GBEST.O;
        
        fileName = sprintf('SO Results Run%d Iteration #%d.mat', repNo, iter);
        GBEST_current = Swarm.GBEST;
        save(fileName, 'GBEST_current');

        et = toc; tt = tt + et;
        fprintf('Iteration Time   : %9.2f sec \n',et);
        fprintf('Total time so far: %9.2f min\n',tt/60);
        fprintf('--------------------------------------------\n');
    end

    GBEST(repNo) = Swarm.GBEST;
    opti = Swarm.GBEST.O;
    BestResults(repNo,:) = [repNo, opti];
end

save('Eff_Mass_History.mat', 'eff_history', 'mass_history');

if Vis == 1
    iterations = 1: IterationsNo;
    semilogy(iterations,cgCurve(repNo,:));
    grid on; hold on;
    title('Convergence curve(s)');
    xlabel('Iterations');
    ylabel('Objective Function Value');
    saveas(gcf,'Convergence_Curves.png')
end

end
