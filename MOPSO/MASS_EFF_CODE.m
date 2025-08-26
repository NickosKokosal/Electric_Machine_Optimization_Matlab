eff_all = zeros(25,1);
mass_all = zeros(25,1);

for i = 1:25
    filename = sprintf('MO Results after iteration # %d.mat', i);  
    data = load(filename);
    X = data.Swarm.GBEST.X;

    try
        [~, eff, M] = MObjFunc_peirama(X);  % Εκτελεί το FEMM
    catch
        fprintf('FEMM failed on iteration %d\n', i);
        eff = 0;
        M = Inf;
    end

    eff_all(i) = eff;
    mass_all(i) = M;

    fprintf('Iter %2d: Eff = %.2f %% | Mass = %.3f kg\n', i, eff*100, M);
end