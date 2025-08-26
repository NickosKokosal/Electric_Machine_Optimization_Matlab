eff_all = zeros(25,1);
mass_all = zeros(25,1);

for i = 1:25
    fname = sprintf('SO Results after iteration # %d.mat', i);
    data = load(fname);
    X = data.Swarm.GBEST.X;
    [~, eff, M] = ObjFunc(X);
    eff_all(i) = eff;
    mass_all(i) = M;
    fprintf('Iteration %2d: Eff = %.2f%% | Mass = %.2f kg\n', i, eff*100, M);
end
