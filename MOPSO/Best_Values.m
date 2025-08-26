X = [0.500, 19.4066, 18.5599, 0.6104, 2.8116]; % Best Values

[f1, f2] = MObjFunc_peirama(X);  % Run FEMM design

disp(['Efficiency (reversed): ', num2str(f1)])
disp(['Mass: ', num2str(f2)])
