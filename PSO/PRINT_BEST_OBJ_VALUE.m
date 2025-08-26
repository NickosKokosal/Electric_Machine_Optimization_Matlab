X_best = [0.5, 16.101585, 15.877706, 0.565706, 3.469779];

% Run of the Objactive Function Value
[opti, eff, M] = ObjFunc_peirama(X_best);

% PRINT OF THE RESULTS
fprintf('BEST POINT\n');
fprintf('Betad (X1): %.4f\n', X_best(1));
fprintf('ThickyokeR (X2): %.4f\n', X_best(2));
fprintf('ThickyokeS (X3): %.4f\n', X_best(3));
fprintf('Beta (X4): %.4f\n', X_best(4));
fprintf('Thickmagnet (X5): %.4f\n', X_best(5));
fprintf('-------------------------\n');
fprintf('Efficiency: %.2f %%\n', eff*100);
fprintf('Mass: %.3f kg\n', M);
fprintf('Objective Function Value: %.6f\n', opti);