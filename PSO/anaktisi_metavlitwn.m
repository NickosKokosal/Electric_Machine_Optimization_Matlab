load('SO_results.mat');  % Φόρτωσε το αρχείο

% Έλεγχος τι υπάρχει
whos  % Δείξε τι περιλαμβάνει το .mat αρχείο

% Αν υπάρχουν πίνακες:
% - BestCost (Nx1) → τιμή objective function κάθε iteration
% - GlobalBestPosition (Nx5) → μεταβλητές κάθε iteration

% Πίνακας αποτελεσμάτων
n_iter = length(BestCost);
tableResults = zeros(n_iter, 6); % [ObjF, X1, X2, X3, X4, X5]

for i = 1:n_iter
    tableResults(i,1) = BestCost(i);                   % ObjF
    tableResults(i,2:6) = GlobalBestPosition(i,:);     % Variables 1-5
end

% Προβολή
disp(array2table(tableResults, ...
    'VariableNames', {'ObjF', 'X1', 'X2', 'X3', 'X4', 'X5'}))