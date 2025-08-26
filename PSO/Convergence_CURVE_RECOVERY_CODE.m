% recover_convergence_curve.m

N = 25;  % Αριθμός iterations (προσαρμόστε ανάλογα με τα αποθηκευμένα .mat αρχεία)
convergence_curve = zeros(N, 1);

for i = 1:N
    filename = sprintf('SO Results after iteration # %d.mat', i);
    data = load(filename);

    if isfield(data, 'Swarm') && isfield(data.Swarm, 'GBEST') && isfield(data.Swarm.GBEST, 'O')
        convergence_curve(i) = data.Swarm.GBEST.O;
    else
        warning('Δεν βρέθηκε objective function στην iteration %d', i);
        convergence_curve(i) = NaN;
    end
end

% Αφαίρεση NaNs για ασφάλεια
valid_idx = ~isnan(convergence_curve);

% Πλοτ
figure;
plot(find(valid_idx), convergence_curve(valid_idx), 'b-o', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Objective Function Value');
title('Recovered Convergence Curve');
grid on;

% Αποθήκευση ως εικόνα PNG
saveas(gcf, 'Recovered_Convergence_Curve.png');

fprintf('✅ Το convergence curve αποθηκεύτηκε ως "Recovered_Convergence_Curve.png"\n');