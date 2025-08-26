i = 4;
filename = sprintf('SO Results after iteration # %d.mat', i);
data = load(filename);
X = data.Swarm.GBEST.X;

disp('X που σπάει στο iteration 4:');
disp(X);