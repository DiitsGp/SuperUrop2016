%% Initialize
close all;
rng(0);

% Generate the Adjacency Matrix
A = [1 1 0 0;...
    0 0 1 1];
[m, n] = size(A);
A_BP = [1 1 0 0 0;...
    1 0 0 1 1;...
    0 0 1 1 0];
[mbp, nbp] = size(A_BP);

% Generate a random initial condition with one unit of flow in total
F0 = rand(m, 1);
F0 = F0/sum(F0);
F = F0;
F_tilde = zeros(m, 1);

% Set the cost matrices
C_coeff = [1 0;...
    0 1;...
    0 1;...
    1 0];
C = zeros(n, 1);
C_BP = zeros(nbp, 1);
C_BPcoeff = [1 0;...
    0 1;...
    0 1;...
    1 0;...
    0 0];
C_BP = zeros(nbp, 1);

% Simulate for T = 10^4 timesteps
T = 10^4;
eta = 1/sqrt(T);
F_evol = zeros(m, T); % for plotting
F_BPevol = zeros(mbp, T);
C_evol = zeros(n, T);
C_BPevol = zeros(nbp, T);

%% Simulate the standard potential function without the zero cost edge
for i = 1:T
   F_evol(:, i) = F;
   G = A'*F;
   for j = 1:n
      c = C_coeff(j, :);
      g = G(j);
      C(j) = c(1)*g + c(2);
   end
   C_evol(:, i) = C;
   grad_phi = A*C;
   F_tilde = F - eta*grad_phi;
   
   F = projsplx(F_tilde);
end

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium Without Zero Cost Edge');
legend('flow along route 1', 'route 2');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes Without Zero Cost Edge');
legend('route 1', 'route 2');

%% Simulate the standard potential function with the zero cost edge
rng(0);
F0 = rand(mbp, 1);
F0 = F0/sum(F0);
F = F0;

for i = 1:T
   F_BPevol(:, i) = F;
   G = A_BP'*F;
   for j = 1:nbp
      c = C_BPcoeff(j, :);
      g = G(j);
      C_BP(j) = c(1)*g + c(2);
   end
   C_BPevol(:, i) = C_BP;
   grad_phi = A_BP*C_BP;
   F_tilde = F - eta*grad_phi;
   
   F = projsplx(F_tilde);
end

figure
plot(1:1:T, F_BPevol);
title('Traffic Flow to Equilibrium With the Zero Cost Edge');
legend('flow along route 1', 'route 2', 'route 3');

figure
plot(1:1:T, A_BP*C_BPevol);
title('Costs Along Routes With the Zero Cost Edge');
legend('route 1', 'route 2', 'route 3');