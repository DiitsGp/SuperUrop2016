%% Initialize
close all;

A = [1 0 1;...
    0 1 1];
[m, n] = size(A);

rng(0);
Costs = rand(3, 2); % linear costs on 3 edges

F0 = rand(m, 1);
F0 = F0/sum(F0);
F = F0;
C = zeros(n, 1);

T = 10^4;
eta = 1/sqrt(T);

F_evol = zeros(m, T);
C_evol = zeros(n, T);

%% Simulate Centralized Social Optimum
F = F0;

for i = 1:T
   F_evol(:, i) = F;
   G = A'*F;
   
   C_evol(:, i) = C;
   grad_edges = 2*Costs(:, 1) .* G + Costs(:, 2);
   grad_phi = A*grad_edges;
   F_tilde = F - eta*grad_phi;
   
   F = projsplx(F_tilde);
end

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium Centralized Social Optimum');
legend('flow along route 1', 'route 2');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes Centralized Social Optimum');
legend('route 1', 'route 2');

%% Distributed (trivially the same in this case?)
