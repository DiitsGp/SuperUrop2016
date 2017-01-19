%% Initialize
close all;
clear;

A = [1 0 1 1 0 0; ...
    0 1 1 1 0 0; ...
    0 0 0 0 1 1];
[m, n] = size(A);

%rng(30);
Costs = rand(n, 2); % linear costs on n edges

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
   if i == 1
       display('grad_phi');
       display(grad_phi);
   end
   F_tilde = F - eta*grad_phi;
   
   F = projsplx(F_tilde);
end

display('Final Flow Vector Centralized');
display(F);

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium Centralized Social Optimum');
legend('flow along route 1', 'route 2', 'route 3');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes Centralized Social Optimum');
legend('route 1', 'route 2', 'route 3');

%% Distributed
F = F0;

Axs = [1 0 1 0 0 0; ...
    0 1 1 0 0 0; ...
    0 0 0 0 1 1]; % chi hops

Axs1 = [1 0 0 0 0 0; ...
    0 1 0 0 0 0; ...
    0 0 0 0 1 0]; % one hop

Axn1 = [1 0 1 1 0 0; ...
    0 1 1 1 0 0; ...
    0 0 0 0 1 0]; % chi hops

Axn11 = [1 0 1 0 0 0; ...
    0 1 1 0 0 0; ...
    0 0 0 0 0 0]; % one hop

Axn2 = [1 0 1 1 0 0; ...
    0 1 1 1 0 0; ...
    0 0 0 0 0 1]; % chi hops

Axn21 = [0 0 1 1 0 0; ...
    0 0 1 1 0 0; ...
    0 0 0 0 0 0]; % one hop

Axn3 = [1 0 0 1 0 0; ...
    0 1 0 1 0 0; ...
    0 0 0 0 1 1]; % chi hops

Axn31 = [0 0 0 0 0 0; ...
    0 0 0 0 0 0; ...
    0 0 0 0 1 1]; % one hop

Axd = [0 0 1 1 0 0; ...
    0 0 1 1 0 0; ...
    0 0 0 0 1 1]; % chi hops

Axd1 = [0 0 0 1 0 0; ...
    0 0 0 1 0 0; ...
    0 0 0 0 0 1]; % one hop

for i = 1:T
   F_evol(:, i) = F;
   G = A'*F;
   
   C_evol(:, i) = C;
   grad_edges = 2*Costs(:, 1) .* G + Costs(:, 2);
   grad_phi = A*grad_edges;
   F_tilde = F - eta*grad_phi;
   
   grad_phi_s = Axs*grad_edges;
   grad_phi_n1 = Axn1*grad_edges;
   grad_phi_n2 = Axn2*grad_edges;
   grad_phi_n3 = Axn3*grad_edges;
   grad_phi_d = Axd*grad_edges;
   
   grad_phi_s1 = Axs1*grad_edges;
   grad_phi_n11 = Axn11*grad_edges;
   grad_phi_n21 = Axn21*grad_edges;
   grad_phi_n31 = Axn31*grad_edges;
   grad_phi_d1 = Axd1*grad_edges;
   
   grad_sum = grad_phi_s + grad_phi_n1 + grad_phi_n2 + grad_phi_n3 + grad_phi_d;
   %grad_sum = grad_phi_s1 + grad_phi_n11 + grad_phi_n21 + grad_phi_n31 + grad_phi_d1;
   
   F = F - eta*grad_sum;
   F = projsplx(F);
end

display('Final Flow Vector Distributed');
display(F);

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium Distributed Social Optimum');
legend('flow along route 1', 'route 2', 'route 3');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes Distributed Social Optimum');
legend('route 1', 'route 2', 'route 3');