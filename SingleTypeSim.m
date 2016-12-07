%% SingleTypeSim.m
% For this code I will assume the following:
% There are m routes from a source S to a destination D, with flows denoted
%    by f_1 ... f_m
% There are n edges from the same source S to the same destination D,
%   with flows denoted by g_1 ... g_n and cost functions c_1(g_1) ...
%   c_n(g_n)
% Let F be the vector [f_1 ... f_m]'
% Let G be the vector [g_1 ... g_n]'
% Let C be the vector [c_1(g_1) ... c_n(g_n)]'
% So that F, G, C are vector-functions of time
% Let A be the adjency matrix of size nxm so that:
%   A_ij = 1 iff edge i belongs to route j and A_ij = 0 otherwise
% I will simulate the dynamics of the system from t=1 to t=T for the
%   example given in Marden et al.'s 'Joint Strategy Fictitious Play With
%   Inertia for Potential Games'

%% Initialize
close all;
rng(0);

% In the Marden example there are 10 edges, each corresponding to a unique
%   route. Each route contains exactly one edge.
A = eye(10);
[m, n] = size(A);

% Generate a random initial condition with one unit of flow in total
F0 = rand(m, 1);
F0 = F0/sum(F0);
F = F0;
F_tilde = zeros(m, 1);

% Generate the coefficients for quadratic cost functions as in the Marden
%   et al. paper
C_coeff = rand(n, 3); %a_i x^2 + b_i x + c_i
C = zeros(n, 1);

% Simulate for T = 10^3 timesteps
T = 10^4;
eta = 1/sqrt(T);
F_evol = zeros(m, T); % for plotting
C_evol = zeros(n, T);

%% Simulate the standard potential function
for i = 1:T
   F_evol(:, i) = F;
   G = A'*F;
   for j = 1:n
      c = C_coeff(j, :);
      g = G(j);
      C(j) = c(1)*g^2 + c(2)*g + c(3);
   end
   C_evol(:, i) = C;
   grad_phi = A*C;
   F_tilde = F - eta*grad_phi;
   
   F = projsplx(F_tilde);
end

F_standard = F;
route_cost_standard = A*C;

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium of the Cost Minimizing Potential Function on Marden Example');
legend('flow along route 1', 'route 2', 'route 3', 'route 4', 'route 5', 'route 6', 'route 7', 'route 8', 'route 9', 'route 10');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes of the Standard Potential Function on Marden Example');
legend('route 1', 'route 2', 'route 3', 'route 4', 'route 5', 'route 6', 'route 7', 'route 8', 'route 9', 'route 10');

%% Simulate the potential function given in Smith
F = F0;

for i = 1:T
   F_evol(:, i) = F;
   G = A'*F;
      for j = 1:n
      c = C_coeff(j, :);
      g = G(j);
      C(j) = c(1)*g^2 + c(2)*g + c(3);
   end
   C_evol(:, i) = C;
   route_cost = A*C;
   
   grad_phi = zeros(m, 1);
   for j = 1:m
       for k = 1:m
           if C(j) > C(k)
               delta_jk = zeros(m, 1);
               delta_jk(j) = -1; delta_jk(k) = 1;
               grad_phi = grad_phi + (F(j)*(C(j)-C(k)))*delta_jk;
           end
       end
   end
   F_tilde = F + eta*grad_phi;
   
   F = projsplx(F_tilde);
end

F_smith = F;
route_cost_smith = route_cost;

figure
plot(1:1:T, F_evol);
title('Traffic Flow to Equilibrium of the Load Balancing Potential Function on Marden Example');
legend('flow along route 1', 'route 2', 'route 3', 'route 4', 'route 5', 'route 6', 'route 7', 'route 8', 'route 9', 'route 10');

figure
plot(1:1:T, A*C_evol);
title('Costs Along Routes of the Smith Potential Function on Marden Example');
legend('route 1', 'route 2', 'route 3', 'route 4', 'route 5', 'route 6', 'route 7', 'route 8', 'route 9', 'route 10');