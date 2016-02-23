% Test find_roots with steady-piston rear boundary condition.
% Case is neutral ZND solution.
close all; clear all; clc

% Physical free parameters.
q = 1.7; theta = 2.35;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% N - grid resolution; lambda_tol - closeness of lambda to equilibrium value;
% M - length of the computational domain.
N = 10000;
lambda_tol = 1.0 / N;
M = compute_x_at_point(1-lambda_tol, params);

grid = linspace(0, M, N);
znd_all = compute_znd(grid, params);

tol = 1e-12; % Tolerance for optimization.
guess = [-9.0e-04; 6.8e-01];
result = find_roots_steady_piston(guess, grid, znd_all, params, tol);

disp(result);

root = result.root;
alpha_c = root(1) + 1j*root(2);
sol = compute_linearized_problem(alpha_c, grid, znd_all, params);

save('results/2016-02-18-steady-piston-bc/neutral.mat');

figure
plot(grid, sol(1, :))
hold on
plot(grid, sol(3, :))
xlabel('x')
ylabel('Real part of perturbations')
legend('Re u', 'Re \lambda', 'Location', 'southwest');
export_fig_in_pdf( ...
    'results/2016-02-18-steady-piston-bc/neutral_pert_real.pdf', [6 3.7]);

figure
plot(grid, sol(2, :))
hold on
plot(grid, sol(4, :))
xlabel('x')
ylabel('Imaginary part of perturbations')
legend('Im u', 'Im \lambda', 'Location', 'southwest');
export_fig_in_pdf( ...
    'results/2016-02-18-steady-piston-bc/neutral_pert_imag.pdf', [6 3.7]);