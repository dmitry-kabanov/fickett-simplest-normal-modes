function [params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess)
%SOLVE_EIGENVALUE_PROBLEM Summary of this function goes here
%   Detailed explanation goes here

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% lambda_tol - closeness of lambda to equilibrium value;
% M - length of the computational domain.
lambda_tol = 1.0 / N;
M = compute_x_at_point(1-lambda_tol, params);

grid = linspace(0, M, N);
znd_all = compute_znd(grid, params);

tol = 1e-12; % Tolerance for optimization.
result = find_roots_steady_piston(guess, grid, znd_all, params, tol);

disp(result);

root = result.root;
alpha_c = root(1) + 1j*root(2);
pert = compute_linearized_problem(alpha_c, grid, znd_all, params);
end

