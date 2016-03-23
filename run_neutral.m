% Test find_roots with steady-piston rear boundary condition.
% Case is neutral ZND solution.
close all; clear all; clc

resultdir = 'results/2016-03-23-new-radiation-condition-tol=1e-4/';
matfile = strcat(resultdir, 'neutral.mat');
picfile = strcat(resultdir, 'neutral.pdf');

% Physical free parameters.
q = 1.7; theta = 2.35;

% Grid resolution.
N = 10000;

% Guess for eigenvalue.
guess.alpha_re = 0.00;
guess.alpha_im = 0.00;

[params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess);
save(matfile);

% Plotting part.
plot_znd_and_perturbations(matfile, picfile);