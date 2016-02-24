% Test find_roots with steady-piston rear boundary condition.
% Case is neutral ZND solution.
close all; clear all; clc

resultdir = 'results/2016-02-18-steady-piston-bc/';
matfile = strcat(resultdir, 'neutral.mat');

% Physical free parameters.
q = 1.7; theta = 2.35;

% Grid resolution.
N = 10000;

% Guess for eigenvalue.
guess = [-9.0e-04; 6.8e-01];

[params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess);
save(matfile);

% Plotting part.
plot_znd_and_perturbations(grid, znd_all, pert);
export_fig_in_pdf(strcat(resultdir, 'neutral.pdf'), [4.5 8.3436]);