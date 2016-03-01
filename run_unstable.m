% Test find_roots with steady-piston rear boundary condition.
% Case is unstable ZND solution.
close all; clear all; clc

resultdir = 'results/2016-03-01-towards-shock/';
matfile = strcat(resultdir, 'unstable.mat');
picfile = strcat(resultdir, 'unstable.pdf');

% Physical free parameters.
q = 1.7; theta = 2.4;

% Grid resolution;
N = 1e4;

% Guess for eigenvalue.
guess = [3.3e-02; 6.7e-01];
pert_guess = [0.0; 0.0; 0.0; 0.0];

[params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess, pert_guess);
save(matfile);

% Plotting part.
plot_znd_and_perturbations(matfile, picfile);