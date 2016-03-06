% Test find_roots with steady-piston rear boundary condition.
% Case is unstable ZND solution.
close all; clear all; clc

resultdir = 'results/2016-03-06-steady-piston-bc-fsolve/';
matfile = strcat(resultdir, 'unstable.mat');
picfile = strcat(resultdir, 'unstable.pdf');

% Physical free parameters.
q = 1.7; theta = 2.4;

% Grid resolution;
N = 10000;

% Guess for eigenvalue.
guess = [3.3e-02; 6.7e-01];

[params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess);
save(matfile);

% Plotting part.
plot_znd_and_perturbations(matfile, picfile);