% Test find_roots with steady-piston rear boundary condition.
% Case is unstable ZND solution.
close all; clear all; clc

resultdir = 'results/2016-03-24-new-radiation-condition-tol=1e-6/';
assert(exist(resultdir, 'dir') == 7);

picfilename = strcat(resultdir, 'unstable');
matfile = strcat(picfilename, '.mat');

% Physical free parameters.
q = 1.7; theta = 2.4;

% Grid resolution;
N = 1000000;

% Guess for eigenvalue.
guess.alpha_re = 0.032711279984620;
guess.alpha_im = 0.663600985237249;

[params, grid, znd_all, result, pert] = solve_eigenvalue_problem(q, theta, N, guess);
save(matfile);

% Plotting part.
plot_znd_and_perturbations(picfilename);