close all
clear all
clc

% Physical free parameters.
q = 1.7; theta = 2.4;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% Resolution.
N = 1e5;

% Lambda tolerance.
lambda_tol = 1.0 / N;

% Left boundary of the computational domain.
M = compute_x_at_point(1-lambda_tol, params);

grid = linspace(0, M, N);

znd_all = compute_znd(grid, params);

alpha_re = 3.2689475021629215e-02;
alpha_im = 6.6361473999774434e-01;
alpha = alpha_re + 1j * alpha_im;
H = compute_boundedness_function(alpha, grid, znd_all, params);
