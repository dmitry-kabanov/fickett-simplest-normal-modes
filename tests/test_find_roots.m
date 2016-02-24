close all
clear all

% Physical free parameters.
q = 1.7; theta = 2.4;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% Resolution.
N = 10000;

% Lambda tolerance.
lambda_tol = 1.0 / N;

M = compute_x_at_point(1-lambda_tol, params);

tol = 1e-12;

grid = linspace(0, M, N);

znd_all = compute_znd(grid, params);

alpha_re_guess = 3.2689475021629215e-02;
alpha_im_guess = 6.6361473999774434e-01;
[root, norm_f] = find_roots([alpha_re_guess; alpha_im_guess], grid, znd_all, params, tol);