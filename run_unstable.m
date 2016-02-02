% Compute growth rate and frequencies for the physical parameters
% that are known from the direct linear simulationsto to give unstable
% ZND solution.
close all
clear all

% Physical free parameters.
q = 1.7; theta = 2.4;

% Domain length.
M = 13;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% Carpet search parameters: lb - lower bound, ub - upper bound, n - resolution.
cp.lb_re = -5e-2;
cp.ub_re = 5e-2;
cp.lb_im = 0;
cp.ub_im = 1;
cp.n_re = 10;
cp.n_im = 10;


[alpha_re, alpha_im, H] = compute_carpet(params, M, cp);
plot_carpet(alpha_re, alpha_im, H);

% Find indices of minimum on the carpet.
[xmin, ymin] = find_carpet_minimum(H);

% Make initial guess.
guess.alpha_re = alpha_re(xmin, ymin);
guess.alpha_im = alpha_im(xmin, ymin);

H(xmin, ymin)
guess

% % Find minimum of $abs(H)$.
% minimize_boundedness_function(guess, M, params);
