close all
clear all

% Physical free parameters.
q = 1.7; theta = 2.4;

% Domain length.
M = 13;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

znd_sol = get_znd_sol(M, params);

alpha = 3.269e-2 + 1j * 6.636e-1;
H = compute_boundedness_function(alpha, M, params, znd_sol);