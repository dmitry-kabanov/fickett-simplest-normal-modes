close all
clear

% Physical free parameters.
q = 1.7; theta = 2.4;

% Domain length. This domain length corresponds to the lambda_tol=1e-6.
M = 12.317732705436374;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

alpha_re = 3.2689475021629215e-02;
alpha_im = 6.6361473999774434e-01;
alpha = alpha_re + 1j * alpha_im;
H = compute_boundedness_function(alpha, M, params);