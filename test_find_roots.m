close all
clear variables

% Physical free parameters.
q = 1.7; theta = 2.4;

% Lambda tolerance.
lambda_tol = 1e-3;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

alpha = [3.269e-2; 6.636e-1];

%alpha_re = 3.2689475021629215e-02;
%alpha_im = 6.6361473999774434e-01;
%alpha = alpha_re + 1j * alpha_im;
H = find_roots(alpha, lambda_tol, params);