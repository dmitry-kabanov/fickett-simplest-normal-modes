% Check that we compute ZND structure correctly when we use function
% `compute_znd_data_at_point`.
close all;
clear;

% Physical free parameters.
q = 1.7; theta = 2.4;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% Lambda tolerance, that is, closeness of lambda to unity.
lambda_tol = 1e-4;

lambda = linspace(0, 1-lambda_tol, 1000);
u = zeros('like', lambda);
du_dx = zeros('like', lambda);
du_dl = zeros('like', lambda);
dl_dx = zeros('like', lambda);
dw_du = zeros('like', lambda);
dw_dl = zeros('like', lambda);
w = zeros('like', lambda);

for i = 1:length(lambda)
    znd = compute_znd_data_at_point(lambda(i), params);

    u(i) = znd.u;
    du_dl(i) = znd.du_dl;
    dw_du(i) = znd.dw_du;
    dw_dl(i) = znd.dw_dl;
    w(i) = znd.w;
end
