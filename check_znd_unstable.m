% Check that we compute ZND structure correctly when we use function
% `compute_znd_data_at_point`.
close all
clear all
clc

% Physical free parameters.
q = 1.7; theta = 2.4;

% Struct with free and dependent parameters.
params = compute_aux_params(q, theta);

% Lambda tolerance, that is, closeness of lambda to unity.
lambda_tol = 1e-6;

N = 10000;

lambda = linspace(0, 1-lambda_tol, N)';
u      = zeros(N, 1);
du_dx  = zeros(N, 1);
du_dl  = zeros(N, 1);
dl_dx  = zeros(N, 1);
dw_du  = zeros(N, 1);
dw_dl  = zeros(N, 1);
w      = zeros(N, 1);
x      = zeros(N, 1);

for i = 1:length(lambda)
    znd = compute_znd_data_at_point(lambda(i), params);

    u(i) = znd.u;
    du_dl(i) = znd.du_dl;
    du_dx(i) = znd.du_dx;
    dl_dx(i) = znd.dl_dx;
    dw_du(i) = znd.dw_du;
    dw_dl(i) = znd.dw_dl;
    w(i) = znd.w;
    x(i) = compute_x_at_point(lambda(i), params);
end
